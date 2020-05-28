#define BOOST_TEST_MODULE SearchTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <iomanip>
#include "spectrum_search.h"
#include "precursor_match.h"
#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "../protein/protein_digest.h"
#include "../protein/protein_ptm.h"
#include "../glycan/glycan_builder.h"
#include <chrono> 

namespace engine{
namespace search {

BOOST_AUTO_TEST_CASE( precusor_match_test ) 
{
    // read spectrum
    std::string path = "/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf";
    std::unique_ptr<util::io::SpectrumParser> parser = 
        std::make_unique<util::io::MGFParser>(path, util::io::SpectrumType::EThcD);
    util::io::SpectrumReader spectrum_reader(path, std::move(parser));
    spectrum_reader.Init();
    BOOST_CHECK(spectrum_reader.GetSpectrum(8030).PrecursorMZ() == 1010.6698);


    // read fasta and build peptides
    util::io::FASTAReader fasta_reader("/home/yu/Documents/MultiGlycan-Cpp/data/haptoglobin.fasta");
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
 
    engine::protein::Digestion digest;
    digest.SetProtease(engine::protein::Proteases::Trypsin);
    std::unordered_set<std::string> seqs = digest.Sequences(proteins.front().Sequence(),
         engine::protein::ProteinPTM::ContainsNGlycanSite);
    digest.SetProtease(engine::protein::Proteases::GluC);
    std::vector<std::string> peptides;
    for(auto& it : seqs)
    {
        std::unordered_set<std::string> seq = digest.Sequences(it,
         engine::protein::ProteinPTM::ContainsNGlycanSite);
        peptides.insert(peptides.end(), seq.begin(), seq.end());
    }
    BOOST_CHECK(std::find(peptides.begin(), peptides.end(), "MVSHHNLTTGATLINE") != peptides.end());

    // // build glycans
    engine::glycan::GlycanBuilder builder(12, 12, 5, 4, 0);
    builder.Build();
    model::glycan::NGlycanComplex glycan;
    std::map<model::glycan::Monosaccharide, int> composite;
    composite[model::glycan::Monosaccharide::GlcNAc] = 5;
    composite[model::glycan::Monosaccharide::Man] = 3;
    composite[model::glycan::Monosaccharide::Gal] = 3;
    composite[model::glycan::Monosaccharide::Fuc] = 1;
    composite[model::glycan::Monosaccharide::NeuAc] = 1;  
    glycan.set_composition(composite);
    std::string glycan_name = glycan.Name();
    BOOST_CHECK(builder.Isomer().QueryMass(glycan_name) == util::mass::GlycanMass::Compute(composite));
    BOOST_CHECK(util::mass::GlycanMass::Compute(composite) == 
        util::mass::GlycanMass::Compute(model::glycan::Glycan::Interpret(glycan_name)));
    std::vector<std::string> collection = builder.Isomer().Collection();
    BOOST_CHECK(std::find(collection.begin(), collection.end(), glycan.Name()) != collection.end());

    // spectrum matching
    PrecursorMatcher precursor_runner(10, algorithm::search::ToleranceBy::PPM, builder.Isomer());
    std::vector<std::string> glycans_str = builder.Isomer().Collection();

    precursor_runner.Init(peptides, glycans_str);
    SpectrumSearcher spectrum_runner(0.01, algorithm::search::ToleranceBy::Dalton, builder.Mass(), builder.Isomer());

    auto special_spec = spectrum_reader.GetSpectrum(8172);
    double special_target = util::mass::SpectrumMass::Compute(special_spec.PrecursorMZ(), special_spec.PrecursorCharge());
    MatchResultStore special_r = precursor_runner.Match(special_target, 2);    
    // std::cout << special_spec.Scan() << " : " << std::endl;
    // for(auto it : special_r.Map())
    // {
    //     std::cout << it.first << std::endl;
    //     for(auto g: it.second)
    //     {
    //         std::cout << g << std::endl;
    //     }
    // }
    BOOST_CHECK(!special_r.Empty());
    

    // special_spec = spectrum_reader.GetSpectrum(8001);
    // special_target = util::mass::SpectrumMass::Compute(special_spec.PrecursorMZ(), special_spec.PrecursorCharge());
    // special_r = precursor_runner.Match(special_target, 2);
    // spectrum_runner.set_candidate(special_r);
    // spectrum_runner.set_spectrum(special_spec);
    // std::vector<SearchResult> special_res = spectrum_runner.Search();
    // std::string special_string = "PEPTIDE";
    // for (int i = 0; i < (int) special_string.length()-1; i++)
    // {
    //     std::cout << util::mass::IonMass::Compute(special_string.substr(0, i+1), 
    //         util::mass::IonType::c) << std::endl;
    // }
    // for (int i = 0; i < (int) special_string.length(); i++)
    // {
    //     std::cout << util::mass::IonMass::Compute(special_string.substr(i, special_string.length() - i), 
    //         util::mass::IonType::z) << std::endl;
    // }   
    // for(auto it : special_res)
    // {
    //     std::cout << it.glycan << std::endl;
    //     std::cout << it.peptide << std::endl;
    //     std::cout << it.score << std::endl;
    // }
    // BOOST_CHECK(!special_res.empty());
    

    std::cout << "Start to scan\n"; 
    auto start = std::chrono::high_resolution_clock::now(); 
    for(auto& spec : spectrum_reader.GetSpectrum())
    {
        double target = util::mass::SpectrumMass::Compute(spec.PrecursorMZ(), spec.PrecursorCharge());
        MatchResultStore r = precursor_runner.Match(target, 2);
        if (r.Empty()) continue;

        spectrum_runner.set_spectrum(spec);
        spectrum_runner.set_candidate(r);
        std::vector<SearchResult> results = spectrum_runner.Search();
        if (results.empty()) continue;
        std::cout << "Scan " << spec.Scan() << ":" << std::endl;

        for(auto it : results)
        {
            std::cout << it.glycan << std::endl;
            std::cout << it.peptide << std::endl;
            std::cout << it.score << std::endl;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

}

} // namespace search
} // namespace engine