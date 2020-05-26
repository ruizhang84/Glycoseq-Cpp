#define BOOST_TEST_MODULE SearchTest
#include <boost/test/unit_test.hpp>

#include <iostream>
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

    // read fasta and build peptides
    util::io::FASTAReader fasta_reader("/home/yu/Documents/MultiGlycan-Cpp/data/haptoglobin.fasta");
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
 
    engine::protein::Digestion digest;
    digest.SetProtease(engine::protein::Proteases::Trypsin);
    std::vector<std::string> seqs = digest.Sequences(proteins.front().Sequence(),
         engine::protein::FindPTM::ContainsNGlycanSite);
    digest.SetProtease(engine::protein::Proteases::GluC);
    std::vector<std::string> peptides;
    for(auto& it : seqs)
    {
        std::vector<std::string> seq = digest.Sequences(it,
         engine::protein::FindPTM::ContainsNGlycanSite);
        peptides.insert(peptides.end(), seq.begin(), seq.end());
    }

    // // build glycans
    engine::glycan::GlycanBuilder builder(12, 12, 5, 4, 0);
    builder.Build();

    engine::glycan::GlycanStore store = builder.Isomer();
    std::vector<std::string> glycans_str = store.Collection();

    auto start = std::chrono::high_resolution_clock::now(); 
    // spectrum matching
    PrecursorMatcher precursor_runner(0.01, algorithm::search::ToleranceBy::Dalton);
    precursor_runner.Init(peptides, glycans_str);
    for(auto& spec : spectrum_reader.GetSpectrum())
    {
        double target = util::mass::SpectrumMass::Compute(spec.PrecursorMZ(), spec.PrecursorCharge());
        std::vector<SearchResult> r = precursor_runner.Match(target);
        if (! r.empty())
        std::cout << spec.Scan() << " : " << std::endl;
        for(auto it : r)
        {
            std::cout << it.glycan << std::endl;
            std::cout << it.peptide << std::endl;
        }
    }
    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << duration.count() << std::endl; 


}

} // namespace search
} // namespace engine