#include <fstream>
#include <iostream>
#include <chrono> 

#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "../../engine/protein/protein_digest.h"
#include "../../engine/protein/protein_ptm.h"
#include "../../engine/glycan/glycan_builder.h"
#include "../../engine/spectrum/normalize.h"
#include "../../engine/search/precursor_match.h"
#include "../../engine/search/spectrum_search.h"
#include "../../engine/search/fdr_filter.h"

int main(int argc, char *argv[])
{
    std::string out = "search.csv";

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

    // // build glycans
    engine::glycan::GlycanBuilder builder(12, 12, 5, 4, 0);
    builder.Build();
    

    // search setup
    engine::search::PrecursorMatcher precursor_runner(10, algorithm::search::ToleranceBy::PPM, builder.Isomer());
    std::vector<std::string> glycans_str = builder.Isomer().Collection();
    precursor_runner.Init(peptides, glycans_str);
    engine::search::SpectrumSearcher spectrum_runner(0.01, algorithm::search::ToleranceBy::Dalton, builder.Mass(), builder.Isomer());

    std::cout << "Start to scan\n"; 
    auto start = std::chrono::high_resolution_clock::now();

    // process spectrum by normalization
    for(auto& spec : spectrum_reader.GetSpectrum())
    {
       engine::spectrum::Normalizer::Transform(spec);
    }

    // seraching targets 
    std::vector<engine::search::SearchResult> targets;
    for(auto& spec : spectrum_reader.GetSpectrum())
    {
        double target = util::mass::SpectrumMass::Compute(spec.PrecursorMZ(), spec.PrecursorCharge());
        engine::search::MatchResultStore r = precursor_runner.Match(target, 2);
        if (r.Empty()) continue;

        spectrum_runner.set_spectrum(spec);
        spectrum_runner.set_candidate(r);
        std::vector<engine::search::SearchResult> res = spectrum_runner.Search();
        if (res.empty()) continue;

        targets.insert(targets.end(), res.begin(), res.end());
    }
    
    // seraching decoys 
    std::vector<engine::search::SearchResult> decoys;
    double pseudo_mass = 50;
    for(auto& spec : spectrum_reader.GetSpectrum())
    {
        double target = util::mass::SpectrumMass::Compute(spec.PrecursorMZ(), spec.PrecursorCharge());
        engine::search::MatchResultStore r = precursor_runner.Match(target + pseudo_mass, 2);
        if (r.Empty()) continue;

        spectrum_runner.set_spectrum(spec);
        spectrum_runner.set_candidate(r);
        std::vector<engine::search::SearchResult> res = spectrum_runner.Search();
        if (res.empty()) continue;

        decoys.insert(decoys.end(), res.begin(), res.end());
    }

    // fdr filtering
    engine::search::FDRFilter fdr_runner(0.01);
    fdr_runner.set_data(targets);
    fdr_runner.set_decoy(decoys);
    fdr_runner.Init();

    std::vector<engine::search::SearchResult> results = fdr_runner.Filter();

    std::ofstream outfile;
    outfile.open (out);
    outfile << "scan#,peptide,glycan,score\n";
    
    for(auto it : results)
    {
        outfile << it.scan << ",";
        outfile << it.peptide << ",";
        outfile << it.glycan << ",";
        outfile << it.score << "\n";
    }
    outfile.close();

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

}