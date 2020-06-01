#include <fstream>
#include <iostream>
#include <thread>  
#include <mutex> 
#include <chrono> 

#include "search_parameter.h"
#include "search_dispatcher.h"

#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "../../engine/protein/protein_digest.h"
#include "../../engine/protein/protein_ptm.h"
#include "../../engine/glycan/glycan_builder.h"
#include "../../engine/spectrum/normalize.h"
#include "../../engine/search/precursor_match.h"
#include "../../engine/search/spectrum_search.h"
#include "../../engine/search/search_result.h"
#include "../../engine/analysis/svm_analyzer.h"
#include "../../engine/search/fdr_filter.h"



void ScoringWorker(
    std::vector<engine::search::SearchResult>& results,
    engine::analysis::SVMAnalyzer& analyzer)
{
    for(auto it : results)
    {
        double prob = analyzer.PredictingProbability(it);
        it.set_score(prob);
    }
}

int main(int argc, char *argv[])
{
    std::string spectra_path = "/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf";
    std::string fasta_path = "/home/yu/Documents/MultiGlycan-Cpp/data/haptoglobin.fasta";
    std::string out_path = "search.csv";
    SearchParameter parameter;

    // read spectrum
    std::unique_ptr<util::io::SpectrumParser> parser = 
        std::make_unique<util::io::MGFParser>(spectra_path, util::io::SpectrumType::EThcD);
    util::io::SpectrumReader spectrum_reader(spectra_path, std::move(parser));
    spectrum_reader.Init();


    // read fasta and build peptides
    util::io::FASTAReader fasta_reader(fasta_path);
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
    std::unique_ptr<engine::glycan::NGlycanBuilder> builder =
        std::make_unique<engine::glycan::NGlycanBuilder>(parameter.hexNAc_upper_bound, 
            parameter.hex_upper_bound, parameter.fuc_upper_bound, 
                parameter.neuAc_upper_bound, parameter.neuGc_upper_bound);
    builder->Build();
    

    // search
    std::cout << "Start to scan\n"; 
    auto start = std::chrono::high_resolution_clock::now();

    // seraching decoys 
    SearchDispatcher decoy_searcher(spectrum_reader.GetSpectrum(4240, 4500), builder.get(), peptides, parameter);
    std::vector<engine::search::SearchResult> decoys = decoy_searcher.Dispatch();

    // seraching targets 
    parameter.pseudo_mass = 0;
    SearchDispatcher target_searcher(spectrum_reader.GetSpectrum(4240, 4500), builder.get(), peptides, parameter);
    std::vector<engine::search::SearchResult> targets = target_searcher.Dispatch();

    // set up scorer
    engine::analysis::SVMAnalyzer analyzer;
    analyzer.Training(targets, decoys);

    std::thread scorer_first(ScoringWorker, std::ref(targets), std::ref(analyzer));
    std::thread scorer_second(ScoringWorker, std::ref(decoys), std::ref(analyzer));   
    scorer_first.join();
    scorer_second.join();

    // fdr filtering
    engine::search::FDRFilter fdr_runner(parameter.fdr_rate);
    fdr_runner.set_data(targets);
    fdr_runner.set_decoy(decoys);
    fdr_runner.Init();

    std::vector<engine::search::SearchResult> results = fdr_runner.Filter();

    // output analysis results
    std::ofstream outfile;
    outfile.open (out_path);
    outfile << "scan#,peptide,glycan,score\n";
    
    for(auto it : results)
    {
        outfile << it.Scan() << ",";
        outfile << it.Sequence() << ",";
        outfile << it.Glycan() << ",";
        outfile << it.Score() << "\n";
    }
    
    outfile.close();

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

}