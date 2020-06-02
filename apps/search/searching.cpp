#include <algorithm>
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
    for(auto& it : results)
    {
        int match = analyzer.Predicting(it);
        std::vector<double> prob = analyzer.PredictingProbability(it);
        if(match > 0)
        {
            it.set_score(std::max(prob[0], prob[1]));
        }
        else
        {
            it.set_score(std::min(prob[0], prob[1]));
        }        
    }
}

// void ScoringWorker(
//     std::vector<engine::search::SearchResult>& results,
//     std::unordered_map<engine::search::SearchType, double> weights)
// {
//     engine::search::SimpleScorer scorer(weights);
//     for(auto& it : results)
//     {
//         double score = scorer.ComputeScore(it);
//         it.set_score(score);
//     }
// }

std::unordered_set<std::string> GeneratePeptdies(const std::string& fasta_path)
{
    util::io::FASTAReader fasta_reader(fasta_path);
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
   
    engine::protein::Digestion digest;
    digest.SetProtease(engine::protein::Proteases::Trypsin);
    std::unordered_set<std::string> seqs = digest.Sequences(proteins.front().Sequence(),
         engine::protein::ProteinPTM::ContainsNGlycanSite);
    digest.SetProtease(engine::protein::Proteases::GluC);
    std::unordered_set<std::string> double_seqs;
    for(auto& it : seqs)
    {
        std::unordered_set<std::string> seq = digest.Sequences(it,
            engine::protein::ProteinPTM::ContainsNGlycanSite);
        double_seqs.insert(seq.begin(), seq.end());
    }
    return double_seqs;
}

int main(int argc, char *argv[])
{
    std::string spectra_path = "/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf";
    std::string fasta_path = "/home/yu/Documents/MultiGlycan-Cpp/data/haptoglobin.fasta";
    std::string out_path = "search3.csv";
    SearchParameter parameter;

    // read spectrum
    std::unique_ptr<util::io::SpectrumParser> parser = 
        std::make_unique<util::io::MGFParser>(spectra_path, util::io::SpectrumType::EThcD);
    util::io::SpectrumReader spectrum_reader(spectra_path, std::move(parser));
    spectrum_reader.Init();

    // read fasta and build peptides
    std::vector<std::string> peptides, decoy_peptides;
    std::unordered_set<std::string> double_seqs = GeneratePeptdies(fasta_path);
    for(const auto& s : double_seqs)
    {
        std::string decoy_s(s);
        peptides.push_back(s);
        std::reverse(decoy_s.begin(), decoy_s.end());
        decoy_peptides.push_back(decoy_s);
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

    // seraching targets 
    parameter.pseudo_mass = 0;
    SearchDispatcher target_searcher(spectrum_reader.GetSpectrum(), builder.get(), peptides, parameter);
    std::vector<engine::search::SearchResult> targets = target_searcher.Dispatch();

    // seraching decoys 
    SearchDispatcher decoy_searcher(spectrum_reader.GetSpectrum(), builder.get(), decoy_peptides, parameter);
    std::vector<engine::search::SearchResult> decoys = decoy_searcher.Dispatch();

    // set up scorer
    engine::analysis::SVMAnalyzer analyzer;
    analyzer.Training(targets, decoys);

    std::thread scorer_first(ScoringWorker, std::ref(targets), std::ref(analyzer));
    std::thread scorer_second(ScoringWorker, std::ref(decoys), std::ref(analyzer));   
    scorer_first.join();
    scorer_second.join();
    // const std::unordered_map<engine::search::SearchType, double> weights 
    // {
    //     {engine::search::SearchType::Core, 1.0}, 
    //     {engine::search::SearchType::Branch, 1.0}, 
    //     {engine::search::SearchType::Terminal, 1.0},
    //     {engine::search::SearchType::Peptide, 1.0}, 
    //     {engine::search::SearchType::Oxonium, 1.0},
    // };
    // std::thread scorer_first(ScoringWorker, std::ref(targets), weights);
    // std::thread scorer_second(ScoringWorker, std::ref(decoys), weights);   
    // scorer_first.join();
    // scorer_second.join();

    // fdr filtering
    engine::search::FDRFilter fdr_runner(parameter.fdr_rate);
    fdr_runner.set_data(targets, decoys);
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
   
    // for(auto it : targets)
    // {
    //     outfile << it.Scan() << ",";
    //     outfile << it.Sequence() << ",";
    //     outfile << it.Glycan() << ",";
    //     outfile << it.Score() << ",target\n";
    // }
   

    // for(auto it : decoys)
    // {
    //     outfile << it.Scan() << ",";
    //     outfile << it.Sequence() << ",";
    //     outfile << it.Glycan() << ",";
    //     outfile << it.Score() << ",decoy\n";
    // }
   
    // for(auto i : targets)
    // {
    //     for(auto j : decoys)
    //     {
    //         if (i.Scan() == j.Scan())
    //         {
    //             std::cout << i.Scan() << std::endl;
    //             std::cout << i.Sequence() << std::endl;
    //             std::cout << i.Glycan() << std::endl;
    //             std::cout << j.Sequence() << std::endl;
    //             std::cout << j.Glycan() << std::endl;
    //         }
               
    //     }
    // }
    
    outfile.close();

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

}