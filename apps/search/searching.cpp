#include <algorithm>
#include <fstream>
#include <iostream>
#include <thread>  
#include <mutex> 
#include <chrono> 

#include "search_parameter.h"
#include "search_dispatcher.h"
#include "search_helper.h"

#include "../../util/io/mgf_parser.h"
#include "../../util/io/fasta_reader.h"
#include "../../engine/protein/protein_digest.h"
#include "../../engine/protein/protein_ptm.h"
#include "../../engine/glycan/glycan_builder.h"
#include "../../engine/spectrum/normalize.h"
#include "../../engine/search/precursor_match.h"
#include "../../engine/search/spectrum_search.h"
#include "../../engine/search/search_result.h"
#include "../../engine/search/fdr_filter.h"


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
    std::unordered_set<std::string> seqs = PeptidesDigestion(fasta_path);
    for(const auto& s : seqs)
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
    SearchDispatcher target_searcher(spectrum_reader.GetSpectrum(), builder.get(), peptides, parameter);
    std::vector<engine::search::SearchResult> targets = target_searcher.Dispatch();

    // seraching decoys 
    SearchDispatcher decoy_searcher(spectrum_reader.GetSpectrum(), builder.get(), decoy_peptides, parameter);
    std::vector<engine::search::SearchResult> decoys = decoy_searcher.Dispatch();

    // set up scorer
    std::thread scorer_first(ScoringWorker, std::ref(targets), parameter.weights);
    std::thread scorer_second(ScoringWorker, std::ref(decoys), parameter.weights);   
    scorer_first.join();
    scorer_second.join();

    // remove the lower score
    targets = ScoreFilter(targets);
    decoys = ScoreFilter(decoys);

    // fdr filtering
    engine::search::FDRFilter fdr_runner(parameter.fdr_rate);
    fdr_runner.set_data(targets, decoys);
    fdr_runner.Init();

    std::vector<engine::search::SearchResult> results = fdr_runner.Filter();

    // output analysis results
    ReportResults(out_path, results);
   
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

    auto stop = std::chrono::high_resolution_clock::now(); 
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start); 
    std::cout << "Total Time: " << duration.count() << std::endl; 

}