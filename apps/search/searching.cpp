#include <fstream>
#include <iostream>
#include <thread>  
#include <mutex> 
#include <chrono> 

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


std::mutex mutex; 

void SearchingWorker(
    std::vector<engine::search::SearchResult>& results, 
    engine::glycan::NGlycanBuilder* builder,
    std::vector<model::spectrum::Spectrum> spectra,
    std::vector<std::string> peptides, 
    double ms1_tol, algorithm::search::ToleranceBy ms1_by,
    double ms2_tol, algorithm::search::ToleranceBy ms2_by,
    int isotopic_count, double pseudo_mass)
{
    engine::search::PrecursorMatcher precursor_runner(ms1_tol, ms1_by, builder->Isomer());
    engine::search::SpectrumSearcher spectrum_runner(ms2_tol, ms2_by, builder);
    std::vector<std::string> glycans_str = builder->Isomer().Collection();
    precursor_runner.Init(peptides, glycans_str);
    spectrum_runner.Init();

    std::vector<engine::search::SearchResult> temp_result;
    for(auto& spec : spectra)
    {
        // precusor
        double target = util::mass::SpectrumMass::Compute(spec.PrecursorMZ(), spec.PrecursorCharge());
        engine::search::MatchResultStore r = precursor_runner.Match(target + pseudo_mass, spec.PrecursorCharge(), isotopic_count);
        if (r.Empty()) continue;

        // process spectrum by normalization
        engine::spectrum::Normalizer::Transform(spec);

        // msms
        spectrum_runner.set_spectrum(spec);
        spectrum_runner.set_candidate(r);
        std::vector<engine::search::SearchResult> res = spectrum_runner.Search();
        if (res.empty()) continue;
        temp_result.insert(temp_result.end(), res.begin(), res.end());
    }

    mutex.lock();
        results.insert(results.end(), temp_result.begin(), temp_result.end());
    mutex.unlock();
}

void ScoringWorker(
    std::vector<engine::search::SearchResult>& results,
    std::unordered_map<engine::search::SearchType, double> parameter)
{
    engine::search::SimpleScorer scorer(parameter);
    for(auto it : results)
    {
        double score = scorer.ComputeScore(it);
        it.set_score(score);
    }
}

int main(int argc, char *argv[])
{
    std::string spectra_path = "/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf";
    std::string fasta_path = "/home/yu/Documents/MultiGlycan-Cpp/data/haptoglobin.fasta";
    std::string out_path = "search2.csv";
    int n_thread = 6;
    int hexNAc_upper_bound = 12;
    int hex_upper_bound = 12;
    int fuc_upper_bound = 5;
    int neuAc_upper_bound = 4;
    int neuGc_upper_bound = 0;
    double ms1_tol = 10;
    algorithm::search::ToleranceBy ms1_by = algorithm::search::ToleranceBy::PPM;
    double ms2_tol = 0.01;
    algorithm::search::ToleranceBy ms2_by = algorithm::search::ToleranceBy::Dalton;
    int isotopic_count = 2;
    double pseudo_mass = 50;
    double fdr_rate = 1.01;
    double core = 1.0;
    double branch = 1.0;
    double terminal = 1.0;
    double peptide = 1.0;
    double oxonium = 1.0;

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
        std::make_unique<engine::glycan::NGlycanBuilder>(hexNAc_upper_bound, 
        hex_upper_bound, fuc_upper_bound, neuAc_upper_bound, neuGc_upper_bound);
    builder->Build();
    

    // search
    std::cout << "Start to scan\n"; 
    auto start = std::chrono::high_resolution_clock::now();

    // seraching targets 
    std::vector<engine::search::SearchResult> targets;
    std::vector<model::spectrum::Spectrum> spectra = spectrum_reader.GetSpectrum();
    int work_load = std::ceil(spectra.size() / n_thread);
    std::vector< std::thread> target_thread_pool;
    for (int i = 0; i < n_thread; i ++)
    {
        int begin = work_load * i;
        int end = begin + work_load > (int) spectra.size() ? 
            (int) spectra.size() : begin + work_load;

        std::thread worker(SearchingWorker, std::ref(targets), builder.get(), 
            std::vector<model::spectrum::Spectrum>(spectra.begin() + begin, spectra.begin() + end),  
                peptides, ms1_tol, ms1_by, ms2_tol, ms2_by, isotopic_count, 0);
        target_thread_pool.push_back(std::move(worker));
    }
    for (auto& worker : target_thread_pool)
    {
        worker.join();
    }

    // seraching decoys 
    std::vector<engine::search::SearchResult> decoys;
    std::vector< std::thread> decoy_thread_pool;
    spectra = spectrum_reader.GetSpectrum();
    for (int i = 0; i < n_thread; i ++)
    {
        int begin = work_load * i;
        int end = begin + work_load > (int) spectra.size() ? 
            (int) spectra.size() : begin + work_load;
        
        std::thread worker(SearchingWorker, std::ref(decoys), builder.get(), 
            std::vector<model::spectrum::Spectrum>(spectra.begin() + begin, spectra.begin() + end), 
                peptides, ms1_tol, ms1_by, ms2_tol, ms2_by, isotopic_count, pseudo_mass);
        decoy_thread_pool.push_back(std::move(worker));
    }
    for (auto& worker : decoy_thread_pool)
    {
        worker.join();
    }


    // set up scorer
    const std::unordered_map<engine::search::SearchType, double> parameter 
    {
        {engine::search::SearchType::Core, core}, 
        {engine::search::SearchType::Branch, branch}, 
        {engine::search::SearchType::Terminal, terminal},
        {engine::search::SearchType::Peptide, peptide}, 
        {engine::search::SearchType::Oxonium, oxonium},
    };
    std::thread scorer_first(ScoringWorker, std::ref(targets), parameter);
    std::thread scorer_second(ScoringWorker, std::ref(decoys), parameter);   
    scorer_first.join();
    scorer_second.join();

    // fdr filtering
    engine::search::FDRFilter fdr_runner(fdr_rate);
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