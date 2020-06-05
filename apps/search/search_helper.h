#include <map>
#include <unordered_map>
#include <fstream>

#include "../../util/io/fasta_reader.h"
#include "../../engine/protein/protein_digest.h"
#include "../../engine/protein/protein_ptm.h"
#include "../../engine/search/search_result.h"
#include "../../engine/analysis/svm_analyzer.h"

// generate peptides by digestion
std::unordered_set<std::string> PeptidesDigestion
    (const std::string& fasta_path, SearchParameter parameter)
{
    util::io::FASTAReader fasta_reader(fasta_path);
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
   
    engine::protein::Digestion digest;
    digest.set_miss_cleavage(parameter.miss_cleavage);
    std::unordered_set<std::string> peptides;

    // digestion
    std::deque<engine::protein::Proteases> proteases(parameter.proteases);
    engine::protein::Proteases enzyme = proteases.front();
    digest.SetProtease(enzyme);
    proteases.pop_front();
    for(const auto& protein : proteins)
    {
        std::unordered_set<std::string> seqs = digest.Sequences(protein.Sequence(),
            engine::protein::ProteinPTM::ContainsNGlycanSite);
        peptides.insert(seqs.begin(), seqs.end());
    }
        
    // double digestion or more
    while (proteases.size() > 0)
    {
        std::unordered_set<std::string> double_seqs;
        enzyme = proteases.front();
        digest.SetProtease(enzyme);
        proteases.pop_front();
        for(const auto& seq : peptides)
        {
            std::unordered_set<std::string> seqs = digest.Sequences(seq,
                engine::protein::ProteinPTM::ContainsNGlycanSite);
            double_seqs.insert(seqs.begin(), seqs.end());
        }
        peptides.insert(double_seqs.begin(), double_seqs.end());
    }   

    return peptides;
}

// assign score to searching results
void ScoringWorker(
    std::vector<engine::search::SearchResult>& results,
    std::map<engine::search::SearchType, double> weights)
{
    engine::search::SimpleScorer scorer(weights);
    for(auto& it : results)
    {
        double score = scorer.ComputeScore(it);
        it.set_score(score);
    }

    // compute coelution of peptide sequence
    std::unordered_map<std::string, int> count;
    for(auto i : results)
    {
        if (count.find(i.Sequence()) == count.end())
        {
            count[i.Sequence()] = 0;
        }
        count[i.Sequence()] += 1;
    }
    std::vector<std::unordered_map<std::string, int>> bucket_count;
    int bucket =  50;
    for(int i = 0; i < bucket; i ++)
    {
        std::unordered_map<std::string, int> temp;
        bucket_count.push_back(temp);
    }
    for(auto i : results)
    {
        int scan = i.Scan();
        int index = scan / 500;
        std::unordered_map<std::string, int>& temp = bucket_count[index];
        
        if (temp.find(i.Sequence()) == temp.end())
        {
            temp[i.Sequence()] = 0;
        }
        temp[i.Sequence()] += 1;
    }
    for(auto& it : results)
    {
        std::string s = it.Sequence();
        int scan = it.Scan();
        int index = scan / 500;
        int counts = bucket_count[index][s];
        if (index > 0)
        {
            std::unordered_map<std::string, int> temp = bucket_count[index-1];
            if (temp.find(s) != temp.end())
                counts += temp[s];
        }

        std::unordered_map<std::string, int> temp = bucket_count[index+1];
        if (temp.find(s) != temp.end())
            counts += temp[s];
          
        
        double score = counts * 1.0 / count[s];
        it.set_score(it.Score() + score);
    }
}

void SVMScoringWorker(
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

// process score, removing lower score
std::vector<engine::search::SearchResult> ScoreFilter
    (const std::vector<engine::search::SearchResult>& results)
{
    std::vector<engine::search::SearchResult> res;
    std::map<int, std::vector<engine::search::SearchResult>> ranker;
    for(auto it : results)
    {
        int scan = it.Scan();
        if (ranker.find(scan) == ranker.end())
        {
            ranker[scan] = std::vector<engine::search::SearchResult>();
        }
        else
        {
            double score = ranker[scan].front().Score();
            if (it.Score() < score) continue;
            else if (it.Score() > score) ranker[scan].clear();
        }
        ranker[scan].push_back(it);
    }
    for(auto it : ranker)
    {
        std::vector<engine::search::SearchResult>& r = it.second;
        res.insert(res.end(), r.begin(), r.end());
    }
    return res;
}

// report glycopeptide identification of spectrum
void ReportResults(const std::string& out_path,
    const std::vector<engine::search::SearchResult>&  results)
{
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
}
