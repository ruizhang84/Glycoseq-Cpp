#include <map>
#include <fstream>

#include "../../util/io/fasta_reader.h"
#include "../../engine/protein/protein_digest.h"
#include "../../engine/protein/protein_ptm.h"
#include "../../engine/search/search_result.h"
#include "../../engine/analysis/svm_analyzer.h"

// generate peptides by digestion
std::unordered_set<std::string> PeptidesDigestion(const std::string& fasta_path)
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

// assign score to searching results
void ScoringWorker(
    std::vector<engine::search::SearchResult>& results,
    std::unordered_map<engine::search::SearchType, double> weights)
{
    engine::search::SimpleScorer scorer(weights);
    for(auto& it : results)
    {
        double score = scorer.ComputeScore(it);
        it.set_score(score);
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
