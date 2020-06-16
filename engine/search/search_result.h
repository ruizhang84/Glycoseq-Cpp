#ifndef ENGINE_SEARCH_SEARCH_SCORE_H
#define ENGINE_SEARCH_SEARCH_SCORE_H


#include <vector>
#include <map>
#include <unordered_map>
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/glycan.h"
#include <iostream>

namespace engine{
namespace search{

enum class SearchType { Core, Branch, Terminal, Oxonium, Peptide };
enum class ScoreType { Precursor };

class SearchResult
{
public:
    int Scan() const { return scan_; }
    int ModifySite() const { return pos_; }
    std::string Sequence() const { return peptide_; }
    std::string Glycan() const { return glycan_; }
    double Score() const { return score_; }
    double ExtraScore(ScoreType type) const 
    {
        const auto& it = extra_.find(type);
        if (it == extra_.end())
            return 0;
        return it->second;
    }
    void set_scan(int scan) { scan_ = scan; }
    void set_site(int pos) { pos_ = pos; }
    void set_sequence(std::string seq) { peptide_ = seq; }
    void set_glycan(std::string glycan) { glycan_ = glycan; }
    void set_score(double score) { score_ = score; }
    void set_extra(double score, ScoreType type) { extra_[type] = score; }

    static double PeakValue(const std::vector<model::spectrum::Peak>& peaks)
    { 
        double sum = 0;
        for(const auto& it : peaks)
        {
            sum += it.Intensity();
        }
        return sum;
    }
    static double PrecursorValue(const std::string peptide, const std::string composite,
        double precursor_mass, double isotopic)
    {
        double mass = util::mass::PeptideMass::Compute(peptide)
            + util::mass::GlycanMass::Compute(model::glycan::Glycan::Interpret(composite));
        
        double ppm = kPPM;
        for (int i = 0; i <= isotopic; i ++)
        {
            double isotopic_mass = util::mass::SpectrumMass::kIon * i + mass;
            double isotopic_ppm = util::mass::SpectrumMass::ComputePPM(isotopic_mass, precursor_mass);
            ppm = (ppm <= isotopic_ppm) ? ppm : isotopic_ppm;
        }
        return 1.0 - ppm / kPPM;
    }
    static constexpr double kPPM = 50.0;

protected:
    int scan_;
    std::string peptide_;
    std::string glycan_;
    int pos_;
    double score_;
    std::map<ScoreType, double> extra_;
    
};


class ResultCollector{
public:
    ResultCollector():  oxonium_(-1){}

    std::vector<SearchResult> Result()
    {
        // remove if hits over 20
        if ((int) results_.size() > max_hits)
        {
            std::sort(results_.begin(), results_.end(), 
                [](const SearchResult& r1, const SearchResult& r2) -> bool { return r1.Score() > r1.Score(); });
            results_.erase(results_.begin() + max_hits, results_.end());
        }

        return results_;
    }
    std::vector<SearchResult> BestResult()
    {
        std::vector<SearchResult> res;
        double max_score = 0;
        for (const auto& it : results_)
        {
            double score = it.Score();
            if (score >= max_score){
                if (score > max_score){
                    res.clear();
                }
                res.push_back(it);
                max_score = score;
            }
        }
        for(auto& it : res)
        {
            it.set_extra(SearchResult::PrecursorValue(it.Sequence(), it.Glycan(), 
                precursor_mass_, isotopic_), ScoreType::Precursor);
        }

        return res;
    }

    void InitCollect()
    {
        oxonium_ = -1;
        peptide_.clear();
        glycan_core_.clear();
        glycan_branch_.clear();
        glycan_terminal_.clear();
    }
    void OxoniumCollect(const std::vector<model::spectrum::Peak>& oxonium_peaks)
    {
        if (oxonium_peaks.empty()) return;
        oxonium_ = SearchResult::PeakValue(oxonium_peaks);
    }
    void PeptideCollect(const std::vector<model::spectrum::Peak>& peptide_peaks, int pos)
    {
        if (!peptide_peaks.empty())
        {  
            peptide_[pos] = SearchResult::PeakValue(peptide_peaks);
        }
    }
    void GlycanCollect(const std::vector<model::spectrum::Peak>& glycan_peaks, 
        std::string isomer, SearchType type)
    {
        switch (type)
        {
            case SearchType::Core:
                glycan_core_[isomer] = SearchResult::PeakValue(glycan_peaks);
                break;
            case SearchType::Branch:
                glycan_branch_[isomer] =  SearchResult::PeakValue(glycan_peaks);
            case SearchType::Terminal:
                glycan_terminal_[isomer] = SearchResult::PeakValue(glycan_peaks);
            default:
                break;
        }
    }
    void PrecursorCollect(double precursor_mass, int isotopic)
    {
        precursor_mass_ = precursor_mass;
        isotopic_ = isotopic;
    }
    bool OxoniumMiss() { return oxonium_ <= 0; }
    bool PeptideMiss()
    {
        return peptide_.empty();
    }
    bool GlycanMiss(const std::string& isomer)
    {
        return (glycan_core_.find(isomer) == glycan_core_.end());
    }
    bool Empty() { return results_.empty(); }
    void Update(int scan, const std::string& peptide, const std::string& composite)
    {
        for(const auto& pos_it : peptide_)
        {
            for(const auto& isomer_it : glycan_core_)
            {
                SearchResult res;
                std::string isomer = isomer_it.first;
                double score = glycan_core_[isomer] + glycan_branch_[isomer] + glycan_terminal_[isomer] 
                    + pos_it.second + oxonium_; 
                res.set_scan(scan);
                res.set_glycan(composite);
                res.set_sequence(peptide);
                res.set_site(pos_it.first);
                res.set_score(score);
                results_.push_back(res);
            }
        }
    }

protected:
    const int max_hits = 20;
    double oxonium_;
    std::unordered_map<int, double> peptide_;
    std::unordered_map<std::string, double> glycan_core_, glycan_branch_, glycan_terminal_;
    double precursor_mass_; 
    int isotopic_;
    std::vector<SearchResult> results_;

};



} // namespace engine
} // namespace search

#endif