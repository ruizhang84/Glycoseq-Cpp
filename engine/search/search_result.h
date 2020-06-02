#ifndef ENGINE_SEARCH_SEARCH_SCORE_H
#define ENGINE_SEARCH_SEARCH_SCORE_H


#include <vector>
#include <unordered_map>
#include "../../model/spectrum/spectrum.h"

#include <iostream>

namespace engine{
namespace search{

enum class SearchType { Core, Branch, Terminal, Peptide, Oxonium};

class SearchResult
{
public:
    int Scan() { return scan_; }
    int ModifySite() { return pos_; }
    std::string Sequence() { return peptide_; }
    std::string Glycan() { return glycan_; }
    double Score() const { return score_; }
    const std::unordered_map<SearchType, double> Match() const 
        { return match_; }

    void set_scan(int scan) { scan_ = scan; }
    void set_site(int pos) { pos_ = pos; }
    void set_sequence(std::string seq) { peptide_ = seq; }
    void set_glycan(std::string glycan) { glycan_ = glycan; }
    void set_score(double score) { score_ = score; }
    void set_match(std::unordered_map<SearchType, double> match)
        { match_ = match; }
    
    void Add(double value, SearchType type)
    {
        match_[type] = value;
    }

    static double PeakValue(const std::vector<model::spectrum::Peak>& peaks)
    { 
        double sum = 0;
        for(const auto& it : peaks)
        {
            sum += it.Intensity();
        }
        return sum;
    }
   
protected:
    int scan_;
    std::string peptide_;
    std::string glycan_;
    int pos_;
    double score_;
    std::unordered_map<SearchType, double> match_;
    
};

class SimpleScorer
{
public:
    SimpleScorer(const std::unordered_map<SearchType, double>& parameter):
        param_(parameter){}

    const std::unordered_map<SearchType, double> Parameter() const 
        { return param_; }
    void set_parameter(const std::unordered_map<SearchType, double>& param)
        { param_ = param; }

    virtual double ComputeScore(const SearchResult& result)
    {
        double score = 0;
        for(const auto& it : result.Match())
        {
            if (param_.find(it.first) != param_.end())
            {
                score += it.second * param_[it.first];
            }
            
        }
        return score;
    }

protected:
    std::unordered_map<SearchType, double> param_;
};

} // namespace engine
} // namespace search

#endif