#ifndef ENGINE_SEARCH_SEARCH_SCORE_H
#define ENGINE_SEARCH_SEARCH_SCORE_H


#include <vector>
#include <map>
#include "../../model/spectrum/spectrum.h"

#include <iostream>

namespace engine{
namespace search{

enum class SearchType { Core, Branch, Terminal, Peptide, Oxonium, Matches, Precursor, Coelution };

class SearchResult
{
public:
    int Scan() const { return scan_; }
    int ModifySite() const { return pos_; }
    std::string Sequence() const { return peptide_; }
    std::string Glycan() const { return glycan_; }
    double Score() const { return score_; }
    const std::map<SearchType, double> Match() const 
        { return match_; }
    double Value(SearchType type) 
    { 
        if (match_.find(type) != match_.end())
            return match_[type];
        return  0;
    }

    void set_scan(int scan) { scan_ = scan; }
    void set_site(int pos) { pos_ = pos; }
    void set_sequence(std::string seq) { peptide_ = seq; }
    void set_glycan(std::string glycan) { glycan_ = glycan; }
    void set_score(double score) { score_ = score; }
    void set_match(std::map<SearchType, double> match)
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
    std::map<SearchType, double> match_;
    
};


} // namespace engine
} // namespace search

#endif