#ifndef ENGINE_SEARCH_FDR_FILTER_H
#define ENGINE_SEARCH_FDR_FILTER_H

#include <vector>
#include <algorithm>
#include "spectrum_search.h"

namespace engine {
namespace search {

class FDRFilter
{
public:
    FDRFilter(double fdr): fdr_(fdr), cutoff_(-1){}

    virtual void Init()
    {
         CutoffInit();
    }
    
    std::vector<SearchResult>& Data() { return data_; }
    std::vector<SearchResult>& Decoy() { return decoy_; }
    double Cutoff() const { return cutoff_; }
    void set_data(std::vector<SearchResult>& data) { data_ = data; }
    void set_decoy(std::vector<SearchResult>& decoy) { decoy_ = decoy;}
    void set_cutoff(double cutoff) { cutoff_ = cutoff; }

    virtual std::vector<SearchResult> Filter()
    {
        std::vector<SearchResult> res;
        if (cutoff_ < 0) return res;
        for(const auto& it : data_)
        {
            if (ComputeScore(it) >= cutoff_)
            {
                res.push_back(it);
            }
        }
        return res;
    }

protected:
    static double ComputeScore(const SearchResult& result)
        { return result.score; }

    static bool ScoreGreater(const SearchResult& r1, const SearchResult& r2) 
        { return (ComputeScore(r1) > ComputeScore(r2)); }

    virtual void CutoffInit()
    {
        // init
        cutoff_ = -1;
        if (decoy_.size() == 0 || data_.size() == 0 || 
            (decoy_.size() * 1.0 / (decoy_.size() + data_.size()) < fdr_))   //trivial case
        {
            return;
        }

        std::sort(data_.begin(), data_.end(), ScoreGreater);
        std::sort(decoy_.begin(), decoy_.end(), ScoreGreater);

        // compare and compute
        int i = 0, j = 0;
        while (i < (int) data_.size() && j < (int) decoy_.size())
        {

            if (ScoreGreater(data_[i], decoy_[j]))
            {
                double rate = j * 1.0 / (i + j + 1);
                if (rate <= fdr_)
                {
                    if (cutoff_ < 0)
                        cutoff_ = ComputeScore(decoy_[j]);
                    else
                        cutoff_ = std::min(cutoff_, ComputeScore(decoy_[j]));
                }
                i++;
            }
            else
            {
                j++;
            }
        }
    }

    double fdr_;
    double cutoff_;
    std::vector<SearchResult> data_;
    std::vector<SearchResult> decoy_;

};

} // namespace search
} // namespace engine 

#endif