#ifndef ENGINE_SCORE_SCORER_H
#define ENGINE_SCORE_SCORER_H

#include "../search/search_result.h"

namespace engine{
namespace score {

class SimpleScorer
{
public:
    SimpleScorer(const std::map<engine::search::SearchType, double>& weight):
        weight_(weight){}

    const std::map<engine::search::SearchType, double> Weight() const 
        { return weight_; }
    void set_weight(const std::map<engine::search::SearchType, double>& weight)
        { weight_ = weight; }

    virtual double ComputeScore(const engine::search::SearchResult& result)
    {
        double score = 0;
        for(const auto& it : result.Match())
        {
            if (weight_.find(it.first) != weight_.end())
            {
                score += it.second * weight_[it.first];
            }
            
        }
        return score;
    }

protected:
    std::map<engine::search::SearchType, double> weight_;
};


}   // namespace score 
}   // namespace engine



#endif