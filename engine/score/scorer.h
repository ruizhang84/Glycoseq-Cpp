#ifndef ENGINE_SCORE_SCORER_H
#define ENGINE_SCORE_SCORER_H

#include "../search/search_result.h"

namespace engine{
namespace score {

class SimpleScorer
{
public:
    SimpleScorer(){}
    virtual double ComputeScore(const engine::search::SearchResult& result)
    {
        return result.Score();
    }
};


}   // namespace score 
}   // namespace engine



#endif