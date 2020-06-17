#ifndef ENGINE_ANALYSIS_MULTI_COMPARISON_H
#define ENGINE_ANALYSIS_MULTI_COMPARISON_H

#include <algorithm>
#include  "../search/search_result.h"
#include <boost/math/distributions/normal.hpp>
using boost::math::normal_distribution;

namespace engine{
namespace analysis{

class MultiComparison
{
public:
    MultiComparison(double fdr): fdr_(fdr){}

    std::vector<engine::search::SearchResult> Tests(
        const std::vector<engine::search::SearchResult>& targets,
        const std::vector<engine::search::SearchResult>& decoys)
    {
        std::vector<engine::search::SearchResult> results;

        std::map<int, std::vector<double>> v;
        for(const auto& it : decoys)
        {
            int i = it.Scan();
            if (v.find(i) == v.end())
            {
                v[i] = std::vector<double>();
            }
            v[i].push_back(it.Score());
        }

        std::map<int, engine::search::SearchResult> s_results;
        std::map<int, double> s_v;
        for(const auto& it : targets)
        {
            if (v.find(it.Scan()) == v.end())
            {
                results.push_back(it);
            }
            else
            {
                std::vector<double> score_list = v[it.Scan()];
                double avg = mean(score_list);
                if (avg > it.Score()) continue;

                double p = pValue(score_list, it.Score());

                s_v[it.Scan()] = p;
                s_results[it.Scan()] = it;
            }
        }
        std::vector<double> p_values;
        for(const auto& it : s_v)
        {
            p_values.push_back(it.second);
        }

        std::sort(p_values.begin(), p_values.end());
        std::map<int, double> r_s_v;
        int size = (int) p_values.size();
        double c = 0;
        for(int i = 1; i <= size; i ++)
        {
            c += 1.0 /i;
        }
        for(const auto& it : s_v)
        {
            double p = it.second;
            
            auto i = std::find(p_values.begin(), p_values.end(), p);
            int d = std::distance(p_values.begin(), i);
            double r_p = d * p * 1.0 / size / c;
            if (r_p < fdr_)
            {
                results.push_back(s_results[it.first]);
            }
        }
    }

protected:
    double fdr_;

    static double mean(std::vector<double> v)
    {
        int i = 0;
        double m = 0;
        for(auto it : v)
        {
            m += it; ++i;
        }
        return m * 1.0 / i;
    }

    static double stdv(std::vector<double> v)
    {
        int i = 0;
        double m = 0;
        double avg = mean(v);
        for(auto it : v)
        {
            m += (it - avg) * (it - avg); ++i;
        }
        return std::sqrt( m * 1.0/ i);
    }

    static double pValue(std::vector<double> v, double q)
    {
        double m = mean(v);
        double s = stdv(v);
        return 0.5 * erfc( (q-m) * 1.0 / (s * std::sqrt(2)) );
    }

};

}
}


#endif