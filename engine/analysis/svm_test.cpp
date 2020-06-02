#define BOOST_TEST_MODULE SVMTest
#include <boost/test/unit_test.hpp>

# include "svm_analyzer.h"

namespace engine {
namespace analysis {

using namespace engine::search;


BOOST_AUTO_TEST_CASE( search_engine_test ) {

SVMAnalyzer analyzer;
std::vector<SearchResult> v1, v2;
for (int i = 0; i < 1000; i++)
{
    SearchResult r;
    r.Add(i, SearchType::Core);
    v1.push_back(r);
}
for (int i = 0; i < 1000; i++)
{
    SearchResult r;
    r.Add(i+2, SearchType::Branch);
    r.Add(i+3, SearchType::Terminal);
    r.Add(i+4, SearchType::Peptide);
    r.Add(i+5, SearchType::Oxonium);
    v2.push_back(r);
}

analyzer.set_problem(v1, v2);
analyzer.set_problem(v1, v2);
analyzer.set_problem(v1, v2);
analyzer.set_problem(v1, v2);

analyzer.Training(v1, v2);
SearchResult r;
int i = 0;
    r.Add(i, SearchType::Core);
// r.Add(i+2, SearchType::Branch);
// r.Add(i+3, SearchType::Terminal);
// r.Add(i+4, SearchType::Peptide);
// r.Add(i+5, SearchType::Oxonium);
std::cout << analyzer.Predicting(r) << std::endl;
for (auto p : analyzer.PredictingProbability(r) )
{
    std::cout << p << std::endl;
}

BOOST_CHECK(analyzer.Predicting(r) == 0);

}

} // namespace analysis 
} // namespace engine



