#define BOOST_TEST_MODULE SVMTest
#include <boost/test/unit_test.hpp>

# include "svm_analyzer.h"

namespace engine {
namespace analysis {

using namespace engine::search;


BOOST_AUTO_TEST_CASE( search_engine_test ) {

SVMAnalyzer analyzer;
std::vector<SearchResult> v1, v2;
for (int i = 0; i < 10; i++)
{
    SearchResult r;
    r.Add(1.0, SearchType::Core);
    v1.push_back(r);
}
for (int i = 0; i < 10; i++)
{
    SearchResult r;
    r.Add(2.0, SearchType::Branch);
    r.Add(2.0, SearchType::Terminal);
    r.Add(2.0, SearchType::Peptide);
    r.Add(2.0, SearchType::Oxonium);
    v2.push_back(r);
}

analyzer.Training(v1, v2);
SearchResult r;
r.Add(1.0, SearchType::Core);
std::cout << analyzer.Predicting(r) << std::endl;


}

} // namespace analysis 
} // namespace engine



