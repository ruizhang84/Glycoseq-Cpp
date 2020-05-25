#define BOOST_TEST_MODULE ProteinTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <iostream>
#include "search.h"

using namespace std;
namespace algorithm {
namespace search {

std::shared_ptr<Point<double>> CreatePoint(double value)
{
    return make_shared<Point<double>>(value, std::vector<double>{ value });
}

BOOST_AUTO_TEST_CASE( Monosaccharide_test ) 
{
    SearchBase<double> searcher(20.0, ToleranceBy::Dalton);
    std::vector<std::shared_ptr<Point<double>>> box; 

    for(int i=0; i<100; i++)
    {
        box.push_back(CreatePoint(i));
    }
    searcher.set_data(box);
    std::vector<double> res = searcher.Search(40);
    BOOST_CHECK(res.size() == 39);

}

} // namespace algorithm
} // namespace search 