#define BOOST_TEST_MODULE LSHClusteringTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>
#include "union_find.h"

namespace algorithm {
namespace base {

BOOST_AUTO_TEST_CASE(  union_find_test ) 
{
    UnionFind finder;

    finder.Union(1, 7);
    finder.Union(2, 7);
    finder.Union(3, 2);

    BOOST_CHECK(finder.IsSameSet(1, 3));
    BOOST_CHECK(finder.Find(3) == 7);

}    

} // namespace base
} // namespace algorithm