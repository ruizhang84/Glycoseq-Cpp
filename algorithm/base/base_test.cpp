#define BOOST_TEST_MODULE LSHClusteringTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>
#include "union_find.h"
#include "binpacking.h"

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

BOOST_AUTO_TEST_CASE(  binpacking_test ) 
{
    BinPacking<int> packer(10, 0, 100);
    std::vector<int> data;
    for(int i = 0; i <= 100; i++)
    {
        data.push_back(i);
    }

    std::vector<std::vector<int>> result = 
        packer.Packing(data);

    for(auto& it : result)
    {   
        std::cout << "bucket: " << std::endl;
        for(auto i : it)
        {
            std::cout << i << std::endl;
        }
    }

}    

} // namespace base
} // namespace algorithm