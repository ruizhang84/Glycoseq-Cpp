#define BOOST_TEST_MODULE LSHTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>
#include <random>

#include "lsh.h"

namespace util {
namespace calc {

BOOST_AUTO_TEST_CASE( calc_test ) 
{
    Calc calculator;
    std::vector<double> v1 {1, 2, 3, 4};
    std::vector<double> v2 {2, 3, 4, 5};
    std::vector<double> v3 {1,-2, 3, -4};
    std::vector<double> v4 {1,-2, 3};
   
    BOOST_CHECK(calculator.DotProduct(v1, v2) == 40.0); 
    BOOST_CHECK(calculator.DotProduct(v1, v3) == -10.0); 
}

BOOST_AUTO_TEST_CASE( LSH_test ) 
{
    LSH mapper(100, 12);
    
    std::random_device generator;
    std::uniform_real_distribution<double> distribution(-10.0, 0);
    std::vector<double> v;

    int nrolls = 100;
    for (int i = 0; i < nrolls; ++i) 
    {
        v.push_back(distribution(generator));
    }

    int key = mapper.Key(v);
    BOOST_CHECK( key >= 0); 
}

} // namespace io
} // namespace calc


