#define BOOST_TEST_MODULE LSHClusteringTest
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>
#include <random>
#include "lsh_clustering.h"
// #include "../../util/io/mgf_parser.h"

namespace algorithm {
namespace clustering {

std::vector<double> GenData()
{
    std::random_device generator;
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);
    std::vector<double> v;

    int nrolls = 100;
    for (int i = 0; i < nrolls; ++i) 
    {
        v.push_back(distribution(generator));
    }
    return v;
}

BOOST_AUTO_TEST_CASE( lsh_clustering_test ) 
{
    LSHClustering cluster_runner(100, 15, 1);
    std::unordered_map<int, std::vector<double>> data;
    for(int i = 0; i < 100; i++)
    {
        data.emplace(i, GenData());
    }
    cluster_runner.set_data(data);
    std::unordered_map<int, std::vector<int>> result = 
        cluster_runner.Clustering();
    for(auto i : result)
    {
        std::cout << "The number: " << i.first 
            << " contains " << std::endl;
        for(auto j : i.second)
        {
            std::cout << j << std::endl;
        }
    }

}

} // namespace clustering
} // namespace algorithm