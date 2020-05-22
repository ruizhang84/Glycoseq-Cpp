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
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-10.0, 10.0);
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
    LSHClustering cluster_runner(100, 15, 100);
    std::unordered_map<int, std::vector<double>> data;
    for(int i = 0; i < 10; i++)
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

// BOOST_AUTO_TEST_CASE( lsh_clustering_test ) 
// {
//     util::io::MGFParser 
//         parser("/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf", 
//             model::spectrum::SpectrumType::EThcD);
//     parser.Init();

//     std::vector<util::io::Peak> pk1 = parser.Peaks(64);
//     std::vector<util::io::Peak> pk2 = parser.Peaks(64);
// }

} // namespace clustering
} // namespace algorithm