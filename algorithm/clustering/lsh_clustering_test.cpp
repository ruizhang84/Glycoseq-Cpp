#include <iostream>
#include <vector>
#include <random>
#include "lsh_clustering.h"

using namespace algorithm::clustering;

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

int main(int argc, char *argv[])
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

