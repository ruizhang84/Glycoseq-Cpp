#include <iostream>
#include <vector>
#include <random>
#include <string>
#include "lsh_clustering.h"
#include "../../util/io/mgf_parser.h"

using namespace algorithm::clustering;

using namespace util::io;

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
    std::string path = "/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf";
    std::unique_ptr<SpectrumParser> parser = 
        std::make_unique<MGFParser>(path, SpectrumType::EThcD);
    std::unique_ptr<SpectrumReader> spectrum_reader = 
        std::make_unique<SpectrumReader>(path, std::move(parser));
    spectrum_reader->Init();

    LSHClustering cluster_runner(100, 15, 2, 0.6, spectrum_reader.get());
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

