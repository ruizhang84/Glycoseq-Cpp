#include <iostream>
#include <vector>
#include <random>
#include <string>
#include "lsh_clustering.h"
#include "../../util/io/mgf_parser.h"
#include "../../engine/spectrum/spectrum_binpacking.h"


int main(int argc, char *argv[])
{
    std::string path = "/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf";
    std::unique_ptr<util::io::MGFParser> parser = 
        std::make_unique<util::io::MGFParser>(path, model::spectrum::SpectrumType::EThcD);
    util::io::SpectrumReader spectrum_reader(path, std::move(parser));
    spectrum_reader.Init();

    algorithm::clustering::LSHUnionFind finder(&spectrum_reader, 0.6); 

    double tol = 0.05, lower = 200, upper = 2000;
    engine::spectrum::SpectrumBinPacking bin_packer(tol, lower, upper);
    int bucket_size = bin_packer.BinSize();

    algorithm::clustering::LSHClustering cluster_runner(bucket_size, 15, 2, 5, &finder);
    // finder.Union(3, 64);
    cluster_runner.union_finder_->Union(3, 64);
    cluster_runner.union_finder_->Union(3, 64);
}

