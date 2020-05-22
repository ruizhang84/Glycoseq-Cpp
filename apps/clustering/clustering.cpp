#include <vector>
#include <string>
#include <iostream>
#include <unistd.h>
#include <cstdio>
#include <unordered_map> 

#include "../../algorithm/clustering/lsh_clustering.h"
#include "../../util/io/spectrum_reader.h"
#include "../../util/io/mgf_parser.h"
#include "../../engine/spectrum/spectrum_binpacking.h"

using namespace util::io;
using namespace engine::spectrum;
using namespace algorithm::clustering;

int main(int argc, char *argv[]){
    extern char *optarg;
    int opt;
    std::string path;
    while ((opt = getopt(argc, argv, ":f:h")) != EOF)
        switch(opt)
        {
            case 'f': 
                path = optarg;
                std::cout <<"Clustering the file located at " << path << std::endl; 
                break;
            case 'h': 
                std::cout <<"clustering mgf spectrum, -f [file]" << std::endl; 
                break;
            case ':':  
                std::cout << "option needs a value" << std::endl; 
                break;  
            case '?':  
                std::cout << "unknown option " << std::endl; 
                break;  
            default: 
                std::cout <<"-h for help" << std::endl;
        }

    path = "/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf";
    std::unique_ptr<SpectrumParser> parser = 
        std::make_unique<MGFParser>(path, SpectrumType::EThcD);
    SpectrumReader spectrum_reader(path, std::move(parser));
    spectrum_reader.Init();

    std::cout << spectrum_reader.GetFirstScan() << std::endl;
    std::cout << spectrum_reader.GetLastScan() << std::endl;

    std::vector<Spectrum> spectra = spectrum_reader.GetSpectrum();

    double tol = 0.1;
    double lower = 200;
    double top = 2000;
    int hash_func_num = 15;
    int cluster_num = 100;
    SpectrumBinPacking bin_packer(tol, lower, top);
    LSHClustering lsh(bin_packer.BinSize(), hash_func_num, cluster_num);
    std::unordered_map<int, std::vector<double>> data_set;
    for(auto spec : spectra)
    {
        std::vector<double> v = bin_packer.Packing(spec);
        data_set[spec.Scan()] = v;
        
    }
    lsh.set_data(data_set);
    std::unordered_map<int, std::vector<int>> clustred = lsh.Clustering();

    return 0;
}