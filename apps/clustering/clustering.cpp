#include <vector>
#include <string>
#include <iostream>
#include <unistd.h>
#include <cstdio>
#include <unordered_map> 
#include <fstream>
#include <chrono> 

#include "../../algorithm/clustering/lsh_clustering.h"
#include "../../util/io/mgf_parser.h"
#include "../../engine/spectrum/spectrum_binpacking.h"
#include "../../util/calc/spectrum_sim.h"

using namespace util::io;
using namespace engine::spectrum;
using namespace algorithm::clustering;
using namespace std::chrono; 


int main(int argc, char *argv[]){
    extern char *optarg;
    int opt;
    std::string path, out_path;
    double tol = 0.05;
    double lower = 200;
    double upper = 2000;
    int hash_func_num = 16;
    int cluster_num = 10;
    int thread = 20;
    double cosine = 0.7;

    while ((opt = getopt(argc, argv, ":t:l:u:k:i:f:o:h:d:c")) != EOF)
        switch(opt)
        {
            case 'f': 
                path = optarg;
                std::cout <<"Clustering the file located at " << path << std::endl; 
                break;
            case 'o': 
                out_path = optarg;
                std::cout <<"Output the file at " << out_path << std::endl; 
                break;
            case 't': 
                tol = atof(optarg);
                std::cout <<"The tolerance " << tol << std::endl; 
                break;
            case 'l': 
                lower = atof(optarg);
                std::cout <<"the lower bound of spectrum " << lower << std::endl; 
                break;
            case 'u': 
                upper = atof(optarg);
                std::cout <<"the upper bound of spectrum " << upper << std::endl; 
                break;
            case 'k': 
                hash_func_num = atoi(optarg);
                std::cout <<"the number of hashing funciton to use " << hash_func_num << std::endl; 
                break;
            case 'i': 
                cluster_num =  atoi(optarg);
                std::cout <<"Clustering iterations " << cluster_num << std::endl; 
                break;
            case 'c': 
                cosine = atof(optarg);
                std::cout <<"default cosine similarity threshold " << cosine << std::endl; 
            case 'd': 
                thread =  atoi(optarg);
                std::cout <<"default thread number " << thread << std::endl; 
                break;
            case 'h': 
                std::cout <<"clustering mgf spectrum, " 
                    << "-f [file] -t [tolerance] -l [lower_bound] -u [upper_bound] "
                    << "-k [function_num] -i [iterations] " << std::endl; 
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

    auto start = high_resolution_clock::now(); 
    path = "/home/yu/Documents/MultiGlycan-Cpp/data/test_EThcD.mgf";
    out_path = "cluster.txt";
    std::unique_ptr<SpectrumParser> parser = 
        std::make_unique<MGFParser>(path, SpectrumType::EThcD);
    SpectrumReader spectrum_reader(path, std::move(parser));
    spectrum_reader.Init();

    std::vector<Spectrum> spectra = spectrum_reader.GetSpectrum();
    std::unordered_map<int, std::vector<double>> data_set;

    SpectrumBinPacking bin_packer(tol, lower, upper);
    int bucket_size = bin_packer.BinSize();

    for(auto spec : spectra)
    {
        std::vector<double> v = bin_packer.Packing(spec);
        data_set[spec.Scan()] = v;   
    }
    
    LSHUnionFind finder (&spectrum_reader, cosine);
    LSHClustering lsh(bucket_size, hash_func_num, cluster_num, &finder);
    lsh.set_data(data_set);
    std::unordered_map<int, std::vector<int>> clustred = lsh.Clustering();

    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

    std::ofstream outfile;
    outfile.open (out_path);
    outfile << "scan num \n";
    for(auto& it : clustred)
    {
        for(auto& scan_num : it.second)
        {
            outfile << scan_num << " ";
        }
        outfile <<"\n";
    }
   
    outfile.close();

    return 0;
}