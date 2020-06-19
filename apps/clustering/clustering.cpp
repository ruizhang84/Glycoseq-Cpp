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
    // set up default 
    double tol = 0.05;
    double lower = 200;
    double upper = 2000;
    int hash_func_num = 15;
    int cluster_num = 100;
    int thread = 5;
    double cosine = 0.8;
    path = "/home/yu/Documents/GlycoSeq-Cpp/data/ZC_20171218_H68_R1.mgf";
    out_path = "cluster.txt";

    // pharser parameter
    while ((opt = getopt(argc, argv, ":f:o:t:l:u:K:L:s:p:h")) != EOF)
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
            case 'K': 
                hash_func_num = atoi(optarg);
                std::cout <<"the number of hashing funciton to use " << hash_func_num << std::endl; 
                break;
            case 'L': 
                cluster_num =  atoi(optarg);
                std::cout <<"Clustering iterations " << cluster_num << std::endl; 
                break;
            case 's': 
                cosine = atof(optarg);
                std::cout <<"default cosine similarity threshold " << cosine << std::endl; 
            case 'p': 
                thread =  atoi(optarg);
                std::cout <<"default thread number " << thread << std::endl; 
                break;
            case 'h': 
                std::cout <<"clustering mgf spectrum, parameter:\n" 
                    << "-f [file] -o [output] \n"
                    << "-t [tolerance] -l [lower_bound] -u [upper_bound] \n"
                    << "-K [hash_function_num] -L [iterations] -s [cosine] \n" 
                    << "-p [thread_num] " << std::endl; 
                return 1;
            case ':':  
                std::cout << "option needs a value" << std::endl;  
                return 1;
            case '?':  
                std::cout << "unknown option " << std::endl; 
                return 1;
            default: 
                std::cout <<"-h for help" << std::endl;
                return 1;
        }

    // read spectrum
    auto start = high_resolution_clock::now(); 
    std::unique_ptr<SpectrumParser> parser = 
        std::make_unique<MGFParser>(path, SpectrumType::EThcD);
    SpectrumReader spectrum_reader(path, std::move(parser));
    spectrum_reader.Init();

    std::vector<Spectrum> spectra = spectrum_reader.GetSpectrum();
    std::unordered_map<int, std::vector<double>> data_set;

    // binpacking spectra
    SpectrumBinPacking bin_packer(tol, lower, upper);
    int bucket_size = bin_packer.BinSize();

    for(auto spec : spectra)
    {
        std::vector<double> v = bin_packer.Packing(spec);
        data_set[spec.Scan()] = v;   
    }
    
    // start clustering
    LSHUnionFind finder (&spectrum_reader, cosine);
    LSHClustering lsh(bucket_size, hash_func_num, cluster_num, thread, &finder);
    lsh.set_data(data_set);
    std::unordered_map<int, std::vector<int>> clustred = lsh.Clustering();

    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
    std::cout << duration.count() << std::endl; 

    // write out result
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