#include <vector>
#include <string>
#include <iostream>
#include <unistd.h>
#include <cstdio>

#include "../../algorithm/clustering/lsh_clustering.h"
#include "../../util/io/spectrum_reader.h"
#include "../../util/io/mgf_parser.h"


int main(int argc, char *argv[]){
    extern char *optarg;
    int opt;
    std::string path;
    char c;
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

    return 0;
}