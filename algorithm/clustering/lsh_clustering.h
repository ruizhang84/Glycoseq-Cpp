#ifndef ALGORITHM_CLUSTERING_LSH_CLUSTERING_H
#define ALGORITHM_CLUSTERING_LSH_CLUSTERING_H

#include <vector>
#include <unordered_map> 
#include <mutex> 
#include <thread>  
#include <iostream>

#include "../base/union_find.h"
#include "../../util/calc/lsh.h"
#include "../../util/io/spectrum_reader.h"
#include "../../util/calc/spectrum_sim.h"
#include "../../util/io/mgf_parser.h"

namespace algorithm {
namespace clustering {

class LSHUnionFind : public algorithm::base::UnionFind
{
public:
    LSHUnionFind(util::io::SpectrumReader* parser, double threshold): 
        threshold_(threshold), parser_(parser) {}

    void set_threshold(double threshold)
        { threshold_ = threshold; }
    void set_spectrum_reader(util::io::SpectrumReader* parser)
        { parser_ = parser; }

    virtual void Union(int i, int j) 
    { 
        if (!IsSameSet(i, j)) 
        {
            // check similarity
            if (distance_.find(EncodeKey(i, j)) != distance_.end())
            {
                double cos = distance_[EncodeKey(i, j)];
                if (cos < threshold_)
                    return;
            }
            else
            {
                std::vector<model::spectrum::Peak> pk1 = parser_->GetSpectrum(i).Peaks();
                std::vector<model::spectrum::Peak> pk2 = parser_->GetSpectrum(j).Peaks();
                double cos = calculator_.ComputeCosine(pk1, pk2);

                distance_[EncodeKey(i, j)] = cos;
                if (cos < threshold_)
                    return ;
            }

            int x = Find(i), y = Find(j);
            // rank is used to keep the tree short
            if (rank_[x] > rank_[y]) 
            { 
                map_[y] = x; 
                size_[x] += size_[y]; 
            }
            else                   
            { 
                map_[x] = y; 
                size_[y] += size_[x];
                if (rank_[x] == rank_[y]) 
                rank_[y]++; 
            } 
        } 
    }

protected:
    std::string EncodeKey(int x, int y){ 
        return x <= y ? std::to_string(x) + "+" + std::to_string(y) :
            std::to_string(y) + "+" + std::to_string(x);
    }
    double threshold_;
    util::io::SpectrumReader* parser_;
    std::unordered_map<std::string, double> distance_;
    util::calc::SpectrumSim calculator_;
};

class LSHClustering
{
public:
    LSHClustering(int bin_size, int hash_func_num, int cluster_num,
        LSHUnionFind* union_finder): bin_size_(bin_size), hash_func_num_(hash_func_num), 
            cluster_num_(cluster_num), union_finder_(union_finder) { Init(); }

    void Init()
    {
        for(int i = 0; i < cluster_num_; i++)
        {
            cluster_.push_back(GenLSH());
        }
    }

    void set_bin_size(int size) { bin_size_ = size; }
    void set_hash_func_num(int num) { hash_func_num_ = num; }
    void set_cluster_num(int num) { cluster_num_ = num; }
    void set_data(std::unordered_map<int, std::vector<double>> data)
        { data_ = data; }

    virtual std::unordered_map<int, std::vector<int>> Clustering();

// protected:
    virtual void ParClustering
        (util::calc::LSH cluster, std::vector<std::tuple<int, int>>& union_set);

    util::calc::LSH GenLSH() 
        { return util::calc::LSH(bin_size_, hash_func_num_); }
    int bin_size_;
    int hash_func_num_;
    int cluster_num_;
    std::vector<util::calc::LSH> cluster_;
    LSHUnionFind* union_finder_;
    std::unordered_map<int, std::vector<double>> data_; 
    std::mutex mutex_;  
};

} // namespace clustering
} // namespace algorithm


#endif