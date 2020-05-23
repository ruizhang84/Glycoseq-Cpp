#ifndef ALGORITHM_CLUSTERING_LSH_CLUSTERING_H
#define ALGORITHM_CLUSTERING_LSH_CLUSTERING_H

#include <vector>
#include <unordered_map> 
#include <mutex> 
#include <thread>  
#include "../base/union_find.h"
#include "../../util/calc/lsh.h"

namespace algorithm {
namespace clustering {

class LSHClustering
{
public:
    LSHClustering(int bin_size, int hash_func_num, int cluster_num): 
            bin_size_(bin_size), hash_func_num_(hash_func_num), 
                cluster_num_(cluster_num) { Init(); }

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

protected:
    virtual void ParClustering
        (util::calc::LSH cluster, std::vector<std::tuple<int, int>>& union_set);

    util::calc::LSH GenLSH() 
        { return util::calc::LSH(bin_size_, hash_func_num_); }
    int bin_size_;
    int hash_func_num_;
    int cluster_num_;
    std::vector<util::calc::LSH> cluster_;
    algorithm::base::UnionFind union_finder_;
    std::unordered_map<int, std::vector<double>> data_; 
    std::mutex mutex_;  
};

} // namespace clustering
} // namespace algorithm


#endif