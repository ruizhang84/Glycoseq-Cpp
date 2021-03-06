#include "lsh_clustering.h"

namespace algorithm {
namespace clustering {

std::unordered_map<int, std::vector<int>> LSHClustering::Clustering()
{
    std::unordered_map<int, std::vector<int>> result;
    std::vector<std::tuple<int, int>> union_set;
    
    // cluster
    if ((int) cluster_.size() <= thread_)
    {
        std::vector<std::thread> thread_pool;
        for(auto& cluster : cluster_)
        {   
            std::thread t(&LSHClustering::ParClustering, this, std::ref(cluster), std::ref(union_set));
            thread_pool.push_back(std::move(t));
        }
        for(auto& t: thread_pool)
        {
            t.join();
        }
    }
    else
    {
        for(int i = 0; i < (int) cluster_.size(); i += thread_)
        {
            std::vector<std::thread> thread_pool;
            int end = (int) cluster_.size() < (i + thread_)
                ? (int) cluster_.size() : (i + thread_);

            std::vector<util::calc::LSH> par_cluster(cluster_.begin() + i, cluster_.begin() + end);
            for(auto& cluster : par_cluster)
            {   
                std::thread t(&LSHClustering::ParClustering, this, std::ref(cluster), std::ref(union_set));
                thread_pool.push_back(std::move(t));
            }
            for(auto& t: thread_pool)
            {
                t.join();
            }
        }
    }



    for(auto x : union_set)
    {
        union_finder_->Union(std::get<0>(x), std::get<1>(x));
    }

    // generate clustering result
    for(auto i : data_)
    {
        int index = union_finder_->Find(i.first);
        if (result.find(index) == result.end())
        {
            result.emplace(index, std::vector<int>());
        }
        result[index].push_back(i.first);
    }
    union_finder_->Clear();

    return result;
}

void LSHClustering::ParClustering
    (util::calc::LSH cluster, std::vector<std::tuple<int, int>>& union_set)
{
    std::unordered_map<int, std::vector<int>> hash_table;
    std::vector<std::tuple<int, int>> temp_set;
    // create hash table
    for (auto it : data_)
    {
        int key = cluster.Key(it.second);
        if (hash_table.find(key) == hash_table.end())
        {
            hash_table.emplace(key, std::vector<int>());
        }
        hash_table[key].push_back(it.first);
    }
    
    // union find
    for(auto const& bucket : hash_table)
    {
        // union
        if (bucket.second.size() > 1)
        {
            for (size_t i = 0; i < bucket.second.size()-1; i++)
            {
                for (size_t j = i + 1; j < bucket.second.size(); j++){
                    int x = bucket.second[i], y = bucket.second[j];
                    temp_set.push_back(std::make_tuple(x, y));
                }
            }
        }
    }
    mutex_.lock();
    union_set.insert(union_set.end(), temp_set.begin(), temp_set.end());
    mutex_.unlock();
}

} // namespace clustering
} // namespace algorithm