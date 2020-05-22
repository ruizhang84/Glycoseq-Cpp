#include "lsh_clustering.h"

namespace algorithm {
namespace clustering {

std::unordered_map<int, std::vector<int>> LSHClustering::Clustering()
{
    std::unordered_map<int, std::vector<int>> result;
    std::unordered_map<int, std::vector<int>> hash_table;
    for(auto cluster : cluster_)
    {   
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
        for(auto const& item : hash_table)
        {
            // union
            if (item.second.size() > 1)
            {
                auto x = item.second.front();
                for(auto y : item.second)
                {
                    union_finder_.Union(x, y);
                }
            }
        }
        // clean up
        hash_table.clear();
    }

    // generate clustering result
    for(auto i : data_)
    {
        int index = union_finder_.Find(i.first);
        if (result.find(index) == result.end())
        {
            result.emplace(index, std::vector<int>());
        }
        result[index].push_back(i.first);
    }
    union_finder_.Clear();

    return result;
}


} // namespace clustering
} // namespace algorithm