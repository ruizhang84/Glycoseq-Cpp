#include "lsh.h"

namespace util {
namespace calc {

int LSH::Key(std::vector<double>& vect)
{
    int key = 0;
    for (int i = 0; i < hash_func_num_; i++)
    {
        key <<= 1; 
        key += RandomProjection(vect, weights_[i]);
    }

    return key;
}

std::vector<double> LSH::GenNormalDistribWeight()
{
    std::vector<double> weight;
    std::random_device generator;
    std::normal_distribution<double> distribution(kMean, kSTD);

    for (int i = 0; i < weight_size_; i++) 
    {
        double number = distribution(generator);
        weight.push_back(number);
    }
    return weight;
}

void LSH::GenHashFunc()
{
    for (int i = 0; i < hash_func_num_; i++)
    {
        weights_.push_back(GenNormalDistribWeight());
    }

}


} // namespace io
} // namespace util
