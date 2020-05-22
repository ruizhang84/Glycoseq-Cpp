#ifndef UTIL_CALC_LSH_H
#define UTIL_CALC_LSH_H

#include <random>
#include <vector>
#include "calc.h"


namespace util {
namespace calc {

class LSH
{
public:
    LSH(int size, int num): 
        weight_size_(size), hash_func_num_(num){ Init(); };
        
    void Init() 
    { 
        if (!weights_.empty())
            weights_.clear(); 
        GenHashFunc(); 
    }
    int Key(std::vector<double>& vect);
    
    int WeightSize() { return weight_size_; }
    int HashFuncNum() { return hash_func_num_; }
    void set_weight_size(int size) { weight_size_ = size; }
    void set_hash_func_num(int num) { hash_func_num_ = num; } 


protected:
    int RandomProjection(std::vector<double>& vect, 
        std::vector<double> weight)
    {
        return calculator_.DotProduct(vect, weight) 
                >= 0 ? 1 : 0;
    }
    std::vector<double> GenNormalDistribWeight();
    void GenHashFunc();

    Calc calculator_;
    std::vector<std::vector<double>> weights_;
    int weight_size_; 
    int hash_func_num_;
    const int kMean = 0;
    const int kSTD = 1; 
};


} // namespace io
} // namespace util


#endif