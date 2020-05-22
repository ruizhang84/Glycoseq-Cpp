#ifndef UTIL_CALC_CALC_H
#define UTIL_CALC_CALC_H

#include <vector>
#include <numeric>

namespace util {
namespace calc {

class Calc
{
public:
    Calc() = default;
    double DotProduct(std::vector<double>&, 
        std::vector<double>&);
};


} // namespace io
} // namespace calc


#endif