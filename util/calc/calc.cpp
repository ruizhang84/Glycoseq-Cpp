#include "calc.h"

namespace util {
namespace calc {

double Calc::DotProduct(
    std::vector<double>& vect1, 
    std::vector<double>& vect2) const
{
    return std::inner_product(vect1.begin(), vect1.end(),
            vect2.begin(), 0);
}

} // namespace calc
} // namespace util


