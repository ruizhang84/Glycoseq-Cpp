#ifndef UTIL_MASS_ION_H
#define UTIL_MASS_ION_H

#include <string>
#include <cmath>
#include "peptide.h"

namespace util {
namespace mass {

class SpectrumMass
{
public:
    static double Compute(const double mz, const int charge)
    {
        return (mz - kIon) * charge;
    }

    static double ComputeMZ(const double mass, const int charge)
    {
        return (mass + kIon * charge) / charge;
    }
    static double ComputePPM(const double expected, const double observed)
    {
        return std::abs(expected - observed) / expected * 1000000.0;
    }
protected:
    static constexpr double kIon = 1.0078;
};

} // namespace mass
} // namespace util

#endif