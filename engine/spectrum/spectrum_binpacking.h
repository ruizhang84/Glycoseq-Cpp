#ifndef ENGINE_SPECTRUM_BINPACKING_H
#define ENGINE_SPECTRUM_BINPACKING_H
#include <vector>
#include <cmath> 
#include "../../model/spectrum/spectrum.h"

namespace engine {
namespace spectrum {

class SpectrumBinPacking
{
public:
    SpectrumBinPacking(double tol, double lower, double upper):
        tolerance_(tol), lower_(lower), upper_(upper) {};

    virtual std::vector<double> Packing(
        model::spectrum::Spectrum spec){};


    double Tolerance() { return tolerance_; }
    void set_lower(double lower) { lower_ = lower; }
    void set_upper(double upper) { upper_ = upper; }
    void set_tolerance(double tol) { tolerance_ = tol; }

protected:
    int Index(int pos)
        { return (int) ceil((pos - lower_) / tolerance_); }
    double tolerance_;
    double lower_;
    double upper_;

};

} // namespace spectrum
} // namespace engine

#endif