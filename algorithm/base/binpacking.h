#ifndef ENGINE_SPECTRUM_BINPACKING_H
#define ENGINE_SPECTRUM_BINPACKING_H
#include <vector>
#include <cmath> 

namespace algorithm {
namespace base {

template <class T>
class BinPacking
{
public:
    BinPacking(double tol, double lower, double upper):
        tolerance_(tol), lower_(lower), upper_(upper) {};

    std::vector<std::vector<T>> Packing(
        const std::vector<T>& vect) const
    {
        std::vector<std::vector<T>> bins;
        bins.assign(Bucket(), std::vector<T>());
        for(auto& it : vect)
        {
            int index = Index(Position(it));
            bins[index].push_back(it);
        }
        return bins;
    }

    double Tolerance() { return tolerance_; }
    void set_lower(double lower) { lower_ = lower; }
    void set_upper(double upper) { upper_ = upper; }
    void set_tolerance(double tol) { tolerance_ = tol; }

protected:
    virtual int Index(int pos) const
        { return (int) floor((pos - lower_) / tolerance_); }
    virtual int Position(const T& elem) const
        { return (int) elem; }
    virtual int Bucket() const
        { return (int) ceil((upper_ - lower_ + 1) / tolerance_); }

    double tolerance_;
    double lower_;
    double upper_;

};

} // namespace base
} // namespace algorithm

#endif