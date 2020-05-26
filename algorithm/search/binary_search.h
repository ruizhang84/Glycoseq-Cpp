#ifndef ALGORITHM_BINARY_SEARCH_H
#define ALGORITHM_BINARY_SEARCH_H

#include <vector>
#include <cstdlib> 
#include <algorithm>
#include "search.h"
#include "../../util/mass/spectrum.h"

namespace algorithm {
namespace search {

// enhanced binary search
class BinarySearch
{
public:
    BinarySearch(double tol, ToleranceBy by):
        tolerance_(tol), by_(by) {};

    virtual void Init() 
    {
        if (!data_.empty())
            std::sort(data_.begin(), data_.end());
    }
    double Tolerance() const { return tolerance_; }
    std::vector<double>& Data() { return data_; }
    ToleranceBy ToleranceType() const { return by_; }
    void set_tolerance(double tol) { tolerance_ = tol; }
    void set_tolerance_by(ToleranceBy by) { by_ = by; }
    void set_data(std::vector<double> data) { data_ = data; Init(); }


    virtual bool Search(const double target)
    {
        if (data_.empty()) 
            return false;

        int start = 0, end = data_.size()-1;
        while (start <= end)
        {
            int mid = (end - start) / 2 + start;
            if (Match(data_[mid], target))
                return true;
            else if (data_[mid] < target)
                start = mid + 1;
            else
                end = mid - 1;
        }
        return false;
    }

protected:
    virtual bool Match(const double p, const double target)
    {
        switch (by_)
        {
        case ToleranceBy::PPM:
            return util::mass::SpectrumMass::ComputePPM(p, target) < tolerance_;
        case ToleranceBy::Dalton:
            return std::abs(p - target) < tolerance_;
        default:
            break;
        }
        return false;
    }

    double tolerance_; 
    ToleranceBy by_;
    std::vector<double> data_;
};

} // namespace algorithm
} // namespace search 

#endif
