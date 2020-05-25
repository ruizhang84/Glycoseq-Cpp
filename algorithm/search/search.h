#ifndef ALGORITHM_SEARCH_H
#define ALGORITHM_SEARCH_H

#include <vector>
#include <memory>
#include <cstdlib> 
#include "point.h"
#include "../../util/mass/spectrum.h"
#include <iostream>
namespace algorithm {
namespace search {

enum class ToleranceBy { PPM, Dalton}; 

template <class T>
class SearchBase
{
typedef std::vector<std::shared_ptr<Point<T>>> Points;
public:
    SearchBase(double tol, ToleranceBy by):
        tolerance_(tol), by_(by) {};

    virtual void Init() {}
    double Tolerance() const { return tolerance_; }
    Points& Data() { return data_; }
    ToleranceBy ToleranceType() const { return by_; }
    void set_tolerance(double tol) { tolerance_ = tol; }
    void set_tolerance_by(ToleranceBy by) { by_ = by; }
    void set_data(std::vector<std::shared_ptr<Point<T>>> data)
        { data_ = data; }

    virtual std::vector<T> Query(const double target)
    {
        std::vector<T> result;
        for(const auto& it: data_)
        {
            if (Match(it.get(), target))
            {
                result.push_back(it->Content());
            }

        }
        return result;
    }

    virtual bool Search(const double target)
    {
        std::vector<T> result;
        for(const auto& it: data_)
        {
            if (Match(it.get(), target))
            {
                return true;
            }

        }
        return false;
    }

protected:
    virtual bool Match(const Point<T>* p, const double target)
    {
        double diff; 
        switch (by_)
        {
        case ToleranceBy::PPM:
            diff = util::mass::SpectrumMass::ComputePPM(p->Value(), target);
            return diff < tolerance_;
        case ToleranceBy::Dalton:
            diff = p->Value() - target;
            return std::abs(diff) < tolerance_;
        default:
            break;
        }
        return false;
    }

    double tolerance_; 
    ToleranceBy by_;
    Points data_;
};

} // namespace algorithm
} // namespace search 

#endif
