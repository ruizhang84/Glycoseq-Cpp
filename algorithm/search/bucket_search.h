#ifndef ALGORITHM_BUCKET_SEARCH_H
#define ALGORITHM_BUCKET_SEARCH_H

#include <vector>
#include <memory>
#include "point.h"

namespace algorithm {
namespace search {

enum class ToleranceBy { PPM, Dalton}; 

template <class T>
class BucketSearch{
public:
    BucketSearch(double tol, ToleranceBy by):
        tolerance_(tol), by_(by) { Init(); };

    void Init()
    {

    }

    // std::vector<std::shared_ptr<Point<T>>> Search(double target); 
    // {
    //     std::vector<std::shared_ptr<Point<T>>> result;
    //     return result;
    // }

protected:
    virtual bool Match(Point<T>* p, double target)
    {
        switch (by_)
        {
        case ToleranceBy::PPM:
            /* code */
            break;
        case ToleranceBy::Dalton:
            double val = p->Value() - target;
            val = val >= 0 ? val : - val;
            return val < tolerance_;

        default:
            break;
        }
        return false;
    }

    std::vector<std::shared_ptr<Point<T>>> data_;
    ToleranceBy by_;
    double tolerance_; 
    double min_;
};

} // namespace algorithm
} // namespace search 

#endif
