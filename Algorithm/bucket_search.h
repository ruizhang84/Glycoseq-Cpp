#ifndef BUCKET_SEARCH_H
#define BUCKET_SEARCH_H
#include <vector>
#include <memory>
#include "point.h"

class BucketSearch{
private:
    std::vector<std::vector<std::shared_ptr<Point>>> points;
    double min_value;
    double tolerance;
    bool Match(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2);
public:
    BucketSearch() = default;
    BucketSearch(std::vector<std::shared_ptr<Point>>& data, double tol);  
    std::vector<std::shared_ptr<Point>> Search(std::shared_ptr<Point> target); 
};

#endif
