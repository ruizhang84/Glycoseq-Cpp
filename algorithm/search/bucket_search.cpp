#include <algorithm>
#include "bucket_search.h"
#include <iostream>
using namespace std;

BucketSearch::BucketSearch(vector<shared_ptr<Point>>& data, double tol) :tolerance(tol == 0 ? 1 : tol){
    if (data.size() == 0)
        return ;
    
    // get lower and upper bound of data
    double lower_bound = INT_MAX;
    double upper_bound = INT_MIN;
    for(int i = 0; i < data.size(); i++){
        lower_bound = min(lower_bound, data[i]->get_value());
        upper_bound = max(upper_bound, data[i]->get_value());
    }
    // intialize 
    this->min_value = lower_bound;
    int bucket_size = (upper_bound - lower_bound) / tolerance + 1;
    for(int i = 0; i < bucket_size; i++){
        points.push_back(vector<shared_ptr<Point>>());
    }

    // fill the bucket
    for(int i = 0; i < data.size(); i++){
       int position = (data[i]->get_value() - lower_bound) / tolerance;
       points[position].push_back(data[i]);
    }
}

bool BucketSearch::Match(std::shared_ptr<Point>& p1, std::shared_ptr<Point>& p2)
{
    double val = p1->get_value() - p2->get_value();
    val = val >= 0 ? val : - val;
    return val < tolerance;
}

vector<shared_ptr<Point>> BucketSearch::Search(shared_ptr<Point> pt)
{
    vector<shared_ptr<Point>> result;
    int position = (pt->get_value() - this->min_value) / this->tolerance;
    for (int i = position - 1; i <= position + 1; i++){
        if (i >= 0 && i < points.size()){
            vector<shared_ptr<Point>> data = points[i];
            for(int j = 0; j < data.size(); j++){
                if (Match(data[j], pt)){
                    result.push_back(data[j]);
                }
            }
        }
    }
    return result;
}




