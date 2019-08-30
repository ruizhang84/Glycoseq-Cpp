// #include "BucketSearch.h"
#include <vector>
#include <iostream>
#include "mass_point.h"
#include "bucket_search.h"

using namespace std;

int main(){
    vector<shared_ptr<Point>> data;
    for (int i = 0; i < 12; i++){
        // shared_ptr<MassPoint> p = make_shared<MassPoint>(i);
        shared_ptr<MassPoint> p (new MassPoint(i));
        data.push_back(p);
    }

    // for (int i = 0; i < 12; i++){
    //     cout << data[i]->get_value() << endl;
    // }

    BucketSearch bucket(data, 0.6);
    shared_ptr<MassPoint> n = make_shared<MassPoint>(3.5);
    vector<shared_ptr<Point>> res = bucket.Search(n);

    for (int i = 0; i < res.size(); i++){
        cout << res[i]->get_value() << endl;
    }
}