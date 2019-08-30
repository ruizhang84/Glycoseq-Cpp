#include "complex_nglycan.h"
#include <iostream>

using namespace std;

int main(){
    ComplexNGlycan nglycan;
    nglycan.set_table(0, 2);
    nglycan.set_table(1, 3);
    nglycan.set_table(2, 2);
    nglycan.set_table(3, 1);
    nglycan.set_table(4, 1);
    nglycan.set_table(5, 1);
    nglycan.set_table(6, 0);
    nglycan.set_table(7, 0);
    nglycan.set_table(8, 1);
    nglycan.set_table(9, 1);
    nglycan.set_table(10, 0);
    nglycan.set_table(11, 0);
    nglycan.set_table(12, 0);
    nglycan.set_table(13, 0);
    nglycan.set_table(14, 0);
    nglycan.set_table(15, 0);
    nglycan.set_table(16, 0);
    nglycan.set_table(17, 0);
    nglycan.set_table(18, 0);
    nglycan.set_table(19, 0);
    nglycan.set_table(20, 0);
    nglycan.set_table(21, 0);
    nglycan.set_table(22, 0);
    nglycan.set_table(23, 0);    

    vector<shared_ptr<Glycan>> glycans = nglycan.Grow(Suger::GlcNAc);
    for (int i = 0; i < glycans.size(); i++){
        cout << glycans[i]->get_name() << endl;
        cout << glycans[i]->get_id() << endl;
    }

    // ComplexNGlycan g = nglycan.CreateByAddGlcNAcCore();
    // cout << g.get_table()[0] << endl;

    // ComplexNGlycan q = g.CreateByAddGlcNAcCore();
    // cout << q.get_table()[0] << endl;

    return 0;
}


