#include "glycan.h"

using namespace std;

vector<int> Glycan::get_table(){
    return table;
}

void Glycan::set_table(std::vector<int>& table){
    this->table = table;
}

void Glycan::set_table(int index, int num){
    if (index >= 0 && index < table.size()){
        table[index] = num;
    }
}