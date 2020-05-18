#include "glycan.h"
#include <iostream>

using namespace std;
using namespace model::glycan;

void dfs(Monosaccharide* root) {
    std::cout << root-> Name() 
            << std::endl;

    if (root->get_child().size() > 0){
        for (auto child : root->get_child()){
            dfs(child);
        }
    }
}

int main(int argc, char* argv[]){
    Glycan glycan;
    GlcNAc a;
    Man b;
    Gal c;

    Monosaccharide* root = &a;

    glycan.set_tree(root);

    std:: cout << glycan.get_tree()->Name()
        << std::endl;

    a.get_child().push_back(&b);
    b.set_parent(&a);
    b.get_child().push_back(&c);

    dfs(glycan.get_tree());

    return 0;
}


