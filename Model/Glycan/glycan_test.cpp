#include "glycan.h"
#include <iostream>

using namespace std;
using namespace model::glycan;

void dfs(Monosaccharide* root) {
    std::cout << root-> Name() 
            << std::endl;

    if (root->Child().size() > 0){
        for (auto child : root->Child()){
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

    std:: cout << glycan.Root()->Name()
        << std::endl;

    a.Child().push_back(&b);
    b.set_parent(&a);
    b.Child().push_back(&c);

    dfs(glycan.Root());

    return 0;
}


