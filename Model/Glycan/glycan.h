#ifndef MODEL_GLYCAN_GLYCAN_H
#define MODEL_GLYCAN_GLYCAN_H

#include <string>
#include <vector>
#include <memory>
#include <utility>

namespace model {
namespace glycan {

struct Monosaccharide {
protected:
    std::string name_;
    double mass_;
    Monosaccharide* parent_ = NULL;
    std::vector<Monosaccharide*> child_;
public:
    Monosaccharide() = default;
    
    std::string Name() { return name_; }
    double Mass() { return mass_; }
    void set_name(std::string name) { name_ = name; }
    void set_mass(double mass) { mass_ = mass; }
    std::vector<Monosaccharide*>& Child() { return child_; }
    Monosaccharide* Parent() { return parent_; }
    void set_parent(Monosaccharide* parent){ 
        parent_ = parent_;
    }

};

class Glycan
{
protected:
    std::string name_;
    std::vector<int> composite_;
    Monosaccharide* root_;

public:
    Glycan() = default;
    // get tree nodes
    Monosaccharide* Root() { return root_; }
    void set_tree(Monosaccharide* root) { root_ = root;}
    // set composition
    void set_composite(std::vector<int>& composite){ 
        composite_ = std::move(composite);
    }
    bool set_composite(int index, int num){
        if (index >= 0 && index < composite_.size() && num >= 0){
            composite_[index] = num;
            return true;
        }
        return false;
    }
    
    virtual std::vector<int>& Composite() { return composite_; }
    virtual std::string Name() { return name_; }
};

struct GlcNAc : public Monosaccharide {
    GlcNAc() { name_ = "GlycNAc", mass_ = 203.0794; }
};

struct Man : public Monosaccharide {
    Man(){ name_ = "Man", mass_ = 162.0528; }
};

struct Gal : public Monosaccharide {
    Gal(){ name_ = "Gal", mass_ = 162.0528; } 
};

struct Fuc : public Monosaccharide {
    Fuc(){ name_ = "Fuc", mass_ = 146.0579; }
};

struct NeuAc : public Monosaccharide {
    NeuAc(){ name_ = "NeuAc", mass_ = 291.0954; } 
};

struct NeuGc : public Monosaccharide {
    NeuGc(){ name_ = "NeuGc", mass_ = 307.0903; }
};

}  //  namespace glycan
}  //  namespace model

#endif

