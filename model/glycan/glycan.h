#ifndef MODEL_GLYCAN_GLYCAN_H
#define MODEL_GLYCAN_GLYCAN_H

#include <string>
#include <vector>
#include <memory>
#include <utility>

namespace model {
namespace glycan {

struct Monosaccharide 
{
protected:
    std::string name_;
    std::string uid_;
    double mass_;

    Monosaccharide* parent_ = NULL;
    std::vector<Monosaccharide*> child_;
public:
    Monosaccharide() = default;
    
    std::string Name() { return name_; }
    std::string ID() { return uid_; }
    double Mass() { return mass_; }
    Monosaccharide* Parent() { return parent_; }
    std::vector<Monosaccharide*>& Child() { return child_; }

    void set_name(std::string name) { name_ = name; }
    void set_mass(double mass) { mass_ = mass; }
    void set_parent(Monosaccharide* parent)
    { 
        parent_ = parent_;
    }

};

class Glycan
{
protected:
    std::string name_;
    std::string uid_;
    double mass_;
    std::vector<int> composite_;
    Monosaccharide* root_ = NULL;
    // types of glycan
    bool isNGlycanComplex_ = false;

public:
    Glycan() = default;

    std::string Name() { return name_; }
    std::string ID() { return uid_; }
    double Mass() { return mass_; }
    std::vector<int>& Composite() { return composite_; }
    Monosaccharide* Root() { return root_; }

    void set_name(std::string name) { name_ = name; }
    void set_mass(double mass) { mass_ = mass; }
    void set_root(Monosaccharide* root) { root_ = root;}
    void set_composite(std::vector<int>& composite){ 
        composite_ = std::move(composite);
    }

    // types of glycan
    bool isNGlycanComplex() { return isNGlycanComplex_; }
    
    // add Monosaccharide
    virtual std::vector<Glycan*> Add(Monosaccharide* suger) 
    { 
        std::vector<Glycan*> result; 
        return result; 
    }
};

struct GlcNAc : public Monosaccharide 
{
    GlcNAc() { name_ = "GlycNAc", mass_ = 203.0794; uid_ = "a"; }
};

struct Man : public Monosaccharide 
{
    Man(){ name_ = "Man", mass_ = 162.0528; uid_ = "b"; }
};

struct Gal : public Monosaccharide 
{
    Gal(){ name_ = "Gal", mass_ = 162.0528; uid_ = "c"; } 
};

struct Fuc : public Monosaccharide 
{
    Fuc(){ name_ = "Fuc", mass_ = 146.0579; uid_ = "d"; }
};

struct NeuAc : public Monosaccharide 
{
    NeuAc(){ name_ = "NeuAc", mass_ = 291.0954; uid_ = "e"; } 
};

struct NeuGc : public Monosaccharide 
{
    NeuGc(){ name_ = "NeuGc", mass_ = 307.0903; uid_ = "f"; }
};

}  //  namespace glycan
}  //  namespace model

#endif

