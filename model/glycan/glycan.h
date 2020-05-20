#ifndef MODEL_GLYCAN_GLYCAN_H
#define MODEL_GLYCAN_GLYCAN_H

#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <unordered_map> 

namespace model {
namespace glycan {

enum class Monosaccharide
{ GlcNAc, Man, Gal, Fuc, NeuAc, NeuGc};

class Moiety
{
public:
    Moiety(Monosaccharide name): name_(name){}
    Moiety(const Moiety& other)
    {
        name_ = other.name_;
        parent_ = other.parent_;
        for(auto& child: other.children_)
        {
            std::unique_ptr<Moiety> new_child = child->Clone();
            new_child->set_parent(this);
            children_.push_back(std::move(new_child));
        }
    }

    Monosaccharide Name() { return name_; }
    Moiety* Parent() { return parent_; }
    void set_parent(Moiety* parent) { parent_ = parent; }
    std::vector<std::unique_ptr<Moiety>>& Children() { return children_; } 

    std::unique_ptr<Moiety> Clone()
    {
        std::unique_ptr<Moiety> root_ = 
            std::make_unique<Moiety>(name_);

        for(auto& child: children_)
        {
            std::unique_ptr<Moiety> new_child = child->Clone();
            new_child->set_parent(root_.get());
            root_->Children().push_back(std::move(new_child));
            
        }
        return root_;
    }

protected:
    Monosaccharide name_;
    Moiety* parent_;
    std::vector<std::unique_ptr<Moiety>> children_;
};

class Glycan
{
public:
    Glycan() = default;
    Glycan(const Glycan& other)
    {
        composite_ = other.composite_;
        isNGlycanComplex_ = other.isNGlycanComplex_;
        root_ = std::move(other.root_->Clone());
    }
    
    Moiety* Root() { return root_.get(); }
    void set_root(std::unique_ptr<Moiety>& root)
    {
        root_ = std::move(root);
    }
    std::unordered_map<Monosaccharide, int>& Composite() { return composite_; } 
    std::unique_ptr<Glycan> Clone()
    {
        return std::make_unique<Glycan>(*this);
    }

    // types of glycan
    bool IsNGlycanComplex() { return isNGlycanComplex_; }

    // add monosaccharide moiety
    virtual std::vector<std::unique_ptr<Glycan>> Add(Monosaccharide suger) 
    { 
        std::vector<std::unique_ptr<Glycan>> result; 
        return result; 
    }

protected:
    std::unique_ptr<Moiety> root_;
    std::unordered_map<Monosaccharide, int> composite_; 
    bool isNGlycanComplex_;
};


}  //  namespace glycan
}  //  namespace model

#endif

