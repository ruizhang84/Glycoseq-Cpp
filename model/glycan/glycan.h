#ifndef MODEL_GLYCAN_GLYCAN_H
#define MODEL_GLYCAN_GLYCAN_H

#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <unordered_map> 
#include <iostream>

namespace model {
namespace glycan {

enum class Monosaccharide
{ GlcNAc, Man, Gal, Fuc, NeuAc, NeuGc};

class Glycan
{
protected:
    std::string uid_;
    std::unordered_map<Monosaccharide, int> composite_; 
    std::vector<Monosaccharide> terminal_;
    bool isNGlycanComplex_;

public:
    Glycan() = default;

    std::string ID() { return uid_; }
    std::vector<Monosaccharide>& Terminal() { return terminal_; }
    std::unordered_map<Monosaccharide, int>& Composite() { return composite_; } 
    void set_uid(std::string uid) { uid_ = uid; }

    // types of glycan
    bool IsNGlycanComplex() { return isNGlycanComplex_; }

    // add monosaccharide moiety
    virtual std::unique_ptr<Glycan> Clone()
    {
        return std::make_unique<Glycan>(*this);
    }
    virtual bool feasilibty_check
        (Monosaccharide place, Monosaccharide suger) const
            { return true; }

    virtual std::vector<std::unique_ptr<Glycan>> Add(Monosaccharide suger) 
    { 
        std::vector<std::unique_ptr<Glycan>> result; 
        return result; 
    }
};


}  //  namespace glycan
}  //  namespace model

#endif

