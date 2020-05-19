#ifndef MODEL_GLYCAN_GLYCAN_H
#define MODEL_GLYCAN_GLYCAN_H

#include <string>
#include <vector>
#include <memory>
#include <utility>

namespace model {
namespace glycan {

enum class SugerType 
{ HexNAc, Hex, Fuc, NeuAc, NeuGc };

struct Monosaccharide 
{
protected:
    std::string name_;
    std::string uid_;
    double mass_;
    SugerType type_;
    std::shared_ptr<Monosaccharide> parent_;
    std::vector<std::shared_ptr<Monosaccharide>> child_;

public:
    Monosaccharide () = default;
    Monosaccharide
    (std::string name, std::string uid, double mass, SugerType type)
    :name_(name), uid_(uid), mass_(mass), type_(type){}
    
    std::string Name() { return name_; }
    std::string ID() { return uid_; }
    double Mass() { return mass_; }
    SugerType Type() { return type_; }
    std::shared_ptr<Monosaccharide>& Parent() { return parent_; }
    std::vector<std::shared_ptr<Monosaccharide>>& Child() { return child_; }

    void set_name(std::string name) { name_ = name; }
    void set_mass(double mass) { mass_ = mass; }
    void set_parent(Monosaccharide& parent) 
    { 
        parent_ = std::make_shared<Monosaccharide>(parent); 
    }
};

class MonosaccharideFactory 
{
public:
    std::shared_ptr<Monosaccharide> GlcNAc() 
    {
        return std::make_shared<Monosaccharide>(Monosaccharide("GlcNAc", "a", 203.0794, SugerType::HexNAc));
    }
 
    std::shared_ptr<Monosaccharide> Man() 
    {
        return std::make_shared<Monosaccharide>(Monosaccharide("Man", "b", 162.0528, SugerType::Hex));
    }

    std::shared_ptr<Monosaccharide> Gal() 
    {
        return std::make_shared<Monosaccharide>(Monosaccharide("Gal", "c", 162.0528, SugerType::Hex));
    }

    std::shared_ptr<Monosaccharide> Fuc() 
    {
        return std::make_shared<Monosaccharide>(Monosaccharide("Fuc", "d", 146.0579, SugerType::Fuc));
    }

    std::shared_ptr<Monosaccharide> NeuAc() 
    {
        return std::make_shared<Monosaccharide>(Monosaccharide("NeuAc", "e", 291.0954, SugerType::NeuAc));
    }

    std::shared_ptr<Monosaccharide> NeuGc() 
    {
        return std::make_shared<Monosaccharide>(Monosaccharide("NeuGc", "f", 307.0903, SugerType::NeuGc));
    }   
};


// class Glycan
// {
// protected:
//     std::string name_ = "";
//     std::string uid_ = "";
//     double mass_ = 0;
//     std::vector<int> composite_ 
//     {0, 0, 0, 0}; // HexAc, Hex, Fuc, NeuAc
//     std::shared_ptr<Monosaccharide> root_;
//     std::vector<std::shared_ptr<Monosaccharide>> terminal_;
//     // types of glycan
//     bool isNGlycanComplex_ = false;

// public:
//     Glycan() = default;

//     std::string Name() { return name_; }
//     std::string ID() { return uid_; }
//     double Mass() { return mass_; }
//     std::shared_ptr<Monosaccharide>& Root() { return root_; }
//     std::vector<std::shared_ptr<Monosaccharide>>& Terminal() { return terminal_; }
//     std::vector<int>& Composite() { return composite_; }

//     void set_name(std::string name) { name_ = name; }
//     void set_mass(double mass) { mass_ = mass; }
//     void set_root(Monosaccharide& root)
//     { 
//         root_ = std::make_shared<Monosaccharide>(root); 
//     }

//     // types of glycan
//     bool isNGlycanComplex() { return isNGlycanComplex_; }
    
//     // add Monosaccharide
//     // virtual bool feasilibty_check
//     // (std::shared_ptr<Monosaccharide>& place, std::shared_ptr<Monosaccharide>& suger)
//     // {
//     //     return true;
//     // }
//     // virtual std::vector<std::shared_ptr<Glycan>> Add(std::shared_ptr<Monosaccharide>& suger) 
//     // { 
//     //     std::vector<std::shared_ptr<Glycan>> result; 
//     //     for(size_t i = 0; i < terminal_.size(); i++)
//     //     {
//     //         if (feasilibty_check(terminal_[i], suger)){
//     //             std::shared_ptr<Glycan> new_glycan = 
//     //                 std::make_shared<Glycan>(*this);
//     //             std::shared_ptr<Monosaccharide> new_suger = 
//     //                 std::make_shared<Monosaccharide>(*suger);

//     //             new_suger->set_parent(new_glycan->terminal_[i]);
//     //             new_glycan->terminal_[i]->Child().push_back(std::make_shared<Monosaccharide>(*suger));
//     //             new_glycan->terminal_[i] = new_glycan->terminal_[i]->Child().back();
//     //             result.push_back(std::move(new_glycan));
//     //         }
//     //     }
//     //     return result; 
//     // }
// };


}  //  namespace glycan
}  //  namespace model

#endif

