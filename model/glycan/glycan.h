#ifndef MODEL_GLYCAN_GLYCAN_H
#define MODEL_GLYCAN_GLYCAN_H

#include <string>
#include <vector>
#include <unordered_map> 
#include <memory>
#include <sstream>
#include <algorithm>
#include <iterator>

namespace model {
namespace glycan {

enum class Monosaccharide
{ GlcNAc, Man, Gal, Fuc, NeuAc, NeuGc};

class Glycan
{
public:
    Glycan() = default;
    virtual ~Glycan(){}
    
    virtual std::string Name() const 
    { 
        std::string name = "";
        for (const auto& it : composite_)
        {
            switch (it.first)
            {
            case Monosaccharide::GlcNAc:
                name += " GlcNAc-" + std::to_string(it.second);
                break;
            case Monosaccharide::Man:
                name += " Man-" + std::to_string(it.second);
                break;
            case Monosaccharide::Gal:
                name += " Gal-" + std::to_string(it.second);
                break;
            case Monosaccharide::Fuc:
                name += " Fuc-" + std::to_string(it.second);
                break;    
            case Monosaccharide::NeuAc:
                name += " NeuAc-" + std::to_string(it.second);
                break;
            case Monosaccharide::NeuGc:
                name += " NeuGc-" + std::to_string(it.second);
                break;        
            default:
                break;
            }
        }
        return name;
    } // for print
    virtual std::string ID() const { return Serialize(); }  // use as key

    void set_name(const std::string& name) 
        { name_ = name; }
    void set_id(const std::string& id) 
        { id_ = id; }

    std::vector<int>& Table() { return table_; }
    void set_table(const std::vector<int>& table) 
        { table_ = table; }
    void set_table(int index, int num)
    {
        if (index >= 0 && index < (int) table_.size())
            table_[index] = num;
    }
    virtual std::string Serialize() const
    {
        std::stringstream result;
        std::copy(table_.begin(), table_.end(), 
            std::ostream_iterator<int>(result, " "));
        return result.str();
    }

    virtual void Deserialize(std::string table_str)
    {
        std::istringstream iss(table_str);
        std::string item;
        std::vector<std::string> tokens 
        {
            std::istream_iterator<std::string>{iss}, 
            std::istream_iterator<std::string>{}
        };
        table_.clear();
        for (auto& s : tokens)
        {
            table_.push_back(std::stoi(s));
        }
    }

    std::unordered_map<Monosaccharide, int>&  Composition()
        { return composite_; }

    std::unordered_map<Monosaccharide, int>  CompositionConst() const
        { return std::unordered_map<Monosaccharide, int>(composite_); }

    virtual std::vector<std::unique_ptr<Glycan>> Grow(Monosaccharide suger)
    {
        std::vector<std::unique_ptr<Glycan>> result;
        return result;
    }

protected:
    std::string name_;
    std::string id_;
    std::vector<int> table_;
    std::unordered_map<Monosaccharide, int> composite_; 

};


}  //  namespace glycan
}  //  namespace model

#endif

