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
    
    virtual std::string Name() const { return name_;} // for print
    void set_name(const std::string& name) 
        { name_ = name; }
    virtual std::string ID() const { return id_; }  // use as key
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

