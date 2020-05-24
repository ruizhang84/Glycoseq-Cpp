#ifndef MODEL_GLYCAN_GLYCAN_H
#define MODEL_GLYCAN_GLYCAN_H

#include <string>
#include <vector>
#include <unordered_map> 
#include <memory>

namespace model {
namespace glycan {

enum class Monosaccharide
{ GlcNAc, Man, Gal, Fuc, NeuAc, NeuGc};

class Glycan
{
public:
    Glycan() = default;
    virtual ~Glycan(){}
    
    virtual std::string Name() const { return name_;}
    void set_name(const std::string& name) 
        { name_ = name; }
    virtual std::string ID() const { return id_; }
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

