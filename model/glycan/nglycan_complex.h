#ifndef MODEL_GLYCAN_NGLYCAN_COMPLEX_H
#define MODEL_GLYCAN_NGLYCAN_COMPLEX_H
// The complex type of n-glycan

#include <sstream>
#include <iterator>
#include <typeinfo>
#include "glycan.h"


namespace model {
namespace glycan {
class NGlycanComplex : public Glycan 
{
public:
    NGlycanComplex()
    { 
        table_.assign(24, 0);
    }

    std::string Name() const override
    {
        std::string name = "NGlycanComplex: ";
        if (table_[2] > 0)
            name += "fucose ";
        if (table_[3] > 0)
            name += "bisect ";
        for (auto& it : composite_)
        {
            name += typeid(it.first).name();
            name +=  "-" + std::to_string(it.second) + " ";
        }
        return name;
    }

    std::string ID() const override
    {
        std::stringstream result;
        std::copy(table_.begin(), table_.end(), 
            std::ostream_iterator<int>(result, " "));
        return result.str();
    }

    std::vector<std::unique_ptr<Glycan>> Grow(Monosaccharide suger) override;

protected:
    bool ValidAddGlcNAcCore();
    bool ValidAddGlcNAc();
    NGlycanComplex CreateByAddGlcNAcCore();

    bool ValidAddGlcNAcBisect();
    NGlycanComplex CreateByAddGlcNAcBisect();

    bool ValidAddGlcNAcBranch();
    std::vector<NGlycanComplex> CreateByAddGlcNAcBranch();

    bool ValidAddMan();
    NGlycanComplex CreateByAddMan();

    bool ValidAddGal();
    std::vector<NGlycanComplex> CreateByAddGal();

    bool ValidAddFucCore();
    NGlycanComplex CreateByAddFucCore();

    bool ValidAddFucTerminal();
    std::vector<NGlycanComplex> CreateByAddFucTerminal();

    bool ValidAddNeuAc();
    std::vector<NGlycanComplex> CreateByAddNeuAc();

    bool ValidAddNeuGc();
    std::vector<NGlycanComplex> CreateByAddNeuGc();

}; 

}  //  namespace glycan
}  //  namespace model


#endif