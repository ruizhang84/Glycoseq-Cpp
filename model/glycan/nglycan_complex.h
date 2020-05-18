#ifndef MODEL_GLYCAN_NGLYCAN_COMPLEX_H
#define MODEL_GLYCAN_NGLYCAN_COMPLEX_H
// The complex type of n-glycan
#include "glycan.h"

namespace model {
namespace glycan {
class NGlycanComplex : public Glycan 
{
public:
    NGlycanComplex() 
    { 
        isNGlycanComplex_ = true;
        for (int i = 0; i < kComposite; i++){
            composite_.push_back(0);
        } 
    };
    std::vector<Glycan*> Add(Monosaccharide* suger) override;

protected:
    void UpdateMass(Monosaccharide* suger);
    void UpdateComposition(Monosaccharide* suger);
    const int kComposite = 4; // HexAc, Hex, Fuc, NeuAc
}; 

}  //  namespace glycan
}  //  namespace model


#endif