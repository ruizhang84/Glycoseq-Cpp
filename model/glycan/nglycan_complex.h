#ifndef MODEL_GLYCAN_NGLYCAN_COMPLEX_H
#define MODEL_GLYCAN_NGLYCAN_COMPLEX_H
// The complex type of n-glycan
#include "glycan.h"

namespace model {
namespace glycan {
class NGlycanComplex : public Glycan 
{
public:
    NGlycanComplex() { isNGlycanComplex_ = true; };

    std::vector<Glycan*> Add(Monosaccharide* suger) override;
}; 

}  //  namespace glycan
}  //  namespace model


#endif