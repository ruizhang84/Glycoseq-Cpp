#ifndef MODEL_GLYCAN_NGLYCAN_COMPLEX_H
#define MODEL_GLYCAN_NGLYCAN_COMPLEX_H
// The complex type of n-glycan
#include "glycan.h"

namespace model {
namespace glycan {
// class NGlycanComplex : public Glycan 
// {
// public:
//     NGlycanComplex() 
//     { 
//         isNGlycanComplex_ = true;
//     };

//     std::unique_ptr<Glycan> Clone() override; 
//     bool feasibility_check 
//         (Monosaccharide place, Monosaccharide suger) const override;
//     std::vector<std::unique_ptr<Glycan>> Add(Monosaccharide suger) override;

// protected:
//     bool feasibility_check_core 
//         (Monosaccharide place, Monosaccharide suger) const;
//     bool feasibility_check_branch 
//         (Monosaccharide place, Monosaccharide suger) const;
//     bool IsCore() const
//     { 
//         auto it = composite_.find(Monosaccharide::GlcNAc);
//         if (it != composite_.end() && it->second > 2)
//             return false;
//         it = composite_.find(Monosaccharide::Gal);
//         if (it != composite_.end() && it->second >= 3)
//             return false;
//         return true;
//     }
// }; 

}  //  namespace glycan
}  //  namespace model


#endif