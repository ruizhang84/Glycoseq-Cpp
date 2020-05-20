#include "nglycan_complex.h"

namespace model {
namespace glycan {

// std::unique_ptr<Glycan> NGlycanComplex::Clone()
// {
//     return std::make_unique<NGlycanComplex>(*this);
// }

// bool NGlycanComplex::feasibility_check_core
//     (Monosaccharide place, Monosaccharide suger) const
// { 
//     switch (suger)
//     {
//     case Monosaccharide::GlcNAc:
//         if (place != Monosaccharide::GlcNAc)
//             break;
//         auto it = composite_.find(Monosaccharide::GlcNAc);
//         if (it == composite_.end() ||  it->second < 2)
//             return true;
//         break;

//     case Monosaccharide::Gal:
//         if (place != Monosaccharide::GlcNAc)
//             break;
//         auto it = composite_.find(Monosaccharide::Gal);
//         if (it == composite_.end() || it -> second < 3)
//             return true;
//         break;

//     case Monosaccharide::Fuc:
//         if (place != Monosaccharide::GlcNAc)
//             break;
//         auto it = composite_.find(Monosaccharide::GlcNAc);
//         if (it != composite_.end() || it -> second == 1)
//             return true;
//         break;

//     default:
//         break;
//     }

//     return false; 
// }

// bool NGlycanComplex::feasibility_check_branch
//     (Monosaccharide place, Monosaccharide suger) const
// { 
//     return true; 
// }


// bool NGlycanComplex::feasibility_check
//     (Monosaccharide place, Monosaccharide suger) const
// { 
//     if (IsCore())
//         return feasibility_check_core(place, suger);
//     return feasibility_check_branch(place, suger);
// }

// std::vector<std::unique_ptr<Glycan>> NGlycanComplex::Add(Monosaccharide suger) 
// { 
//     std::vector<std::unique_ptr<Glycan>> result; 
//     return result; 
// }

}  //  namespace glycan
}  //  namespace glycan
