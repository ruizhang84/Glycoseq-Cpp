#ifndef UTIL_MASS_GLYCAN_H
#define UTIL_MASS_GLYCAN_H

#include "../../model/glycan/glycan.h"

namespace util {
namespace mass {

class GlycanMass
{
public:
    static double Compute(model::glycan::Glycan& glycan) 
    {
        double mass = 0;
        for(auto& it : glycan.Composite())
        {
            switch (it.first)
            {
            case model::glycan::Monosaccharide::GlcNAc:
                mass += kHexNAc * it.second;
                break;
            case model::glycan::Monosaccharide::Gal:
                mass += kHex * it.second;
                break;     
            case model::glycan::Monosaccharide::Man:
                mass += kHex * it.second;
                break;     
            case model::glycan::Monosaccharide::Fuc:
                mass += kFuc * it.second;
                break;     
            case model::glycan::Monosaccharide::NeuAc:
                mass += kNeuAc * it.second;
                break;           
            case model::glycan::Monosaccharide::NeuGc:
                mass += kNeuGc * it.second;
                break;    
            default:
                break;
            }
        }
        return mass;
    }

protected:
    static const double kHexNAc = 203.0794;
    static const double kHex = 162.0528;
    static const double kFuc = 146.0579;
    static const double kNeuAc = 291.0954;
    static const double kNeuGc = 307.0903;
};



} // namespace mass
} // namespace util

#endif