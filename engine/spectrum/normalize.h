#ifndef ENGINE_SPECTRUM_NORMALIZE_H
#define ENGINE_SPECTRUM_NORMALIZE_H

#include <vector>
#include <algorithm>
#include <numeric> 
#include "../../model/spectrum/spectrum.h"

namespace engine {
namespace spectrum {
 
class Normalizer
{
public:
    static void Transform(model::spectrum::Spectrum& spec)
    {   
        Transform(spec.Peaks());
    }   

    //normalization on total ion intensity sums 
    static void Transform(std::vector<model::spectrum::Peak>& peaks)
    {
        // double max_intensity = std::max_element(peaks.begin(), peaks.end(), IntensityCmp)->Intensity();
        // for(auto& it : peaks)
        // {
        //     it.set_intensity(it.Intensity() / max_intensity);
        // }
        double sum = std::accumulate(peaks.begin(), peaks.end(), 0, IntensitySum);
        for(auto& it : peaks)
        {
            it.set_intensity(it.Intensity() / sum);
        }

    }

protected:
    static bool IntensityCmp (model::spectrum::Peak& i, model::spectrum::Peak& j) 
        { return (i.Intensity() < j.Intensity()); }
    static double IntensitySum (model::spectrum::Peak& i, model::spectrum::Peak& j) 
        { return i.Intensity() + j.Intensity(); }

};

} // namespace spectrum
} // namespace engine

#endif