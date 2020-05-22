#include "spectrum_binpacking.h"

namespace engine {
namespace spectrum {

using namespace model::spectrum;

std::vector<double> SpectrumBinPacking::Packing
    (model::spectrum::Spectrum& spec)
{
    std::vector<double> result;
    std::vector<Peak> peaks = spec.Peaks();
    std::vector<std::vector<Peak>> peak_bins = BinPacking::Packing(peaks); 

    return result;
}

} // namespace spectrum
} // namespace engine