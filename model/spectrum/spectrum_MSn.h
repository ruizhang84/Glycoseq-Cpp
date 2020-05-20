#ifndef MODEL_SPECTRUM_SPECTRUM_MSN_H_
#define MODEL_SPECTRUM_SPECTRUM_MSN_H_

#include "spectrum.h"

namespace model {
namespace spectrum {
class SpectrumMSn : public Spectrum{
public:
    SpectrumMSn() = default;
    SpectrumMSn(const SpectrumMSn&);
    SpectrumMSn& operator=(const SpectrumMSn&);

    double PrecursorMZ() { return precursor_mz_; }
    double PrecursorCharge() { return precursor_charge_; }

    void set_parent_mz(double mz) { precursor_mz_ = mz;}
    void set_parent_charge(int charge) { precursor_charge_ = charge; }
private:
    double precursor_mz_;
    int precursor_charge_;
};

} //  namespace spectrum
} //  namespace model




#endif