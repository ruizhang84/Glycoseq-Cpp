#ifndef MODEL_SPECTRUM_SPECTRUM_MSN_H_
#define MODEL_SPECTRUM_SPECTRUM_MSN_H_

#include "spectrum.h"

enum TypeOfMSActivation { CID, MPD, ECD, PQD, ETD, HCD, Any, SA, PTR, NETD, NPTR };

class SpectrumMSn : public Spectrum{
public:
    SpectrumMSn() = default;
    SpectrumMSn(const SpectrumMSn&);
    SpectrumMSn& operator=(const SpectrumMSn&);
    TypeOfMSActivation get_activation();
    void set_activation(TypeOfMSActivation activation_type);
    double get_parent_mz();
    void set_parent_mz(double mz);
    int get_parent_charge();
    void set_parent_charge(int charge);
private:
    TypeOfMSActivation activation_type;
    double parent_mz;
    int parent_charge;
};

#endif