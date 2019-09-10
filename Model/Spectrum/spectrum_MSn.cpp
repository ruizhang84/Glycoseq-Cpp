#include "spectrum_MSn.h"

using namespace std;

SpectrumMSn::SpectrumMSn(const SpectrumMSn& spec) : Spectrum(spec){
    activation_type = spec.activation_type;
    parent_mz = spec.parent_mz;
    parent_charge = spec.parent_charge;
}

SpectrumMSn& SpectrumMSn::operator=(const SpectrumMSn& spec){
    Spectrum::operator=(spec);
    this->activation_type = spec.activation_type;
    this->parent_mz = spec.parent_mz;
    this->parent_charge = spec.parent_charge;
    return *this;
}

TypeOfMSActivation SpectrumMSn::get_activation(){
    return activation_type;
}

void SpectrumMSn::set_activation(TypeOfMSActivation type){
    activation_type = type;
}

double SpectrumMSn::get_parent_mz(){
    return parent_mz;
}

void SpectrumMSn::set_parent_mz(double mz){
    parent_mz = mz;
}

int SpectrumMSn::get_parent_charge(){
    return parent_charge;
}

void SpectrumMSn::set_parent_charge(int charge){
    parent_charge = charge;
}

