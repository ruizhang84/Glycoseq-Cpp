#include "peak.h"

Peak::Peak(double mz, double intensity): mz(mz), intensity(intensity){}

Peak::Peak(const Peak& peak){
    mz = peak.mz;
    intensity = peak.intensity;
}

Peak& Peak::operator=(const Peak& peak){
    this->mz = peak.mz;
    this->intensity = peak.intensity;
    return *this;
}

double Peak::get_mz(){
    return mz;
}

void Peak::set_mz(double mz){
    this->mz = mz;
}

double Peak::get_intensity(){
    return intensity;
}

void Peak::set_intensity(double intensity){
    this->intensity = intensity;
}