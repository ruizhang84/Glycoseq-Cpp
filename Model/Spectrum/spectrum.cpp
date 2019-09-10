#include "spectrum.h"

using namespace std;

Spectrum::Spectrum(const Spectrum& spectrum){
    peaks = spectrum.peaks;
    MSn_order = spectrum.MSn_order;
    scan_num = spectrum.scan_num;
}

Spectrum& Spectrum::operator=(const Spectrum& spectrum){
    this->peaks = spectrum.peaks;
    this->MSn_order = spectrum.MSn_order;
    this->scan_num = spectrum.scan_num;
    return *this;
}

vector<Peak> Spectrum::get_peaks(){
    return peaks;
}

void Spectrum::Add(Peak peak){
    peaks.push_back(peak);
}

void Spectrum::Clear(){
    peaks.clear();
}

int Spectrum::get_MSn_order(){
    return MSn_order;
}

void Spectrum::set_MSn_order(int order){
    MSn_order = order;
}

int Spectrum::get_scan_num(){
    return scan_num;
}

void Spectrum::set_scan_num(int num){
    scan_num = num;
}

