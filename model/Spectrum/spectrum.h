#ifndef MODEL_SPECTRUM_H_
#define MODEL_SPECTRUM_H_

#include <vector>
#include "peak.h"


class Spectrum{
public:
    Spectrum() = default;
    Spectrum(const Spectrum&);
    Spectrum& operator=(const Spectrum&);
    std::vector<Peak> get_peaks();
    void Add(Peak peak);
    void Clear();
    int get_MSn_order();
    void set_MSn_order(int MSn_order);
    int get_scan_num();
    void set_scan_num(int scan_num);

protected:
    std::vector<Peak> peaks;
    int MSn_order;
    int scan_num;
};

#endif