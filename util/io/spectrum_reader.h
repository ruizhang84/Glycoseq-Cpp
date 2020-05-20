#ifndef UTIL_IO_READ_MGF_H_
#define UTIL_IO_READ_MGF_H_

#include <string>
#include <vector>
#include "../../model/spectrum/peak.h"
#include "../../model/spectrum/spectrum_MSn.h"

namespace util {
namespace io {

using namespace model::spectrum;

class SpectrumParser
{
public:
    double ParentMZ(int scan_num);
    int ParentCharge(int scan_num);
    std::vector<Peak> Peaks(int scan_num);
    
    std::string Path() { return path_; }
    void set_path(std::string path) { path_ = path; }

protected:
    std::string path_;
};


class SpectrumReader
{
public:
    SpectrumReader(std::string path, SpectrumParser parser);

    std::string Path() { return path_; }
    void set_path(std::string path) { path_ = path; }
    void set_parser(SpectrumParser parser) { parser_ = parser_; }

    int GetFirstScan();
    int GetLastScan();
    SpectrumType GetSpectrumType(int scan_num);
    Spectrum GetSpectrum(int scan_num);

protected:
    SpectrumParser parser_;
    std::string path_;

};

} // namespace io
} // namespace util


#endif