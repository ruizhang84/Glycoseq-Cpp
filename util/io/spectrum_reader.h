#ifndef UTIL_IO_SPECTRUM_READER_H_
#define UTIL_IO_SPECTRUM_READER_H_

#include <memory>
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
    virtual double ParentMZ(int scan_num) = 0;
    virtual int ParentCharge(int scan_num) = 0;
    virtual int GetFirstScan() = 0;
    virtual int GetLastScan() = 0;
    virtual SpectrumType GetSpectrumType(int scan_num) = 0;
    virtual std::vector<Peak> Peaks(int scan_num) = 0;
    
    std::string Path() { return path_; }
    void set_path(std::string path) { path_ = path; }

protected:
    std::string path_;
};


class SpectrumReader
{
public:
    SpectrumReader(std::string path, SpectrumParser* parser):
        path_(path), parser_(std::move(parser)){}

    std::string Path() { return path_; }
    void set_path(std::string path) { path_ = path; }
    void set_parser(SpectrumParser* parser) { parser_ = std::move(parser_); }

    int GetFirstScan() { return parser_->GetFirstScan(); }
    int GetLastScan() { return parser_->GetLastScan(); }
    SpectrumType GetSpectrumType(int scan_num) 
        { return parser_->GetSpectrumType(scan_num); }
    virtual Spectrum GetSpectrum(int scan_num);
    virtual std::vector<Spectrum> GetSpectrum();

protected:
    std::unique_ptr<SpectrumParser> parser_;
    std::string path_;

};

} // namespace io
} // namespace util


#endif