#ifndef UTIL_IO_SPECTRUM_READER_H_
#define UTIL_IO_SPECTRUM_READER_H_

#include <memory>
#include <string>
#include <vector>
#include "../../model/spectrum/spectrum_MSn.h"

namespace util {
namespace io {

using namespace model::spectrum;

class SpectrumParser
{
public:
    SpectrumParser(std::string path): path_(path){}
    virtual ~SpectrumParser(){}

    virtual double ParentMZ(int scan_num) = 0;
    virtual int ParentCharge(int scan_num) = 0;
    virtual int GetFirstScan() = 0;
    virtual int GetLastScan() = 0;
    virtual SpectrumType GetSpectrumType(int scan_num) = 0;
    virtual std::vector<Peak> Peaks(int scan_num) = 0;
    virtual std::string GetScanInfo(int scan_num) = 0;
    virtual double RTFromScanNum(int scan_num) = 0;
    virtual bool Exist(int scan_num) = 0;
    virtual void Init() = 0;
    
    std::string Path() { return path_; }
    void set_path(std::string path) { path_ = path; Init(); }

protected:
    std::string path_;
};


class SpectrumReader
{
public:
    SpectrumReader(std::string path,
        std::unique_ptr<SpectrumParser> parser):
            path_(path), parser_(std::move(parser)){}
    virtual void Init() {  parser_->Init(); }

    std::string Path() { return path_; }
    void set_path(std::string path) 
        { path_ = path; parser_->set_path(path); }
    void set_parser(std::unique_ptr<SpectrumParser> parser) 
        { parser_ = std::move(parser); }

    virtual int GetFirstScan() { return parser_->GetFirstScan(); }
    virtual int GetLastScan() { return parser_->GetLastScan(); }


    virtual std::string GetScanInfo(int scan_num) 
    {   
        if (! parser_->Exist(scan_num))
            return "Not Exist!";
        return parser_->GetScanInfo(scan_num); 
    }
    
    virtual SpectrumType GetSpectrumType(int scan_num) 
    { 
        if (! parser_->Exist(scan_num))
            return SpectrumType::NONE;
        return parser_->GetSpectrumType(scan_num); 
    }
    virtual Spectrum GetSpectrum(int scan_num)
    {
        SpectrumMSn spectrum;
        if (parser_->Exist(scan_num))
        {
            std::vector<Peak> peaks = parser_->Peaks(scan_num);
            SpectrumType type = GetSpectrumType(scan_num);
            double mz = parser_->ParentMZ(scan_num);
            int charge = parser_->ParentCharge(scan_num);
            spectrum.set_peaks(peaks);
            spectrum.set_scan(scan_num);
            spectrum.set_type(type);
            spectrum.set_parent_mz(mz);
            spectrum.set_parent_charge(charge);
        }
        return spectrum;
    }
    virtual double RTFromScanNum(int scan_num) 
    { 
        if (! parser_->Exist(scan_num))
            return -1;
        return parser_->RTFromScanNum(scan_num); 
    }
    virtual std::vector<Spectrum> GetSpectrum()
    {
        std::vector<Spectrum> result;
        int start = GetFirstScan();
        int last = GetLastScan();
        for (int scan_num = start; scan_num <= last; scan_num++)
        {
            if (parser_->Exist(scan_num))
                result.push_back(GetSpectrum(scan_num));
        }
        return result;
    }

protected:
    std::string path_;
    std::unique_ptr<SpectrumParser> parser_;
    
};

} // namespace io
} // namespace util


#endif