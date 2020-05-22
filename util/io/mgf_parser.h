#ifndef UTIL_IO_MGF_PARSER_H_
#define UTIL_IO_MGF_PARSER_H_

#include <string>
#include <map> 
#include "spectrum_reader.h"

namespace util {
namespace io {

class MGFParser : public SpectrumParser
{   
public:
    MGFParser(std::string path, SpectrumType type): 
        SpectrumParser(path), type_(type){}

    void Init() override;

    double ParentMZ(int scan_num) override 
    { 
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.pep_mass;
        }
        return 0;
    }
    int ParentCharge(int scan_num) override
    { 
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.charge;
        }
        return 0;
    }
    int GetFirstScan() override 
    {
        auto it = data_set_.begin();
        if (it != data_set_.end())
        {
            return it->first;
        }
        return -1;
    }
    int GetLastScan() override
    {
        auto it = data_set_.rbegin();
        if (it != data_set_.rend())
        {
            return it->first;
        }
        return -1;
    }
    std::vector<Peak> Peaks(int scan_num) override
    { 
        std::vector<Peak> peaks;
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            for(size_t i = 0; 
                i < it->second.intensity.size(); i++)
            {
                Peak pk(it->second.mz[i], 
                        it->second.intensity[i]);
                peaks.push_back(pk);
            }
        }
        return peaks;
    }
    SpectrumType GetSpectrumType(int scan_num) override 
        { return type_; };
    double RTFromScanNum(int scan_num) override
    {
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.rt_seconds;
        }
        return -1;
    }
    std::string GetScanInfo(int scan_num) override
    {
        auto it = data_set_.find(scan_num); 
        if (it != data_set_.end())
        {
            return it->second.title;
        }
        return "";
    }
    bool Exist(int scan_num) override
    {
        return data_set_.find(scan_num) != data_set_.end();
    }
    
private:
    class MGFData
    {
    public:
        MGFData() = default;

        std::vector<double> mz;
        std::vector<double> intensity;
        double pep_mass;
        int charge;
        double rt_seconds;
        int scans;
        std::string title;
    };
    SpectrumType type_;
    std::map<int, MGFData> data_set_;
};


} // namespace io
} // namespace util


#endif