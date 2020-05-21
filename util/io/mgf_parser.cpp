#include <fstream>
#include <regex>
#include "mgf_parser.h"

namespace util {
namespace io {

void MGFParser::Init()
{
    MGFData data;
    int scan_num = -1;

    std::ifstream file(path_);
    std::string line;

    std::smatch result;
    std::regex start("BEGIN\\s+IONS");
    std::regex end("END\\s+IONS");
    std::regex title("TITLE=(.*)");
    std::regex pepmass("PEPMASS=(\\d+\\.?\\d*)");
    std::regex charge("CHARGE=(\\d+)");
    std::regex rt_second("RTINSECONDS=(\\d+\\.?\\d*)");
    std::regex scan("SCANS=(\\d+)");
    std::regex mz_intensity("^(\\d+\\.?\\d*)\\s+(\\d+\\.?\\d*)");

    if (file.is_open()){
        while(std::getline(file, line)){
            if (regex_search(line, result, start))
            {
                data = MGFData();
                scan_num++;
            }else if (regex_search(line, result, mz_intensity))
            {
                data.mz.push_back(stod(result[1]));
                data.intensity.push_back( stod(result[2]));
            }
            else if (regex_search(line, result, pepmass))
            {
                data.pep_mass = stod(result[1]);
            }
            else if (regex_search(line, result, charge)){
                data.charge = stoi(result[1]);
            }
            else if (regex_search(line, result, scan))
            {
                scan_num = stoi(result[1]);
                data.scans = scan_num;
            }
            else if (regex_search(line, result, title))
            {
                data.title = std::string(result[1]);
            }
            else if (regex_search(line, result, rt_second))
            {
                data.rt_seconds = stod(result[1]);
            }
            else if (regex_search(line, result, end))
            {
                data_set_.emplace(scan_num, data);
            } 
        }
    }
}

} // namespace io
} // namespace util

