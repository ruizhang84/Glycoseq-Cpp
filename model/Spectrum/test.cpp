#include "spectrum_MSn.h"
#include "peak.h"
#include <iostream>
#include <fstream>
#include <string>
#include <regex>

using namespace std;

int main(){
    vector<SpectrumMSn> spectra;
    SpectrumMSn spec; 

    ifstream file("test_CID.mgf");
    string line;
    
    regex start("BEGIN\\s+IONS");
    regex end("END\\s+IONS");
    regex title("TITLE=File:");
    regex ms("PEPMASS=(\\d+\\.?\\d*)\\s+(\\d+\\.?\\d*)");
    regex chr("CHARGE=(\\d+)");
    regex rt("RTINSECONDS=(\\d+)");
    regex sc("SCANS=(\\d+)");
    regex pk("^(\\d+\\.?\\d*)\\s+(\\d+\\.?\\d*)");
    smatch result;

    if (file.is_open()){
        while(getline(file, line)){
            if (regex_search(line, result, start)){
                spec = SpectrumMSn(); 
                spec.set_MSn_order(2);
                spec.set_activation(TypeOfMSActivation::CID);
            }else if (regex_search(line, result, pk)){
                Peak pk(stod(result[1]), stod(result[2]));
                spec.Add(pk);
            }else if (regex_search(line, result, ms)){
                spec.set_parent_mz(stod(result[1]));
            }else if (regex_search(line, result, chr)){
                spec.set_parent_charge(stoi(result[1]));
            }else if (regex_search(line, result, sc)){
                spec.set_scan_num(stoi(result[1]));
            }else if (regex_search(line, result, end)){
                spectra.push_back(spec);
            } 
        }
    }

    for(vector<SpectrumMSn>::iterator it = spectra.begin();
        it != spectra.end(); it++){
            cout << it->get_scan_num() << endl;
            vector<Peak> peaks = it->get_peaks();
            for(int i = 0; i < peaks.size(); i++){
                cout << peaks[i].get_mz() << endl;
            }
        }

    return 0;
}