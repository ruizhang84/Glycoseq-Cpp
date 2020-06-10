#ifndef ENGINE_SEARCH_DECOY_SEARCH_H
#define ENGINE_SEARCH_DECOY_SEARCH_H

#include "spectrum_search.h"

namespace engine{
namespace search{

class DecoySearcher : public SpectrumSearcher
{
public:
    DecoySearcher(const double tol, const algorithm::search::ToleranceBy by, int isotope,
        engine::glycan::NGlycanBuilder* builder): SpectrumSearcher(tol, by, isotope, builder){};

    std::vector<SearchResult> Search() override
    {
        SearchInit();

        std::vector<SearchResult> res;
        std::unordered_map<std::string, SearchResult> res_map;
        
        int peak_size = (int) spectrum_.Peaks().size();
        for(const auto& peptide : candidate_.Peptides())
        {
            std::vector<model::spectrum::Peak> result_oxonium = SearchOxonium();
            if (result_oxonium.empty()) continue;
            double oxonium = SearchResult::PeakValue(result_oxonium);
            int oxonium_matched = 0;
            oxonium_matched += result_oxonium.size();

            for(const auto& composite: candidate_.Glycans(peptide))
            {
                std::vector<int> positions = engine::protein::ProteinPTM::FindNGlycanSite(peptide);
                std::unordered_map<int, double> result_position;
                std::unordered_map<int, int> matched_position;
                for (const auto& pos : positions)
                {
                    std::vector<model::spectrum::Peak> result_temp = SearchPeptides(peptide, composite, pos);
                    if (!result_temp.empty())
                    {  
                        result_position[pos] = SearchResult::PeakValue(result_temp);
                        matched_position[pos] = (int) result_temp.size();
                    }
                }
                if (result_position.empty()) continue;

                std::unordered_set<std::string> glycan_ids = glycan_isomer_.Query(composite);
                std::unordered_map<std::string, double> result_core, result_branch, result_terminal;
                std::unordered_map<std::string, int> matched_isomer;

                for(const auto & isomer : glycan_ids)
                {
                    std::vector<model::spectrum::Peak> result_temp;
                    result_temp = SearchGlycans(peptide, isomer, glycan_core_);
                    if (result_temp.empty()) continue;
                    
                    result_core[isomer] = SearchResult::PeakValue(result_temp);
                    matched_isomer[isomer] = (int) result_temp.size();
                    result_temp = SearchGlycans(peptide, isomer, glycan_branch_);
                    result_branch[isomer] = SearchResult::PeakValue(result_temp);
                    matched_isomer[isomer] += (int) result_temp.size();
                    result_temp = SearchGlycans(peptide, isomer, glycan_terminal_);
                    result_terminal[isomer] = SearchResult::PeakValue(result_temp);
                    matched_isomer[isomer] += (int) result_temp.size(); 
                }


                for(const auto& pos_it : result_position)
                {
                    for(const auto& isomer_it : result_core)
                    {
                        std::string isomer = isomer_it.first;
                        double score = result_core[isomer] + result_branch[isomer] + result_terminal[isomer];
                        if (score == 0) continue;
                        score  += pos_it.second + oxonium; 

                       
                        std::string glycopeptide = peptide + composite;
                        if (res_map.find(glycopeptide) == res_map.end())
                        {
                            SearchResult best;
                            best.set_scan(spectrum_.Scan());
                            best.set_glycan(composite);
                            best.set_sequence(peptide);
                            best.set_score(score);
                            res_map[glycopeptide] = best;
                        }
                        
                        SearchResult& best = res_map[glycopeptide];
                        if (score <= best.Score()) continue;

                        best.set_site(pos_it.first);
                        best.Add(oxonium, SearchType::Oxonium);
                        best.Add(pos_it.second, SearchType::Peptide);
                        best.Add(result_core[isomer], SearchType::Core);
                        best.Add(result_branch[isomer], SearchType::Branch);
                        best.Add(result_terminal[isomer], SearchType::Terminal);
                        double matched = (oxonium_matched + matched_position[pos_it.first]
                            +  matched_isomer[isomer] + 0.0) / peak_size;
                        best.Add(matched, SearchType::Matches);
                    }
                    
                }
            }
        }
        if (res_map.empty())
            return res;
        

        // compute precursor differ
        double precursor_mass = 
            util::mass::SpectrumMass::Compute(spectrum_.PrecursorMZ(), spectrum_.PrecursorCharge());
        for(auto& it : res_map)
        {
            SearchResult best = it.second;
            double mass = util::mass::PeptideMass::Compute(best.Sequence())
                + util::mass::GlycanMass::Compute(model::glycan::Glycan::Interpret(best.Glycan()));
            
            double ppm = kPPM;
            for (int i = 0; i <= isotopic_; i ++)
            {
                double isotopic_mass = util::mass::SpectrumMass::kIon * i + mass;
                double isotopic_ppm = util::mass::SpectrumMass::ComputePPM(isotopic_mass, precursor_mass);
                ppm = (ppm <= isotopic_ppm) ? ppm : isotopic_ppm;
            }
            double ratio = 1.0 - ppm / kPPM;
            it.second.Add(ratio, SearchType::Precursor);
        }
        
        // save 
        for(const auto& it : res_map)
        {
            res.push_back(it.second);
        }

        // remove if hits over 20
        if (res.size() > 20)
        {
            std::sort(res.begin(), res.end(), 
                [](const SearchResult& r1, const SearchResult& r2) -> bool { return r1.Score() > r1.Score(); });
            res.erase(res.begin() + 20);
        }
  
        return res;   
    }

};

}  // namespace search
}  // namespace engine

#endif