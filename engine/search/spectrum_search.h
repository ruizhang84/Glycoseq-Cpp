#ifndef ENGINE_SEARCH_SPECTRUM_MATCH_H
#define ENGINE_SEARCH_SPECTRUM_MATCH_H

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include "precursor_match.h"

#include "../../algorithm/search/bucket_search.h"
#include "../../algorithm/search/binary_search.h"
#include "../../util/mass/peptide.h"
#include "../../model/glycan/glycan.h"
#include "../../model/spectrum/spectrum.h"
#include "../../util/mass/glycan.h"
#include "../../util/mass/spectrum.h"
#include "../../util/mass/ion.h"
#include "../../engine/glycan/glycan_builder.h"
#include "../../engine/protein/protein_ptm.h"

#include <iostream>

namespace engine{
namespace search{

struct SearchResult
{
    int scan;
    std::string peptide;
    std::string glycan;
    int posit;
    double score;
    
};

class SpectrumSearcher
{
public:
    SpectrumSearcher(double tol, algorithm::search::ToleranceBy by, 
        const engine::glycan::GlycanMassStore& subset, const engine::glycan::GlycanStore& isomer):
            tolerance_(tol), by_(by), glycan_mass_(subset), glycan_isomer_(isomer),
                searcher_(algorithm::search::BucketSearch<model::spectrum::Peak>(tol, by)),
                binary_(algorithm::search::BinarySearch(tol, by)){}

    model::spectrum::Spectrum& Spectrum() { return spectrum_; }
    MatchResultStore& Candidate() { return candidate_; }
    void set_spectrum(const model::spectrum::Spectrum& spectrum) { spectrum_ = spectrum; }
    void set_candidate(const MatchResultStore& candidate) { candidate_ = candidate; }

    double Tolerance() const { return tolerance_; }
    algorithm::search::ToleranceBy ToleranceType() const { return by_; }
    void set_tolerance(double tol) 
        { tolerance_ = tol; searcher_.set_tolerance(tol); searcher_.Init(); }
    void set_tolerance_by(algorithm::search::ToleranceBy by) 
        { by_ = by; searcher_.set_tolerance_by(by); searcher_.Init(); }

    virtual std::vector<SearchResult> Search()
    {
        SearchInit();

        std::vector<SearchResult> res;
        
        double max_score = 0;
        SearchResult best;

        for(const auto& peptide : candidate_.Peptides())
        {
            std::vector<model::spectrum::Peak> result_oxonium = SearchOxonium();
            if (result_oxonium.empty()) continue;

            for(const auto& composite: candidate_.Glycans(peptide))
            {
                std::vector<int> positions = engine::protein::ProteinPTM::FindNGlycanSite(peptide);
                std::unordered_map<int, std::vector<model::spectrum::Peak>> result_position;
                for (const auto& pos : positions)
                {
                    std::vector<model::spectrum::Peak> result_temp = SearchPeptides(peptide, composite, pos);
                    if (!result_temp.empty())
                    {
                        result_position[pos] = result_temp;
                    }
                }

                std::unordered_set<std::string> glycan_ids = glycan_isomer_.Query(composite);
                std::unordered_map<std::string, std::vector<model::spectrum::Peak>> result_isomer;
                for(const auto & isomer : glycan_ids)
                {
                    std::vector<model::spectrum::Peak> result_temp = SearchGlycans(peptide, isomer);
                    if (!result_temp.empty())
                    {
                        result_isomer[isomer] = result_temp;
                    }
                }

                if (result_position.empty())
                {
                    for(const auto& isomer_it : result_isomer)
                    {
                        double score = IntensitySum(result_oxonium) + 
                            + IntensitySum(isomer_it.second);
                         
                        if (score > max_score)
                        {
                            max_score = score;
                            best.glycan = composite;
                            best.peptide = peptide;
                            best.score = score;
                        }
                    }
                }

                for(const auto& pos_it : result_position)
                {
                    for(const auto& isomer_it : result_isomer)
                    {
                        double score = IntensitySum(result_oxonium) + 
                            IntensitySum(pos_it.second) + IntensitySum(isomer_it.second);
                         
                        if (score > max_score)
                        {
                            max_score = score;
                            best.glycan = composite;
                            best.peptide = peptide;
                            best.score = score;
                        }
                    }
                }
            }
        }
        if (max_score > 0)
        {
            best.scan = spectrum_.Scan();
            res.push_back(best);
        }
        return res;   
    }

protected:
    virtual void SearchInit()
    {
        std::vector<std::shared_ptr<algorithm::search::Point<model::spectrum::Peak>>> points;
        for(const auto& it : spectrum_.Peaks())
        {
            for (int charge = 1; charge <= spectrum_.PrecursorCharge(); charge++)
            {
                double mass = util::mass::SpectrumMass::Compute(it.MZ(), charge); 
                std::shared_ptr<algorithm::search::Point<model::spectrum::Peak>> p = 
                    std::make_shared<algorithm::search::Point<model::spectrum::Peak>>(mass, it);
                points.push_back(std::move(p));
            }
        }
        searcher_.set_data(std::move(points));
        searcher_.Init();
    }

    virtual std::vector<model::spectrum::Peak> SearchOxonium()
    {
        std::vector<model::spectrum::Peak> res;
        for (const auto& mass : oxonium_)
        {
            std::vector<model::spectrum::Peak> p = searcher_.Query(mass);
            if (! p.empty())
            {
                res.push_back(*std::max_element(p.begin(), p.end(), IntensityCmp));
            }
        }
        return res;
    }

    virtual std::vector<model::spectrum::Peak> SearchPeptides
        (const std::string& seq, const std::string& composite, const int pos)
    {
        std::vector<model::spectrum::Peak> res;
        std::vector<double> peptides_mz;
        std::vector<double> peptides_mass = ComputePTMPeptideMass(seq, pos);
        double extra = util::mass::GlycanMass::Compute(model::glycan::Glycan::Interpret(composite));
        for (const auto& mass : peptides_mass)
        {
            for(int charge = 1; charge <= spectrum_.PrecursorCharge(); charge++)
            {
                double mz = util::mass::SpectrumMass::ComputeMZ(mass + extra, charge);
                peptides_mz.push_back(mz);
            }
        }
        // peptides_mass = ComputeNonePTMPeptideMass(seq, pos);
        // for (const auto& mass : peptides_mass)
        // {
        //     for(int charge = 1; charge <= spectrum_.PrecursorCharge(); charge++)
        //     {
        //         double mz = util::mass::SpectrumMass::ComputeMZ(mass, charge);
        //         peptides_mz.push_back(mz);
        //     }
        // }

        binary_.set_data(peptides_mz);
        binary_.Init();

        for(const auto& pk : spectrum_.Peaks())
        {
            if (binary_.Search(pk.MZ()))
            {
                res.push_back(pk);
            }
        }
        return res;
    }

    virtual std::vector<model::spectrum::Peak> SearchGlycans
        (const std::string& seq, const std::string& id)
    {
        std::vector<model::spectrum::Peak> res;
        std::unordered_set<double> subset = glycan_mass_.Query(id);
        std::vector<double> subset_mass;
        subset_mass.insert(subset_mass.end(), subset.begin(), subset.end());
        binary_.set_data(subset_mass);
        binary_.Init();

        double extra = util::mass::PeptideMass::Compute(seq);
        for(const auto& pk : spectrum_.Peaks())
        {
            for(int charge = 1; charge <= spectrum_.PrecursorCharge(); charge++)
            {
                double mass = util::mass::SpectrumMass::Compute(pk.MZ(), charge) - extra;
                if (binary_.Search(mass))
                {
                    res.push_back(pk);
                    break;
                }
            }
        }
        return res;
    }

    // for computing the peptide ions
    static std::vector<double> ComputePTMPeptideMass(const std::string& seq, const int pos)
    {
        std::vector<double> mass_list;
        for (int i = pos; i < (int) seq.length() - 1; i++) // seldom at n
        {
            double mass = util::mass::IonMass::Compute(seq.substr(0, i+1), util::mass::IonType::c);
            mass_list.push_back(mass);
        }
        for (int i = 1; i <= pos; i++)
        {
            double mass = util::mass::IonMass::Compute(seq.substr(i, seq.length()-i), util::mass::IonType::z);
            mass_list.push_back(mass);
        }
        return mass_list;
    }

    static std::vector<double> ComputeNonePTMPeptideMass(const std::string& seq, const int pos)
    {
        std::vector<double> mass_list;
        for (int i = 0; i < pos; i++) // seldom at n
        {
            double mass = util::mass::IonMass::Compute(seq.substr(0, i+1), util::mass::IonType::c);
            mass_list.push_back(mass);
        }
        for (int i = pos + 1; i < (int) seq.length(); i++)
        {
            double mass = util::mass::IonMass::Compute(seq.substr(i, seq.length()-i), util::mass::IonType::z);
            mass_list.push_back(mass);
        }
        return mass_list;
    }


    static bool IntensityCmp(const model::spectrum::Peak& p1, const model::spectrum::Peak& p2)
        { return (p1.Intensity() < p2.Intensity()); }
    
    static double IntensitySum(const std::vector<model::spectrum::Peak>& peaks)
    { 
        double sum = 0;
        for(const auto& it : peaks)
        {
            sum += it.Intensity();
        }
        return sum;
    }

    double tolerance_;
    algorithm::search::ToleranceBy by_;
    engine::glycan::GlycanMassStore glycan_mass_;
    engine::glycan::GlycanStore glycan_isomer_;
    algorithm::search::BucketSearch<model::spectrum::Peak> searcher_;
    algorithm::search::BinarySearch binary_;
    MatchResultStore candidate_;
    model::spectrum::Spectrum spectrum_;

    const std::vector<double> oxonium_ 
    {
        util::mass::GlycanMass::kHexNAc,
        util::mass::GlycanMass::kHexNAc - util::mass::GlycanMass::kWater,
        util::mass::GlycanMass::kHexNAc - util::mass::GlycanMass::kWater * 2,
        util::mass::GlycanMass::kHexNAc + util::mass::GlycanMass::kHex
    };
}; 

} // namespace engine
} // namespace search

#endif