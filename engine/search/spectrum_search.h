#ifndef ENGINE_SEARCH_SPECTRUM_MATCH_H
#define ENGINE_SEARCH_SPECTRUM_MATCH_H

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
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

namespace engine{
namespace search{

struct SearchResult
{
    std::string peptide;
    std::string glycan;
    int posit;
    double score;
};

class GlycanMassStore
{
public:
    GlycanMassStore(engine::glycan::GlycanStore subset_store){ Init(subset_store); }
    std::unordered_set<double> Mass(const std::string& name)
    {
        if(map_.find(name) != map_.end())
            return map_[name];
        return std::unordered_set<double>();
    }
    std::unordered_map<std::string, std::unordered_set<double>>& Map()
        { return map_; }

protected:
    virtual void Init(engine::glycan::GlycanStore& subset_store)
    {
        for(const auto& it : subset_store.Map())
        {
            map_[it.first] = std::unordered_set<double>();
            for(const auto& s : it.second)
            {
                if (memory_.find(s) == memory_.end())
                {
                    memory_[s] = util::mass::GlycanMass::Compute(
                        model::glycan::NGlycanComplex::InterpretID(s));
                }
                map_[it.first].insert(memory_[s]);
            }
        }

    }
    std::unordered_map<std::string, double> memory_;
    std::unordered_map<std::string, std::unordered_set<double>> map_;
};


class SpectrumSearcher
{
public:
    SpectrumSearcher(double tol, algorithm::search::ToleranceBy by, 
        const engine::glycan::GlycanStore subset):
            tolerance_(tol), by_(by), glycan_mass_(GlycanMassStore(subset)), 
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
        for(const auto& pep : candidate_.Peptides())
        {
            std::vector<model::spectrum::Peak> p1 = SearchOxonium(pep);
            if (p1.empty()) continue;
            for(const auto& composite: candidate_.Glycans(pep))
            {
                std::vector<int> positions = engine::protein::ProteinPTM::FindNGlycanSite(pep);
                for (const auto& pos : positions)
                {
                    std::vector<model::spectrum::Peak> p2 = SearchPeptides(pep, composite, pos);
                    if (p2.empty()) continue;
                }
            }
        }
        return res;   
    }



protected:
    virtual void SearchInit()
    {
        std::vector<std::shared_ptr<algorithm::search::Point<model::spectrum::Peak>>> points;
        for(const auto& it : spectrum_.Peaks())
        {
            std::shared_ptr<algorithm::search::Point<model::spectrum::Peak>> p = 
                std::make_shared<algorithm::search::Point<model::spectrum::Peak>>(it.MZ(), it);
            points.push_back(std::move(p));
        }
        searcher_.set_data(std::move(points));
        searcher_.Init();
    }

    virtual std::vector<model::spectrum::Peak> SearchOxonium(const std::string& seq)
    {
        double mass = util::mass::PeptideMass::Compute(seq);
        std::vector<model::spectrum::Peak> res;
        for (int i = 1; i <= 2; i++)
        {
            std::vector<model::spectrum::Peak> p =
                searcher_.Query(mass + util::mass::GlycanMass::kHexNAc * i);
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
        std::vector<double> peptides_mass = ComputePTMPeptideMass(seq, pos);
        binary_.set_data(peptides_mass);
        for(const auto& pk : spectrum_.Peaks())
        {
            if (binary_.Search(pk.MZ()))
            {
                res.push_back(pk);
            }
        }
        return res;
    }

    // for computing the peptide ions
    static std::vector<double> ComputePTMPeptideMass(const std::string& seq, const int pos)
    {
        std::vector<double> mass_list;
        for (int i = pos; i < seq.length() - 1; i++) // seldom at n
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

    static bool IntensityCmp(const model::spectrum::Peak& p1, const model::spectrum::Peak& p2)
        { return (p1.Intensity() < p2.Intensity()); }
    
    // virtual std::vector<model::spectrum::Peak> SearchGlycans(double target);

    double tolerance_;
    algorithm::search::ToleranceBy by_;
    GlycanMassStore glycan_mass_;
    algorithm::search::BucketSearch<model::spectrum::Peak> searcher_;
    algorithm::search::BinarySearch binary_;
    MatchResultStore candidate_;
    model::spectrum::Spectrum spectrum_;
}; 

} // namespace engine
} // namespace search

#endif