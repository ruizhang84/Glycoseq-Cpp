#ifndef ENGINE_SEARCH_PRECURSOR_MATCH_H
#define ENGINE_SEARCH_PRECURSOR_MATCH_H

#include <string>
#include <vector>
#include "../../algorithm/search/bucket_search.h"
#include "../../util/mass/peptide.h"
#include "../../model/glycan/glycan.h"
#include "../../util/mass/glycan.h"
#include "../../util/mass/spectrum.h"

namespace engine{
namespace search{

struct SearchResult
{
    std::string peptide;
    std::string glycan;
};

class PrecursorMatcher
{
public:
    PrecursorMatcher(double tol, algorithm::search::ToleranceBy by):
        tolerance_(tol), by_(by), 
            searcher_(algorithm::search::BucketSearch<std::string>(tol, by)){}

    void Init(const std::vector<std::string>& peptides, const std::vector<std::string>& glycans)
    {
        // set up glycans
        set_glycans(glycans);
        // set up search box
        set_peptides(peptides);
    }

    std::vector<std::string>& Glycans() { return glycans_; }
    void set_glycans(const std::vector<std::string>& glycans) { glycans_ = glycans; }

    void set_peptides(const std::vector<std::string>& peptides)
    {
        std::vector<std::shared_ptr<algorithm::search::Point<std::string>>> points;
        for(const auto& it : peptides)
        {
            double mass = util::mass::PeptideMass::Compute(it);
            std::shared_ptr<algorithm::search::Point<std::string>> p = 
                std::make_shared<algorithm::search::Point<std::string>>(mass, it);
            points.push_back(std::move(p));
        }
        searcher_.set_data(std::move(points));
        searcher_.Init();
    }

    double Tolerance() const { return tolerance_; }
    algorithm::search::ToleranceBy ToleranceType() const { return by_; }
    void set_tolerance(double tol) 
        { tolerance_ = tol; searcher_.set_tolerance(tol); searcher_.Init(); }
    void set_tolerance_by(algorithm::search::ToleranceBy by) 
        { by_ = by; searcher_.set_tolerance_by(by); searcher_.Init(); }

    std::vector<SearchResult> Match(const double mz, const int charge)
    {
        std::vector<SearchResult> res;
        double mass = util::mass::SpectrumMass::Compute(mz, charge);
        for(const auto& it : glycans_)
        {
            double delta = mass -
                 util::mass::GlycanMass::Compute(model::glycan::Glycan::Interpret(it));
            if (delta <= 0 ) continue;

            std::vector<std::string> peptides = searcher_.Query(delta);
            for(auto& p : peptides)
            {
                SearchResult r;
                r.glycan = it;
                r.peptide = p;
                res.push_back(r);
            }
        }
        return res;
    }

protected:
    double tolerance_;
    algorithm::search::ToleranceBy by_;
    algorithm::search::BucketSearch<std::string> searcher_;
    std::vector<std::string> glycans_;

}; 

} // namespace engine
} // namespace search

#endif