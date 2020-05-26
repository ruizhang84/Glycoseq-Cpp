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

class MatchResultStore
{
public:
    void Add(const std::string& peptide, const std::string& glycan)
    {
        if (map_.find(peptide) == map_.end())
        {
            peptides_.push_back(peptide);
            map_[peptide] = std::vector<std::string>();
        }
        map_[peptide].push_back(glycan);
    }
    std::unordered_map<std::string, 
        std::vector<std::string>> Map() { return map_; }
    bool Empty() { return peptides_.size() == 0; }
    std::vector<std::string> Peptides() { return peptides_; }
    std::vector<std::string> Glycans(const std::string& peptide)
    {
        if (map_.find(peptide) != map_.end())
        {
            return map_[peptide];
        }
        return std::vector<std::string>();
    }

protected:
    std::vector<std::string> peptides_;
    std::unordered_map<std::string, 
        std::vector<std::string>> map_;
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
    std::vector<std::string>& Peptides() { return peptides_; }
    virtual void set_glycans(const std::vector<std::string>& glycans) { glycans_ = glycans; }
    virtual void set_peptides(const std::vector<std::string>& peptides)
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

    virtual MatchResultStore Match(const double target)
    {
        MatchResultStore res;
        for(const auto& it : glycans_)
        {
            double delta = target -
                 util::mass::GlycanMass::Compute(model::glycan::Glycan::Interpret(it));
            if (delta <= 0 ) continue;

            std::vector<std::string> peptides = searcher_.Query(delta);
            for(auto& p : peptides)
            {
                res.Add(p, it);
            }
        }
        return res;
    }

protected:
    double tolerance_;
    algorithm::search::ToleranceBy by_;
    algorithm::search::BucketSearch<std::string> searcher_;
    std::vector<std::string> glycans_;
    std::vector<std::string> peptides_;

}; 

} // namespace engine
} // namespace search

#endif