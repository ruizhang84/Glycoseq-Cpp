#ifndef APP_SEARCH_SEARCH_PARAMETER_H
#define APP_SEARCH_SEARCH_PARAMETER_H

#include <deque>
#include <map>

#include "../../algorithm/search/search.h"
#include "../../engine/protein/protein_digest.h"
#include "../../engine/search/search_result.h"

struct SearchParameter
{

    // upper bound of glycan seaerch
    int n_thread = 6;
    int hexNAc_upper_bound = 12;
    int hex_upper_bound = 12;
    int fuc_upper_bound = 5;
    int neuAc_upper_bound = 4;
    int neuGc_upper_bound = 0;
    // searching precision
    double ms1_tol = 10;
    algorithm::search::ToleranceBy ms1_by =
        algorithm::search::ToleranceBy::PPM;
    double ms2_tol = 0.01;
    algorithm::search::ToleranceBy ms2_by = 
        algorithm::search::ToleranceBy::Dalton;
    // isotopic effects on precursor
    int isotopic_count = 2;
    // fdr
    double fdr_rate = 0.01;
    // protease
    std::deque<engine::protein::Proteases> proteases
    {
        engine::protein::Proteases::Trypsin,
        engine::protein::Proteases::GluC
    };
    int miss_cleavage = 2;
    // scoring weight
    std::map<engine::search::SearchType, double> weights 
    {
        {engine::search::SearchType::Core, 1.0}, 
        {engine::search::SearchType::Branch, 1.0}, 
        {engine::search::SearchType::Terminal, 1.0},
        {engine::search::SearchType::Peptide, 1.0}, 
        {engine::search::SearchType::Oxonium, 1.0},
        {engine::search::SearchType::Precursor, 1.0},
    };
    void set_search_weight(double w, engine::search::SearchType type)
        { weights[type] = w; }
};



#endif