#ifndef ENGINE_GLYCAN_GLYCAN_BUILDER_H
#define ENGINE_GLYCAN_GLYCAN_BUILDER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <deque>
#include <memory>
#include "../../model/glycan/nglycan_complex.h"

namespace engine{
namespace glycan {

struct GlycanStore
{
    std::string composite_str;
    std::vector<std::string> tables_str;
};

struct GlycanMassStore
{
    std::string table_str;
    std::vector<double> subset_mass;
};

using namespace model::glycan;

class GlycanBuilder
{

public:
    GlycanBuilder(int hexNAc, int hex, int fuc, int neuAc, int neuGc):
        hexNAc_(hexNAc), hex_(hex), fuc_(fuc), neuAc_(neuAc), neuGc_(neuGc),
            candidates_({Monosaccharide::GlcNAc, Monosaccharide::Man, Monosaccharide::Gal,
                Monosaccharide::Fuc, Monosaccharide::NeuAc}){}

    
    void Build()
    {
        std::unique_ptr<NGlycanComplex> root = 
            std::make_unique<NGlycanComplex>();
    
        std::deque<std::unique_ptr<Glycan>> queue;
        std::unordered_map<std::string, NGlycanComplex*> visited;

        queue.push_back(std::move(root));
        while (!queue.empty())
        {
            std::unique_ptr<Glycan> node = std::move(queue.front());
            queue.pop_front();
            for(const auto& it : candidates_)
            {
                std::vector<std::unique_ptr<Glycan>> res = node->Grow(it);
                for(auto& g : res)
                {
                    queue.push_back(std::move(g));
                }
            }
        }

    }

    std::vector<Monosaccharide> Candidates() { return candidates_; }
    int HexNAc() { return hexNAc_; }
    int Hex() { return hex_; }
    int Fuc() { return fuc_; }
    int NeuAc() { return neuAc_; }
    int NeuGc() { return neuGc_; }
    void set_candidates(std::vector<Monosaccharide> sugers) { candidates_ = sugers; }
    void set_HexNAc(int num) { hexNAc_ = num; }
    void set_Hex(int num) { hex_ = num; }
    void set_Fuc(int num) { fuc_ = num; }
    void set_NeuAc(int num) { neuAc_ = num; }
    void set_NeuGc(int num) { neuGc_ = num; }

protected:
    int hexNAc_;
    int hex_;
    int fuc_;
    int neuAc_;
    int neuGc_;
    std::vector<Monosaccharide> candidates_;

};


} // namespace engine
} // namespace glycan




#endif