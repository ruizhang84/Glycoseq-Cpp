#ifndef ENGINE_GLYCAN_GLYCAN_BUILDER_H
#define ENGINE_GLYCAN_GLYCAN_BUILDER_H

#include <deque>
#include <memory>
#include "glycan_store.h"
#include "../../model/glycan/nglycan_complex.h"
#include "../../util/mass/glycan.h"

namespace engine{
namespace glycan {

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
        queue.push_back(std::move(root));

        while (!queue.empty())
        {
            std::unique_ptr<Glycan> node = std::move(queue.front());
            isomer_store_.Add(node->Name(), node->ID());
            isomer_store_.Add(node->Name(), 
                util::mass::GlycanMass::Compute(node->Composition()));

            queue.pop_front();
            for(const auto& it : candidates_)
            {
                std::vector<std::unique_ptr<Glycan>> res = node->Grow(it);
                for(auto& g : res)
                {
                    if (SatisfyCriteria(g.get()))
                    {
                        std::string id = g->ID();
                        if (!mass_store_.Contains(id))
                        {
                            queue.push_back(std::move(g));
                        }
                        mass_store_.AddSubset(id, node->ID(), 
                            util::mass::GlycanMass::Compute(node->Composition()));
                    }
                }
            }
        }
    }

    GlycanStore Isomer() { return isomer_store_; }
    GlycanMassStore Mass() { return mass_store_; }
    void Clear() 
    {
        isomer_store_.Clear();
        mass_store_.Clear();
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
    bool SatisfyCriteria(const Glycan* glycan) const
    {
        int hexNAc = 0, hex = 0, fuc = 0, neuAc = 0, neuGc = 0;
        for(auto& it : glycan->CompositionConst())
        {
            switch (it.first)
            {
            case Monosaccharide::GlcNAc:
                hexNAc += it.second;
                break;
            case Monosaccharide::Gal:
                hex += it.second;
                break;
            case Monosaccharide::Man:
                hex += it.second;
                break;    
            case Monosaccharide::Fuc:
                fuc += it.second;
                break;   
            case Monosaccharide::NeuAc:
                neuAc += it.second;
                break;   
            case Monosaccharide::NeuGc:
                neuGc += it.second;
                break;           
            default:
                break;
            }
        }
        return (hexNAc <= hexNAc_ && hex <= hex_ && fuc <= fuc_
                && neuAc <= neuAc_ && neuGc <= neuGc_);
    }

    int hexNAc_;
    int hex_;
    int fuc_;
    int neuAc_;
    int neuGc_;
    std::vector<Monosaccharide> candidates_;
    GlycanStore isomer_store_;
    GlycanMassStore mass_store_;

};


} // namespace engine
} // namespace glycan




#endif