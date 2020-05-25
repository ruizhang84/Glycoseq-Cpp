#ifndef ENGINE_GLYCAN_GLYCAN_BUILDER_H
#define ENGINE_GLYCAN_GLYCAN_BUILDER_H

#include <string>
#include <vector>
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

class GlycanBuilder
{
public:
    GlycanBuilder(int hexNAc, int hex, int fuc, int neuAc, int neuGc):
        hexNAc_(hexNAc), hex_(hex), fuc_(fuc), neuAC_(neuAc), neuGc_(neuGc){};

    
    void Build()
    {
        model::glycan::NGlycanComplex root;
    }

protected:
    int hexNAc_;
    int hex_;
    int fuc_;
    int neuAC_;
    int neuGc_;
};


} // namespace engine
} // namespace glycan




#endif