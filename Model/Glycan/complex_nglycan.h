#ifndef COMPLEX_N_GLYCAN_H
#define COMPLEX_N_GLYCAN_H

#include "glycan.h"

class ComplexNGlycan : public Glycan
{
private:
    bool ValidAddGlcNAcCore();
    bool ValidAddGlcNAc();
    ComplexNGlycan CreateByAddGlcNAcCore();

    bool ValidAddGlcNAcBisect();
    ComplexNGlycan CreateByAddGlcNAcBisect();

    bool ValidAddGlcNAcBranch();
    std::vector<ComplexNGlycan> CreateByAddGlcNAcBranch();

    bool ValidAddMan();
    ComplexNGlycan CreateByAddMan();

    bool ValidAddGal();
    std::vector<ComplexNGlycan> CreateByAddGal();

    bool ValidAddFucCore();
    ComplexNGlycan CreateByAddFucCore();

    bool ValidAddFucTerminal();
    std::vector<ComplexNGlycan> CreateByAddFucTerminal();

    bool ValidAddNeuAc();
    std::vector<ComplexNGlycan> CreateByAddNeuAc();

    bool ValidAddNeuGc();
    std::vector<ComplexNGlycan> CreateByAddNeuGc();

public:
    ComplexNGlycan();
    ComplexNGlycan(const ComplexNGlycan&);
    ComplexNGlycan& operator= (const ComplexNGlycan&);
    std::string get_name() override;
    std::string get_id() override;
    std::vector<int> get_composition() override;
    std::vector<std::shared_ptr<Glycan>> Grow(Suger) override;
};


#endif