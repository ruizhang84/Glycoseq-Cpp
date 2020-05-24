#include "nglycan_complex.h"

namespace model {
namespace glycan {


std::vector<std::unique_ptr<Glycan>> NGlycanComplex::Grow(Monosaccharide suger){
   std::vector<std::unique_ptr<Glycan>>  glycans;
    switch (suger)
    {
    case Monosaccharide::GlcNAc:
        if (ValidAddGlcNAcCore()){
            NGlycanComplex g = CreateByAddGlcNAcCore();
            std::unique_ptr<NGlycanComplex> ptr = std::make_unique<NGlycanComplex> (g);
            glycans.push_back(std::move(ptr));
        }else if (ValidAddGlcNAc()){
            if (ValidAddGlcNAcBisect()){
                NGlycanComplex g = CreateByAddGlcNAcBisect();
                std::unique_ptr<NGlycanComplex> ptr = std::make_unique<NGlycanComplex> (g);
                glycans.push_back(std::move(ptr));
            }
            if (ValidAddGlcNAcBranch()){
                std::vector<NGlycanComplex> gs = CreateByAddGlcNAcBranch();
                for (int i = 0; i < gs.size(); i ++){
                    std::unique_ptr<NGlycanComplex> ptr = std::make_unique<NGlycanComplex> (gs[i]);
                    glycans.push_back(std::move(ptr));
                }
            }
        }
        break;

    case Monosaccharide::Man:
        if (ValidAddMan()){
            NGlycanComplex g = CreateByAddMan();
            std::unique_ptr<NGlycanComplex> ptr = std::make_unique<NGlycanComplex> (g);
            glycans.push_back(std::move(ptr));
        }
        break;

    case Monosaccharide::Gal:
        if (ValidAddGal()){
            std::vector<NGlycanComplex> gs = CreateByAddGal();
            for (int i = 0; i < gs.size(); i ++){
                std::unique_ptr<NGlycanComplex> ptr = std::make_unique<NGlycanComplex> (gs[i]);
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::Fuc:
        if (ValidAddFucCore()){
            NGlycanComplex g = CreateByAddFucCore();
            std::unique_ptr<NGlycanComplex> ptr = std::make_unique<NGlycanComplex> (g);
            glycans.push_back(std::move(ptr));
        }
        else if (ValidAddFucTerminal()){
            std::vector<NGlycanComplex> gs = CreateByAddFucTerminal();
            for (int i = 0; i < gs.size(); i ++){
                std::unique_ptr<NGlycanComplex> ptr = std::make_unique<NGlycanComplex> (gs[i]);
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::NeuAc:
        if (ValidAddNeuAc()){
            std::vector<NGlycanComplex> gs = CreateByAddNeuAc();
            for (int i = 0; i < gs.size(); i ++){
                std::unique_ptr<NGlycanComplex> ptr = std::make_unique<NGlycanComplex> (gs[i]);
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    case Monosaccharide::NeuGc:
        if (ValidAddNeuGc()){
            std::vector<NGlycanComplex> gs = CreateByAddNeuGc();
            for (int i = 0; i < gs.size(); i ++){
                std::unique_ptr<NGlycanComplex> ptr = std::make_unique<NGlycanComplex> (gs[i]);
                glycans.push_back(std::move(ptr));
            }
        }
        break;

    default:
        break;
    }
    return glycans;
}


bool NGlycanComplex::ValidAddGlcNAcCore(){
    if (table_[0] < 2)
        return true;
    return false;
}

bool NGlycanComplex::ValidAddGlcNAc(){
    if (table_[0] == 2 && table_[1] == 3)
        return true;
    return false;
}

NGlycanComplex NGlycanComplex::CreateByAddGlcNAcCore(){
    NGlycanComplex g = *this;
    g.set_table(0, table_[0]+1);
    return g;
}

bool NGlycanComplex::ValidAddGlcNAcBisect(){
    if (table_[1] == 3 && table_[3] == 0 && table_[4] == 0) //bisect 0, not extanding on GlcNAc
        return true;
    return false;

}

NGlycanComplex NGlycanComplex::CreateByAddGlcNAcBisect(){
    NGlycanComplex g = *this;
    g.set_table(3, 1);
    return g;
}

bool NGlycanComplex::ValidAddGlcNAcBranch(){
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 4] < table_[i + 3]) // make it order
        {
            if (table_[i + 4] == table_[i + 8] && table_[i + 12] == 0 && table_[i + 16] == 0 && table_[i + 20] == 0)
            //equal GlcNAc Gal, no Fucose attached at terminal, no terminal NeuAc, NeuGc
            {
                return true;
            }
        }
    }
    return false;

}

std::vector<NGlycanComplex> NGlycanComplex::CreateByAddGlcNAcBranch(){
    std::vector<NGlycanComplex> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 4] < table_[i + 3]) // make it order
        {
            if (table_[i + 4] == table_[i + 8] && table_[i + 12] == 0 && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                NGlycanComplex g = *this;
                g.set_table(i + 4, table_[i + 4] + 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}

bool NGlycanComplex::ValidAddMan(){
    if (table_[0] == 2 && table_[1] < 3)
        return true;
    return false;
}

NGlycanComplex NGlycanComplex::CreateByAddMan(){
    NGlycanComplex g = *this;
    g.set_table(1, table_[1] + 1);
    return g;
}

bool NGlycanComplex::ValidAddGal(){
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
        {
            if (table_[i + 4] == table_[i + 8] + 1)
            {
                return true;
            }
        }
    }
    return false;
}
    
std::vector<NGlycanComplex> NGlycanComplex::CreateByAddGal(){
    std::vector<NGlycanComplex> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 8] < table_[i + 7]) // make it order
        {
            if (table_[i + 4] == table_[i + 8] + 1)
            {
                NGlycanComplex g = *this;
                g.set_table(i + 8, table_[i + 8] + 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}

bool NGlycanComplex::ValidAddFucCore()
{
    if (table_[0] == 1 && table_[1] == 0 && table_[2] == 0){   //core
        return true;
    }
    return false;
}

NGlycanComplex NGlycanComplex::CreateByAddFucCore(){
    NGlycanComplex g = *this;
    g.set_table(2, 1);
    return g;
}

bool NGlycanComplex::ValidAddFucTerminal()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 12] < table_[i + 11]) // make it order
        {
            if (table_[i + 12] == 0 && table_[i + 4] > 0)
            {
                return true;
            }
        }
    }
    return false;
}
        
std::vector<NGlycanComplex> NGlycanComplex::CreateByAddFucTerminal()
{
    std::vector<NGlycanComplex> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 12] < table_[i + 11]) // make it order
        {
            if (table_[i + 12] == 0 && table_[i + 4] > 0)
            {
                NGlycanComplex g = *this;
                g.set_table(i + 12, 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}

bool NGlycanComplex::ValidAddNeuAc()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 16] < table_[i + 15]) // make it order
        {
            if (table_[i + 4] > 0 && table_[i + 4] == table_[i + 8] && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                return true;
            }
        }
    }
    return false;
}

std::vector<NGlycanComplex> NGlycanComplex::CreateByAddNeuAc()
{
    std::vector<NGlycanComplex> glycans;
     for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 16] < table_[i + 15]) // make it order
        {
            if (table_[i + 4] > 0 && table_[i + 4] == table_[i + 8] && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                NGlycanComplex g = *this;
                g.set_table(i + 16, 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}

bool NGlycanComplex::ValidAddNeuGc()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 20] < table_[i + 19]) // make it order
        {
            if (table_[i + 4] > 0 && table_[i + 4] == table_[i + 8] && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                return true;
            }
        }
    }
    return false;
}

std::vector<NGlycanComplex> NGlycanComplex::CreateByAddNeuGc()
{
    std::vector<NGlycanComplex> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table_[i + 20] < table_[i + 19]) // make it order
        {
            if (table_[i + 4] > 0 && table_[i + 4] == table_[i + 8] && table_[i + 16] == 0 && table_[i + 20] == 0)
            {
                NGlycanComplex g = *this;
                g.set_table(i + 20, 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}


}  //  namespace glycan
}  //  namespace model
