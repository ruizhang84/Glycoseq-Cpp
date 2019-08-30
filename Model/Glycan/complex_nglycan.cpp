#include "complex_nglycan.h"
#include "numeric"
#include <algorithm>
#include <sstream>
#include <iostream>

using namespace std;

ComplexNGlycan::ComplexNGlycan() : Glycan(){
    table.assign(24, 0);
}

ComplexNGlycan::ComplexNGlycan(const ComplexNGlycan& glycan ){
    table = glycan.table;
}

ComplexNGlycan& ComplexNGlycan::operator= (const ComplexNGlycan& glycan){
    table = glycan.table;
    return *this;
}

string ComplexNGlycan::get_name(){
    vector<int> composition = get_composition();
    string name = "ComplexNGlycan: ";
    if (table[2] > 0)
        name += "fucose ";
    if (table[3] > 0)
        name += "bisect ";
    name += "[" + to_string(composition[0]) + "," + to_string(composition[1]) + 
    "," + to_string(composition[2]) + "," + to_string(composition[3]) +"]";
    return name;
}

string ComplexNGlycan::get_id(){
    stringstream result;
    copy(table.begin(), table.end(), ostream_iterator<int>(result, " "));
    return result.str();
}

vector<int> ComplexNGlycan::get_composition(){
    vector<int> composition;
    composition.push_back(table[0] + table[3] + table[4] + table[5] + table[6] + table[7]);
    composition.push_back(table[1] + table[8] + table[9] + table[10] + table[11]);
    composition.push_back(table[2] + table[12] + table[13] + table[14] + table[15]);
    composition.push_back(table[16] + table[17] + table[18] + table[19]);
    composition.push_back(table[20] + table[21] + table[22] + table[23]);
    return composition;
}

vector<shared_ptr<Glycan>> ComplexNGlycan::Grow(Suger suger){
    vector<shared_ptr<Glycan>> glycans;
    switch (suger)
    {
    case GlcNAc:
        if (ValidAddGlcNAcCore()){
            ComplexNGlycan g = CreateByAddGlcNAcCore();
            shared_ptr<ComplexNGlycan> ptr = make_shared<ComplexNGlycan> (g);
            glycans.push_back(ptr);
        }else{
            if (ValidAddGlcNAcBisect()){
                ComplexNGlycan g = CreateByAddGlcNAcBisect();
                shared_ptr<ComplexNGlycan> ptr = make_shared<ComplexNGlycan> (g);
                glycans.push_back(ptr);
            }
            if (ValidAddGlcNAcBranch()){
                vector<ComplexNGlycan> gs = CreateByAddGlcNAcBranch();
                for (int i = 0; i < gs.size(); i ++){
                    shared_ptr<ComplexNGlycan> ptr = make_shared<ComplexNGlycan> (gs[i]);
                    glycans.push_back(ptr);
                }
            }
        }
        break;

    case Man:
        if (ValidAddMan()){
            ComplexNGlycan g = CreateByAddMan();
            shared_ptr<ComplexNGlycan> ptr = make_shared<ComplexNGlycan> (g);
            glycans.push_back(ptr);
        }
        break;

    case Gal:
        if (ValidAddGal()){
            vector<ComplexNGlycan> gs = CreateByAddGal();
            for (int i = 0; i < gs.size(); i ++){
                shared_ptr<ComplexNGlycan> ptr = make_shared<ComplexNGlycan> (gs[i]);
                glycans.push_back(ptr);
            }
        }
        break;

    case Fuc:
        if (ValidAddFucCore()){
            ComplexNGlycan g = CreateByAddFucCore();
            shared_ptr<ComplexNGlycan> ptr = make_shared<ComplexNGlycan> (g);
            glycans.push_back(ptr);
        }
        else if (ValidAddFucTerminal()){
            vector<ComplexNGlycan> gs = CreateByAddFucTerminal();
            for (int i = 0; i < gs.size(); i ++){
                shared_ptr<ComplexNGlycan> ptr = make_shared<ComplexNGlycan> (gs[i]);
                glycans.push_back(ptr);
            }
        }
        break;

    case NeuAc:
        if (ValidAddNeuAc()){
            vector<ComplexNGlycan> gs = CreateByAddNeuAc();
            for (int i = 0; i < gs.size(); i ++){
                shared_ptr<ComplexNGlycan> ptr = make_shared<ComplexNGlycan> (gs[i]);
                glycans.push_back(ptr);
            }
        }
        break;

    case NeuGc:
        if (ValidAddNeuGc()){
            vector<ComplexNGlycan> gs = CreateByAddNeuGc();
            for (int i = 0; i < gs.size(); i ++){
                shared_ptr<ComplexNGlycan> ptr = make_shared<ComplexNGlycan> (gs[i]);
                glycans.push_back(ptr);
            }
        }
        break;

    default:
        break;
    }
    return glycans;
}


bool ComplexNGlycan::ValidAddGlcNAcCore(){
    if (table[0] < 2)
        return true;
    return false;
}

ComplexNGlycan ComplexNGlycan::CreateByAddGlcNAcCore(){
    ComplexNGlycan g = *this;
    g.set_table(0, table[0]+1);
    return g;
}

bool ComplexNGlycan::ValidAddGlcNAcBisect(){
    if (table[1] == 3 && table[3] == 0 && table[4] == 0) //bisect 0, not extanding on GlcNAc
        return true;
    return false;

}

ComplexNGlycan ComplexNGlycan::CreateByAddGlcNAcBisect(){
    ComplexNGlycan g = *this;
    g.set_table(3, 1);
    return g;
}

bool ComplexNGlycan::ValidAddGlcNAcBranch(){
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 4] < table[i + 3]) // make it order
        {
            if (table[i + 4] == table[i + 8] && table[i + 12] == 0 && table[i + 16] == 0 && table[i + 20] == 0)
            //equal GlcNAc Gal, no Fucose attached at terminal, no terminal NeuAc, NeuGc
            {
                return true;
            }
        }
    }
    return false;

}

vector<ComplexNGlycan> ComplexNGlycan::CreateByAddGlcNAcBranch(){
    vector<ComplexNGlycan> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 4] < table[i + 3]) // make it order
        {
            if (table[i + 4] == table[i + 8] && table[i + 12] == 0 && table[i + 16] == 0 && table[i + 20] == 0)
            {
                ComplexNGlycan g = *this;
                g.set_table(i + 4, table[i + 4] + 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}

bool ComplexNGlycan::ValidAddMan(){
    if (table[0] == 2 && table[1] < 3)
        return true;
    return false;
}

ComplexNGlycan ComplexNGlycan::CreateByAddMan(){
    ComplexNGlycan g = *this;
    g.set_table(1, table[1] + 1);
    return g;
}

bool ComplexNGlycan::ValidAddGal(){
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 8] < table[i + 7]) // make it order
        {
            if (table[i + 4] == table[i + 8] + 1)
            {
                return true;
            }
        }
    }
    return false;
}
    
vector<ComplexNGlycan> ComplexNGlycan::CreateByAddGal(){
    vector<ComplexNGlycan> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 8] < table[i + 7]) // make it order
        {
            if (table[i + 4] == table[i + 8] + 1)
            {
                ComplexNGlycan g = *this;
                g.set_table(i + 8, table[i + 8] + 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}

bool ComplexNGlycan::ValidAddFucCore()
{
    if (table[0] == 1 && table[1] == 0 && table[2] == 0){   //core
        return true;
    }
    return false;
}

ComplexNGlycan ComplexNGlycan::CreateByAddFucCore(){
    ComplexNGlycan g = *this;
    g.set_table(2, 1);
    return g;
}

bool ComplexNGlycan::ValidAddFucTerminal()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 12] < table[i + 11]) // make it order
        {
            if (table[i + 12] == 0 && table[i + 4] > 0)
            {
                return true;
            }
        }
    }
    return false;
}
        
vector<ComplexNGlycan> ComplexNGlycan::CreateByAddFucTerminal()
{
    vector<ComplexNGlycan> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 12] < table[i + 11]) // make it order
        {
            if (table[i + 12] == 0 && table[i + 4] > 0)
            {
                ComplexNGlycan g = *this;
                g.set_table(i + 12, 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}

bool ComplexNGlycan::ValidAddNeuAc()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 16] < table[i + 15]) // make it order
        {
            if (table[i + 4] > 0 && table[i + 4] == table[i + 8] && table[i + 16] == 0 && table[i + 20] == 0)
            {
                return true;
            }
        }
    }
    return false;
}

vector<ComplexNGlycan> ComplexNGlycan::CreateByAddNeuAc()
{
    vector<ComplexNGlycan> glycans;
     for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 16] < table[i + 15]) // make it order
        {
            if (table[i + 4] > 0 && table[i + 4] == table[i + 8] && table[i + 16] == 0 && table[i + 20] == 0)
            {
                ComplexNGlycan g = *this;
                g.set_table(i + 16, 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}

bool ComplexNGlycan::ValidAddNeuGc()
{
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 20] < table[i + 19]) // make it order
        {
            if (table[i + 4] > 0 && table[i + 4] == table[i + 8] && table[i + 16] == 0 && table[i + 20] == 0)
            {
                return true;
            }
        }
    }
    return false;
}

vector<ComplexNGlycan> ComplexNGlycan::CreateByAddNeuGc()
{
    vector<ComplexNGlycan> glycans;
    for (int i = 0; i < 4; i++)
    {
        if (i == 0 || table[i + 20] < table[i + 19]) // make it order
        {
            if (table[i + 4] > 0 && table[i + 4] == table[i + 8] && table[i + 16] == 0 && table[i + 20] == 0)
            {
                ComplexNGlycan g = *this;
                g.set_table(i + 20, 1);
                glycans.push_back(g);
            }
        }
    }
    return glycans;
}


