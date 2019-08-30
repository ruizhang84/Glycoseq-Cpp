#ifndef GLYCAN_H
#define GLYCAN_H

#include <string>
#include <vector>
#include <memory>

enum Suger { GlcNAc, Man, Gal, Fuc, NeuAc, NeuGc};

class Glycan
{
protected:
    std::vector<int> table;
   
public:
    Glycan() = default;
    
    virtual std::string get_name() = 0;
    virtual std::string get_id() = 0;
    std::vector<int> get_table();
    void set_table(std::vector<int>& table);
    void set_table(int index, int num);
    virtual std::vector<int> get_composition() = 0;
    virtual std::vector<std::shared_ptr<Glycan>> Grow(Suger) = 0;
};

#endif

