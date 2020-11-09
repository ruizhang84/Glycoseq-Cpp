#define BOOST_TEST_MODULE ProteinTest
#include <boost/test/unit_test.hpp>
#include <vector>
#include <string>
#include <iostream>

#include "../../apps/search/search_parameter.h"
#include "protein_digest.h"
#include "protein_ptm.h"

#include <fstream>
#include <iostream>
#include "../../engine/protein/protein_digest.h"
#include "../../util/io/fasta_reader.h"

std::unordered_set<std::string> PeptidesDigestion
    (const std::string& fasta_path, SearchParameter parameter)
{
    util::io::FASTAReader fasta_reader(fasta_path);
    std::vector<model::protein::Protein> proteins = fasta_reader.Read();
   
    engine::protein::Digestion digest;
    digest.set_miss_cleavage(parameter.miss_cleavage);
    std::unordered_set<std::string> peptides;

    // digestion
    std::deque<engine::protein::Proteases> proteases(parameter.proteases);
    engine::protein::Proteases enzyme = proteases.front();
    digest.SetProtease(enzyme);
    proteases.pop_front();
    for(const auto& protein : proteins)
    {
        std::unordered_set<std::string> seqs = digest.Sequences(protein.Sequence(),
            engine::protein::ProteinPTM::ContainsNGlycanSite);
        peptides.insert(seqs.begin(), seqs.end());
    }
        
    // double digestion or more
    while (proteases.size() > 0)
    {
        std::unordered_set<std::string> double_seqs;
        enzyme = proteases.front();
        digest.SetProtease(enzyme);
        proteases.pop_front();
        for(const auto& seq : peptides)
        {
            std::unordered_set<std::string> seqs = digest.Sequences(seq,
                engine::protein::ProteinPTM::ContainsNGlycanSite);
            double_seqs.insert(seqs.begin(), seqs.end());
        }
        peptides.insert(double_seqs.begin(), double_seqs.end());
    }   

    return peptides;
}

BOOST_AUTO_TEST_CASE( Monosaccharide_test ) 
{

    // std::string seq = "MSALGAVIALLLWGQLFAVDSGNDSVTDIADDGCP"
    //                 "KPPEIAHGYVEHSVRYQCKNYYKLRTEGDGVYTLND";
    // engine::protein::Digestion digest;
    // std::unordered_set<std::string> seqs = digest.Sequences(seq, engine::protein::ProteinPTM::ContainsNGlycanSite);
    SearchParameter parameter;
    std::unordered_set<std::string> seqs = PeptidesDigestion("/home/ruiz/Documents/Glycoseq-Cpp/data/haptoglobin.fasta", parameter);
    for(auto& s: seqs)
    {
        std::cout << s << std::endl;
    }
}










