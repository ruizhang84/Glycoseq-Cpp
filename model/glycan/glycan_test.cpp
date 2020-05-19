#define BOOST_TEST_MODULE GlycanTest
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "glycan.h"
#include "nglycan_complex.h"

using namespace std;
using namespace model::glycan;

class GlycanTemp : public Glycan
{
public:
    std::unique_ptr<Glycan> Clone() override
    {
        return std::make_unique<GlycanTemp>(*this);
    }

    std::vector<std::unique_ptr<Glycan>> Add(Monosaccharide suger) override
    { 
        std::vector<std::unique_ptr<Glycan>> result; 
        std::unique_ptr<Glycan> new_glycan = Clone();
        new_glycan->Terminal().push_back(suger);
        result.push_back(std::move(new_glycan));
        return result; 
    }
};

BOOST_AUTO_TEST_CASE( glycan_test ) {
    std::unique_ptr<GlycanTemp> glycan = std::make_unique<GlycanTemp>();
    std::vector<std::unique_ptr<Glycan>> vect = glycan->Add(Monosaccharide::GlcNAc);
    BOOST_CHECK( vect.front()->Terminal().front() == Monosaccharide::GlcNAc);   

    GlycanTemp* new_glycan = (GlycanTemp*) vect.front().get();
    std::vector<std::unique_ptr<Glycan>> vect2 = new_glycan->Add(Monosaccharide::Fuc);
    BOOST_CHECK( vect2.front()->Terminal().front() == Monosaccharide::GlcNAc);   
    BOOST_CHECK( vect2.front()->Terminal().back() == Monosaccharide::Fuc);   

}


// BOOST_AUTO_TEST_CASE( glycan_add_test )
// {
//     Glycan glycan;
//     GlcNAc a; Man b; Gal c; Fuc d; NeuAc e;

//     glycan.set_root(&a);
//     glycan.Terminal().push_back(glycan.Root());

//     Monosaccharide* curr = glycan.Terminal().front();
//     BOOST_CHECK( curr->Name() == "GlcNAc"); 

//     // Glycan* new_glycan = glycan.Add(&b).front();
//     // curr = new_glycan->Terminal().front();
//     // BOOST_CHECK( curr->Name() == "Man"); 
//     Glycan new_glycan = glycan;
//     curr = new_glycan.Terminal().front();
//     BOOST_CHECK( curr->Name() == "GlcNAc"); 

//     // BOOST_CHECK_MESSAGE( curr->Name() == "Man", 
//     // "The terminal is " << new_glycan->Terminal().front()->Name()); 

// }

// BOOST_AUTO_TEST_CASE( nglycan_complex_test ) {
//     NGlycanComplex nglycan;
//     GlcNAc a; Man b; Gal c; Fuc d; NeuAc e;
//     nglycan.set_root(&a);

//     BOOST_CHECK( nglycan.Root()->Name() == "GlcNAc"); 
//     BOOST_CHECK( nglycan.Root()->Type() == SugerType::HexNAc); 

// }

// void dfs(Monosaccharide* root, std::vector<string>& result)
// {
//     result.push_back(root->Name());
//     if (root->Child().size() > 0){
//         for (auto child : root.Child()){
//             dfs(child, result);
//         }
//     }
// }

// int add( int i, int j ) { return i+j; }

// BOOST_AUTO_TEST_CASE( my_test )
// {
//     // seven ways to detect and report the same error:
//     BOOST_CHECK( add( 2,2 ) == 4 );        // #1 continues on error

//     BOOST_REQUIRE( add( 2,2 ) == 4 );      // #2 throws on error

//     if( add( 2,2 ) != 4 )
//       BOOST_ERROR( "Ouch..." );            // #3 continues on error

//     if( add( 2,2 ) != 4 )
//       BOOST_FAIL( "Ouch..." );             // #4 throws on error

//     if( add( 2,2 ) != 4 ) throw "Ouch..."; // #5 throws on error

//     BOOST_CHECK_MESSAGE( add( 2,2 ) == 4,  // #6 continues on error
//                          "add(..) result: " << add( 2,2 ) );

//     BOOST_CHECK_EQUAL( add( 2,2 ), 4 );	  // #7 continues on error
// }








