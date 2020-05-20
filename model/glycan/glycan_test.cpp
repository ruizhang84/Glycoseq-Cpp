#define BOOST_TEST_MODULE GlycanTest
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "glycan.h"
#include "nglycan_complex.h"

using namespace std;
using namespace model::glycan;

std::unique_ptr<Moiety> CreateRoot()
{
    std::unique_ptr<Moiety> a = 
        std::make_unique<Moiety>(Monosaccharide::GlcNAc);
    std::unique_ptr<Moiety> b = 
        std::make_unique<Moiety>(Monosaccharide::Gal);
    std::unique_ptr<Moiety> c = 
        std::make_unique<Moiety>(Monosaccharide::Man);
    std::unique_ptr<Moiety> d = 
        std::make_unique<Moiety>(Monosaccharide::Fuc);

    c->Children().push_back(std::move(d));
    b->Children().push_back(std::move(c));
    a->Children().push_back(std::move(b));
    
    std::unique_ptr<Moiety> e = a->Clone();
    return e;
}

BOOST_AUTO_TEST_CASE( Monosaccharide_test ) 
{

    std::unique_ptr<Moiety> e = CreateRoot();
    BOOST_CHECK( e->Name() ==  Monosaccharide::GlcNAc); 
    BOOST_CHECK( e->Children().front()->Name() ==  Monosaccharide::Gal); 
    BOOST_CHECK( e->Children().front()->Children().front()
                    ->Children().front()->Name() ==  Monosaccharide::Fuc); 
    BOOST_CHECK( e->Children().front()->Children().front()
                    ->Parent()->Name() ==  Monosaccharide::Gal); 
}

BOOST_AUTO_TEST_CASE( glycan_test ) 
{
    Glycan nglycan;
    std::unique_ptr<Moiety> e = CreateRoot();

    nglycan.set_root(e);
    BOOST_CHECK( nglycan.Root()->Children().front()->Children().front()
                ->Parent()->Name() ==  Monosaccharide::Gal); 

    std::unique_ptr<Glycan> nglycan_ptr = nglycan.Clone();
    std::unique_ptr<Moiety> f = CreateRoot();
    nglycan_ptr->set_root(f);
    BOOST_CHECK( nglycan_ptr->Root()->Children().front()->Children().front()
                    ->Parent()->Name() ==  Monosaccharide::Gal); 
}

// BOOST_AUTO_TEST_CASE( nglycan_complex_test ) {
//     std::unique_ptr<NGlycanComplex> nglycan = std::make_unique<NGlycanComplex>();


// }

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








