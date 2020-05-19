#define BOOST_TEST_MODULE GlycanTest
#include <boost/test/unit_test.hpp>
#include <iostream>

#include "glycan.h"
#include "nglycan_complex.h"

using namespace std;
using namespace model::glycan;

void dfs(Monosaccharide& root, std::vector<string>& result);

BOOST_AUTO_TEST_CASE( glycan_test ) {
    // Glycan glycan;
    MonosaccharideFactory factory;
    std::shared_ptr<Monosaccharide> a = factory.GlcNAc(); 
    std::shared_ptr<Monosaccharide> b = factory.Man();
    std::shared_ptr<Monosaccharide> c = factory.Gal();

    // glycan.set_root(*a);
    // BOOST_CHECK( glycan.Root()->Name() == "GlcNAc"); 

    a->Child().push_back(b);
    b->set_parent(*a);
    b->Child().push_back(c);
    c->set_parent(*b);

    // std::vector<string> result;
    // dfs(*glycan.Root(), result);
    // Monosaccharide& curr = *glycan.Root();

    // BOOST_REQUIRE( !curr.Child().empty() );
    // curr = *curr.Child().front();
    // BOOST_CHECK( curr.Name() == "Man"); 

    // BOOST_REQUIRE( !curr.Child().empty() );
    // curr = curr.Child().front();
    // BOOST_CHECK( curr.Name() == "Gal"); 

    // BOOST_REQUIRE( curr.Parent() != NULL );
    // curr = curr.Parent();
    // BOOST_CHECK( curr.Name() == "Man"); 
    // NGlycanComplex test;

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








