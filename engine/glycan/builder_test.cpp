#define BOOST_TEST_MODULE BUilderTest
#include <boost/test/unit_test.hpp>

#include "glycan_builder.h"


namespace engine{
namespace glycan {

BOOST_AUTO_TEST_CASE( glycan_builder_test ) 
{
    GlycanBuilder builder(2, 3, 1, 0, 0);
    builder.Build();

}

}
}