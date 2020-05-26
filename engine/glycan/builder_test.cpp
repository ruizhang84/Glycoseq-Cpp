#define BOOST_TEST_MODULE BuilderTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include "glycan_builder.h"


namespace engine{
namespace glycan {

BOOST_AUTO_TEST_CASE( glycan_builder_test ) 
{
    GlycanBuilder builder(2, 3, 1, 0, 0);
    builder.Build();
    std::cout <<  builder.Isomer().Map().size() << std::endl;
    for(auto& it :  builder.Isomer().Map()){
        std::cout << it.first << std::endl;
        for (auto& j: it.second)
        {
            std::cout << j << std::endl;
        }
        std::cout << std::endl;
    }
    BOOST_CHECK(builder.Isomer().Map().size() > 10);

}

}
}