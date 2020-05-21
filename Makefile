#Created by Rui 5/17/20

CC = g++
CPPFLAGS =-g -Wall -std=c++14 
INCLUDES = -I/usr/local/include -L/usr/local/lib -lboost_unit_test_framework -static

TEST_CASES := glycan_test mgf_parser_test

mgf_parser_test:
	$(CC) $(CPPFLAGS) $(INCLUDES) \
	-o test/mgf_parser_test util/io/mgf_parser_test.cpp util/io/mgf_parser.cpp

glycan_test:
	$(CC) $(CPPFLAGS) $(INCLUDES) \
	-o test/glycan_test model/glycan/glycan_test.cpp model/glycan/nglycan_complex.cpp 

# test
test: ${TEST_CASES}

# clean up
clean:
	rm -f core test/*