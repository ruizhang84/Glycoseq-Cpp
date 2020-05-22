#Created by Rui 5/17/20

CC = g++
CPPFLAGS =-g -Wall -std=c++14 
INCLUDES = -I/usr/local/include -L/usr/local/lib -lboost_unit_test_framework -static

TEST_CASES := glycan_test mgf_parser_test lsh_test lsh_clustering_test

lsh_clustering_test:
	$(CC) $(CPPFLAGS) $(INCLUDES) \
	-o test/lsh_clustering_test algorithm/clustering/lsh_clustering_test.cpp \
	 algorithm/clustering/lsh_clustering.cpp util/calc/lsh.cpp util/calc/calc.cpp

lsh_test:
	$(CC) $(CPPFLAGS) $(INCLUDES) \
	-o test/lsh_test util/calc/lsh_test.cpp util/calc/lsh.cpp util/calc/calc.cpp

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