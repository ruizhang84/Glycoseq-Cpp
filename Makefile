#Created by Rui 5/17/20

CC = g++
CPPFLAGS =-g -Wall -std=c++14 
INCLUDES = -I/usr/local/include -L/usr/local/lib -lboost_unit_test_framework -static -lpthread
LIB = -I/usr/local/include -L/usr/local/lib -lpthread

TEST_CASES := algorithm_base_test glycan_test mgf_parser_test lsh_test lsh_clustering_test
OBJS = lsh_clustering.o


lsh_clustering_test:
	$(CC) $(CPPFLAGS) -o test/lsh_clustering_test \
	 algorithm/clustering/lsh_clustering_test.cpp algorithm/clustering/lsh_clustering.cpp \
	 util/calc/lsh.cpp util/calc/calc.cpp util/io/mgf_parser.cpp $(INCLUDES)

clustering:
	$(CC) $(CPPFLAGS)  -o clustering \
	apps/clustering/clustering.cpp util/io/mgf_parser.cpp engine/spectrum/spectrum_binpacking.cpp \
	util/calc/lsh.cpp util/calc/calc.cpp algorithm/clustering/lsh_clustering.cpp $(LIB)


#  test
sim_test:
	$(CC) $(CPPFLAGS) -o test/sim_test \
	util/calc/spectrum_sim_test.cpp util/calc/calc.cpp util/io/mgf_parser.cpp $(INCLUDES)

algorithm_base_test:
	$(CC) $(CPPFLAGS) -o test/algorithm_base_test \
	algorithm/base/base_test.cpp $(INCLUDES)

lsh_test:
	$(CC) $(CPPFLAGS) -o test/lsh_test \
	util/calc/lsh_test.cpp util/calc/lsh.cpp util/calc/calc.cpp $(INCLUDES)

lsh_clustering_test:
	$(CC) $(CPPFLAGS) -o test/lsh_clustering_test \
	 algorithm/clustering/lsh_clustering_test.cpp \
	 algorithm/clustering/lsh_clustering.cpp util/calc/lsh.cpp util/calc/calc.cpp $(LIB)

mgf_parser_test:
	$(CC) $(CPPFLAGS) -o test/mgf_parser_test \
	util/io/mgf_parser_test.cpp util/io/mgf_parser.cpp  $(INCLUDES)

glycan_test:
	$(CC) $(CPPFLAGS) -o test/glycan_test \
	 model/glycan/glycan_test.cpp model/glycan/nglycan_complex.cpp $(INCLUDES)

# test
test: ${TEST_CASES}

# clean up
clean:
	rm -f core test/* *.o clustering