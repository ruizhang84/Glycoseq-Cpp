#Created by Rui 5/17/20

CC = g++
CPPFLAGS =-g -Wall -std=c++14 
INCLUDES = -I/usr/local/include -L/usr/local/lib -lboost_unit_test_framework -static

TEST = glycan_test

glycan_test:
	$(CC) $(CPPFLAGS) $(INCLUDES) \
	-o glycan_test model/glycan/glycan_test.cpp  model/glycan/nglycan_complex.cpp 

# clean up
clean:
	rm -f core $(TEST)