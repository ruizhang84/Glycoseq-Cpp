#Created by Rui 5/17/20

CC = g++
CPPFLAGS =-g -Wall
GLYCAN = glycan_test nglycan_complex.o

glycan_test: nglycan_complex.o
	$(CC) $(CPPFLAGS) -o glycan_test model/glycan/glycan_test.cpp nglycan_complex.o
nglycan_complex.o:
	$(CC) $(CPPFLAGS) -c model/glycan/nglycan_complex.cpp 



# clean up
clean:
	rm -f core $(GLYCAN)