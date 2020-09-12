

# make the subdir first
CXX = g++
FLAG = -O3 -fopenmp 
SHELL :=/bin/bash

# compile directory
all:
	$(CXX) DeepCNF_Pred.cpp DeepCNF_Misc.cpp DeepCNF.cpp Chebyshev.cpp getopt.cpp $(FLAG) -o DeepCNF_Pred
	mv DeepCNF_Pred ../
	
# clean all the temp files
clean:
