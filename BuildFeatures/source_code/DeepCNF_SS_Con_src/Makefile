
# make the subdir first
CXX = g++
#DEFINE = -D_MPI
FLAG = -O3 
#FLAG = -O3
SHELL :=/bin/bash

# compile directory
all:
	$(CXX) DeepCNF_SS_Con.cpp DeepCNF.cpp Chebyshev.cpp DeepCNF_Misc.cpp getopt.cpp profile.cpp seq.cpp $(FLAG) $(DEFINE) -o DeepCNF_SS_Con
	mv DeepCNF_SS_Con ../
	
# clean all the temp files
clean:
