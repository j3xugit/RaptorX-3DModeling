

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:
	# make AcconPred
	$(CXX) profile.cpp seq.cpp template.cpp CNF_calc.cpp AcconPred.cpp $(FLAG) -o AcconPred
	mv AcconPred ../

# clean all the temp files
clean:
