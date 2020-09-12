
# make the subdir first
CXX = g++
FLAG = -O3 
SHELL :=/bin/bash

all:
	#-> util
	$(CXX) MSA_To_SEQ.cpp $(FLAG) -o MSA_To_SEQ
	$(CXX) A3M_To_A2M.cpp $(FLAG) -o A3M_To_A2M
	$(CXX) A3M_Seq_Refine.cpp $(FLAG) -o A3M_Seq_Refine
	#-> move
	mkdir -p ../util
	mv MSA_To_SEQ A3M_To_A2M A3M_Seq_Refine ../util

# clean all the temp files
clean:

