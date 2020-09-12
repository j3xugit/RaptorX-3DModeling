
# make the subdir first
CC = gcc
CXX = g++
FLAG = -O3
SHELL :=/bin/bash

all:
	cd GCNNzsco_To_TBL_src&&make
	$(CXX) EstiDist_To_TBL.cpp $(FLAG) -o EstiDist_To_TBL
	$(CXX) SEQ_To_Input.cpp $(FLAG) -o SEQ_To_Input
	$(CXX) SSE_To_TBL.cpp $(FLAG) -o SSE_To_TBL
	$(CXX) DGSA_File_Mod.cpp $(FLAG) -o DGSA_File_Mod
	mv GCNNzsco_To_TBL EstiDist_To_TBL SEQ_To_Input SSE_To_TBL DGSA_File_Mod ../bin

# clean all the temp files
clean:
