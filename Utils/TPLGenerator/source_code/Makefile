

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all: util main

util:
	# ------ utils -------#
	$(CXX) util_src/PDB_File_Cut.cpp $(FLAG) -o PDB_File_Cut
	mv PDB_File_Cut ../util

main:
	$(CXX) TGT_To_TPL2.cpp $(FLAG) -o TGT_To_TPL2
	$(CXX) TGTpkl_To_TPL2.cpp $(FLAG) -o TGTpkl_To_TPL2
	mv TGT_To_TPL2 TGTpkl_To_TPL2 ../

# clean all the temp files
clean:
