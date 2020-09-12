

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:

	# ------ A3M to TGT -------#
	$(CXX) Verify_FASTA.cpp $(FLAG) -o Verify_FASTA
	$(CXX) A3M_To_PSI.cpp $(FLAG) -o A3M_To_PSI
	$(CXX) SS8_To_SS3.cpp $(FLAG) -o SS8_To_SS3
	$(CXX) ACC_ReFormat.cpp $(FLAG) -o ACC_ReFormat
	# move to ../ directory
	mv Verify_FASTA A3M_To_PSI SS8_To_SS3 ACC_ReFormat ../

# clean all the temp files
clean:
