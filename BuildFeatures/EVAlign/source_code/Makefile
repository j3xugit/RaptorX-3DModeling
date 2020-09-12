

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:
	#------ util executable -----#
	$(CXX) Gen_EVfold_config.cpp $(FLAG) -o Gen_EVfold_config
	$(CXX) Gen_EVfold_config2.cpp $(FLAG) -o Gen_EVfold_config2
	$(CXX) Gen_EVfold_config3.cpp $(FLAG) -o Gen_EVfold_config3
	$(CXX) EVfold_MSA_Trans.cpp $(FLAG) -o EVfold_MSA_Trans
	$(CXX) Verify_FASTA.cpp $(FLAG) -o Verify_FASTA
	# move to ../util directory
	mv Gen_EVfold_config Gen_EVfold_config2 Gen_EVfold_config3 EVfold_MSA_Trans Verify_FASTA ../util

# clean all the temp files
clean:
