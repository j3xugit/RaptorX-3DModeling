
# make the subdir first
CXX = g++
FLAG = -O3 
PARAL = -fopenmp
SHELL :=/bin/bash

all: tgt util sse acc


tgt:
	# ---  MSA_To_PSSM ---#
	cd MSA_To_PSSM_cdhit_src&&make
	mv MSA_To_PSSM_cdhit_src/MSA_To_PSSM ../util
	# --- printTGT ---- #
	cd printTGT_src&&make
	mv printTGT_src/printTGT ../bin
	# ---  A3M_To_TGT ----#
	$(CXX) A3M_To_TGT.cpp $(FLAG) -o A3M_To_TGT
	mv A3M_To_TGT ../util

util:
	# --- tgt related ---#
	$(CXX) util_src/ACC_ReFormat.cpp $(FLAG) -o ACC_ReFormat
	$(CXX) util_src/SS8_To_SS3.cpp $(FLAG) -o SS8_To_SS3
	$(CXX) util_src/Verify_FASTA.cpp $(FLAG) -o Verify_FASTA
	mv ACC_ReFormat SS8_To_SS3 Verify_FASTA ../util
	# --- contact related --#
	$(CC) util_src/alnstats.c $(FLAG) -lm -o alnstats
	$(CC) util_src/cov21stats.c $(FLAG) -o alnstats
	$(CXX) util_src/alnstats_omp.cpp $(FLAG) $(PARAL) -lm -o alnstats_omp
	$(CXX) util_src/A3M_To_PSI.cpp $(FLAG) -o A3M_To_PSI
	mv cov21stats ../bin
	mv alnstats ../bin
	mv alnstats_omp ../bin
	mv A3M_To_PSI ../util

sse:
	# ---- DeepCNF_SS_Con ---- #
	cd DeepCNF_SS_Con_src&&make
	mv DeepCNF_SS_Con ../bin

acc:
	# ---- AcconPred ----- #
	cd AcconPred_src&&make
	mv AcconPred ../bin


# clean all the temp files
clean:
