CXX = g++
CCFLAGS = -O3 -Wall -ffast-math -fopenmp
#CCFLAGS = -O3 -Wall -ffast-math -DNO_OPENMP
OBJ = MSA_To_PSSM.o getopt.o blast_mat.o blast_util.o blast_stat.o blast_pssm.o cdhit_util.o cdhit_scomat.o cdhit_option.o cdhit_param.o cdhit_seq.o cdhit_table.o cdhit_buffer.o cdhit_seqdb.o

make: $(OBJ)
	${CXX} $(CCFLAGS)  -o MSA_To_PSSM $(OBJ)
	rm -f $(OBJ)

MSA_To_PSSM.o  : 
	${CXX} ${CCFLAGS} -c MSA_To_PSSM.cpp
getopt.o : 
	${CXX} ${CCFLAGS} -c getopt.cpp
blast_mat.o :
	${CXX} ${CCFLAGS} -c blast_mat.cpp
blast_util.o :
	${CXX} ${CCFLAGS} -c blast_util.cpp
blast_stat.o :
	${CXX} ${CCFLAGS} -c blast_stat.cpp
blast_pssm.o :
	${CXX} ${CCFLAGS} -c blast_pssm.cpp
cdhit_util.o :
	${CXX} ${CCFLAGS} -c cdhit_util.cpp
cdhit_scomat.o :
	${CXX} ${CCFLAGS} -c cdhit_scomat.cpp
cdhit_option.o :
	${CXX} ${CCFLAGS} -c cdhit_option.cpp
cdhit_param.o :
	${CXX} ${CCFLAGS} -c cdhit_param.cpp
cdhit_seq.o :
	${CXX} ${CCFLAGS} -c cdhit_seq.cpp
cdhit_table.o :
	${CXX} ${CCFLAGS} -c cdhit_table.cpp
cdhit_buffer.o :
	${CXX} ${CCFLAGS} -c cdhit_buffer.cpp
cdhit_seqdb.o :
	${CXX} ${CCFLAGS} -c cdhit_seqdb.cpp

clean:
	rm -f MSA_To_PSSM $(OBJ)

