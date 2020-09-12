CXX = g++
CCFLAGS = -O3 -Wall -ffast-math -fopenmp
#CCFLAGS = -O3 -Wall -ffast-math -DNO_OPENMP
OBJ = meff_cdhit.o cdhit_util.o cdhit_scomat.o cdhit_option.o cdhit_param.o cdhit_seq.o cdhit_table.o cdhit_buffer.o cdhit_seqdb.o getopt.o

make: $(OBJ)
	${CXX} $(CCFLAGS)  -o meff_cdhit $(OBJ)
	rm -f $(OBJ)

meff_cdhit.o  : 
	${CXX} ${CCFLAGS} -c meff_cdhit.cpp
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
getopt.o :
	${CXX} ${CCFLAGS} -c getopt.cpp

clean:
	rm -f meff_cdhit $(OBJ)

