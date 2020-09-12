CXX = g++
CCFLAGS = -O3 -Wall -ffast-math -fopenmp
#CCFLAGS = -O3 -Wall -ffast-math -DNO_OPENMP
OBJ = self_filter.o cdhit_util.o cdhit_scomat.o cdhit_option.o cdhit_param.o cdhit_seq.o cdhit_table.o cdhit_buffer.o cdhit_seqdb.o

make: $(OBJ)
	${CXX} $(CCFLAGS)  -o self_filter $(OBJ)
	rm -f $(OBJ)

self_filter.o  : 
	${CXX} ${CCFLAGS} -c self_filter.cpp
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
	rm -f self_filter $(OBJ)

