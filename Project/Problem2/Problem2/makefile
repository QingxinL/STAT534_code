MPI_LSTDFLG = -lstdc++ -llapack -lcblas -lm
MPI_INCLUDE = -I/usr/include/
MPI_LIB = -L/usr/lib/atlas/
MPI_OBJS = parallelR2

all:	${MPI_OBJS}
	rm -f *.o

parallelR2.o: parallelR2.cpp
	mpic++ -g -c parallelR2.cpp -o parallelR2.o ${MPI_INCLUDE}

parallelR2: parallelR2.o
	mpic++ parallelR2.o -o parallelR2 ${MPI_LIB} ${MPI_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${MPI_OBJS}