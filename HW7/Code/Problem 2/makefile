LSTDFLG = -lstdc++ -lm -lgsl -lgslcblas
INCLUDE = -I/usr/include/
LIB = -L/usr/lib/
OBJS = randbin

all:	${OBJS}
	rm -f *.o

vectors.o: vectors.cpp vectors.h
	gcc -g -c vectors.cpp -o vectors.o

main.o: main.cpp vectors.cpp
	mpic++ -g -c main.cpp -o main.o ${INCLUDE} ${LIB}

randbin: main.o vectors.o
	mpic++ main.o vectors.o -o randbin ${LIB} ${LSTDFLG}

clean:
	rm -f *.o
	rm -f ${OBJS}