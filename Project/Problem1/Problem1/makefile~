VECTORS_LSTDFLG = -lstdc++ -lm
VECTORS_INCLUDE = -I/usr/include/
VECTORS_LIB = -L/usr/lib/
VECTORS_OBJS = vectors

all:	${VECTORS_OBJS}
	rm -f *.o

vectors.o: vectors.cpp vectors.h
	gcc -g -c vectors.cpp -o vectors.o

main.o: main.cpp vectors.h
	gcc -g -c main.cpp -o main.o ${VECTORS_INCLUDE} ${VECTORS_LIB}

vectors: main.o vectors.o
	gcc main.o vectors.o -o vectors ${VECTORS_LIB} ${VECTORS_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${VECTORS_OBJS}