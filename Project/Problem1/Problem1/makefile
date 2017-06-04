BINTREE_LSTDFLG = -lstdc++ -lm
BINTREE_INCLUDE = -I/usr/include/
BINTREE_LIB = -L/usr/lib/
BINTREE_OBJS = bintree

all:	${BINTREE_OBJS}
	rm -f *.o

vectors.o: vectors.cpp vectors.h
	gcc -g -c vectors.cpp -o vectors.o

tree.o:	tree.cpp tree.h
	gcc -g -c tree.cpp -o tree.o

main.o: main.cpp vectors.h
	gcc -g -c main.cpp -o main.o ${BINTREE_INCLUDE} ${BINTREE_LIB}

bintree: main.o tree.o vectors.o
	gcc main.o tree.o vectors.o -o bintree ${BINTREE_LIB} ${BINTREE_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${BINTREE_OBJS}