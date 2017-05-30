LSTDFLG = -lstdc++ -lm -lgsl -lgslcblas
INCLUDE = -I/usr/include/
LIB = -L/usr/lib/
OBJS = randmvn
CFLAGS = -ffast-math -ftree-vectorize -ftree-loop-linear -mfpmath=sse

all:	${OBJS}
	rm -f *.o

main.o: main.cpp
	gcc -m64 -g -c -O3 $(CFLAGS) main.cpp -o main.o ${INCLUDE} ${LIB}

randmvn: main.o
	gcc -m64 main.o -o randmvn ${LIB} ${LSTDFLG}

clean:
	rm -f *.o
	rm -f ${OBJS}