REGRESSIONS_LSTDFLG = -lstdc++ -lgsl -lgslcblas -lm
REGRESSIONS_INCLUDE = -I/usr/include/
REGRESSIONS_LIB = -L/usr/lib/
REGRESSIONS_OBJS = regressions

all:	${REGRESSIONS_OBJS}
	rm -f *.o

regmodels.o: regmodels.cpp regmodels.h
	gcc -g -c regmodels.cpp -o regmodels.o ${REGRESSIONS_INCLUDE}

main.o: main.cpp
	gcc -g -c main.cpp -o main.o

regressions: main.o regmodels.o
	gcc main.o regmodels.o -o regressions ${REGRESSIONS_LIB} ${REGRESSIONS_LSTDFLG}

clean:
	rm -f *.o