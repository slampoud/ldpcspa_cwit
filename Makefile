CC = gcc -Wall
LIBS = -lm
CFLAGS = -g $(INC)
BINARY=ldpcspa

simulation: elements.h elements.o simulation.h simulation.c SPA.h SPA.o voutput.o norm.o norm.h
	${CC} ${CFLAGS} -o ${BINARY} simulation.c elements.o SPA.o voutput.o norm.o $(LIBS)

elements.o: elements.c elements.h
	${CC} ${CFLAGS} -c elements.c 

voutput.o: voutput.c voutput.h
	${CC} ${CFLAGS} -c voutput.c 

SPA.o: SPA.h SPA.c
	${CC} ${CFLAGS} -c SPA.c 

norm.o: norm.h norm.c
	${CC} ${CFLAGS} -c norm.c 

flip: flip.c
	${CC} -o flip flip.c
clean:
	rm -rf *.o ${BINARY} *dSYM *~
