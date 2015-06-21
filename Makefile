#Makefile for tsp solving by genetic algorithm

CC = gcc
MPICC = mpicc
DEBUG = -g
OPT = 
CFLAGS = $(OPT) -I./ $(DEBUG)
LFLAGS = -lm

all: tsp ptsp

tsp: tsp.o parser.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

ptsp: ptsp.o parser.o
	$(MPICC) $(CFLAGS) $^ -o $@ $(LFLAGS)

tsp.o: tsp.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(LFLAGS)

ptsp.o: ptsp.c
	$(MPICC) $(CFLAGS) -c $^ -o $@ $(LFLAGS)

parser.o: parser.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(LFLAGS)

clean:
	rm *.o tsp
