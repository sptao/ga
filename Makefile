#Makefile for tsp solving by genetic algorithm

CC = gcc
DEBUG = -g
OPT = 
CFLAGS = $(OPT) -I./ $(DEBUG)
LFLAGS = -lm
TARGET = tsp

all = $(TARGET)

$(TARGET): tsp.o parser.o
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

tsp.o: tsp.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(LFLAGS)

parser.o: parser.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(LFLAGS)

clean:
	rm *.o tsp
