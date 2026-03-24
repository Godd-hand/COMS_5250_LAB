CC = gcc
CFLAGS = -O3 -Wall -std=c99 -fopenmp

.PHONY: clean

main.exe: main.o
	$(CC) $(CFLAGS) main.o -o main.exe

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

clean:
	rm -f *.o main.exe *.dat *.png