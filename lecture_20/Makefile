CC = gcc
CFLAGS = -O3 -Wall -std=c99 -fopenmp

OBJECTS = main.o adaptive_quad.o

.PHONY: clean

main.exe: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o main.exe -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o main.exe *.dat