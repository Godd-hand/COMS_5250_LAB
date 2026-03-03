CC = gcc
FFLAGS = -O3 -Wall -std=c99
OBJECTS = main.o matrix.o noise.o

.PHONY: clean

main.exe: $(OBJECTS)
	$(CC) $(FFLAGS) $(OBJECTS) -o main.exe -lm

%.o: %.c
	$(CC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) main.exe