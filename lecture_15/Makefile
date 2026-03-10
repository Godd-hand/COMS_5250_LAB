CC = gcc
CFLAGS = -O3 -Wall -std=c99

OBJECTS = main.o matrix.o trimatrix.o Hessenberg.o FindZero.o GetShift.o QRdecomp.o QRA.o

.PHONY: clean

main.exe: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o main.exe -lm

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJECTS) main.exe