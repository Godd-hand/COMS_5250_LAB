CC = mpicc
CFLAGS = -Wall -O3
LDFLAGS = -lm

all: serial_heat heat2d_mpi

# Serial version
serial_heat: serial_heat2d.c
	gcc $(CFLAGS) -o serial_heat serial_heat2d.c $(LDFLAGS)

# MPI version
heat2d_mpi: main.o grid.o
	$(CC) $(CFLAGS) -o heat2d_mpi main.o grid.o $(LDFLAGS)

main.o: main.c grid.h
	$(CC) $(CFLAGS) -c main.c

grid.o: grid.c grid.h
	$(CC) $(CFLAGS) -c grid.c

clean:
	rm -f *.o serial_heat heat2d_mpi output_*.dat *.png