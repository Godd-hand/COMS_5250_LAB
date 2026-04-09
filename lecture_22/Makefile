CC_OMP = gcc
CC_MPI = mpicc
CFLAGS = -O3 -Wall -std=c99

all: pi_openmp.exe pi_mpi.exe

pi_openmp.exe: pi_openmp.c
	$(CC_OMP) $(CFLAGS) -fopenmp pi_openmp.c -o pi_openmp.exe

pi_mpi.exe: pi_mpi.c
	$(CC_MPI) $(CFLAGS) pi_mpi.c -o pi_mpi.exe

clean:
	rm -f *.exe *.txt *.out
