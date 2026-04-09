#!/bin/bash

#SBATCH --nodes=2   # Number of nodes to use
#SBATCH --ntasks-per-node=32   # Use 32 processor cores per node 
#SBATCH --time=0-0:10:0   # Walltime limit (DD-HH:MM:SS)
#SBATCH --qos=instruction   # Quality of service
#SBATCH --job-name="monte_carlo_pi"   # Job name to display in squeue
#SBATCH --mail-user=georgeth@iastate.edu   # Email address
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output="results_pi.out"   # Job standard output file
#SBATCH --error="errors_pi.out"   # Job standard error file
#SBATCH --partition=instruction

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load intel

# Compile the codes
make clean
make

# Number of points
N=1000000000

echo "=========================================="
echo " OpenMP PI"
echo "=========================================="

./pi_openmp.exe $N

echo ""
echo "=========================================="
echo " MPI PI"
echo "=========================================="

# Clear the MPI results text file before starting so appending is clean
rm -f mpi_pi_results.txt

mpirun -np 2 ./pi_mpi.exe $N
mpirun -np 4 ./pi_mpi.exe $N
mpirun -np 8 ./pi_mpi.exe $N
mpirun -np 16 ./pi_mpi.exe $N
mpirun -np 32 ./pi_mpi.exe $N
