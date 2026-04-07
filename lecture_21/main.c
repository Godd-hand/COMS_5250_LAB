/*
* Lab Assignment 21: Pass a Baton
*/

#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char* argv[])
{
    int comm_sz;
    int my_rank;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (comm_sz < 2) 
    {
        if (my_rank == 0) 
        {
            printf("Error: Need at least 2 processors\n");
        }
        MPI_Finalize();
        return 1;
    }

    int token;

    if (my_rank == 0)
    {
        // Rank 0 initializes the token
        token = 1000;
        int next_rank = 1;
        int prev_rank = comm_sz - 1;

        // Rank 0 sends the initial token to Rank 1
        MPI_Send(&token, 1, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
        printf("Rank %d sent token %d to Rank %d\n", my_rank, token, next_rank);

        // Rank 0 waits to receive the token back from the last rank
        MPI_Recv(&token, 1, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank %d received token %d from Rank %d\n", my_rank, token, prev_rank);
        
        // Verification
        printf("Rank 0 final token value is %d (Baton successfully passed)\n", token);
    }
    else
    {
        int prev_rank = my_rank - 1;
        int next_rank = (my_rank + 1) % comm_sz;

        // Receive token from the previous rank
        MPI_Recv(&token, 1, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Rank %d received token %d from Rank %d\n", my_rank, token, prev_rank);

        // Increment the token
        token += 1;

        // Send the token to the next rank or back to rank 0 if this is the last rank
        MPI_Send(&token, 1, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
        printf("Rank %d sent token %d to Rank %d\n", my_rank, token, next_rank);
    }

    MPI_Finalize();
    return 0;
}