#include <mpi.h>
#include <fstream>
#include <cstdio>


int main(int argc, char** argv) {
    // MPI initialization
    int proc_rank, proc_size;
    int master_rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    double start_time = MPI_Wtime();

    double end_time = MPI_Wtime();
    double delta = end_time - start_time;
    double max_time;
    MPI_Reduce(&delta, &max_time, 1, MPI_DOUBLE, MPI_MAX, master_rank, MPI_COMM_WORLD);

    if (proc_rank == master_rank) {
    }
    MPI_Finalize();
    return 0;
}
