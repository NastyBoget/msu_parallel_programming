#include <mpi.h>
#include <fstream>
#include <cmath>
#include "include/grid.h"
#include "include/solver_mpi.h"


int main(int argc, char** argv) {
    // MPI initialization
    int proc_rank, proc_size;
    int master_rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    int timeSteps = 20;
    int K = 1000;
    double T = 1.0;

    // parameters parsing
    // N (x,y,z points), L (Lx=Ly=Lz=pi), out_file
    int N = atoi(argv[1]);
    double L = (argc == 4) ? strtod(argv[2], NULL) : M_PI;
    char* filename = (argc == 4) ? argv[3] : argv[2];

    double start_time = MPI_Wtime();

    // equation solving
    Grid grid = Grid(L, L, L, N, T, K);
    SolverMPI solver = SolverMPI(grid);
    double error = solver.Solve(timeSteps);

    // elapsed time computing
    double end_time = MPI_Wtime();
    double delta = end_time - start_time;
    double max_time;
    MPI_Reduce(&delta, &max_time, 1, MPI_DOUBLE, MPI_MAX, master_rank, MPI_COMM_WORLD);

    if (proc_rank == master_rank) {
        std::ofstream fout(filename, std::ios_base::app);
        fout << N << " " << proc_size << " " << error << " " << max_time << std::endl;
        fout.close();
    }
    MPI_Finalize();
    return 0;
}
