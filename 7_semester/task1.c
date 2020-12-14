#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define SIZE 4

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, tasks;
    MPI_Comm comm;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);

    int size[2] = {SIZE, SIZE};
    int periodic[2] = {0};
    // создание транспьютерной матрицы
    MPI_Cart_create(MPI_COMM_WORLD, 2, size, periodic, 0, &comm);
    int coords[2];
    MPI_Cart_coords(comm, rank, 2, coords);
    srand(rank + 4);
    int a = rand() % 1000;
    printf("Coordinates for process %d: (%d, %d)\n", rank, coords[0], coords[1]);
    printf("a[%d][%d] = %d\n", coords[0], coords[1], a);

    // обмен значениями за 6 шагов
    // шаг 1
    int result = 0;
    MPI_Status status;
    if (rank < 4 || rank > 11) {
        if (rank < 4) {
            MPI_Send(&a, 1, MPI_INT, rank + 4, 0, MPI_COMM_WORLD);
        } else {
            MPI_Send(&a, 1, MPI_INT, rank - 4, 0, MPI_COMM_WORLD);
        }
    } else {
        if (rank < 8) {
            MPI_Recv(&result, 1, MPI_INT, rank - 4, 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(&result, 1, MPI_INT, rank + 4, 0, MPI_COMM_WORLD, &status);
        }
        if (result > a) {
            a = result;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // шаг 2
    if (rank > 7 && rank < 12) {
        MPI_Send(&a, 1, MPI_INT, rank - 4, 0, MPI_COMM_WORLD);
    }
    if (rank > 3 && rank < 8) {
        MPI_Recv(&result, 1, MPI_INT, rank + 4, 0, MPI_COMM_WORLD, &status);
        if (result > a) {
            a = result;
        }
    }
    // шаг 3
    if (rank == 4 || rank == 6) {
        MPI_Send(&a, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
    if (rank == 5 || rank == 7) {
        MPI_Recv(&result, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &status);
        if (result > a) {
            a = result;
        }
    }
    // шаг 4
    if (rank == 7) {
        MPI_Send(&a, 1, MPI_INT, 5, 0, MPI_COMM_WORLD);
    }
    if (rank == 5) {
        MPI_Recv(&result, 1, MPI_INT, 7, 0, MPI_COMM_WORLD, &status);
        if (result > a) {
            a = result;
        }
    }
    // шаг 5
    if (rank == 5) {
        MPI_Send(&a, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    if (rank == 1) {
        MPI_Recv(&result, 1, MPI_INT, 5, 0, MPI_COMM_WORLD, &status);
        if (result > a) {
            a = result;
        }
    }
    // шаг 6
    if (rank == 1) {
        MPI_Send(&a, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        MPI_Recv(&result, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
        if (result > a) {
            a = result;
        }
    }
    if (rank == 0) {
        printf("Max result: %d\n", a);
    }
    MPI_Finalize();
    return 0;
}
