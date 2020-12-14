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

    int result = 0;
    int other_coords[2];
    int other_rank = 0;
    MPI_Status status;
    // обмен значениями за 6 шагов
    // шаг 1
    other_coords[1] = coords[1];
    if (coords[0] == 0 || coords[0] == 3) {
        if (coords[0] == 0) {
            other_coords[0] = coords[0] + 1;
            MPI_Cart_rank(comm, other_coords, &other_rank);
            MPI_Send(&a, 1, MPI_INT, other_rank, 0, comm);
        } else {
            other_coords[0] = coords[0] - 1;
            MPI_Cart_rank(comm, other_coords, &other_rank);
            MPI_Send(&a, 1, MPI_INT, other_rank, 0, comm);
        }
    } else {
        if (coords[0] == 1) {
            other_coords[0] = coords[0] - 1;
            MPI_Cart_rank(comm, other_coords, &other_rank);
            MPI_Recv(&result, 1, MPI_INT, other_rank, 0, comm, &status);
        } else {
            other_coords[0] = coords[0] + 1;
            MPI_Cart_rank(comm, other_coords, &other_rank);
            MPI_Recv(&result, 1, MPI_INT, other_rank, 0, comm, &status);
        }
        if (result > a) {
            a = result;
        }
    }
    MPI_Barrier(comm);
    // шаг 2
    if (coords[0] == 2) {
        other_coords[0] = coords[0] - 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&a, 1, MPI_INT, other_rank, 0, comm);
    }
    if (coords[0] == 1) {
        other_coords[0] = coords[0] + 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&result, 1, MPI_INT, other_rank, 0, comm, &status);
        if (result > a) {
            a = result;
        }
    }
    MPI_Barrier(comm);
    // шаг 3
    if (coords[0] == 1 && (coords[1] == 0 || coords[1] == 3)) {
        other_coords[0] = coords[0];
        if (coords[1] == 0) {
            other_coords[1] = coords[1] + 1;
        } else {
            other_coords[1] = coords[1] - 1;
        }
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&a, 1, MPI_INT, other_rank, 0, comm);
    }
    if (coords[0] == 1 && (coords[1] == 1 || coords[1] == 2)) {
        other_coords[0] = coords[0];
        if (coords[1] == 1) {
            other_coords[1] = coords[1] - 1;
        } else {
            other_coords[1] = coords[1] + 1;
        }
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&result, 1, MPI_INT, other_rank, 0, comm, &status);
        if (result > a) {
            a = result;
        }
    }
    MPI_Barrier(comm);
    // шаг 4
    if (coords[0] == 1 && coords[1] == 3) {
        other_coords[0] = coords[0];
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&a, 1, MPI_INT, other_rank, 0, comm);
    }
    if (coords[0] == 1 && coords[1] == 1) {
        other_coords[0] = coords[0];
        other_coords[1] = 3;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&result, 1, MPI_INT, other_rank, 0, comm, &status);
        if (result > a) {
            a = result;
        }
    }
    MPI_Barrier(comm);
    // шаг 5
    if (coords[0] == 1 && coords[1] == 1) {
        other_coords[0] = 0;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&a, 1, MPI_INT, other_rank, 0, comm);
    }
    if (coords[0] == 0 && coords[1] == 1) {
        other_coords[0] = 1;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&result, 1, MPI_INT, other_rank, 0, comm, &status);
        if (result > a) {
            a = result;
        }
    }
    MPI_Barrier(comm);
    // шаг 6
    if (coords[0] == 0 && coords[1] == 1) {
        other_coords[0] = 0;
        other_coords[1] = 0;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Send(&a, 1, MPI_INT, other_rank, 0, comm);
    }
    if (coords[0] == 0 && coords[1] == 0) {
        other_coords[0] = 0;
        other_coords[1] = 1;
        MPI_Cart_rank(comm, other_coords, &other_rank);
        MPI_Recv(&result, 1, MPI_INT, other_rank, 0, comm, &status);
        if (result > a) {
            a = result;
        }
    }
    if (coords[0] == 0 && coords[1] == 0) {
        printf("Max result: %d\n", a);
    }
    MPI_Finalize();
    return 0;
}
