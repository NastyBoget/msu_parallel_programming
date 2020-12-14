#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void swap(int *x, int *y)
{
    int tmp = *x;
    *x = *y;
    *y = tmp;
}

void sort(int size, int* arr)
{
    int i, j;
    for (i = size - 2; i >= 0; i--)
        for (j = 0; j <= i; j++)
            if (arr[j] > arr[j + 1])
                swap (&arr[j], &arr[j + 1]);
}

void merge(const int *in_a, int len_a, const int *in_b, int len_b, int *out) {
    int i, j;
    int out_count = 0;
    for (i = 0, j = 0; i < len_a; i++) {
        while ((in_b[j] < in_a[i]) && j < len_b) {
            out[out_count++] = in_b[j++];
        }
        out[out_count++] = in_a[i];
    }
    while (j < len_b) {
        out[out_count++] = in_b[j++];
    }
}

void pairwise_exchange(int local_n, int *local_a, int send_rank, int recv_rank)
{
// процесс-отправитель отправляет свой массив и ждет результата
// процесс-получатель сортирует массивы и возвращает нужную половину
    int rank;
    int remote[local_n];
    int buf_all[2 * local_n];
    const int merge_tag = 1;
    const int sorted_tag = 2;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == send_rank) {
        MPI_Send(local_a, local_n, MPI_INT, recv_rank, merge_tag, MPI_COMM_WORLD);
        MPI_Recv(local_a, local_n, MPI_INT, recv_rank, sorted_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(remote, local_n, MPI_INT, send_rank, merge_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        merge(local_a, local_n, remote, local_n, buf_all);
        int their_start = 0, my_start = local_n;
        if (send_rank > rank) {
            their_start = local_n;
            my_start = 0;
        }
        MPI_Send(&(buf_all[their_start]), local_n, MPI_INT, send_rank, sorted_tag, MPI_COMM_WORLD);
        for (int i = my_start; i < my_start + local_n; i++)
            local_a[i - my_start] = buf_all[i];
    }
}

void parallel_odd_even_sort(int n, int *a, int rank, int tasks)
{
// предполагается, что размер массива делится нацело на число процессов
    int i;
    int *local_a;
// получаем номер процесса
    local_a = malloc(n / tasks * sizeof(int));
// распределение массива по процессам
    MPI_Scatter(a, n / tasks, MPI_INT, local_a, n / tasks, MPI_INT, 0, MPI_COMM_WORLD);
// сортируем часть массива, которая досталась процессу
    sort(n / tasks, local_a);
// четно-нечетная сортировка
    for (i = 1; i <= tasks; i++) {
        if ((i + rank) % 2 == 0) {  // у номера процесса и номера итерации одинаковая четность
            if (rank < tasks - 1) {
                pairwise_exchange(n / tasks, local_a, rank, rank + 1);
            }
        } else if (rank > 0) {
            pairwise_exchange(n / tasks, local_a, rank - 1, rank);
        }
    }
// собираем части массива в один
    MPI_Gather(local_a, n / tasks, MPI_INT, a, n / tasks, MPI_INT, 0, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, tasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);
    int size = atoi(argv[1]);
    int *a = NULL;
    if (rank == 0) {
        a = malloc(size * sizeof(int));
// Заполнение массива случайными числами
        for (int i = 0; i < size; i++)
            a[i] = rand() % size;
        double start_time = MPI_Wtime();
        parallel_odd_even_sort(size, a, rank, tasks);
        double end_time = MPI_Wtime();
        printf("Parallel time: %f\n", end_time - start_time);
        //for (int i = 0; i < size; ++i) {
        //    printf("%d ", a[i]);
        //}
        //printf("\n");
        free(a);
    } else {
        parallel_odd_even_sort(size, a, rank, tasks);
    }
    MPI_Finalize();
    return 0;
}