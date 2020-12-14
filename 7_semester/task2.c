#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>

#define KILLED_PROCESS 2

int rank, tasks;
char filename[10];
unsigned error_occured = 0;
MPI_Comm main_comm;

void itoa(int n, char s[])
{
    int i = 0;
    do {
        s[i++] = n % 10 + '0';
    } while ((n /= 10) > 0);
    s[i] = '\0';
    int j, k;
    char c;
    for (j = 0, k = strlen(s)-1; j < k; j++, k--) {
        c = s[j];
        s[j] = s[k];
        s[k] = c;
    }
}

static void err_handler(MPI_Comm *pcomm, int *perr, ...) {
    error_occured = 1;
    int err = *perr;
    char errstr[MPI_MAX_ERROR_STRING];
    int size, nf, len;
    MPI_Group group_f;

    MPI_Comm_size(main_comm, &size);
    MPIX_Comm_failure_ack(main_comm);
    MPIX_Comm_failure_get_acked(main_comm, &group_f);
    MPI_Group_size(group_f, &nf);
    MPI_Error_string(err, errstr, &len);
    printf("\nRank %d / %d: Notified of error %s. %d found dead\n", rank, size, errstr, nf);

    // создаем новый коммуникатор без вышедшего из строя процесса
    MPIX_Comm_shrink(main_comm, &main_comm);
    MPI_Comm_rank(main_comm, &rank);
    itoa(rank, filename);
    strcat(filename, ".txt");
}

void print_array_to_file(int *a, int len_a, char *filename) {
    FILE *fp = fopen(filename, "w");
    for (int i = 0; i < len_a; ++i) {
        fprintf(fp, "%d ", a[i]);
    }
    fclose(fp);
}

void read_array_from_file(int *a, int len_a, char *filename) {
    FILE *fp = fopen(filename, "r");
    for (int i = 0; i < len_a; ++i) {
        fscanf(fp, "%d", &a[i]);
    }
    fclose(fp);
}

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
                swap(&arr[j], &arr[j + 1]);
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
    int remote[local_n];
    int buf_all[2 * local_n];
    const int merge_tag = 1;
    const int sorted_tag = 2;
    if (rank == send_rank) {
        MPI_Send(local_a, local_n, MPI_INT, recv_rank, merge_tag, main_comm);
        if (error_occured == 1) {
            return;
        }
        MPI_Recv(local_a, local_n, MPI_INT, recv_rank, sorted_tag, main_comm, MPI_STATUS_IGNORE);
    } else {
        MPI_Recv(remote, local_n, MPI_INT, send_rank, merge_tag, main_comm, MPI_STATUS_IGNORE);
        if (error_occured == 1) {
            return;
        }
        merge(local_a, local_n, remote, local_n, buf_all);
        int their_start = 0, my_start = local_n;
        if (send_rank > rank) {
            their_start = local_n;
            my_start = 0;
        }
        MPI_Send(&(buf_all[their_start]), local_n, MPI_INT, send_rank, sorted_tag, main_comm);
        for (int i = my_start; i < my_start + local_n; i++)
            local_a[i - my_start] = buf_all[i];
    }
}

void parallel_odd_even_sort(int n, int *a)
{
    // размер массива делится нацело на число процессов
    int *local_a = malloc(n / (tasks - 1) * sizeof(int));
    int *sendcounts = malloc(tasks * sizeof(int));
    int *displs = malloc(tasks * sizeof(int));
    displs[0] = 0;
    sendcounts[0] = n / (tasks - 1);
    for (int i = 1; i < tasks - 1; ++i) {
        displs[i] = displs[i - 1] + n / (tasks - 1);
        sendcounts[i] = n / (tasks - 1);
    }
    displs[tasks - 1] = displs[tasks - 2] + n / (tasks - 1);
    sendcounts[tasks - 1] = 0;
    // распределение массива по процессам, последнему процессу ничего не даем
    if (rank == tasks - 1) {
        MPI_Scatterv(a, sendcounts, displs, MPI_INT, local_a, 0, MPI_INT, 0, main_comm);
    } else {
        MPI_Scatterv(a, sendcounts, displs, MPI_INT, local_a, n / (tasks - 1), MPI_INT, 0, main_comm);
        // сортируем часть массива, которая досталась процессу
        sort(n / (tasks - 1), local_a);
        print_array_to_file(local_a, n / (tasks - 1), filename);
    }
    // убиваем один из процессов
    if (rank == KILLED_PROCESS) {
        raise(SIGKILL);
    }
    // четно-нечетная сортировка
    for (int i = 1; i < tasks - 1; i++) {
        // если случилась ошибка в каком-то процессе, соседние процессы заново читают данные из файла и делают итерацию
        checkpoint:
        MPI_Barrier(main_comm);
        read_array_from_file(local_a, n / (tasks - 1), filename);
        if ((i + rank) % 2 == 0) {  // у номера процесса и номера итерации одинаковая четность
            if (rank < tasks - 2) {
                pairwise_exchange(n / (tasks - 1), local_a, rank, rank + 1);
                if (error_occured == 1) {
                    error_occured = 0;
                    goto checkpoint;
                }
            }
        } else if (rank > 0) {
            pairwise_exchange(n / (tasks - 1), local_a, rank - 1, rank);
            if (error_occured == 1) {
                error_occured = 0;
                goto checkpoint;
            }
        }
        MPI_Barrier(main_comm);
        print_array_to_file(local_a, n / (tasks - 1), filename);
    }
    MPI_Barrier(main_comm);
    MPI_Gather(local_a, n / (tasks - 1), MPI_INT, a, n / (tasks - 1), MPI_INT, 0, main_comm);
    free(local_a);
}

// Задание: при запуске программы на счет сразу запустить некоторое дополнительное количество MPI-процессов,
// которые использовать в случае сбоя (в нашем случае последний процесс считается резервным).
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &tasks);
    main_comm = MPI_COMM_WORLD;

    // устанавливаем обработчик ошибок
    MPI_Errhandler errh;
    MPI_Comm_create_errhandler(err_handler, &errh);
    MPI_Comm_set_errhandler(main_comm, errh);
    MPI_Barrier(main_comm);
    // для каждого процесса формируем имя файла для записи данных контрольных точек
    itoa(rank, filename);
    strcat(filename, ".txt");

    int size = atoi(argv[1]);
    // Длина массива должна делиться нацело на число работающих процессов
    if (size % (tasks - 1) != 0) {
        size += tasks - 1 - size % (tasks - 1);
    }
    int *a = NULL;

    if (rank == 0) {
        a = malloc(size * sizeof(int));
        // Заполнение массива случайными числами
        for (int i = 0; i < size; i++)
            a[i] = rand() % size;
        printf("%s\n", "Before sorting");
        for (int i = 0; i < size; ++i) {
            printf("%d ", a[i]);
        }
        double start_time = MPI_Wtime();
        parallel_odd_even_sort(size, a);
        double end_time = MPI_Wtime();
        printf("\nParallel time: %f\n", end_time - start_time);
        printf("%s\n", "After sorting");
        for (int i = 0; i < size; ++i) {
            printf("%d ", a[i]);
        }
        printf("\n");
        free(a);
    } else {
        parallel_odd_even_sort(size, a);
    }
    MPI_Barrier(main_comm);
    remove(filename);
    MPI_Finalize();
    return 0;
}