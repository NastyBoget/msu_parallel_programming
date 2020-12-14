#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
void compare_exchange(int *i, int *j)
{
    if (*i > *j) {
        int tmp = *i;
        *i = *j;
        *j = tmp;
    }
}
//Функция для алгоритма чет-нечетной сортировки
void parallel_odd_even_sort(int* arr, int n, int num_threads)
{
    int upper_bound;
    if (n % 2 == 0)
        upper_bound = n / 2 - 1;
    else
        upper_bound = n / 2;
    for (int i = 0; i < n; i++) {
        if (i % 2 == 0) {
//четная итерация
#pragma omp parallel for num_threads(num_threads)
        for (int j = 0; j < n / 2; j++)
            compare_exchange(&arr[2 * j], &arr[2 * j + 1]);
        } else {
//нечетная итерация
#pragma omp parallel for num_threads(num_threads)
        for (int k = 0; k < upper_bound; k++)
            compare_exchange(&arr[2 * k + 1], &arr[2 * k + 2]);
        }
    }
}
int main(int argc, char **argv)
{
	int num_threads = argv[1][0] - '0';
	size_t size = atoi(argv[2]);
    int a[size];
    // Заполнение массива случайными числами
    for (int i = 0; i < size; i++)
        a[i] = rand() % size;
    double start = omp_get_wtime();
    parallel_odd_even_sort(a, size, num_threads); // вызов функции сортировки
    double end = omp_get_wtime();
    printf("\n%f\n", end - start);
    return 0;
}
