#include <mpi.h>
#include <fstream>
#include <cstdlib>
#include <cmath>


bool in_G(double x, double y, double z) {
    if ((x < 0) or (y < 0) or (z < 0))
        return false;
    return x*x + y*y + z*z <= 1;
}


double F(double x, double y, double z) {
    if (in_G(x, y, z))
        return sin(x*x + z*z) * y;
    return 0;
}


void fill_array_with_points(double *arr, size_t arr_len) {
    for (size_t i = 0; i < arr_len * 3; i += 3) {
        // we generate points in the [0, 1] segment
        arr[i] = rand() * 1.0 / RAND_MAX;
        arr[i + 1] = rand() * 1.0 / RAND_MAX;
        arr[i + 2] = rand() * 1.0 / RAND_MAX;
    }
}


int main(int argc, char** argv) {
    // MPI initialization
    int proc_rank, proc_size;
    int master_rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    // Run program as follows: task2 eps out_file
    double eps = strtod(argv[1], NULL);
    char* filename = argv[2];

    double start_time = MPI_Wtime();

    // the value of the integral, computed manually
    double actual_result = M_PI / 8 * (1 - sin(1));
    // the approximate value of the integral
    double result, error;

    // 0 <= x <= 1; 0 <= y <= 1; 0 <= z <= 1; parallelepiped_volume = (1 - 0) * (1 - 0) * (1 - 0)
    double parallelepiped_volume = 1.0;
    srand(pow(2, proc_rank));

    // this number may be changed
    size_t number_of_generated_points = 1000;

    double *points = new double[3 * number_of_generated_points];

    int generate_more = 1;
    double total_sum = 0;
    size_t total_points = 0;

    while (generate_more) {
        fill_array_with_points(points, number_of_generated_points);

        double local_sum = 0;
        for (size_t i = 0; i < 3 * number_of_generated_points; i += 3) {
            local_sum += F(points[i], points[i + 1], points[i + 2]);
        }
        double sum = 0;
        MPI_Reduce(&local_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, master_rank, MPI_COMM_WORLD);

        // check if we should stop
        if (proc_rank == master_rank) {
            total_sum += sum;
            total_points += number_of_generated_points * proc_size;
            result = parallelepiped_volume * total_sum / total_points;
            error = fabs(actual_result - result);
            generate_more = error >= eps ? 1 : 0;
        }
        MPI_Bcast(&generate_more, 1, MPI_INT, master_rank, MPI_COMM_WORLD);
    }

    delete[] points;

    double end_time = MPI_Wtime();
    double delta = end_time - start_time;
    double max_time;
    MPI_Reduce(&delta, &max_time, 1, MPI_DOUBLE, MPI_MAX, master_rank, MPI_COMM_WORLD);

    if (proc_rank == master_rank) {
        std::ofstream fout(filename, std::ios_base::app);
        fout << result << " " << error << " " << total_points << " " << max_time << std::endl;
        fout.close();
    }
    MPI_Finalize();
    return 0;
}
