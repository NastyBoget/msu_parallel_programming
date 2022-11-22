#include <chrono>
#include <iostream>
#include <cmath>
#include "include/solver.h"
#include "include/grid.h"


int main(int argc, char** argv) {

    bool save = false;
    int timeSteps = 20;
    int K = 1000;
    double T = 1.0;
    // parameters parsing
    // N (x,y,z points), L (Lx=Ly=Lz=pi), out_file
    int N = atoi(argv[1]);
    double L = (argc == 4) ? strtod(argv[2], NULL) : M_PI;
    char* filename = (argc == 4) ? argv[3] : argv[2];

    auto start_time = std::chrono::high_resolution_clock::now();

    Grid grid = Grid(L, L, L, N, T, K);
    Solver solver = Solver(grid);
    double error = solver.Solve(timeSteps, save);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double delta = duration.count() / 1000000.0;

    std::ofstream fout(filename, std::ios_base::app);
    fout << N << " " << 0 << " " << error << " " << delta << std::endl;
    fout.close();

    return 0;
}
