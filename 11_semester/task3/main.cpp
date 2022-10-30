#include <iostream>
#include <cmath>
#include "include/solver.h"
#include "include/grid.h"


int main(int argc, char** argv) {

    int timeSteps = 20;
    // parameters parsing
    // N (x,y,z points), K (t points), T (time), L (Lx=Ly=Lz=pi), out_file
    int N = atoi(argv[1]);
    int K = atoi(argv[2]);
    double T = strtod(argv[3], NULL);
    double L = (argc == 6) ? strtod(argv[4], NULL) : M_PI;
    char* filename = (argc == 6) ? argv[5] : argv[4];

    Grid grid = Grid(L, L, L, N, T, K);
    Solver solver = Solver(grid);
    double error = solver.Solve(timeSteps);
    std::cout << error << std::endl;
    return 0;
}
