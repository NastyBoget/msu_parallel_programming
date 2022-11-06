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

    Grid grid = Grid(L, L, L, N, T, K);
    Solver solver = Solver(grid);
    double error = solver.Solve(timeSteps, save);

    std::ofstream fout(filename, std::ios_base::app);
    fout << "L=" << L << "; N=" << N << "; error=" << error << std::endl;
    fout.close();
    std::cout << error << std::endl;
    return 0;
}
