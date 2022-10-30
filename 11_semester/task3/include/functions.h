#ifndef TASK3_FUNCTIONS_H
#define TASK3_FUNCTIONS_H

#include <cmath>
#include "grid.h"
#include "index.h"


class Functions {

    double a; // a_t in the analytical solution
    Grid g; // grid with the main data related to the solving
    Index ind; // for getting flattened indexes in the 3-d array

public:

    explicit Functions(Grid g) : g(g), ind(g) {
        a = M_PI * sqrt(1.0 / (g.L_x * g.L_x) + 4.0 / (g.L_y * g.L_y) + 9.0 / (g.L_z * g.L_z));
    }

    double AnalyticalSolution(double x, double y, double z, double t) const {
        return sin(M_PI * x / g.L_x) * sin(2 * M_PI * y / g.L_y) * sin(3 * M_PI * z / g.L_z) * cos(a * t);
    }

    double Phi(double x, double y, double z) const {
        return AnalyticalSolution(x, y, z, 0);
    }

    double LaplaceOperator(const double *u, int i, int j, int k) const {
        return (u[ind(i - 1, j, k)] - 2 * u[ind(i, j, k)] + u[ind(i + 1, j, k)]) / (g.h_x * g.h_x) +
               (u[ind(i, j - 1, k)] - 2 * u[ind(i, j, k)] + u[ind(i, j + 1, k)]) / (g.h_y * g.h_y) +
               (u[ind(i, j, k - 1)] - 2 * u[ind(i, j, k)] + u[ind(i, j, k + 1)]) / (g.h_z * g.h_z);
    }

};




#endif //TASK3_FUNCTIONS_H
