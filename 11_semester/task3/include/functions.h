#ifndef TASK3_FUNCTIONS_H
#define TASK3_FUNCTIONS_H

#include <cmath>
#include "grid.h"


class Functions {

    double a; // a_t in the analytical solution
    Grid g; // grid with the main data related to the solving

public:

    explicit Functions(Grid g) : g(g) {
        a = M_PI * sqrt(1.0 / (g.L_x * g.L_x) + 4.0 / (g.L_y * g.L_y) + 9.0 / (g.L_z * g.L_z));
    }

    double AnalyticalSolution(double x, double y, double z, double t) const {
        return sin(M_PI * x / g.L_x) * sin(2 * M_PI * y / g.L_y) * sin(3 * M_PI * z / g.L_z) * cos(a * t);
    }

    double Phi(double x, double y, double z) const {
        return AnalyticalSolution(x, y, z, 0);
    }

};


#endif //TASK3_FUNCTIONS_H
