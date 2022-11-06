#ifndef TASK3_SOLVER_H
#define TASK3_SOLVER_H

#include "grid.h"
#include "index.h"
#include "functions.h"
#include "saver.h"


class Solver {

    Functions f; // the main functions for the task
    Grid g; // grid with the main data related to the solving
    Index ind; // for getting flattened indexes in the 3-d array
    Saver saver; // for saving results
    double **u; // function in the equation for 3 steps in the algorithm
    
public:

    explicit Solver(Grid g) : g(g), f(g), ind(g), saver(g) {
        int layerSize = (g.N + 1) * (g.N + 1) * (g.N + 1);
        u = new double*[3];
        u[0] = new double[layerSize];
        u[1] = new double[layerSize];
        u[2] = new double[layerSize];
    }

    ~Solver() {
        delete[] u[0];
        delete[] u[1];
        delete[] u[2];
        delete[] u;
    }

    double LaplaceOperator(const double *func, int i, int j, int k) const {
        return (func[ind(i - 1, j, k)] - 2 * func[ind(i, j, k)] + func[ind(i + 1, j, k)]) / (g.h_x * g.h_x) +
               (func[ind(i, j - 1, k)] - 2 * func[ind(i, j, k)] + func[ind(i, j + 1, k)]) / (g.h_y * g.h_y) +
               (func[ind(i, j, k - 1)] - 2 * func[ind(i, j, k)] + func[ind(i, j, k + 1)]) / (g.h_z * g.h_z);
    }

    void FillBoundaryValues(int uInd, double t) {
        // Variant 3 -> first kind for x, periodic for y, first kind for z
        for (int i = 0; i <= g.N; i++) {
            for (int j = 0; j <= g.N; j++) {
                // for x
                u[uInd][ind(0, i, j)] = 0;
                u[uInd][ind(g.N, i, j)] = 0;
                // for y
                u[uInd][ind(i, 0, j)] = f.AnalyticalSolution(i * g.h_x, 0, j * g.h_z, t);
                u[uInd][ind(i, g.N, j)] = f.AnalyticalSolution(i * g.h_x, g.L_y, j * g.h_z, t);
                // for z
                u[uInd][ind(i, j, 0)] = 0;
                u[uInd][ind(i, j, g.N)] = 0;
            }
        }
    }

    void InitValues() {
        // boundary (i = 0,N or j = 0,N or k = 0,N)
        FillBoundaryValues(0, 0);
        FillBoundaryValues(1, g.tau);

        // initial values for inner points in u_0
        #pragma omp parallel for collapse(3)
        for (int i = 1; i < g.N; i++)
            for (int j = 1; j < g.N; j++)
                for (int k = 1; k < g.N; k++)
                    u[0][ind(i, j, k)] = f.Phi(i * g.h_x, j * g.h_y, k * g.h_z);
        // initial values for inner points in u_1
        #pragma omp parallel for collapse(3)
        for (int i = 1; i < g.N; i++)
            for (int j = 1; j < g.N; j++)
                for (int k = 1; k < g.N; k++)
                    u[1][ind(i, j, k)] = u[0][ind(i, j, k)] + g.tau * g.tau / 2 * LaplaceOperator(u[0], i, j, k);
    }

    double ComputeLayerError(int uInd, double t) {
        double error = 0;
        // maximum difference between values of u analytical and u computed
        #pragma omp parallel for collapse(3) reduction(max: error)
        for (int i = 0; i <= g.N; i++)
            for (int j = 0; j <= g.N; j++)
                for (int k = 0; k <= g.N; k++)
                    error = std::max(error, fabs(u[uInd][ind(i, j, k)] -
                                                 f.AnalyticalSolution(i * g.h_x, j * g.h_y, k * g.h_z, t)));
        return error;
    }

    double Solve(int steps, bool save) {
        // init u_0 and u_1
        InitValues();

        // calculate the next time layers for u
        for (int step = 2; step <= steps; step++) {
            #pragma omp parallel for collapse(3)
            for (int i = 1; i < g.N; i++)
                for (int j = 1; j < g.N; j++)
                    for (int k = 1; k < g.N; k++)
                        // calculate u_n+1 inside the area
                        u[step % 3][ind(i, j, k)] = 2 * u[(step + 2) % 3][ind(i, j, k)] - u[(step + 1) % 3][ind(i, j, k)] +
                                                    g.tau * g.tau * LaplaceOperator(u[(step + 2) % 3], i, j, k);

            FillBoundaryValues(step % 3, step * g.tau);
        }

        double error = ComputeLayerError(steps % 3, steps * g.tau);
        if (save) {
            saver.SaveLayer(u[steps % 3], steps * g.tau, "numerical.json");
            saver.SaveDifferenceValues(u[steps % 3], steps * g.tau, "difference.json");
            saver.SaveAnalyticalValues(steps * g.tau, "analytical.json");
        }

        return error;
    }

};

#endif //TASK3_SOLVER_H
