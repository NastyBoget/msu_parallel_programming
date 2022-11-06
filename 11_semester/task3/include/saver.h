#ifndef TASK3_SAVER_H
#define TASK3_SAVER_H

#include <fstream>
#include "grid.h"
#include "index.h"
#include "functions.h"


class Saver {

    Grid g; // grid with the main data related to the solving
    Index ind;
    Functions f;

public:

    explicit Saver(Grid g) : g(g), ind(g), f(g) {}

    void SaveLayer(const double *layer, double t, const char *filename) const {
        std::ofstream fout(filename);

        fout << "{" << std::endl;
        fout << "    \"Lx\": " << g.L_x << ", " << std::endl;
        fout << "    \"Ly\": " << g.L_y << ", " << std::endl;
        fout << "    \"Lz\": " << g.L_z << ", " << std::endl;
        fout << "    \"N\": " << g.N << ", " << std::endl;
        fout << "    \"t\": " << t << ", " << std::endl;
        fout << "    \"u\": [" << std::endl;

        bool wasPrinted = false;

        for (int i = 0; i <= g.N; i++) {
            for (int j = 0; j <= g.N; j++) {
                for (int k = 0; k <= g.N; k++) {
                    if (wasPrinted) {
                        fout << ", " << std::endl;
                    }
                    else {
                        wasPrinted = true;
                    }

                    fout << "    " << layer[ind(i, j, k)];
                }
            }
        }

        fout << std::endl;
        fout << "    ]" << std::endl;
        fout << "}" << std::endl;

        fout.close();
    }

    void SaveAnalyticalValues(double t, const char *filename) const {
        double *u_copy = new double[(g.N + 1) * (g.N + 1) * (g.N + 1)];
        #pragma omp parallel for collapse(3)
        for (int i = 0; i <= g.N; i++)
            for (int j = 0; j <= g.N; j++)
                for (int k = 0; k <= g.N; k++)
                    u_copy[ind(i, j, k)] = f.AnalyticalSolution(i * g.h_x, j * g.h_y, k * g.h_z, t);

        SaveLayer(u_copy, t, filename);
    }

    void SaveDifferenceValues(const double *u, double t, const char *filename) const {
        double *u_copy = new double[(g.N + 1) * (g.N + 1) * (g.N + 1)];
        #pragma omp parallel for collapse(3)
        for (int i = 0; i <= g.N; i++)
            for (int j = 0; j <= g.N; j++)
                for (int k = 0; k <= g.N; k++)
                    u_copy[ind(i, j, k)] = fabs(u[ind(i, j, k)] - f.AnalyticalSolution(i * g.h_x, j * g.h_y, k * g.h_z, t));

        SaveLayer(u_copy, t, filename);
    }
};

#endif //TASK3_SAVER_H
