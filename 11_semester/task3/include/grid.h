#ifndef TASK3_GRID_H
#define TASK3_GRID_H

struct Grid {

    double L_x, L_y, L_z; // spatial sizes
    double h_x, h_y, h_z; // spatial steps
    int N; // number of spatial points

    double T; // time range
    double tau; // time step
    int K; // number of time points

    Grid(double L_x, double L_y, double L_z, int N, double T, int K) {
        this->L_x = L_x;
        this->L_y = L_y;
        this->L_z = L_z;
        this->N = N;

        this->T = T;
        this->K = K;

        this->h_x = L_x / N;
        this->h_y = L_y / N;
        this->h_z = L_z / N;

        this->tau = T / K;
    }
};

#endif //TASK3_GRID_H
