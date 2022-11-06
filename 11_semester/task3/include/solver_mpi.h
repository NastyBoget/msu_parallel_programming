#ifndef TASK3_SOLVER_MPI_H
#define TASK3_SOLVER_MPI_H


#include <mpi.h>
#include <iostream>
#include <vector>
#include "grid.h"
#include "index.h"
#include "functions.h"
#include "splitter.h"


class SolverMPI {

    Functions f; // the main functions for the task
    Grid g; // grid with the main data related to the solving
    Index ind; // for getting flattened indexes in the 3-d array
    std::vector< std::vector<double> > u; // function in the equation for 3 steps in the algorithm
    std::vector< std::pair<int, Block> > blocksToSend; // neighbours with proc_rank whom this process should send the data
    std::vector< std::pair<int, Block> > blocksToReceive; // neighbours with proc_rank from whom this process should receive the data
    int proc_rank, proc_size; // information connected with MPI processes

public:

    explicit SolverMPI(Grid g) : g(g), f(g), ind(g) {
        proc_rank = 0; proc_size = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
    }

    std::vector<double> GetSendData(int uInd, const Block block, const Block otherBlock) const {
        std::vector<double> dataToSend(otherBlock.size);

        #pragma omp parallel for collapse(3)
        for (int i = otherBlock.x_min; i <= otherBlock.x_max; i++)
            for (int j = otherBlock.y_min; j <= otherBlock.y_max; j++)
                for (int k = otherBlock.z_min; k <= otherBlock.z_max; k++)
                    dataToSend[ind(i, j, k, otherBlock)] = u[uInd][ind(i, j, k, block)];

        return dataToSend;
    }

    std::vector< std::vector<double> > Exchange(int uInd, const Block block) const {
        std::vector< std::vector<double> > dataToReceive(blocksToReceive.size());
        std::vector<MPI_Request> requests(2);
        std::vector<MPI_Status> statuses(2);

        for (int i = 0; i < blocksToReceive.size(); i++) {
            std::vector<double> dataToSend = GetSendData(uInd, block, blocksToSend[i].second);
            dataToReceive[i] = std::vector<double>(blocksToReceive[i].second.size);
            MPI_Isend(dataToSend.data(), blocksToSend[i].second.size, MPI_DOUBLE, blocksToSend[i].first, 0, MPI_COMM_WORLD, &requests[0]);
            MPI_Irecv(dataToReceive[i].data(), blocksToReceive[i].second.size, MPI_DOUBLE, blocksToReceive[i].first, 0, MPI_COMM_WORLD, &requests[1]);
            MPI_Waitall(2, requests.data(), statuses.data());
        }
        return dataToReceive;
    }

    double FindU(int uInd, int i, int j, int k, const Block b, const std::vector< std::vector<double> > &recieved) const {

        if (b.x_min <= i and i <= b.x_max and b.y_min <= j and j <= b.y_max and b.z_min <= k and k <= b.z_max) {
            return u[uInd][ind(i, j, k, b)];
        }

        for (int r_i = 0; r_i < blocksToReceive.size(); i++) {
            Block otherB = blocksToReceive[r_i].second;

            if (i < otherB.x_min or i > otherB.x_max or
                j < otherB.y_min or j > otherB.y_max or
                k < otherB.z_min or k > otherB.z_max)
                continue;

            return recieved[r_i][ind(i, j, k, otherB)];
        }
        throw std::runtime_error("u value non found");
    }

    double LaplaceOperator(int uInd, int i, int j, int k, const Block b, const std::vector< std::vector<double> > &recieved) const {
        double dx = (FindU(uInd, i, j - 1, k, b, recieved) - 2 * u[uInd][ind(i, j, k, b)] + FindU(uInd, i, j + 1, k, b, recieved)) / (g.h_y * g.h_y);
        double dy = (FindU(uInd, i - 1, j, k, b, recieved) - 2 * u[uInd][ind(i, j, k, b)] + FindU(uInd, i + 1, j, k, b, recieved)) / (g.h_x * g.h_x);
        double dz = (FindU(uInd,i, j, k - 1, b, recieved) - 2 * u[uInd][ind(i, j, k, b)] + FindU(uInd, i, j, k + 1, b, recieved)) / (g.h_z * g.h_z);
        return dx + dy + dz;
    }

    double ComputeLayerError(int uInd, double t, const Block b) const {
        double errorLocal = 0;
        // maximum difference between values of u analytical and u computed
        #pragma omp parallel for collapse(3) reduction(max: error)
        for (int i = b.x_min; i <= b.x_max; i++)
            for (int j = b.y_min; j <= b.y_max; j++)
                for (int k = b.z_min; k <= b.z_max; k++)
                    errorLocal = std::max(errorLocal, fabs(u[uInd][ind(i, j, k, b)] -
                                                           f.AnalyticalSolution(i * g.h_x, j * g.h_y, k * g.h_z, t)));
        double error;
        MPI_Reduce(&errorLocal, &error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        return errorLocal;
    }

    void FillBoundaryValues(int uInd, double t, const Block b) {
        // Variant 3 -> first kind for x, periodic for y, first kind for z
        if (b.x_min == 0) {
            #pragma omp parallel for collapse(2)
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u[uInd][ind(b.x_min, i, j, b)] = 0;
        }

        if (b.x_max == g.N) {
            #pragma omp parallel for collapse(2)
            for (int i = b.y_min; i <= b.y_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u[uInd][ind(b.x_max, i, j, b)] = 0;
        }

        if (b.y_min == 0) {
            #pragma omp parallel for collapse(2)
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u[uInd][ind(i, b.y_min, j, b)] = f.AnalyticalSolution(i * g.h_x, 0, j * g.h_z, t);
        }

        if (b.y_max == g.N) {
            #pragma omp parallel for collapse(2)
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.z_min; j <= b.z_max; j++)
                    u[uInd][ind(i, b.y_max, j, b)] = f.AnalyticalSolution(i * g.h_x, g.L_y, j * g.h_z, t);
        }

        if (b.z_min == 0) {
            #pragma omp parallel for collapse(2)
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u[uInd][ind(i, j, b.z_min, b)] = 0;
        }

        if (b.z_max == g.N) {
            #pragma omp parallel for collapse(2)
            for (int i = b.x_min; i <= b.x_max; i++)
                for (int j = b.y_min; j <= b.y_max; j++)
                    u[uInd][ind(i, j, b.z_max, b)] = 0;
        }
    }

    void InitValues(const Block b) {
        // boundary (i = 0,N or j = 0,N or k = 0,N)
        FillBoundaryValues(0, 0, b);
        FillBoundaryValues(1, g.tau, b);

        // compute the boundaries of the current block
        int x1 = std::max(b.x_min, 1); int x2 = std::min(b.x_max, g.N - 1);
        int y1 = std::max(b.y_min, 1); int y2 = std::min(b.y_max, g.N - 1);
        int z1 = std::max(b.z_min, 1); int z2 = std::min(b.z_max, g.N - 1);

        // initial values for inner points in u_0
        #pragma omp parallel for collapse(3)
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                for (int k = z1; k <= z2; k++)
                    u[0][ind(i, j, k, b)] = f.Phi(i * g.h_x, j * g.h_y, k * g.h_z);

        std::vector< std::vector<double> > recieved = Exchange(0, b);
        // initial values for inner points in u_1
        #pragma omp parallel for collapse(3)
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                for (int k = z1; k <= z2; k++)
                    u[1][ind(i, j, k, b)] = u[0][ind(i, j, k, b)] + g.tau * g.tau / 2 * LaplaceOperator(0, i, j, k, b, recieved);
    }

    void GetNextU(int step, const Block b) {
        // compute the boundaries of the current block
        int x1 = std::max(b.x_min, 1); int x2 = std::min(b.x_max, g.N - 1);
        int y1 = std::max(b.y_min, 1); int y2 = std::min(b.y_max, g.N - 1);
        int z1 = std::max(b.z_min, 1); int z2 = std::min(b.z_max, g.N - 1);

        std::vector< std::vector<double> > received = Exchange((step + 2) % 3, b);
        // calculate u_n+1 inside the area
        #pragma omp parallel for collapse(3)
        for (int i = x1; i <= x2; i++)
            for (int j = y1; j <= y2; j++)
                for (int k = z1; k <= z2; k++)
                    u[step % 3][ind(i, j, k, b)] = 2 * u[(step + 2) % 3][ind(i, j, k, b)] -
                                                   u[(step + 1) % 3][ind(i, j, k, b)] +
                                                   g.tau * g.tau * LaplaceOperator((step + 2) % 3, i, j, k, b, received);
        FillBoundaryValues(step % 3, step * g.tau, b);
    }

    void GetNeighbours(const std::vector<Block> &blocks) {
        Block block = blocks[proc_rank];

        for (int i = 0; i < proc_size; i++) {
            if (i == proc_rank)
                continue;

            Block otherBlock = blocks[i];
            if (block.x_min == otherBlock.x_max + 1 or otherBlock.x_min == block.x_max + 1) {
                int xSend = block.x_min == otherBlock.x_max + 1 ? block.x_min : block.x_max;
                int xRecv = otherBlock.x_min == block.x_max + 1 ? otherBlock.x_min : otherBlock.x_max;
                int y_min = std::max(block.y_min, otherBlock.y_min); int y_max = std::min(block.y_max, otherBlock.y_max);
                int z_min = std::max(block.z_min, otherBlock.z_min); int z_max = std::min(block.z_max, otherBlock.z_max);
                // add block as a rectangle (it's a border between two processes)
                blocksToSend.emplace_back(i, Block(xSend, xSend, y_min, y_max, z_min, z_max));
                blocksToReceive.emplace_back(i, Block(xRecv, xRecv, y_min, y_max, z_min, z_max));
                continue;
            }
            if (block.y_min == otherBlock.y_max + 1 or otherBlock.y_min == block.y_max + 1) {
                int ySend = block.y_min == otherBlock.y_max + 1 ? block.y_min : block.y_max;
                int yRecv = otherBlock.y_min == block.y_max + 1 ? otherBlock.y_min : otherBlock.y_max;
                int x_min = std::max(block.x_min, otherBlock.x_min); int x_max = std::min(block.x_max, otherBlock.x_max);
                int z_min = std::max(block.z_min, otherBlock.z_min); int z_max = std::min(block.z_max, otherBlock.z_max);
                blocksToSend.emplace_back(i, Block(x_min, x_max, ySend, ySend, z_min, z_max));
                blocksToReceive.emplace_back(i, Block(x_min, x_max, yRecv, yRecv, z_min, z_max));
                continue;
            }
            if (block.z_min == otherBlock.z_max + 1 or otherBlock.z_min == block.z_max + 1) {
                int zSend = block.z_min == otherBlock.z_max + 1 ? block.z_min : block.z_max;
                int zRecv = otherBlock.z_min == block.z_max + 1 ? otherBlock.z_min : otherBlock.z_max;
                int x_min = std::max(block.x_min, otherBlock.x_min); int x_max = std::min(block.x_max, otherBlock.x_max);
                int y_min = std::max(block.y_min, otherBlock.y_min); int y_max = std::min(block.y_max, otherBlock.y_max);
                blocksToSend.emplace_back(i, Block(x_min, x_max, y_min, y_max, zSend, zSend));
                blocksToReceive.emplace_back(i, Block(x_min, x_max, y_min, y_max, zRecv, zRecv));
                continue;
            }
        }
    }

    double Solve(int steps) {
        // split grid between processes
        std::vector<Block> blocks;
        split(0, g.N, 0, g.N, 0, g.N, proc_size, X, blocks);
        Block block = blocks[proc_rank];

        // allocate spase for u
        u.resize(3);
        for (int i = 0; i < 3; i++)
            u[i].resize(block.size);

        // fill blocksToSend and blocksToReceive vectors
        GetNeighbours(blocks);

        // init u_0 and u_1
        InitValues(block);

        // calculate the next time layers for u
        for (int step = 2; step <= steps; step++) {
            GetNextU(step, block);
        }

        return ComputeLayerError(steps % 3, steps * g.tau, block);
    }

};


#endif //TASK3_SOLVER_MPI_H
