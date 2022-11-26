#ifndef TASK4_SOLVER_H
#define TASK4_SOLVER_H

#include <mpi.h>
#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "grid.h"
#include "splitter.h"

const int num_threads = 4;


struct Neighbour {
    int number;
    Block send;
    Block receive;
    Neighbour(int num, Block b1, Block b2) : number(num), send(b1), receive(b2) {}
};


__device__ double analyticalSolution(double x, double y, double z, double t, double a, const Grid g) {
    return sin(M_PI * x / g.L_x) * sin(2 * M_PI * y / g.L_y) * sin(3 * M_PI * z / g.L_z) * cos(a * t);
}


__device__ double phi(double x, double y, double z, double a, const Grid g) {
    return analyticalSolution(x, y, z, 0, a, g);
}


__host__ __device__ inline int ind(int i, int j, int k, const Block b) {
    // get the linear index inside the array of the given grid block
    return (i - b.x_min) * b.y_size * b.z_size + (j - b.y_min) * b.z_size + (k - b.z_min);
}


thrust::host_vector<double> getSendData(thrust::host_vector<double> &u, const Block block, const Block otherBlock) {
    thrust::host_vector<double> dataToSend(otherBlock.size);

    for (int i = otherBlock.x_min; i <= otherBlock.x_max; i++)
        for (int j = otherBlock.y_min; j <= otherBlock.y_max; j++)
            for (int k = otherBlock.z_min; k <= otherBlock.z_max; k++)
                dataToSend[ind(i, j, k, otherBlock)] = u[ind(i, j, k, block)];

    return dataToSend;
}

thrust::host_vector<double> exchange(thrust::host_vector<double> &u, const Block block,
                                     Block *send, Block *receive, int *numbers, int vsize) {
    thrust::host_vector<double> dataToReceive;
    int offset = 0;
    thrust::host_vector<MPI_Request> requests(2);
    thrust::host_vector<MPI_Status> statuses(2);

    for (int i = 0; i < vsize; i++) {
        thrust::host_vector<double> dataToSend = getSendData(u, block, static_cast<Block>(send[i]));
        dataToReceive.insert(dataToReceive.end(), static_cast<Block>(receive[i]).size, 0);

        MPI_Isend(dataToSend.data(), static_cast<Block>(send[i]).size, MPI_DOUBLE, numbers[i], 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(dataToReceive.data() + offset, static_cast<Block>(receive[i]).size, MPI_DOUBLE, numbers[i], 0, MPI_COMM_WORLD, &requests[1]);
        MPI_Waitall(2, requests.data(), statuses.data());
        offset += static_cast<Block>(receive[i]).size;
    }
    return dataToReceive;
}

__device__ double findU(double *u, int i, int j, int k, const Block b, double *recieved,
                        Block *send, Block *receive, int *numbers, int vsize) {

    if (b.x_min <= i and i <= b.x_max and b.y_min <= j and j <= b.y_max and b.z_min <= k and k <= b.z_max) {
        return u[ind(i, j, k, b)];
    }

    int offset = 0;

    for (int r_i = 0; r_i < vsize; r_i++) {
        Block otherB = receive[r_i];

        if (i < otherB.x_min or i > otherB.x_max or j < otherB.y_min or j > otherB.y_max or k < otherB.z_min or k > otherB.z_max) {
            offset += receive[r_i].size;
            continue;
        }
        return recieved[offset + ind(i, j, k, otherB)];
    }
    return 1;
}

__device__ double laplaceOperator(double *u, int i, int j, int k, const Block b, const Grid g, double *recieved,
                                  Block *send, Block *receive, int *numbers, int vsize) {
    double dx = (findU(u, i, j - 1, k, b, recieved, send, receive, numbers, vsize) - 2 * u[ind(i, j, k, b)] +
            findU(u, i, j + 1, k, b, recieved, send, receive, numbers, vsize)) / (g.h_y * g.h_y);
    double dy = (findU(u, i - 1, j, k, b, recieved, send, receive, numbers, vsize) - 2 * u[ind(i, j, k, b)] +
            findU(u, i + 1, j, k, b, recieved, send, receive, numbers, vsize)) / (g.h_x * g.h_x);
    double dz = (findU(u,i, j, k - 1, b, recieved, send, receive, numbers, vsize) - 2 * u[ind(i, j, k, b)] +
            findU(u, i, j, k + 1, b, recieved, send, receive, numbers, vsize)) / (g.h_z * g.h_z);
    return dx + dy + dz;
}

__global__ void computeLayerErrorKernel(double *u, double t, double a, const Block b, const Grid g) {
    int idx = threadIdx.x;

    for (int index = idx; index < b.size; index += num_threads) {
        int i = b.x_min + index / (b.y_size * b.z_size);
        int j = b.y_min + index % (b.y_size * b.z_size) / b.z_size;
        int k = b.z_min + index % b.z_size;
        u[ind(i, j, k, b)] = fabs(u[ind(i, j, k, b)] - analyticalSolution(i * g.h_x, j * g.h_y, k * g.h_z, t, a, g));
    }
}

double computeLayerError(thrust::host_vector<double> &u, double t, const Block b, const Grid g, double a) {
    thrust::device_vector<double> uDevice(u);
    computeLayerErrorKernel<<<1, num_threads>>>(thrust::raw_pointer_cast(&uDevice[0]), t, a, b, g);
    thrust::device_vector<double>::iterator iter = thrust::max_element(uDevice.begin(), uDevice.end());
    double errorLocal = uDevice[iter - uDevice.begin()];
    double error = 0;
    MPI_Reduce(&errorLocal, &error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return error;
}


__global__ void fillFirstKindBoundaryKernel(double *u, const Block b, Axis axis, int i, int v1, int v2, int v1_size, int v2_size) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= v1_size * v2_size)
        return;

    int i1 = v1 + index / v2_size;
    int i2 = v2 + index % v2_size;

    switch (axis) {
        case X:
            u[ind(i, i1, i2, b)] = 0;
            break;
        case Y:
            u[ind(i1, i, i2, b)] = 0;
            break;
        case Z:
            u[ind(i1, i2, i, b)] = 0;
            break;
    }
}

__global__ void fillPeriodicBoundaryKernel(double *u, const Block b, const Grid g, double a, double t,
                                            Axis axis, int i, int v1, int v2, int v1_size, int v2_size) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= v1_size * v2_size)
        return;

    int i1 = v1 + index / v2_size;
    int i2 = v2 + index % v2_size;

    switch (axis) {
        case X:
            u[ind(i, i1, i2, b)] = analyticalSolution(i * g.h_x, i1 * g.h_y, i2 * g.h_z, t, a, g);
            break;
        case Y:
            u[ind(i1, i, i2, b)] = analyticalSolution(i1 * g.h_x, i * g.h_y, i2 * g.h_z, t, a, g);
            break;
        case Z:
            u[ind(i1, i2, i, b)] = analyticalSolution(i1 * g.h_x, i2 * g.h_y, i * g.h_z, t, a, g);
            break;
    }
}

void fillBoundaryValues(thrust::host_vector<double> &u, double t, const Block b, const Grid g, double a) {
    thrust::device_vector<double> uDevice(u);
    // Variant 3 -> first kind for x, periodic for y, first kind for z
    if (b.x_min == 0) {
        fillFirstKindBoundaryKernel<<<((b.y_size * b.z_size + num_threads - 1) / num_threads), num_threads>>>(
                thrust::raw_pointer_cast(&uDevice[0]), b, X, 0, b.y_min, b.z_min, b.y_size, b.z_size);
    }

    if (b.x_max == g.N) {
        fillFirstKindBoundaryKernel<<<((b.y_size * b.z_size + num_threads - 1) / num_threads), num_threads>>>(
                thrust::raw_pointer_cast(&uDevice[0]), b, X, g.N, b.y_min, b.z_min, b.y_size, b.z_size);
    }

    if (b.y_min == 0) {
        fillPeriodicBoundaryKernel<<<((b.x_size * b.z_size + num_threads - 1) / num_threads), num_threads>>>(
                thrust::raw_pointer_cast(&uDevice[0]), b, g, a, t, Y, 0, b.x_min, b.z_min, b.x_size, b.z_size);
    }

    if (b.y_max == g.N) {
        fillPeriodicBoundaryKernel<<<((b.x_size * b.z_size + num_threads - 1) / num_threads), num_threads>>>(
                thrust::raw_pointer_cast(&uDevice[0]), b, g, a, t, Y, g.N, b.x_min, b.z_min, b.x_size, b.z_size);
    }

    if (b.z_min == 0) {
        fillFirstKindBoundaryKernel<<<((b.x_size * b.y_size + num_threads - 1) / num_threads), num_threads>>>(
                thrust::raw_pointer_cast(&uDevice[0]), b, Z, 0, b.x_min, b.y_min, b.x_size, b.y_size);
    }

    if (b.z_max == g.N) {
        fillFirstKindBoundaryKernel<<<((b.x_size * b.y_size + num_threads - 1) / num_threads), num_threads>>>(
                thrust::raw_pointer_cast(&uDevice[0]), b, Z, g.N, b.x_min, b.y_min, b.x_size, b.y_size);
    }
    u = uDevice;
}


__global__ void initZeroLayerKernel(double *u0, double a, int size, int x1, int y1, int z1, int y_size, int z_size, const Grid g, const Block b) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= size)
        return;

    int i = x1 + index / (y_size * z_size);
    int j = y1 + index % (y_size * z_size) / z_size;
    int k = z1 + index % z_size;

    u0[ind(i, j, k, b)] = phi(i * g.h_x, j * g.h_y, k * g.h_z, a, g);
}

__global__ void initFirstLayerKernel(double *u0, double *u1, double *received,
                                     Block *send, Block *receive, int *numbers, int vsize,
                                     int size, int x1, int y1, int z1, int y_size, int z_size, const Grid g, const Block b) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= size)
        return;

    int i = x1 + index / (y_size * z_size);
    int j = y1 + index % (y_size * z_size) / z_size;
    int k = z1 + index % z_size;

    u1[ind(i, j, k, b)] = u0[ind(i, j, k, b)] +
            g.tau * g.tau / 2 * laplaceOperator(u0, i, j, k, b, g, received, send, receive, numbers, vsize);
}

void initValues(const Block b, const Grid g, double a,
                thrust::host_vector<double> &u0, thrust::host_vector<double> &u1,
                Block *send, Block *receive, int *numbers, int vsize) {
    // boundary (i = 0,N or j = 0,N or k = 0,N)
    fillBoundaryValues(u0, 0, b, g, a);
    fillBoundaryValues(u1, g.tau, b, g, a);

    // compute the boundaries of the current block
    int x1 = std::max(b.x_min, 1); int x2 = std::min(b.x_max, g.N - 1);
    int y1 = std::max(b.y_min, 1); int y2 = std::min(b.y_max, g.N - 1);
    int z1 = std::max(b.z_min, 1); int z2 = std::min(b.z_max, g.N - 1);
    // compute length of linearized data
    int x_size = x2 - x1 + 1;
    int y_size = y2 - y1 + 1;
    int z_size = z2 - z1 + 1;
    int size = x_size * y_size * z_size;

    thrust::device_vector<double> u0Device(u0);
    initZeroLayerKernel<<<(size + num_threads - 1) / num_threads, num_threads>>>(
            thrust::raw_pointer_cast(&u0Device[0]), a, size, x1, y1, z1, y_size, z_size, g, b);

    u0 = u0Device;
    thrust::device_vector<double> u1Device(u1);
    thrust::host_vector<double> received = exchange(u0, b, send, receive, numbers, vsize);
    thrust::device_vector<double> receivedDevice(received);

    initFirstLayerKernel<<<(size + num_threads - 1) / num_threads, num_threads>>>(
            thrust::raw_pointer_cast(&u0Device[0]), thrust::raw_pointer_cast(&u1Device[0]), thrust::raw_pointer_cast(&receivedDevice[0]),
            send, receive, numbers, vsize, size, x1, y1, z1, y_size, z_size, g, b);
    u1 = u1Device;
}

__global__ void getNextUKernel(double *u, double *u0, double *u1, double *received,
                               Block *send, Block *receive, int *numbers, int vsize,
                               int size, int x1, int y1, int z1, int y_size, int z_size, const Grid g, const Block b) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= size)
        return;

    int i = x1 + index / (y_size * z_size);
    int j = y1 + index % (y_size * z_size) / z_size;
    int k = z1 + index % z_size;

    u[ind(i, j, k, b)] = 2 * u1[ind(i, j, k, b)] - u0[ind(i, j, k, b)] +
            g.tau * g.tau * laplaceOperator(u1, i, j, k, b, g, received, send, receive, numbers, vsize);
}

void getNextU(int step, const Block b, std::vector< thrust::host_vector<double> > &u, const Grid g, double a,
              Block *send, Block *receive, int *numbers, int vsize) {
    // compute the boundaries of the current block
    int x1 = std::max(b.x_min, 1); int x2 = std::min(b.x_max, g.N - 1);
    int y1 = std::max(b.y_min, 1); int y2 = std::min(b.y_max, g.N - 1);
    int z1 = std::max(b.z_min, 1); int z2 = std::min(b.z_max, g.N - 1);
    // compute length of linearized data
    int x_size = x2 - x1 + 1;
    int y_size = y2 - y1 + 1;
    int z_size = z2 - z1 + 1;
    int size = x_size * y_size * z_size;

    thrust::host_vector<double> received = exchange(u[(step + 2) % 3], b, send, receive, numbers, vsize);
    thrust::device_vector<double> receivedDevice(received);
    thrust::device_vector<double> uDevice(u[step % 3]);
    thrust::device_vector<double> u0(u[(step + 1) % 3]);
    thrust::device_vector<double> u1(u[(step + 2) % 3]);

    getNextUKernel<<<(size + num_threads - 1) / num_threads, num_threads>>>(
            thrust::raw_pointer_cast(&uDevice[0]), thrust::raw_pointer_cast(&u0[0]), thrust::raw_pointer_cast(&u1[0]),
                    thrust::raw_pointer_cast(&receivedDevice[0]), send, receive, numbers, vsize, size, x1, y1, z1, y_size, z_size, g, b);
    u[step % 3] = uDevice;
    fillBoundaryValues(u[step % 3], step * g.tau, b, g, a);
}

inline bool isInside(int xmin1, int xmax1, int ymin1, int ymax1, int xmin2, int xmax2, int ymin2, int ymax2) {
    return xmin2 <= xmin1 && xmax1 <= xmax2 && ymin2 <= ymin1 && ymax1 <= ymax2;
}

void getNeighbours(const std::vector<Block> &blocks, thrust::host_vector<Block> &send, thrust::host_vector<Block> &receive,
                   thrust::host_vector<int> &numbers, int proc_rank, int proc_size) {
    Block block = blocks[proc_rank];

    for (int i = 0; i < proc_size; i++) {
        if (i == proc_rank)
            continue;

        Block otherBlock = blocks[i];
        if (block.x_min == otherBlock.x_max + 1 or otherBlock.x_min == block.x_max + 1) {
            int xSend = block.x_min == otherBlock.x_max + 1 ? block.x_min : block.x_max;
            int xRecv = otherBlock.x_min == block.x_max + 1 ? otherBlock.x_min : otherBlock.x_max;
            int y_min, y_max, z_min, z_max;
            if (isInside(block.y_min, block.y_max, block.z_min, block.z_max,
                         otherBlock.y_min, otherBlock.y_max, otherBlock.z_min, otherBlock.z_max)) {
                y_min = block.y_min; y_max = block.y_max; z_min = block.z_min; z_max = block.z_max;
            } else if (isInside(otherBlock.y_min, otherBlock.y_max, otherBlock.z_min, otherBlock.z_max,
                                block.y_min, block.y_max, block.z_min, block.z_max)) {
                y_min = otherBlock.y_min; y_max = otherBlock.y_max; z_min = otherBlock.z_min; z_max = otherBlock.z_max;
            } else {
                continue;
            }
            // add block as a rectangle (it's a border between two processes)
            send.push_back(Block(xSend, xSend, y_min, y_max, z_min, z_max));
            receive.push_back(Block(xRecv, xRecv, y_min, y_max, z_min, z_max));
            numbers.push_back(i);
            continue;
        }
        if (block.y_min == otherBlock.y_max + 1 or otherBlock.y_min == block.y_max + 1) {
            int ySend = block.y_min == otherBlock.y_max + 1 ? block.y_min : block.y_max;
            int yRecv = otherBlock.y_min == block.y_max + 1 ? otherBlock.y_min : otherBlock.y_max;
            int x_min, x_max, z_min, z_max;
            if (isInside(block.x_min, block.x_max, block.z_min, block.z_max,
                         otherBlock.x_min, otherBlock.x_max, otherBlock.z_min, otherBlock.z_max)) {
                x_min = block.x_min; x_max = block.x_max; z_min = block.z_min; z_max = block.z_max;
            } else if (isInside(otherBlock.x_min, otherBlock.x_max, otherBlock.z_min, otherBlock.z_max,
                                block.x_min, block.x_max, block.z_min, block.z_max)) {
                x_min = otherBlock.x_min; x_max = otherBlock.x_max; z_min = otherBlock.z_min; z_max = otherBlock.z_max;
            } else {
                continue;
            }
            send.push_back(Block(x_min, x_max, ySend, ySend, z_min, z_max));
            receive.push_back(Block(x_min, x_max, yRecv, yRecv, z_min, z_max));
            numbers.push_back(i);
            continue;
        }
        if (block.z_min == otherBlock.z_max + 1 or otherBlock.z_min == block.z_max + 1) {
            int zSend = block.z_min == otherBlock.z_max + 1 ? block.z_min : block.z_max;
            int zRecv = otherBlock.z_min == block.z_max + 1 ? otherBlock.z_min : otherBlock.z_max;
            int x_min, x_max, y_min, y_max;
            if (isInside(block.x_min, block.x_max, block.y_min, block.y_max,
                         otherBlock.x_min, otherBlock.x_max, otherBlock.y_min, otherBlock.y_max)) {
                x_min = block.x_min; x_max = block.x_max; y_min = block.y_min; y_max = block.y_max;
            } else if (isInside(otherBlock.x_min, otherBlock.x_max, otherBlock.y_min, otherBlock.y_max,
                                block.x_min, block.x_max, block.y_min, block.y_max)) {
                x_min = otherBlock.x_min; x_max = otherBlock.x_max; y_min = otherBlock.y_min; y_max = otherBlock.y_max;
            } else {
                continue;
            }
            send.push_back(Block(x_min, x_max, y_min, y_max, zSend, zSend));
            receive.push_back(Block(x_min, x_max, y_min, y_max, zRecv, zRecv));
            numbers.push_back(i);
            continue;
        }
    }
}

double solve(int steps, Grid g) {
    double a = M_PI * sqrt(1.0 / (g.L_x * g.L_x) + 4.0 / (g.L_y * g.L_y) + 9.0 / (g.L_z * g.L_z));
    int proc_rank, proc_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);

    // split grid between processes
    std::vector<Block> blocks;
    split(0, g.N, 0, g.N, 0, g.N, proc_size, X, blocks);
    Block block = blocks[proc_rank];

    // allocate space for u
    std::vector< thrust::host_vector<double> > u(3);
    for (int i = 0; i < 3; i++)
        u[i].resize(block.size);

    // fill blocksToSend and blocksToReceive vectors
    thrust::host_vector<Block> blocksToSend, blocksToReceive;
    thrust::host_vector<int> neighboursRanks;
    getNeighbours(blocks, blocksToSend, blocksToReceive, neighboursRanks, proc_rank, proc_size);
    thrust::device_vector<Block> blocksToSendDevice(blocksToSend);
    thrust::device_vector<Block> blocksToReceiveDevice(blocksToReceive);
    thrust::device_vector<int> neighboursRanksDevice(neighboursRanks);

    // init u_0 and u_1
    initValues(block, g, a, u[0], u[1], thrust::raw_pointer_cast(&blocksToSendDevice[0]),
               thrust::raw_pointer_cast(&blocksToReceiveDevice[0]),
               thrust::raw_pointer_cast(&neighboursRanksDevice[0]), neighboursRanksDevice.size());

    // calculate the next time layers for u
    for (int step = 2; step <= steps; step++) {
        getNextU(step, block, u, g, a, thrust::raw_pointer_cast(&blocksToSendDevice[0]),
                 thrust::raw_pointer_cast(&blocksToReceiveDevice[0]),
                 thrust::raw_pointer_cast(&neighboursRanksDevice[0]), neighboursRanksDevice.size());
    }

    return computeLayerError(u[steps % 3], steps * g.tau, block, g, a);
}


#endif //TASK4_SOLVER_H
