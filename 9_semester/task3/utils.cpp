#include <cmath>
#include <cfloat>
#include <fstream>
#include <iostream>
#include "data_structures.h"


float x(int i, int j) {
    return 10 * i;
}

float y(int i, int j) {
    return 10 * j;
}

vector<Point> generatePoints(int n1, int n2, int procRank, int partArraySize, int procSize) {
    int actualArraySize = n1 * n2;
    int procNumberWithoutFElem = actualArraySize % procSize;
    int minPartArraySize = actualArraySize / procSize;
    // fictive elements are distributed uniformly between processes in the end
    int nonFictiveElemNumber = procRank < procNumberWithoutFElem ? partArraySize : minPartArraySize;
    int startIndex = procRank < procNumberWithoutFElem ? procRank * partArraySize :
            procNumberWithoutFElem * partArraySize + (procRank - procNumberWithoutFElem) * nonFictiveElemNumber;

    vector<Point> result;
    result.reserve(partArraySize);
    for(int index = startIndex;index < startIndex + partArraySize - procNumberWithoutFElem; index++) {
        int i = index / n2;
        int j = index % n2;
        Point newPoint;
        newPoint.coord[0] = x(i, j);
        newPoint.coord[1] = y(i, j);
        newPoint.index = index;
        result.push_back(newPoint);
    }
    if (procNumberWithoutFElem != 0) {
        Point newPoint;
        newPoint.coord[0] = FLT_MAX;
        newPoint.coord[1] = FLT_MAX;
        newPoint.index = -1;
        result.push_back(newPoint);
    }
    return result;
}

void buildDerivedType(MPI_Datatype* message_type_ptr) {
    int block_lengths[2];
    MPI_Datatype typelist[2];
    MPI_Aint displacements[2];

    typelist[0] = MPI_FLOAT;
    typelist[1] = MPI_INT;

    block_lengths[0] = 2;
    block_lengths[1] = 1;

    displacements[0] = offsetof(Point, coord);
    displacements[1] = offsetof(Point, index);

    MPI_Type_create_struct(2, block_lengths, displacements, typelist, message_type_ptr);
    MPI_Type_commit(message_type_ptr);
}

void printResults(const vector<Point> &partArray, const vector<int> &domains, const string &fileName, int n1, int n2) {
    int procRank, procSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    int rank = 0;
    std::ofstream out;
    while (rank < procSize) {
        if (procRank == rank) {
            if (procRank == 0) {
                out.open(fileName);
            } else {
                out.open(fileName, ios::app);
            }
            for(int i = 0; i < partArray.size(); i++) {
                if(partArray[i].index == -1)
                    continue;
                out << partArray[i].index / n2 << " " << partArray[i].index % n2 << " ";
                out << partArray[i].coord[0] << " " << partArray[i].coord[1] << " " << domains[i] << endl;
            }
            out.close();
        }
        rank++;
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void print(const vector<Point>& pointsArray) {
    int procRank, procSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    int rank = 0;
    while (rank < procSize) {
        if (procRank == rank) {
            cout << "Processor #" << procRank << endl;
            for(int i = 0; i < pointsArray.size(); i++) {
                cout << "(" << pointsArray[i].coord[0] << "," << pointsArray[i].coord[1] << ") ";
            }
            cout << endl;
        }
        rank++;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    cout << endl;
}