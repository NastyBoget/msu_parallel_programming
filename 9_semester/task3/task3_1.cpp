#include <iostream>
#include <cfloat>
#include <cmath>
#include "data_structures.h"

#define LOOPS 1

int coordSorted = 0;
int width; // n2
int edgesPerProcess = 0;
MPI_Datatype messageType;


void init_values(int k, int n, int& k1, int& k2, int& n1, int& n2) {
    k1 = (k + 1) / 2;
    k2 = k - k1;
    n1 = n * (k1 / (double) k);
    n2 = n - n1;
}

int getGroup(int procRank, int delimiterProcRank) {
    int color;
    if (delimiterProcRank == 0) {
        color = procRank > delimiterProcRank ? 0 : 1;
    } else {
        color = procRank >= delimiterProcRank ? 0 : 1;
    }
    return color;
}

vector<Point> getElemsToSend(vector<Point>::iterator arrayIt, int elemsNumber, int procSize) {
    int elemToSendSize = ceil(elemsNumber / (double) procSize);
    int newSize = elemToSendSize * procSize;
    vector<Point> newArray(newSize);

    for(int i = 0; i < elemsNumber; i++, arrayIt++){
        newArray[i] = *arrayIt;
        arrayIt->index = -1;
        arrayIt->coord[0] = FLT_MAX;
        arrayIt->coord[1] = FLT_MAX;
    }

    for(int i = elemsNumber; i < elemToSendSize * procSize; i++){
        newArray[i].index = -1;
        newArray[i].coord[0] = FLT_MAX;
        newArray[i].coord[1] = FLT_MAX;
    }
    return newArray;
}

void elemRedistribution(MPI_Comm comm, int procRank, int procSize, int delimiterProcRank,
                        int delimiterPartArray, vector<Point>& array, int& elemsRightAdded) {
    MPI_Status status;
    int partArraySize = array.size();
    if (delimiterProcRank == 0) {
        int elemToSendSize = ceil((partArraySize - delimiterPartArray) / (double) (procSize - 1));
        if (procRank == delimiterProcRank) {
            vector<Point> elemsToSend = getElemsToSend(array.begin() + delimiterPartArray, partArraySize - delimiterPartArray,procSize - 1);
            elemsRightAdded = elemToSendSize;
            for(int i = 1, j = 0; i < procSize; i++, j++)
                MPI_Send(elemsToSend.data() + j * elemToSendSize, elemToSendSize, messageType, i, i, comm);
        } else {
            array.resize(partArraySize + elemToSendSize);
            MPI_Recv(array.data() + partArraySize, elemToSendSize, messageType, delimiterProcRank, procRank, comm, &status);
        }
        return;
    }
    // delimiterProcRank > 0
    if (procRank <= delimiterProcRank) {
        int elemToSendSize = ceil(delimiterPartArray / (double) delimiterProcRank);
        if (procRank == delimiterProcRank) {
            vector<Point> elemsToSend = getElemsToSend(array.begin(), delimiterPartArray, delimiterProcRank);
            for(int i = 0; i < delimiterProcRank; i++)
                MPI_Send(elemsToSend.data() + i * elemToSendSize, elemToSendSize, messageType, i, i, comm);
        } else {
            array.resize(partArraySize + elemToSendSize);
            MPI_Recv(array.data() + partArraySize, elemToSendSize, messageType, delimiterProcRank, procRank, comm, &status);
        }
    }
}

void findBisectCoord(vector<Point> &array, int arrayStart, int n, int n1) {
    sort(array.begin() + arrayStart, array.begin() + arrayStart + n, comp);
    coordSorted = 1 - coordSorted;
}

void findBisectCoordParallel(vector<Point> &array,  MPI_Comm comm, int delimiterProcRank, int delimiterPartArray) {
    int procRank, procSize;
    MPI_Comm_rank(comm, &procRank);
    MPI_Comm_size(comm, &procSize);
    int elemsAdded = 0;
    runSortParallel(array, comm);
    coordSorted = 1 - coordSorted;
    elemsAdded = 0;
    elemRedistribution(comm, procRank, procSize, delimiterProcRank, delimiterPartArray, array, elemsAdded);
}

void localBisect(vector<Point> &array, vector<int>& domains, int domainStart, int k, int arrayStart, int n) {
    if (k == 1) {
        for(int i = 0; i < n; i++) {
            domains[arrayStart + i] = domainStart;
        }
        return;
    }
    int k1, k2, n1, n2;
    init_values(k, n, k1, k2, n1, n2);
    findBisectCoord(array, arrayStart, n, n1);

    localBisect(array, domains, domainStart, k1, arrayStart, n1);
    localBisect(array, domains, domainStart + k1, k2, arrayStart + n1, n2);
}

void removeFictive(vector<Point>& array) {
    vector<Point> newArray;
    for(int i = 0; i < array.size(); i++) {
        if (array[i].index == -1)
            continue;
        newArray.push_back(array[i]);
    }
    array = newArray;
}

void recursiveBisect(vector<Point>& array, vector<int>& domains, int domainStartValue, int k, int n, MPI_Comm comm) {
    int procRank, procSize;
    MPI_Comm_rank(comm, &procRank);
    MPI_Comm_size(comm, &procSize);

    // recursion base
    if (procSize == 1) {
        removeFictive(array);
        int arraySize = array.size();
        domains.resize(arraySize);
        localBisect(array, domains, domainStartValue, k, 0, arraySize);
        return;
    }
    if (k == 1) {
        removeFictive(array);
        int arraySize = array.size();
        domains.resize(arraySize);
        for(int i = 0; i < arraySize; i++)
            domains[i] = domainStartValue;
        return;
    }

    int k1, k2, n1, n2;
    init_values(k, n, k1, k2, n1, n2);
    int delimiterProcRank = n1 / array.size();
    int delimiterPartArray = n1 % array.size();
    findBisectCoordParallel(array, comm, delimiterProcRank, delimiterPartArray);
    MPI_Comm newComm;
    int group = getGroup(procRank, delimiterProcRank);
    MPI_Comm_split(comm, group, procRank, &newComm);
    if ((delimiterProcRank == 0 and procRank == 0) or procRank < delimiterProcRank) {
        recursiveBisect(array, domains, domainStartValue, k1, n1, newComm);
    } else {
        recursiveBisect(array, domains, domainStartValue + k1, k2, n2, newComm);
    }
}


int main(int argc, char **argv) {
    int procRank, procSize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    buildDerivedType(&messageType);

    int k = atoi(argv[1]);
    int n1 = atoi(argv[2]);
    int n2 = atoi(argv[3]);
    int netNum = atoi(argv[4]);
    char* filename = argv[5];
    width = n2;
    int partArraySize = ceil(n1 * n2 / (double) procSize);
    double sumTime = 0;
    vector<Point> array;
    vector<int> domains;

    for(int loop = 0; loop < LOOPS; loop++) {
        array.clear();
        domains.clear();
        edgesPerProcess = 0;
        array = generatePoints(n1, n2, procRank, partArraySize, procSize, netNum);
        MPI_Barrier(MPI_COMM_WORLD);
        double startTime, endTime, maxTime, delta;
        startTime = MPI_Wtime();
        recursiveBisect(array, domains, 0, k, n1 * n2, MPI_COMM_WORLD);
        endTime = MPI_Wtime();
        delta = endTime - startTime;
        MPI_Reduce(&delta, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        sumTime += maxTime;
    }
    int edgesSum;
    MPI_Reduce(&edgesPerProcess, &edgesSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (procRank == 0) {
        cout << sumTime / LOOPS << endl;
        cout << edgesSum << endl;
    }
    printResults(array, domains, filename, n1, n2);
    MPI_Finalize();
    return 0;
}
