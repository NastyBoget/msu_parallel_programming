#include <iostream>
#include <cfloat>
#include <cmath>
#include "data_structures.h"


int coordSorted = 0;
MPI_Datatype messageType;

void init_values(int k, int n, int& k1, int& k2, int& n1, int& n2) {
    k1 = (k + 1) / 2;
    k2 = k - k1;
    n1 = n * (k1 / (double) k);
    n2 = n - n1;
}

void localBisect(vector<Point> &a, vector<int>& domains, int domainStart, int k, int arrayStart, int n) {
    if (k == 1) {
        for(int i = 0; i < n; i++) {
            domains[arrayStart + i] = domainStart;
        }
        return;
    }

    coordSorted = coordSorted == 0 ? 1 : 0;  // TODO
    sort(a.begin() + arrayStart, a.begin() + n, comp);

    int k1, k2, n1, n2;
    init_values(k, n, k1, k2, n1, n2);

    localBisect(a, domains, domainStart, k1, arrayStart, n1);
    localBisect(a, domains, domainStart + k1, k2, arrayStart + n1, n2);
}

void removeFictive(vector<Point>& array) {
    vector<Point> newArray;
    for(auto & i : array) {
        if (i.index == -1)
            continue;
        newArray.push_back(i);
    }
    array = newArray;
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

void recursiveBisect(vector<Point>& array, vector<int>& domains, int domainStartValue, int k, int n,
                     int& partArraySize, MPI_Comm comm) {
    int procRank, procSize;
    MPI_Comm_rank(comm, &procRank);
    MPI_Comm_size(comm, &procSize);

    // recursion base
    if (procSize == 1) {
        removeFictive(array);
        int arraySize = array.size();
        domains.resize(arraySize);
        localBisect(array, domains, domainStartValue, k, 0, arraySize);
        partArraySize = arraySize;
        return;
    }
    if (k == 1) {
        removeFictive(array);
        int arraySize = array.size();
        domains.resize(arraySize);
        for(int i = 0; i < arraySize; i++)
            domains[i] = domainStartValue;
        partArraySize = arraySize;
        return;
    }

    MPI_Status status;
    int k1, k2, n1, n2;
    init_values(k, n, k1, k2, n1, n2);
    int delimiterProcRank = n1 / partArraySize;
    int delimiterPartArray = n1 % partArraySize;
    int color;

    if (delimiterProcRank == 0) {
        color = procRank > delimiterProcRank ? 0 : 1;
    } else {
        color = procRank >= delimiterProcRank ? 0 : 1;
    }
    coordSorted = coordSorted == 0 ? 1 : 0;  // TODO
    runSortParallel(array, comm);
    MPI_Comm newComm;
    MPI_Comm_split(comm, color, procRank, &newComm);

    if (delimiterProcRank == 0) {
        int elemToSendSize = ceil((partArraySize - delimiterPartArray) / (double) procSize);
        if (procRank == delimiterProcRank) {
            vector<Point> elemsToSend = getElemsToSend(array.begin() + delimiterPartArray,
                    partArraySize - delimiterPartArray,procSize - delimiterProcRank - 1);
            for(int i = delimiterProcRank + 1, j = 0; i < procSize; i++, j++)
                MPI_Send(elemsToSend.data() + j * elemToSendSize, elemToSendSize, messageType, i, 0, comm);
            recursiveBisect(array, domains, domainStartValue, k1, n1, partArraySize, newComm);
        } else {
            vector<Point> extendedArray(array);
            extendedArray.resize(partArraySize + elemToSendSize);
            MPI_Recv(extendedArray.data() + partArraySize, elemToSendSize, messageType, delimiterProcRank, 0,
                     comm, &status);
            array = extendedArray;
            partArraySize += elemToSendSize;
            recursiveBisect(array, domains, domainStartValue + k1, k2, n2, partArraySize,
                            newComm);
        }
        return;
    }
    // delimiterProcRank > 0
    if (procRank <= delimiterProcRank) {
        int elemToSendSize = ceil(delimiterPartArray / (double) delimiterProcRank);
        if (procRank == delimiterProcRank) {
            vector<Point> elemsToSend = getElemsToSend(array.begin(), delimiterPartArray, delimiterProcRank);
            for(int i = 0; i < delimiterProcRank; i++)
                MPI_Send(elemsToSend.data() + i * elemToSendSize, elemToSendSize, messageType, i, 0, comm);
        } else {
            vector<Point> extendedArray(array);
            extendedArray.resize(partArraySize + elemToSendSize);
            MPI_Recv(extendedArray.data() + partArraySize, elemToSendSize, messageType, delimiterProcRank, 0,
                     comm, &status);
            array = extendedArray;
            partArraySize += elemToSendSize;
        }
    }
    if (procRank < delimiterProcRank) {
        recursiveBisect(array, domains, domainStartValue, k1, n1, partArraySize, newComm);
    } else {
        recursiveBisect(array, domains, domainStartValue + k1, k2, n2, partArraySize, newComm);
    }
}

int main(int argc, char **argv)
{
    int procRank, procSize;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procSize);
    buildDerivedType(&messageType);

    int k = atoi(argv[1]);
    int n1 = atoi(argv[2]);
    int n2 = atoi(argv[3]);
    int partArraySize = ceil(n1 * n2 / (double) procSize);

    vector<Point> array = generatePoints(n1, n2, procRank, partArraySize, procSize);
    vector<int> domains;

    double startTime, endTime, maxTime, delta;
    startTime = MPI_Wtime();
    recursiveBisect(array, domains, 0, k, n1 * n2, partArraySize, MPI_COMM_WORLD);
    endTime = MPI_Wtime();
    delta = endTime - startTime;
    MPI_Reduce(&delta, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (procRank == 0) {
        cout << maxTime << endl;
    }
    printResults(array, domains, "out.txt", n1, n2);
    // TODO count wedges
    
    MPI_Finalize();
    return 0;
}
