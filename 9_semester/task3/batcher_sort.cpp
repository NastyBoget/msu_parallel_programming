#include <cmath>
#include "data_structures.h"


bool comp(Point i, Point j) {
    if (i.coord[coordSorted] < j.coord[coordSorted]) {
        return true;
    } else {
        if (i.coord[coordSorted] > j.coord[coordSorted]) {
            return false;
        }
    }
    return i.coord[1 - coordSorted] < j.coord[1 - coordSorted];
}

void merge(vector<Point>& array1, const vector<Point>& array2, int proc1, int proc2, int procRank) {
    int arraySize = array1.size();
    vector<Point> tmpArray(2 * arraySize);
    vector<Point> left, right;
    if (proc1 == procRank) {
        left = array1;
        right = array2;
    } else {
        left = array2;
        right = array1;
    }
    for (int i = 0, l = 0, r = 0; i < 2 * arraySize;i++) {
        if (l >= arraySize) {
            tmpArray[i] = right[r++];
        } else {
            if (r >= arraySize) {
                tmpArray[i] = left[l++];
            } else {
                if (comp(left[l], right[r])) {
                    tmpArray[i] = left[l++];
                } else {
                    tmpArray[i] = right[r++];
                }
            }
        }
    }

    if (procRank == proc1) {
        array1.assign(tmpArray.begin(), tmpArray.begin() + arraySize);
    } else {
        array1.assign(tmpArray.begin() + arraySize, tmpArray.end());
    }
}

void compareExchangeParallel(vector<Point>& currentArray, int proc1, int proc2, MPI_Comm comm, int procRank) {
    int arraySize = currentArray.size();
    int otherProc;
    MPI_Request request;
    MPI_Status status;
    vector<Point> otherArray(arraySize);

    if (procRank == proc1) {
        otherProc = proc2;
    } else {
        otherProc = proc1;
    }

    MPI_Isend(currentArray.data(), arraySize, messageType, otherProc, procRank, comm, &request);
    MPI_Recv(otherArray.data(), arraySize, messageType, otherProc, otherProc, comm, &status);
    MPI_Wait(&request, &status);

    merge(currentArray, otherArray, proc1, proc2, procRank);
}

void bSortParallel(vector<Point>& pointsArrayPart, MPI_Comm comm) {
    int procRank, procSize;
    MPI_Comm_rank(comm, &procRank);
    MPI_Comm_size(comm, &procSize);
    int N = procSize;
    int t = int(log2(N));
    int q, r, d;
    int const_q = pow((float) 2, t);
    for (int p = const_q; p > 0; p >>= 1) {
        r = 0;
        d = p;
        q = const_q;
        bool check;
        do {
            check = false;
            for (int i = 0; i < N - d; i++) {
                if ((i & p) == r and (procRank == i or procRank == i + d)) {
                    compareExchangeParallel(pointsArrayPart, i, i + d, comm, procRank);
                }
            }
            if (q != p) {
                d = q - p;
                q >>= 1;
                r = p;
                check = true;
            }
        } while (check);
    }
}

void runSortParallel(vector<Point>& partArray, MPI_Comm comm) {
    sort(partArray.begin(), partArray.end(), comp);
    bSortParallel(partArray, comm);
    MPI_Barrier(comm);
}
