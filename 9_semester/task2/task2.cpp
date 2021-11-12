#include <iostream>
#include <limits>
#include <vector>
#include <functional>
#include <cmath>
#include "mpi.h"

using namespace std;

int procSize, procRank;
bool verbose = false;
int coordSorted = 0;
MPI_Datatype messageType;

struct Point {
    float coord[2];
    int index;
};

float x(int i, int j) {
    return (sin(i) + cos(j)) * (i + j);
}

float y(int i, int j) {
    return (cos(i) + sin(j)) * (i + j);
}

bool comp(Point i, Point j) {
    return i.coord[coordSorted] < j.coord[coordSorted];
}

vector<Point> generatePoints(int n1, int n2, int startIndex, int arraySize) {
    vector<Point> result;
    result.reserve(arraySize);
    for(int index = startIndex; index < startIndex + arraySize; index++) {
        Point newPoint;
        int i = index / n2;
        int j = index % n2;
        newPoint.coord[0] = x(i, j);
        newPoint.coord[1] = y(i, j);
        newPoint.index = index;
        result.push_back(newPoint);
    }
    return result;
}

void print(const vector<Point>& pointsArray, const string& msg) {
	if (!verbose)
		return;
    int rank = 0;
    while (rank < procSize) {
        if (procRank == rank) {
            cout << "Processor #" << procRank << ". " << msg << endl;
            for(int i = 0; i < pointsArray.size(); i++) {
                cout << pointsArray[i].coord[coordSorted] << " ";
            }
            cout << endl;
            fflush(stdout);
        }
        rank ++;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    cout << endl;
}

bool checkCorrectness(const vector<Point>& pointsArray) {
    for(int i = 0; i < pointsArray.size() - 1; i++) {
        if (pointsArray[i].coord[coordSorted] > pointsArray[i + 1].coord[coordSorted]) {
            return false;
        }
    }
    return true;
}

inline void compareExchange(vector<Point>& pointsArray, int i, int j) {
    if (pointsArray[i].coord[coordSorted] > pointsArray[j].coord[coordSorted]) {
        std::swap(pointsArray[i].coord[coordSorted], pointsArray[j].coord[coordSorted]);
    }
}

void buildDerivedType(vector<Point>& indata, MPI_Datatype* message_type_ptr)
{
    int block_lengths[2];
    MPI_Aint displacements[2];
    MPI_Aint addresses[3];
    MPI_Datatype typelist[2];

    typelist[0] = MPI_FLOAT;
    typelist[1] = MPI_INT;

    block_lengths[0] = 2;
    block_lengths[1] = 1;

    MPI_Get_address(&(indata[0]), &addresses[0]);
    MPI_Get_address(&(indata[0].coord), &addresses[1]);
    MPI_Get_address(&(indata[0].index), &addresses[2]);

    displacements[0] = addresses[1] - addresses[0];
    displacements[1] = addresses[2] - addresses[0];

    MPI_Type_create_struct(2, block_lengths, displacements, typelist, message_type_ptr);
    MPI_Type_commit(message_type_ptr);
}

void oddEvenMergeSort(vector<Point>& left, const vector<Point>& right, int proc1, int proc2) {
    int arraySize = left.size();
    int tmpSize = 2 * arraySize;
    vector<Point> tmpArray(tmpSize);

    for(size_t idx = 0; idx < 2; idx++) {
		for(size_t i = idx, l = idx, r = idx; i < tmpSize; i += 2) {
			if (l < arraySize and (r >= arraySize or left[l].coord[coordSorted] <= right[r].coord[coordSorted])) {
				tmpArray[i] = left[l];
				l += 2;
			} else {
				tmpArray[i] = right[r];
				r += 2;
			}
		}
	}
	if (arraySize % 2 == 1) {
	    if (left[arraySize - 1].coord[coordSorted] > right[arraySize - 1].coord[coordSorted]) {
            tmpArray[tmpSize - 1] = left[arraySize - 1];
	    } else {
            tmpArray[tmpSize - 1] = right[arraySize - 1];
	    }
	}

	for(int i = 1; i < tmpSize - 1; i += 2) {
		compareExchange(tmpArray, i, i + 1);
	}
    if (procRank == proc1) {
        std::copy(tmpArray.begin(), tmpArray.begin() + tmpSize / 2, left.begin());
    } else if (procRank == proc2) {
        std::copy(tmpArray.begin() + tmpSize / 2, tmpArray.end(), left.begin());
    }
}

void compareExchangeParallel(vector<Point>& currentArray, int proc1, int proc2) {
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
	
	MPI_Isend(currentArray.data(), arraySize, messageType, otherProc, procRank, MPI_COMM_WORLD, &request);
	MPI_Recv(otherArray.data(), arraySize, messageType, otherProc, otherProc, MPI_COMM_WORLD, &status);
	MPI_Wait(&request, &status);

	oddEvenMergeSort(currentArray, otherArray, proc1, proc2);
}

void bSortParallel(vector<Point>& pointsArrayPart) {
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
                    compareExchangeParallel(pointsArrayPart, i, i + d);
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

void runSort(vector<Point>& pointsArray) {
	double startTime = MPI_Wtime();
	sort(pointsArray.begin(), pointsArray.end(), comp);
	double endTime = MPI_Wtime();
	if (!checkCorrectness(pointsArray)) {
		cout << "Serial: incorrect results" << endl;
	} 	
	print(pointsArray, "Serial: array after sorting: ");
	cout << "Serial time = " << endTime - startTime << endl;
}

void runSortParallel(int n1, int n2) {
    int partArraySize = 0;
    for(; partArraySize * procSize < n1 * n2; partArraySize++) {}

    vector<Point> partArray = generatePoints(n1, n2, procRank * partArraySize, partArraySize);
    print(partArray, "Before sorting: ");
    buildDerivedType(partArray, &messageType);
	double startTime = 0.0, endTime;

	MPI_Barrier(MPI_COMM_WORLD);
    if (procRank == 0) {
        startTime = MPI_Wtime();
    }

    sort(partArray.begin(), partArray.end(), comp);
    bSortParallel(partArray);

    MPI_Barrier(MPI_COMM_WORLD);
	if (procRank == 0) {
		endTime = MPI_Wtime();
		cout << "Parallel time = " << endTime - startTime << endl;
	}
    if (!checkCorrectness(partArray)) {
        cout << "Parallel: incorrect results" << endl;
    }
    print(partArray, "After sorting: ");
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);

    if (procRank == 0) {
        cout << "n1 = " << n1 << ", n2 = " << n2 << ", proc number = " << procSize << endl;
    }

	if (procSize == 1) {
        vector<Point> pointsArray = generatePoints(n1, n2, 0, n1 * n2);
	    runSort(pointsArray);
	} else {
		runSortParallel(n1, n2);
	}
	MPI_Finalize();
	return 0;
}
