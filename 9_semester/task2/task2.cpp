#include <iostream>
#include <limits>
#include <random>
#include <functional>
#include <cmath>
#include "mpi.h"

using namespace std;

int procSize, procRank;
double serialTime, parallelTime;
bool verbose = false;
int coordSorted = 0;
MPI_Datatype messageType;

struct Point {
    float coord[2];
    int index;
};

float x(int i, int j) {
    static std::default_random_engine e;
    static std::uniform_real_distribution<> dis(-1000, 1000);
    return dis(e);
}

float y(int i, int j) {
    return x(i, j);
}

Point* generatePoints(int n1, int n2) {
    auto result = new Point[n1 * n2];
    for(int i = 0; i < n1; i++)
        for(int j = 0; j < n2; j++) {
            result[i * n2 + j].coord[0] = x(i, j);
            result[i * n2 + j].coord[1] = y(i, j);
            result[i * n2 + j].index = i * n2 + j;
        }
    return result;
}

void print(const Point* pointsArray, size_t arraySize, const string& msg) {
	if (!verbose)
		return;
	cout << msg << endl;
    for(size_t i = 0; i < arraySize; i++) {
		cout << pointsArray[i].coord[coordSorted] << " ";
	}
	cout << endl;
}

bool checkCorrectness(const Point* pointsArray, size_t arraySize) {
    for(size_t i = 0; i < arraySize - 1; i++) {
        if (pointsArray[i].coord[coordSorted] > pointsArray[i + 1].coord[coordSorted]) {
            return false;
        }
    }
    return true;
}

inline void compareExchange(Point* pointsArray, int i, int j) {
    if (pointsArray[i].coord[coordSorted] > pointsArray[j].coord[coordSorted]) {
        std::swap(pointsArray[i].coord[coordSorted], pointsArray[j].coord[coordSorted]);
    }
}

void buildDerivedType(Point* indata, MPI_Datatype* message_type_ptr)
{
    int block_lengths[2];
    MPI_Aint displacements[2];
    MPI_Aint addresses[3];
    MPI_Datatype typelist[2];

    typelist[0] = MPI_FLOAT;
    typelist[1] = MPI_INT;

    block_lengths[0] = 2;
    block_lengths[1] = 1;

    MPI_Get_address(indata, &addresses[0]);
    MPI_Get_address(&(indata->coord), &addresses[1]);
    MPI_Get_address(&(indata->index), &addresses[2]);

    displacements[0] = addresses[1] - addresses[0];
    displacements[1] = addresses[2] - addresses[0];

    MPI_Type_create_struct(2, block_lengths, displacements, typelist, message_type_ptr);
    MPI_Type_commit(message_type_ptr);
}

void oddEvenMergeSort(Point* left, const Point* right, int arraySize, int proc1, int proc2) {
    auto tmpSize = 2 * arraySize;
	auto tmpArray = new Point[tmpSize];
	
	for(size_t idx = 0; idx < 2; idx += 1) {
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
	    for (int i = 0; i < arraySize; i++) {
	        left[i] = tmpArray[i];
	    }
	} else if (procRank == proc2) {
        for (int i = 0; i < arraySize; i++) {
            left[i] = tmpArray[i + arraySize];
        }
    }
	delete[] tmpArray;
}

void compareExchangeParallel(Point* currentArray, int arraySize, int proc1, int proc2) {
	int otherProc;
	MPI_Request request;
	MPI_Status status;
	auto otherArray = new Point[arraySize];

	if (procRank == proc1) {
		otherProc = proc2;
	} else {
		otherProc = proc1;
	}
	
	MPI_Isend(currentArray, arraySize, messageType, otherProc, procRank, MPI_COMM_WORLD, &request);
	MPI_Recv(otherArray, arraySize, messageType, otherProc, otherProc, MPI_COMM_WORLD, &status);
	MPI_Wait(&request, &status);

	oddEvenMergeSort(currentArray, otherArray, arraySize, proc1, proc2);
	delete[] otherArray;
}

void bSortParallel(Point* pointsArrayPart, int arraySize) {
    int N = procSize;
    int t = int(log2(N));
    int q, r, d;
    auto const_q = pow(2, t);
    for (int p = pow(2, t); p > 0; p >>= 1) {
        r = 0;
        d = p;
        q = const_q;
        bool check;
        do {
        	check = false;
            for (int i = 0; i < N - d; i += 1) {
            	if ((i & p) == r and (procRank == i or procRank == i + d)) {
                    compareExchangeParallel(pointsArrayPart, arraySize, i, i + d);
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

void bSort(Point* pointsArray, int arraySize) { // TODO use more effective algorithm
    int N = arraySize;
    int t = int(log2(N));
    int q, r, d;
    auto const_q = pow(2, t);
    for (int p = pow(2, t); p > 0; p >>= 1) {
        r = 0;
        d = p;
        q = const_q;
        bool check;
        do {
        	check = false;
            for (int i = 0; i < N - d; i += 1) {
            	if ((i & p) == r) {
					compareExchange(pointsArray, i, i + d);
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

void runSort(Point* pointsArray, int arraySize) {
	double startTime = MPI_Wtime();
    bSort(pointsArray, arraySize);
	double endTime = MPI_Wtime();
	if (!checkCorrectness(pointsArray, arraySize)) {
		cout << "Serial: incorrect results" << endl;
	} 	
	print(pointsArray, arraySize, "Serial: array after sorting: ");
	serialTime = endTime - startTime;
	cout << "Serial time = " << serialTime << endl;
}

Point* resize(Point* array, int oldSize, int newSize) {
    auto result = new Point[newSize];
    for(int i = 0; i < newSize; i++) {
        if (i < min(oldSize, newSize)) {
            result[i] = array[i];
        } else {
            result[i] = Point({numeric_limits<float>::max(), numeric_limits<float>::max(), i});
        }
    }
    delete[] array;
    return result;
}

void runSortParallel(Point* pointsArray, int n1, int n2) {
    int partArraySize = 0;
    for(; partArraySize * procSize < n1 * n2; partArraySize++) {}

    auto readArray = new Point[partArraySize];
	double startTime = 0.0, endTime;
	if (procRank == 0) {
		startTime = MPI_Wtime();
		if (n1 * n2 != procSize * partArraySize) {
            pointsArray = resize(pointsArray, n1 * n2, procSize * partArraySize);
		}
		print(pointsArray, n1 * n2, "Before sorting: ");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	buildDerivedType(readArray, &messageType);
	MPI_Scatter(pointsArray, partArraySize, messageType, readArray, partArraySize, messageType, 0, MPI_COMM_WORLD);

    bSort(readArray, partArraySize);
    bSortParallel(readArray, partArraySize);
	
	MPI_Gather(readArray, partArraySize, messageType, pointsArray, partArraySize, messageType, 0, MPI_COMM_WORLD);
    delete[] readArray;
	if (procRank == 0) {
        if (n1 * n2 != procSize * partArraySize) {
            pointsArray = resize(pointsArray, procSize * partArraySize, n1 * n2);
        }
		endTime = MPI_Wtime();
		print(pointsArray, n1 * n2, "Parallel: after sorting: ");
		parallelTime = endTime - startTime;
		cout << "Parallel time = " << parallelTime << endl;
	}
}

Point* copyArray(const Point* pointsArray, int arraySize) {
    auto copyPointsArray = new Point[arraySize];
    for (int i = 0; i < arraySize; i++) {
        copyPointsArray[i] = pointsArray[i];
    }
    return copyPointsArray;
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);
    Point* pointsArray = nullptr;
    Point* copyPointsArray = nullptr;
    if (procRank == 0) {
        pointsArray = generatePoints(n1, n2);
    }

	if (procSize == 1) {
	    runSort(pointsArray, n1 * n2);
	    delete[] pointsArray;
	} else {
	    if (procRank == 0) {
            copyPointsArray = copyArray(pointsArray, n1 * n2);
	    }

		runSortParallel(pointsArray, n1, n2);
        if (procRank == 0) {
            print(copyPointsArray, n1 * n2, "Before sorting: ");
            runSort(copyPointsArray, n1 * n2);
            delete[] pointsArray;
            delete[] copyPointsArray;
            double effectiveness = serialTime / parallelTime / procSize;
            cout << "Effectiveness: " << effectiveness << endl;
        }
	}
	MPI_Finalize();
	return 0;
}
