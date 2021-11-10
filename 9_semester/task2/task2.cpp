#include <iostream>
#include <limits>
#include <random>
#include <vector>
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

bool comp(Point i, Point j) {
    return i.coord[coordSorted] < j.coord[coordSorted];
}

vector<Point> generatePoints(int n1, int n2) {
    vector<Point> result;
    result.reserve(n1 * n2);
    for(int i = 0; i < n1; i++)
        for(int j = 0; j < n2; j++) {
            result.emplace_back(Point({x(i, j), y(i, j), i * n2 + j}));
        }
    return result;
}

void print(const vector<Point>& pointsArray, const string& msg) {
	if (!verbose)
		return;
	cout << msg << endl;
    for(auto& i : pointsArray) {
		cout << i.coord[coordSorted] << " ";
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

void buildDerivedType(const vector<Point>& indata, MPI_Datatype* message_type_ptr)
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
    auto arraySize = left.size();
    auto tmpSize = 2 * arraySize;
    vector<Point> tmpArray(tmpSize);
	
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
        std::copy(tmpArray.begin(), tmpArray.begin() + tmpSize / 2, left.begin());
    } else if (procRank == proc2) {
        std::copy(tmpArray.begin() + tmpSize / 2, tmpArray.end(), left.begin());
    }
}

void compareExchangeParallel(vector<Point>& currentArray, int proc1, int proc2) {
    auto arraySize = currentArray.size();
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
	serialTime = endTime - startTime;
	cout << "Serial time = " << serialTime << endl;
}

void runSortParallel(vector<Point> pointsArray, int arraySize) {
    int partArraySize = 0;
    for(; partArraySize * procSize < arraySize; partArraySize++) {}

    vector<Point> readArray;
	double startTime = 0.0, endTime;
	if (procRank == 0) {
		startTime = MPI_Wtime();
        pointsArray.resize(procSize * partArraySize,
                Point({numeric_limits<float>::max(), numeric_limits<float>::max(), procSize * partArraySize}));
		print(pointsArray, "Before sorting: ");
	}
    buildDerivedType(readArray, &messageType);
	MPI_Barrier(MPI_COMM_WORLD);
    readArray.resize(partArraySize);
	MPI_Scatter(pointsArray.data(), partArraySize, messageType, readArray.data(),
	        partArraySize, messageType, 0, MPI_COMM_WORLD);

    sort(readArray.begin(), readArray.end(), comp);
    bSortParallel(readArray);
	
	MPI_Gather(readArray.data(), partArraySize, messageType, pointsArray.data(),
	        partArraySize, messageType, 0, MPI_COMM_WORLD);
	if (procRank == 0) {
        pointsArray.resize(arraySize);
		endTime = MPI_Wtime();
		print(pointsArray, "Parallel: after sorting: ");
		parallelTime = endTime - startTime;
		cout << "Parallel time = " << parallelTime << endl;
	}
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	int n1 = atoi(argv[1]);
    int n2 = atoi(argv[2]);
    vector<Point> pointsArray, copyPointsArray;
    if (procRank == 0) {
        pointsArray = generatePoints(n1, n2);
    }

	if (procSize == 1) {
	    runSort(pointsArray);
	} else {
	    if (procRank == 0) {
            std::copy(pointsArray.begin(), pointsArray.end(), std::back_inserter(copyPointsArray));
	    }

		runSortParallel(pointsArray, n1 * n2);
        if (procRank == 0) {
            print(copyPointsArray, "Before sorting: ");
            runSort(copyPointsArray);
            double effectiveness = serialTime / parallelTime / procSize;
            cout << "Effectiveness: " << effectiveness << endl;
        }
	}
	MPI_Finalize();
	return 0;
}
