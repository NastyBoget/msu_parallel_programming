#include <iostream>
#include <vector>
#include <limits>
#include <random>
#include <algorithm>
#include <functional>
#include "mpi.h"

using namespace std;

const long RANGE = 10000000;
int procSize, procRank;
double serialTime, parTime;
bool verbose = false;

void print(const vector<int>& vec, string msg) {
	if (!verbose)
		return;
	cout << msg << endl;
	for(auto i : vec) {
		cout << i << " ";
	}
	cout << endl;
}

int findMult(int mult, int size) {
	int count = 0;
	while(count * mult < size) {
		count += 1;
	}
	return count;
}

vector<int> generateArray(int arraySize) {
	vector<int> res;
	res.reserve(arraySize);
	mt19937 eng;
	uniform_int_distribution<> d(-RANGE, RANGE); 
	auto rand = bind(d, eng);
	for(int i = 0; i < arraySize; i += 1) {
		res.push_back(rand());		
	}
	return res;
}

inline void compareExchange(vector<int>& vec, size_t i, size_t j) {
	if (vec[i] > vec[j]) {
		std::swap(vec[i], vec[j]);
	}
}

inline bool isOdd(const int& n) {
	return n % 2 == 1;
}

bool checkCorrectness(const vector<int>& vec) {
	for(size_t i = 0; i < vec.size() - 1; i += 1) {
		if (vec[i] > vec[i + 1]) {
			return false;
		}
	}
	return true;
}

vector<int> oddEvenMergeSort(const vector<int>& left, const vector<int>& right, int proc1, int proc2) {
	vector<int> tmp(left.size() + right.size());
	
	for(size_t idx = 0; idx < 2; idx += 1) {
		for(size_t i = idx, l = idx, r = idx; i < tmp.size(); i += 2) {
			if (l < left.size() and (r >= right.size() or left[l] <= right[r])) {
				tmp[i] = left[l];
				l += 2;
			} else {
				tmp[i] = right[r];
				r += 2;
			}
		}
	}
	if (isOdd(left.size()) and isOdd(right.size())) {
		tmp[tmp.size() - 1] = max(left[left.size() - 1], right[right.size() - 1]);
	}
	
	for(size_t i = 1; i < tmp.size() - 1; i += 2) {
		compareExchange(tmp, i, i + 1);
	}

	if (procRank == proc1) {
		return std::vector<int>(tmp.begin(), tmp.begin() + tmp.size() / 2);
	} else if (procRank == proc2) {
		return std::vector<int>(tmp.begin() + tmp.size() / 2, tmp.end());
	} else {
		throw "WTF?!";
	}
}

inline bool onThisProcs(int proc1, int proc2) {
	return procRank == proc1 or procRank == proc2;
}

void compareExchangePar(vector<int>& array, int proc1, int proc2) {
	int thisProc, otherProc;
	MPI_Request request;
	MPI_Status status;
	vector<int> arrayFromOtherProc(array.size());

	if (procRank == proc1) {
		thisProc = proc1;
		otherProc = proc2;
	} else {
		thisProc = proc2;
		otherProc = proc1;
	}
	
	MPI_Isend(&array[0], (int) array.size(), MPI_INT, otherProc, thisProc, MPI_COMM_WORLD, &request);
	MPI_Recv(&arrayFromOtherProc[0], arrayFromOtherProc.size(), MPI_INT, otherProc, otherProc, MPI_COMM_WORLD, &status);
	MPI_Wait(&request, &status);

	array = oddEvenMergeSort(array, arrayFromOtherProc, proc1, proc2);	
}

void batcherSortPar(vector<int>& chunk) {
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
            	if ((i & p) == r and onThisProcs(i, i + d)) {
					compareExchangePar(chunk, i, i + d);
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

void batcherSort(vector<int>& array) {
    int N = array.size();
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
					compareExchange(array, i, i + d);
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

void runSort(vector<int>& vec) {
	double startTime = MPI_Wtime();
	batcherSort(vec);
	double endTime = MPI_Wtime();
	if (!checkCorrectness(vec)) {
		cout << "Serial: incorrect results" << endl;
	} 	
	print(vec, "After sorting: ");
	serialTime = endTime - startTime;
	cout << "Serial time = " << serialTime << endl;
}

bool compareResults(vector<int> vec1, vector<int> vec2) {
	if (vec1.size() != vec2.size()) {
		return false;	
	}
	for(size_t i = 0; i < vec1.size(); i += 1) {
		if (vec1[i] != vec2[i]) {
			return false;
		}
	}
	return true;
}

void runSortPar(int arraySize) {
	int actualSizeByProc = findMult(procSize, arraySize);	
	vector<int> readBuf;
	vector<int> inputArray;
	vector<int> inArraySerial;
	double startTime = 0.0, endTime = 0.0;
	if (procRank == 0) {
		inputArray = generateArray(arraySize);		
		inArraySerial = inputArray;
		runSort(inArraySerial);
		startTime = MPI_Wtime();
		inputArray.resize(procSize * actualSizeByProc, std::numeric_limits<int>::max());
		print(inputArray, "Before sorting: ");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	readBuf.resize(actualSizeByProc);
	MPI_Scatter(&(inputArray[0]), actualSizeByProc, MPI_INT, 
				&(readBuf[0]), actualSizeByProc, MPI_INT, 0, MPI_COMM_WORLD);
	
	batcherSort(readBuf);
	
	batcherSortPar(readBuf);
	
	MPI_Gather(&(readBuf[0]), actualSizeByProc, MPI_INT, 
			   &inputArray[0], actualSizeByProc, MPI_INT, 0, MPI_COMM_WORLD);
	if (procRank == 0) {
		inputArray.resize(arraySize);
		endTime = MPI_Wtime();
		print(inputArray, "After sorting: ");
		parTime = endTime - startTime;
		cout << "Parallel time = " << parTime << endl;
		if (!compareResults(inArraySerial, inputArray)) {
			cout << "Incorrect results" << endl;
		}
	}
}


int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

	int arraySize = atoi(argv[1]);

	if (procSize == 1) {
		auto inArray = generateArray(arraySize);
		runSort(inArray);
	} else {
		runSortPar(arraySize);
        if (procRank == 0) {
            double effectiveness = serialTime / parTime / procSize;
            cout << "Effectiveness: " << effectiveness << endl;
        }
	}
	MPI_Finalize();
	return 0;
}
