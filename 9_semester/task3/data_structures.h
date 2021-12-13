#ifndef TASK3_DATA_STRUCTURES_H
#define TASK3_DATA_STRUCTURES_H

#include <vector>
#include <mpi.h>

using namespace std;
extern int coordSorted;
extern MPI_Datatype messageType;

struct Point {
    float coord[2];
    int index;
};

bool comp(Point i, Point j);
void runSortParallel(vector<Point>& partArray, MPI_Comm comm);
vector<Point> generatePoints(int n1, int n2, int procRank, int partArraySize, int procSize);
void buildDerivedType(MPI_Datatype* message_type_ptr);
void printResults(const vector<Point> &partArray, const vector<int> &domains, const string &fileName, int n1, int n2);
void print(const vector<Point>& pointsArray);

#endif //TASK3_DATA_STRUCTURES_H
