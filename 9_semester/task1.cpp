#include <iostream>
#include <cstdlib>
#include <vector>

using std::cout;
using std::vector;
using std::pair;

int nComp = 0;
int nTact = 0;

void printArray(vector<pair<size_t, int> > &arr)
{
    for (auto &it : arr) {
        cout << it.second << " ";
    }
    cout << std::endl;
}

void compare_exchange(pair<size_t, int> &firstElem, pair<size_t, int> &secondElem)
{
    cout << firstElem.first << " " << secondElem.first << std::endl;
    if (firstElem.second > secondElem.second) {
        int tmp = firstElem.second;
        firstElem.second = secondElem.second;
        secondElem.second = tmp;
    }
    nComp++;
}

vector<pair<size_t, int> > sortTwoArrays(vector<pair<size_t, int> > &firstVector,
                                         vector<pair<size_t, int> > &secondVector, int &localTacts)
{
    auto firstSize = firstVector.size();
    auto secondSize = secondVector.size();

    if (firstSize == 0) return secondVector;

    if (secondSize == 0) return firstVector;

    if (firstSize == 1 and secondSize == 1) {
        localTacts = 1;
        compare_exchange(firstVector[0], secondVector[0]);
        firstVector.insert(firstVector.end(), secondVector.begin(), secondVector.end());
        return firstVector;
    }

    vector<pair<size_t, int> > newFirstVector;
    vector<pair<size_t, int> > newSecondVector;
    for(size_t i = 0; i < firstSize; i += 2) {
        newFirstVector.push_back(firstVector[i]);
    }
    for(size_t i = 0; i < secondSize; i += 2) {
        newSecondVector.push_back(secondVector[i]);
    }

    int firstPartTacts = 0;
    auto firstPart = sortTwoArrays(newFirstVector, newSecondVector, firstPartTacts);

    newFirstVector.clear();
    newSecondVector.clear();
    for(size_t i = 1; i < firstSize; i += 2) {
        newFirstVector.push_back(firstVector[i]);
    }
    for(size_t i = 1; i < secondSize; i += 2) {
        newSecondVector.push_back(secondVector[i]);
    }
    int secondPartTacts = 0;
    auto secondPart = sortTwoArrays(newFirstVector, newSecondVector, secondPartTacts);

    localTacts += std::max(firstPartTacts, secondPartTacts);
    localTacts++;
    for(size_t i = 0; i < std::min(firstPart.size() - 1, secondPart.size()); i++) {
        compare_exchange(secondPart[i], firstPart[i + 1]);
    }

    int i = 1;
    for(auto &it : secondPart) {
        vector<pair<size_t, int> >::const_iterator first = firstPart.begin() + i;
        firstPart.insert(first, it);
        i += 2;
    }
    return firstPart;
}

void bsort(vector<pair<size_t, int> > &arr)
{
    // there recursion ends
    auto arrSize = arr.size();
    if (arrSize < 2) return;

    // make two parts
    vector<pair<size_t, int> >::const_iterator first = arr.begin();
    vector<pair<size_t, int> >::const_iterator middle = arr.begin() + (arrSize / 2);
    vector<pair<size_t, int> >::const_iterator last = arr.end();
    vector<pair<size_t, int> > firstArray(first, middle);
    vector<pair<size_t, int> > secondArray(middle, last);
    // sort each part
    bsort(firstArray);
    bsort(secondArray);
    arr = sortTwoArrays(firstArray, secondArray, nTact);
}

int main(int argc, char **argv)
{
    size_t arrSize = atoi(argv[1]);
    cout << arrSize << " " << "0 0" << std::endl;
    vector<pair<size_t, int> > arr(arrSize);
    for (size_t i = 0; i < arrSize; i++)
        arr[i] = std::make_pair(i, rand() % arrSize);

    bsort(arr);
    cout << nComp << std::endl;
    cout << nTact << std::endl;
    return 0;
}