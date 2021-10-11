#include <iostream>
#include <cstdlib>
#include <cmath>

using std::cout;
using std::pair;

int nComp = 0;
int nTact = 0;

void printArray(pair<size_t, int> *arr, size_t arrSize)
{
    for (size_t i = 0; i < arrSize; i++) {
        cout << arr[i].second << " ";
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

void sortTwoArrays(pair<size_t, int> *array, size_t firstSize, size_t secondSize, int &localTacts)
{
    if (firstSize == 0 or secondSize == 0) return;

    if (firstSize == 1 and secondSize == 1) {
        localTacts = 1;
        compare_exchange(array[0], array[1]);
        return;
    }

    // make two arrays with even elements
    size_t evenFirstSize1 = round((float)firstSize / 2);
    size_t evenSecondSize1 = round((float)secondSize / 2);
    auto evenArray = new pair<size_t, int>[evenFirstSize1 + evenSecondSize1];
    for(size_t i = 0; i < evenFirstSize1; i++) {
        evenArray[i] = array[2 * i];
    }
    for(size_t i = 0; i < evenSecondSize1; i++) {
        evenArray[evenFirstSize1 + i] = array[firstSize + 2 * i];
    }
    int firstPartTacts = 0;
    sortTwoArrays(evenArray, evenFirstSize1, evenSecondSize1, firstPartTacts);

    // make two arrays with odd elements
    size_t oddFirstSize2 = firstSize - evenFirstSize1;
    size_t oddSecondSize2 = secondSize - evenSecondSize1;
    auto oddArray = new pair<size_t, int>[oddFirstSize2 + oddSecondSize2];
    for(size_t i = 0; i < oddFirstSize2; i++) {
        oddArray[i] = array[2 * i + 1];
    }
    for(size_t i = 0; i < oddSecondSize2; i++) {
        oddArray[oddFirstSize2 + i] = array[firstSize + 2 * i + 1];
    }
    int secondPartTacts = 0;
    sortTwoArrays(oddArray, oddFirstSize2, oddSecondSize2, secondPartTacts);

    localTacts += std::max(firstPartTacts, secondPartTacts);

    // move sorted odd and even parts to array
    for(size_t i = 0; i < evenFirstSize1; i++) {
        array[2 * i] = evenArray[i];
    }
    for(size_t i = 0; i < evenSecondSize1; i++) {
        array[firstSize + 2 * i] = evenArray[evenFirstSize1 + i];
    }
    for(size_t i = 0; i < oddFirstSize2; i++) {
        array[2 * i + 1] = oddArray[i];
    }
    for(size_t i = 0; i < oddSecondSize2; i++) {
        array[firstSize + 2 * i + 1] = oddArray[oddFirstSize2 + i];
    }
    delete[] evenArray;
    delete[] oddArray;

    localTacts++;
    for(size_t i = 1; i < firstSize + secondSize - 1; i += 2) {
        compare_exchange(array[i], array[i + 1]);
    }
}

void bsort(pair<size_t, int> *arr, size_t arrSize)
{
    // there recursion ends
    if (arrSize < 2) return;
    // sort each part
    bsort(arr, arrSize / 2);
    bsort(arr + arrSize / 2, arrSize - arrSize / 2);
    sortTwoArrays(arr, arrSize / 2, arrSize - arrSize / 2, nTact);
}

int main(int argc, char **argv)
{
    size_t arrSize = atoi(argv[1]);
    cout << arrSize << " " << "0 0" << std::endl;
    auto arr = new pair<size_t, int>[arrSize];
    for (size_t i = 0; i < arrSize; i++)
        arr[i] = std::make_pair(i, rand() % arrSize);
    bsort(arr, arrSize);
    delete[] arr;
    cout << nComp << std::endl;
    cout << nTact << std::endl;
    return 0;
}