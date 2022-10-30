#ifndef TASK3_INDEX_H
#define TASK3_INDEX_H

#include "grid.h"


class Index {

    Grid g;

public:

    explicit Index(Grid g) : g(g) {}

    inline int operator()(int i, int j, int k) const {
        // 3d-array of points is flattened, get the linear index
        return (i * (g.N + 1) + j) * (g.N + 1) + k;
    }
};


#endif //TASK3_INDEX_H
