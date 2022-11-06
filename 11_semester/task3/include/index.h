#ifndef TASK3_INDEX_H
#define TASK3_INDEX_H

#include "grid.h"
#include "splitter.h"


class Index {

    Grid g;

public:

    explicit Index(Grid g) : g(g) {}

    inline int operator()(int i, int j, int k) const {
        // 3d-array of points is flattened, get the linear index
        return (i * (g.N + 1) + j) * (g.N + 1) + k;
    }

    inline int operator()(int i, int j, int k, Block b) const {
        // get the linear index inside the array of the given grid block
        return (i - b.x_min) * b.y_size * b.z_size + (j - b.y_min) * b.z_size + (k - b.z_min);
    }

};


#endif //TASK3_INDEX_H
