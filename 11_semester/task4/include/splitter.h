#ifndef TASK4_SPLIT_H
#define TASK4_SPLIT_H

#include <vector>


enum Axis {
    X, Y, Z
};


struct Block {
    int x_min;
    int x_max;
    int y_min;
    int y_max;
    int z_min;
    int z_max;
    int x_size;
    int y_size;
    int z_size;
    int size;

    Block(int x_min, int x_max, int y_min, int y_max, int z_min, int z_max) : x_min(x_min), x_max(x_max),
        y_min(y_min), y_max(y_max), z_min(z_min), z_max(z_max) {
        x_size = x_max - x_min + 1;
        y_size = y_max - y_min + 1;
        z_size = z_max - z_min + 1;
        size = x_size * y_size * z_size;
    }

    Block() {}
};


void split(int x_min, int x_max, int y_min, int y_max, int z_min, int z_max, int size, Axis axis, std::vector<Block> &blocks) {
    // split the grid recursively into blocks
    if (size == 1) {
        blocks.emplace_back(x_min, x_max, y_min, y_max, z_min, z_max);
        return;
    }

    if (size % 2 == 1) { // if size is odd we make it even
        if (axis == X) {
            int x = x_min + (x_max - x_min) / size;
            blocks.emplace_back(x_min, x, y_min, y_max, z_min, z_max);
            x_min = x + 1;
            axis = Y;
        }
        else if (axis == Y) {
            int y = y_min + (y_max - y_min) / size;
            blocks.emplace_back(x_min, x_max, y_min, y, z_min, z_max);
            y_min = y + 1;
            axis = Z;
        }
        else { // axis == Z
            int z = z_min + (z_max - z_min) / size;
            blocks.emplace_back(x_min, x_max, y_min, y_max, z_min, z);
            z_min = z + 1;
            axis = X;
        }

        size--;
    }

    // now the size is even
    if (axis == X) {
        int x = (x_min + x_max) / 2;
        split(x_min, x, y_min, y_max, z_min, z_max, size / 2, Y, blocks);
        split(x + 1, x_max, y_min, y_max, z_min, z_max, size / 2, Y, blocks);
    }
    else if (axis == Y) {
        int y = (y_min + y_max) / 2;
        split(x_min, x_max, y_min, y, z_min, z_max, size / 2, Z, blocks);
        split(x_min, x_max, y + 1, y_max, z_min, z_max, size / 2, Z, blocks);
    }
    else {
        int z = (z_min + z_max) / 2;
        split(x_min, x_max, y_min, y_max, z_min, z, size / 2, X, blocks);
        split(x_min, x_max, y_min, y_max, z + 1, z_max, size / 2, X, blocks);
    }
}


#endif //TASK4_SPLIT_H
