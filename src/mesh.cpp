//
// Created by Linus on 14.07.2024.
//
#include "mesh.h"

Mesh Mesh::createEquidistantMesh(int intervals, double tf) {
    std::vector<double> grid(intervals + 1);
    std::vector<double> deltaT(intervals);
    double h = tf / intervals;

    for (int i = 0; i <= intervals; ++i) {
        grid[i] = i * h;
    }

    for (int i = 0; i < intervals; ++i) {
        deltaT[i] = h;
    }

    return {intervals, tf, std::move(grid), std::move(deltaT)};
}
