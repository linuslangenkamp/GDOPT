//
// Created by Linus on 14.07.2024.
//

#ifndef IPOPT_DO_MESH_H
#define IPOPT_DO_MESH_H

#include <utility>
#include <vector>

class Mesh {
public:
    Mesh(int intervals,
    double tf,
    std::vector<double> grid,
    std::vector<double> deltaT)
    : intervals(intervals), tf(tf), grid(std::move(grid)), deltaT(std::move(deltaT)) {}

    int intervals;
    double tf;
    std::vector<double> grid;
    std::vector<double> deltaT;
};


Mesh createEquidistantMesh(int intervals, double tf) {
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

#endif //IPOPT_DO_MESH_H
