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

    static Mesh createEquidistantMesh(int intervals, double tf);
};

#endif // IPOPT_DO_MESH_H
