#ifndef IPOPT_DO_MESH_H
#define IPOPT_DO_MESH_H

#include <utility>
#include <vector>

class Mesh {
public:
    Mesh() = default;
    Mesh(int intervals, double tf, std::vector<double> grid, std::vector<double> deltaT) :
         intervals(intervals), tf(tf), grid(std::move(grid)), deltaT(std::move(deltaT)) {}

    int intervals;
    double tf;
    std::vector<double> grid;
    std::vector<double> deltaT;

    static Mesh createEquidistantMesh(int intervals, double tf);
    void update(std::vector<int>&);
};

#endif // IPOPT_DO_MESH_H
