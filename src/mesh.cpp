#include "mesh.h"

#include "util.h"

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

void Mesh::update(std::vector<int>& markedIntervals) {
    int index = 0;
    std::vector<double> newGrid;
    for (int i = 0; i < intervals; i++) {
        newGrid.push_back(grid[i]);
        if (markedIntervals[index] == i) {
            newGrid.push_back(grid[i] + deltaT[i] / 2);
            index++;
        }
    }
    newGrid.push_back(tf);
    grid = newGrid;
    deltaT = {};
    for (int k = 0; k < sz(grid) - 1; k++) {
        deltaT.push_back(grid[k + 1] - grid[k]);
    }
    intervals = sz(deltaT);
}
