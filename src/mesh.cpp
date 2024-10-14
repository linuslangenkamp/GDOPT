/**
 * GDOPT - General Dynamic Optimizer
 * Copyright (C) 2024  Linus Langenkamp
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/

#include "mesh.h"

#include "util.h"

Mesh Mesh::createEquidistantMesh(int intervals, double tf) {
    std::vector<double> grid(intervals + 1);
    std::vector<double> deltaT(intervals);
    double h = tf / intervals;
    for (int i = 0; i <= intervals; i++) {
        grid[i] = i * h;
    }
    for (int i = 0; i < intervals; i++) {
        deltaT[i] = h;
    }
    return {intervals, tf, std::move(grid), std::move(deltaT)};
}

void Mesh::update(std::vector<int>& markedIntervals) {
    int index = 0;
    std::vector<double> newGrid{};
    for (int i = 0; i < intervals; i++) {
        newGrid.push_back(grid[i]);
        if (markedIntervals[index] == i && index < sz(markedIntervals)) {
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
