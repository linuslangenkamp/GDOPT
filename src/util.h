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

#ifndef GDOPT_UTIL_H
#define GDOPT_UTIL_H

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

template <typename T>
inline int sz(const std::vector<T>& vec) {
    return int(vec.size());
}

struct n2hash {
    std::size_t operator()(const std::tuple<int, int>& v) const {
        auto const [t1, t2] = v;
        auto h1 = std::hash<int>{}(t1);
        auto h2 = std::hash<int>{}(t2);
        return h1 ^ (h2 << 1);
    }
};

inline int nnzMatrix(const std::vector<std::vector<int>>& matrix) {
    int sum = 0;
    for (const auto& row : matrix) {
        for (const auto val : row) {
            sum += val;
        }
    }
    return sum;
}

inline void exportSparsity(const int* iRow, const int* jCol, const int L, const std::tuple<int, int> dim, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    auto const [dimRow, dimCol] = dim;
    outFile << "dimRow"
            << ","
            << "dimCol"
            << "\n";
    outFile << dimRow << "," << dimCol << "\n";

    outFile << "row"
            << ","
            << "col"
            << "\n";
    for (int i = 0; i < L; i++) {
        outFile << iRow[i] << "," << jCol[i] << "\n";
    }
    outFile.close();
}

inline std::string double2Str(double value, int prec) {
    std::stringstream out;
    out << std::scientific << std::setprecision(prec) << value;
    return out.str();
}

inline std::string printTime(double value) {
    std::stringstream out;
    out << std::fixed << std::setprecision(5) << value;
    return out.str();
}

inline double calculateMean(const std::vector<double>& vec) {
    double sum = 0;
    for (double num : vec) {
        sum += num;
    }
    return sum / sz(vec);
}

inline double calculateStdDev(const std::vector<double>& vec, double mean) {
    double sum = 0;
    for (double num : vec) {
        sum += (num - mean) * (num - mean);
    }
    return std::sqrt(sum / sz(vec));
}

// functions to read in the initial state guesses
inline std::vector<double> split(const std::string& line, char delimiter) {
    std::vector<double> result;
    std::stringstream stream(line);
    std::string item;

    while (std::getline(stream, item, delimiter)) {
        result.push_back(std::stod(item));
    }
    return result;
}

inline std::vector<std::vector<double>> readInitialValues(const std::string& filePath) {
    std::ifstream file(filePath);
    std::string line;
    std::vector<std::vector<double>> trajectories;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        return trajectories;
    }

    while (std::getline(file, line)) {
        if (!line.empty()) {
            std::vector<double> row = split(line, ',');
            trajectories.push_back(row);
        }
    }

    file.close();
    return trajectories;
}

inline double checkNominalValue(const double value) {
    double nom = std::abs(value);
    if (nom >= 1e-18 and 1e18 >= nom) {
        return 1 / nom;
    }
    else {
        std::cerr << "The program will be terminated, since the nominal value abs(" << value << ") is not in the range [1e-18, 1e18]." << std::endl;
        abort();
    }
}

#endif  // GDOPT_UTIL_H
