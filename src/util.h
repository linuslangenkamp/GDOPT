#ifndef IPOPT_DO_UTIL_H
#define IPOPT_DO_UTIL_H

#include <vector>
#include <cstdlib>
#include <tuple>
#include <fstream>
#include <iostream>

template<typename T>
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

inline int fullSum(const std::vector<std::vector<int>>& matrix) {
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
    outFile << "dimRow" << "," << "dimCol" << "\n";
    outFile << dimRow << "," << dimCol << "\n";

    outFile << "row" << "," << "col" << "\n";
    for (int i = 0; i < L; i++) {
        outFile << iRow[i] << "," << jCol[i] << "\n";
    }
    outFile.close();
}

#endif //IPOPT_DO_UTIL_H
