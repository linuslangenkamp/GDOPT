#ifndef IPOPT_DO_UTIL_H
#define IPOPT_DO_UTIL_H

#include <vector>
#include <cstdlib>
#include <tuple>

template<typename T>
inline int sz(const std::vector<T>& vec) {
    return int(vec.size());
}

struct n2hash {
    // hashing of N x N
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

#endif //IPOPT_DO_UTIL_H
