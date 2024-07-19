#ifndef IPOPT_DO_UTIL_H
#define IPOPT_DO_UTIL_H

#include <vector>
#include <cstdlib>
#include <unordered_map>
#include <array>

template<typename T>
inline int sz(const std::vector<T>& vec) {
    return int(vec.size());
}

struct n2hash {
    // hashing of N x N
    std::size_t operator()(const std::array<int, 2>& v) const {
        auto h1 = std::hash<int>{}(v[0]);
        auto h2 = std::hash<int>{}(v[1]);
        return h1 ^ (h2 << 1);
    }
};

#endif //IPOPT_DO_UTIL_H
