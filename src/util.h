#ifndef IPOPT_DO_UTIL_H
#define IPOPT_DO_UTIL_H

#include <vector>
#include <cstdlib>

template<typename T>
int sz(const std::vector<T>& vec) {
    return static_cast<int>(vec.size());
}

#endif //IPOPT_DO_UTIL_H
