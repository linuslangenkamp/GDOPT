#ifndef IPOPT_DO_UTIL_H
#define IPOPT_DO_UTIL_H

#include <vector>
#include <cstdlib>

template<typename T>
inline int sz(const std::vector<T>& vec) {
    return int(vec.size());
}

#endif //IPOPT_DO_UTIL_H
