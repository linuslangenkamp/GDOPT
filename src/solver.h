#ifndef IPOPT_DO_SOLVER_H
#define IPOPT_DO_SOLVER_H

#include "gdop.h"

class Solver {
public:
    explicit Solver(const SmartPtr<GDOP>& gdop, const int maxMeshIterations);

    SmartPtr<GDOP> gdop;
    const int maxMeshIterations;
    const double tolerance = 1e-14;

    int solve() const;
};

/* add interpolation type
 * add flags: mesh iterations, type of mesh refinement
 */
#endif //IPOPT_DO_SOLVER_H
