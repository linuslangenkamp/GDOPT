#ifndef IPOPT_DO_SOLVER_H
#define IPOPT_DO_SOLVER_H

#include "gdop.h"

enum class LinearSolver {
    MUMPS,
    MA27,
    MA57,
    MA77,
    MA86,
    MA97,
    PARDISO
};

class Solver {
public:
    explicit Solver(const SmartPtr<GDOP>& gdop, int maxMeshIterations, LinearSolver linearSolver);

    SmartPtr<GDOP> gdop;
    LinearSolver linearSolver;
    const int maxMeshIterations;
    const double tolerance = 1e-14;

    int solve() const;
};

/* add interpolation type
 * add flags: mesh iterations, type of mesh refinement
 */
#endif //IPOPT_DO_SOLVER_H
