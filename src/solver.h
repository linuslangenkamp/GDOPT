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
    Solver(const SmartPtr<GDOP>& gdop, int maxMeshIterations, LinearSolver linearSolver);

    SmartPtr<GDOP> gdop;
    LinearSolver linearSolver;
    int meshIteration = 0;
    const int maxMeshIterations;
    const double tolerance = 1e-14;
    std::vector<double> cbValues;           // starting values after refinement

    // important methods
    int solve();
    std::vector<int> basicStochasticStrategy(double) const;
    void refine(std::vector<int>&);
    void finalizeOptimization() const;


    // additional / optional flags, ...
    void postOptimization();
    std::vector<double> objectiveHistory;   // history of objectives in refinement process
    std::string exportOptimumPath;
    bool exportOptimum = false;
    void setExportOptimumPath(const std::string&);

};

#endif //IPOPT_DO_SOLVER_H
