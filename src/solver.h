#ifndef IPOPT_DO_SOLVER_H
#define IPOPT_DO_SOLVER_H

#include <chrono>
#include <unordered_map>
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

enum class MeshAlgorithm {
    NONE,
    BASIC,
    L2_BOUNDARY_NORM
};

struct SolverPrivate;
namespace Ipopt {
    class IpoptApplication;
}


class Solver {
public:
    Solver(GDOP* gdop, int maxMeshIterations, LinearSolver linearSolver, MeshAlgorithm meshAlgorithm);
    ~Solver();

    LinearSolver linearSolver;
    MeshAlgorithm meshAlgorithm;
    int meshIteration = 0;
    const int maxMeshIterations;
    double tolerance = 1e-12;
    std::vector<double> cbValues;           // starting values after refinement
    std::unordered_map<std::string, double> meshParameters;

    // important methods
    int solve();
    std::vector<int> detect() const;
    void refine(std::vector<int>& markedIntervals);
    void finalizeOptimization();

    // detection methods
    std::vector<int> basicStrategy() const;
    std::vector<int> l2BoundaryNorm() const;

    // additional / optional flags, printouts, ...
    std::vector<double> objectiveHistory;
    int initialIntervals = -1;
    std::chrono::_V2::system_clock::time_point solveStartTime;
    std::chrono::duration<double> timedeltaIO{0}; // time in IO operations
    void postOptimization();
    std::string exportOptimumPath;
    bool exportOptimum = false;
    std::string exportHessianPath;
    bool exportHessian = false;
    std::string exportJacobianPath;
    bool exportJacobian = false;
    void printObjectiveHistory();
    void printMeshStats() const;
    void createModelInfo() const;
    void setExportOptimumPath(const std::string& exportPath);
    void setExportJacobianPath(const std::string& exportPath);
    void setExportHessianPath(const std::string& exportPath);
    void initSolvingProcess();
    void setTolerance(double tol);
    void setRefinementParameters();
    void setl2BoundaryNorm();
    void setBasicStrategy();
    void setMeshParameter(const std::string& field, double value);
    void setSolverFlags(Ipopt::IpoptApplication& app);

    // basicStrategy parameters
    double basicStrategySigma = 2.5;

    // L2 norm / boundary parameters
    double L2Level = 0;
    double L2CornerTol = 0.2;
private:
    std::unique_ptr<SolverPrivate> _priv;
};

#endif //IPOPT_DO_SOLVER_H
