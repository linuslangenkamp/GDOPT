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

#ifndef IPOPT_DO_SOLVER_H
#define IPOPT_DO_SOLVER_H

#include <chrono>
#include <set>
#include <unordered_map>

#include "gdop.h"

enum class LinearSolver { MUMPS, MA27, MA57, MA77, MA86, MA97, PARDISO };

enum class RefinementMethod { POLYNOMIAL, LINEAR_SPLINE };

/*
 * NONE: no mesh algorithm
 * BASIC: basic strategy based on deviation to the mean of control vars
 * L2_BOUNDARY_NORM: estimates the steepness by evaluating the L2 norm of the 1st and 2nd diff of the interpolating poly for each control on each sub-interval
 * and calculates the steepness difference between adjacent intervals. If these identifier are above some eps -> insert the interval + the adjacent intervals
 */

// TODO: add L2BN variants: quadsection, non bisective -> blurry split, remove consts!
enum class MeshAlgorithm { NONE, BASIC, L2_BOUNDARY_NORM };

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
    RefinementMethod refinementMethod = RefinementMethod::LINEAR_SPLINE;
    int meshIteration = 0;
    const int maxMeshIterations;
    double tolerance = 1e-12;
    std::vector<double> cbValues;  // starting values after refinement
    std::unordered_map<std::string, double> meshParameters;

    // important methods
    int solve();
    std::vector<int> detect() const;
    void refine(std::vector<int>& markedIntervals);
    void finalizeOptimization();
    void postOptimization();

    // detection methods
    std::vector<int> basicStrategy() const;
    std::vector<int> L2BoundaryNorm() const;

    // interpolation type for the new initial solution
    void refinePolynomial(std::vector<int>& markedIntervals);
    void refineLinear(std::vector<int>& markedIntervals);

    // additional / optional flags, printouts, ...
    std::vector<double> objectiveHistory;
    int initialIntervals;
    std::chrono::_V2::system_clock::time_point solveStartTime;
    std::chrono::duration<double> solveTotalTimeTaken{};   // total time in solver
    std::chrono::duration<double> solveActualTimeTaken{};  // total time in solver - IO
    std::chrono::duration<double> timedeltaIO{0};          // time in IO operations
    std::string exportOptimumPath;
    bool exportOptimum = false;
    std::string exportHessianPath;
    bool exportHessian = false;
    std::string exportJacobianPath;
    bool exportJacobian = false;
    bool userScaling = true;  // use nlp scaling provided by the user, otherwise use gradient-based

    void printObjectiveHistory();
    void printMeshStats() const;
    void createModelInfo() const;
    void setRefinementMethod(RefinementMethod method);
    void setExportOptimumPath(const std::string& exportPath);
    void setExportJacobianPath(const std::string& exportPath);
    void setExportHessianPath(const std::string& exportPath);
    void initSolvingProcess();
    void setTolerance(double tol);
    void setRefinementParameters();
    void setL2BoundaryNorm();
    void setBasicStrategy();
    void setMeshParameter(const std::string& field, double value);
    void setSolverFlags(Ipopt::IpoptApplication& app);
    void setStandardSolverFlags(Ipopt::IpoptApplication& app);

    // basicStrategy parameters
    double basicStrategySigma = 2.5;

    // L2 norm / boundary parameters
    double L2Level = 0;
    double L2CornerTol = 0.2;

private:
    std::unique_ptr<SolverPrivate> _priv;
};

#endif  // IPOPT_DO_SOLVER_H
