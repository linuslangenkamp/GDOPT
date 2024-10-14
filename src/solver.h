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

#ifndef GDOPT_SOLVER_H
#define GDOPT_SOLVER_H

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
    Solver(GDOP* gdop);
    ~Solver();

    // important methods
    int solve();
    std::vector<int> detect() const;
    void refine(std::vector<int>& markedIntervals);

    // detection methods
    std::vector<int> basicStrategy() const;
    std::vector<int> L2BoundaryNorm() const;

    // interpolation type for the new initial solution
    void refinePolynomial(std::vector<int>& markedIntervals);
    void refineLinear(std::vector<int>& markedIntervals);

    // mesh and mesh refinement
    int initialIntervals;          // intervals in the first mesh
    int meshIteration = 0;         // mesh iteration counter
    std::vector<double> cbValues;  // starting values after refinement

    // timing, printing, ...
    std::vector<int> numberOfIntervalsHistory{};
    std::vector<double> ipoptObjectiveHistory{};
    std::vector<int> ipoptIterationHistory{};
    std::vector<double> ipoptIterationTotalTime{};
    std::vector<double> ipoptIterationFuncEvalTime{};
    std::vector<double> ipoptIterationNonfuncEvalTime{};
    std::chrono::high_resolution_clock::time_point solveStartTime;
    std::chrono::duration<double> solveTotalTimeTaken{};   // total time in solver
    std::chrono::duration<double> solveActualTimeTaken{};  // total time in solver - IO
    std::chrono::duration<double> timedeltaIO{0};          // time in IO operations
    double ipoptFuncTime{};                                // sum of time in ipopt func evals
    double ipoptTotalTime{};                               // sum of total time in ipopt
    double ipoptActualTime{};                              // sum of total time - func evals in ipopt

    void initSolvingProcess();
    void finalizeOptimization(Ipopt::IpoptApplication& app);
    void postOptimization(Ipopt::IpoptApplication& app);
    void setSolverFlags(Ipopt::IpoptApplication& app);
    void setStandardSolverFlags(Ipopt::IpoptApplication& app);
    void createModelInfo(Ipopt::IpoptApplication& app) const;

    void printASCIIArt() const;
    void printMeshIterationHistory();
    void setExportSparsityPath();

private:
    std::unique_ptr<SolverPrivate> _priv;
};

#endif  // GDOPT_SOLVER_H
