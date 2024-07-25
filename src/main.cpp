#include <IpIpoptApplication.hpp>
#include <chrono>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "hypersensitive.h"
#include "parameterSweep.h"
#include "batchReactor.h"
#include "solver.h"

using namespace Ipopt;

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_batchReactor());
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps6);
    Mesh mesh = Mesh::createEquidistantMesh(15, 1);
    LinearSolver linearSolver = LinearSolver::MUMPS;
    int meshIterations = 5;

    Solver solver = Solver(new GDOP(problem, mesh, rk, initVars), meshIterations, linearSolver);

    // set solver flags
    solver.setExportOptimumPath("/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData/batchReactorRefinement");

    // optimize
    auto solveStartTime = std::chrono::high_resolution_clock::now();
    int status = solver.solve();
    auto timeTaken = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - solveStartTime).count();
    std::cout << "\nSolving took: " << timeTaken << " seconds" << std::endl;

    return status;
}
