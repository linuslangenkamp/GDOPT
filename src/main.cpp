#include <IpIpoptApplication.hpp>
#include <chrono>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "hypersensitive.h"
#include "parameterSweep.h"
#include "solver.h"

using namespace Ipopt;

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_hypersensitive());
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps7);
    Mesh mesh = Mesh::createEquidistantMesh(50, 10000);
    LinearSolver linearSolver = LinearSolver::MA57;
    int meshIterations = 25;

    Solver solver = Solver(new GDOP(problem, mesh, rk, initVars), meshIterations, linearSolver);

    auto solveStartTime = std::chrono::high_resolution_clock::now();
    int status = solver.solve();
    auto timeTaken = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - solveStartTime).count();
    std::cout << "\nSolving took: " << timeTaken << " seconds" << std::endl;

    return status;
}
