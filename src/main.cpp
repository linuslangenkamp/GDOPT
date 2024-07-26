#include <IpIpoptApplication.hpp>
#include <chrono>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "hypersensitive.h"
#include "parameterSweep.h"
#include "rocketTrajectory.h"
#include "batchReactor.h"
#include "solver.h"

using namespace Ipopt;

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_batchReactor());
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps2);
    Mesh mesh = Mesh::createEquidistantMesh(50, 1);
    LinearSolver linearSolver = LinearSolver::MUMPS;
    int meshIterations = 25;

    Solver solver = Solver(new GDOP(problem, mesh, rk, initVars), meshIterations, linearSolver);

    // set solver flags
    solver.setExportOptimumPath("/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData/rocketTrajectory");
    solver.setTolerance(1e-14);

    // optimize
    int status = solver.solve();
    return status;
}
