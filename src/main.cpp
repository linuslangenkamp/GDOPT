#include <IpIpoptApplication.hpp>
#include <chrono>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "hypersensitive.h"
#include "parameterSweep.h"
#include "rocketTrajectory.h"
#include "trivialBangBang.h"
#include "batchReactor.h"
#include "solver.h"

using namespace Ipopt;

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_batchReactor());
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps2);
    Mesh mesh = Mesh::createEquidistantMesh(50, 1);
    LinearSolver linearSolver = LinearSolver::MA57;
    MeshAlgorithm meshAlgorithm = MeshAlgorithm::BASIC;
    int meshIterations = 10;

    Solver solver = Solver(new GDOP(problem, mesh, rk, initVars), meshIterations, linearSolver, meshAlgorithm);

    // set solver flags
    solver.setExportOptimumPath("/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData/batchReactorRefinement");
    solver.setTolerance(1e-14);

    // optimize
    int status = solver.solve();
    return status;
}
