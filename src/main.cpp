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
    LinearSolver linearSolver = LinearSolver::MUMPS;
    MeshAlgorithm meshAlgorithm = MeshAlgorithm::L2_BOUNDARY_NORM;
    int meshIterations = 5;
    std::unordered_map<std::string, double> meshParameters;
    meshParameters.emplace("level", 0);
    meshParameters.emplace("clevel", 0);

    Solver solver = Solver(new GDOP(problem, mesh, rk, initVars), meshIterations, linearSolver, meshAlgorithm, meshParameters);

    // set solver flags
    solver.setExportOptimumPath("/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData");
    solver.setTolerance(1e-14);

    // optimize
    int status = solver.solve();
    return status;
}
