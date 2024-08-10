#include <IpIpoptApplication.hpp>
#include <chrono>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "hypersensitive.h"
#include "parameterSweep.h"
#include "rocketTrajectory.h"
#include "analyticHypersensitive.h"
#include "trivialBangBang.h"
#include "satellite.h"
#include "dieselMotor.h"
#include "batchReactor.h"
#include "solver.h"

using namespace Ipopt;

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_dieselMotor());
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps7);
    Mesh mesh = Mesh::createEquidistantMesh(500, 0.5);
    LinearSolver linearSolver = LinearSolver::MA57;
    MeshAlgorithm meshAlgorithm = MeshAlgorithm::L2_BOUNDARY_NORM;
    int meshIterations = 0;

    Solver solver = Solver(new GDOP(problem, mesh, rk, initVars), meshIterations, linearSolver, meshAlgorithm);

    // set solver flags
    solver.setExportOptimumPath("/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData");
    solver.setTolerance(1e-13);

    // set solver mesh parameters
    solver.setMeshParameter("level", 0);
    solver.setMeshParameter("ctol", 0.1);

    // optimize
    int status = solver.solve();
    return status;
}
