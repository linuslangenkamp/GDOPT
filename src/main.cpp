#include <IpIpoptApplication.hpp>
#include <chrono>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"
#include "../examples/hypersensitive.h"
#include "../examples/rocketTrajectory.h"
#include "../examples/analyticHypersensitive.h"
#include "../examples/trivialBangBang.h"
#include "../examples/satellite.h"
#include "../examples/dieselMotor.h"
#include "../codegen/examples/simpleParameterGenerated.h"
#include "../examples/batchReactor.h"


using namespace Ipopt;

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_satellite());
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps12);
    Mesh mesh = Mesh::createEquidistantMesh(3, 100);
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
