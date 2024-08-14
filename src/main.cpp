#include <IpIpoptApplication.hpp>
#include <chrono>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"
#include "../examples/hypersensitive.h"
#include "../examples/analyticHypersensitive.h"
#include "../examples/trivialBangBang.h"
#include "../codegen/examples/satelliteGenerated.h"
#include "../codegen/examples/dieselMotorGenerated.h"
#include "../codegen/examples/simpleParameterGenerated.h"
#include "../codegen/examples/pureParameterGenerated.h"
#include "../codegen/examples/invertedPendulumGenerated.h"
#include "../examples/batchReactor.h"

using namespace Ipopt;

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_invertedPendulum());
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps1);
    Mesh mesh = Mesh::createEquidistantMesh(10000, 12);
    LinearSolver linearSolver = LinearSolver::MA57;
    MeshAlgorithm meshAlgorithm = MeshAlgorithm::L2_BOUNDARY_NORM;
    int meshIterations = 5;

    Solver solver = Solver(new GDOP(problem, mesh, rk, initVars), meshIterations, linearSolver, meshAlgorithm);

    // set solver flags
    // "/home/linus/Documents/outputsGDOP"
    solver.setExportOptimumPath("/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData");
    // solver.setExportHessianPath("/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/Sparsity/hessianSparsity.csv");
    // solver.setExportJacobianPath("/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/Sparsity/jacobianSparsity.csv");
    solver.setTolerance(1e-13);

    // set solver mesh parameters
    solver.setMeshParameter("level", 0);
    solver.setMeshParameter("ctol", 0.1);

    // optimize
    int status = solver.solve();
    return status;
}
