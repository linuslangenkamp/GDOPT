#include <chrono>

#include "gdop.h"
#include "integrator.h"
#include "mesh.h"
#include "solver.h"

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_batchReactor());
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps3);
    Mesh mesh = Mesh::createEquidistantMesh(100, 0.5);
    LinearSolver linearSolver = LinearSolver::MA57;
    MeshAlgorithm meshAlgorithm = MeshAlgorithm::L2_BOUNDARY_NORM;
    int meshIterations = 5;

    Solver solver = Solver(create_gdop(problem, mesh, rk, initVars), meshIterations, linearSolver, meshAlgorithm);
    
    // set solver mesh parameters
    solver.setMeshParameter("level", 0);
    solver.setMeshParameter("ctol", 0.1);

    // optimize
    int status = solver.solve();
    return status;
}
