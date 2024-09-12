#include "gdop.h"
#include "integrator.h"
#include "mesh.h"
#include "solver.h"
#include <chrono>

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_batchReactor());
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps3);
    Mesh mesh = Mesh::createEquidistantMesh(100, 0.5);
    LinearSolver linearSolver = LinearSolver::MA57;
    MeshAlgorithm meshAlgorithm = MeshAlgorithm::L2_BOUNDARY_NORM;
    int meshIterations = 5;

    Solver solver = Solver(create_gdop(problem, mesh, rk, initVars), meshIterations, linearSolver, meshAlgorithm);

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
