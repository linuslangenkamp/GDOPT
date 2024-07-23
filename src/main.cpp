#include <IpIpoptApplication.hpp>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "hypersensitive.h"
#include "parameterSweep.h"
#include "solver.h"

using namespace Ipopt;

int main() {
    Problem problem = createProblem_hypersensitive();
    InitVars initVars = InitVars::CONST;
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps3);
    Mesh mesh = Mesh::createEquidistantMesh(100, 100);
    LinearSolver linearSolver = LinearSolver::MA57;
    int meshIterations = 0;

    Solver solver = Solver(new GDOP(std::move(problem), mesh, rk, initVars), meshIterations, linearSolver);
    int status = solver.solve();

    return status;
}
