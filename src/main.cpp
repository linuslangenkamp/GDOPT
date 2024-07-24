#include <IpIpoptApplication.hpp>
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
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps3);
    Mesh mesh = Mesh::createEquidistantMesh(200, 10000);
    LinearSolver linearSolver = LinearSolver::MA57;
    int meshIterations = 25;

    Solver solver = Solver(new GDOP(problem, mesh, rk, initVars), meshIterations, linearSolver);
    int status = solver.solve();

    return status;
}
