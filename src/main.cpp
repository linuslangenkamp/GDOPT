#include <IpIpoptApplication.hpp>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "parameterSweep.h"
#include "solver.h"

using namespace Ipopt;

int main() {
    Problem problem = createProblem_parameterSweep();
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps3);
    Mesh mesh = Mesh::createEquidistantMesh(4000, 5);

    Solver solver = Solver(new GDOP(std::move(problem), mesh, rk, InitVars::CONST), 0);
    int status = solver.solve();

    return status;
}
