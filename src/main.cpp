#include <iostream>
#include <IpIpoptApplication.hpp>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "hypersensitive.h"

using namespace Ipopt;

int main() {
    Problem problem = createProblem_hypersensitive();
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps3);
    Mesh mesh = Mesh::createEquidistantMesh(10000, 100);

    SmartPtr<GDOP> DOP{new GDOP(std::move(problem), mesh, rk, InitVars::CONST)};
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    // numeric jacobian and hessian
    // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    // app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");

    app->Options()->SetNumericValue("tol", 1e-14);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("timing_statistics", "yes");
    app->Options()->SetStringValue("output_file", "ipopt.out");

    ApplicationReturnStatus status = app->Initialize();

    status = app->OptimizeTNLP(DOP);

    return status;
}
