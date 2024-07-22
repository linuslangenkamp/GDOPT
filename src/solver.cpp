#include "solver.h"
#include "IpIpoptApplication.hpp"

// TODO: Output solution to csv, refine

Solver::Solver(const SmartPtr<GDOP>& gdop, const int maxMeshIterations) : gdop(gdop),
                                                                            maxMeshIterations(maxMeshIterations) {}

int Solver::solve() const {
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    // numeric jacobian and hessian
    // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    // app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");

    app->Options()->SetNumericValue("tol", tolerance);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("timing_statistics", "yes");
    app->Options()->SetStringValue("linear_solver", "ma97");
    app->Options()->SetStringValue("hsllib", "/usr/local/lib/libcoinhsl.so.2.2.5");

    app->Options()->SetStringValue("output_file", "ipopt.out");

    ApplicationReturnStatus status = app->Initialize();

    int iteration = 0;
    bool stopRefinement = false;
    while (iteration <= maxMeshIterations || stopRefinement){
        status = app->OptimizeTNLP(gdop);
        // get intervals that ve to be refined -> returns interval indices (i)
        // interpolate x and u
        // rebuild gdop -> run OptimizeTNLP again
        iteration++;
    }

    return status;
}
