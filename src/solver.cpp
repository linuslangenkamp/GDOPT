#include "solver.h"
#include "IpIpoptApplication.hpp"

// TODO: Output solution to csv, refine

Solver::Solver(const SmartPtr<GDOP>& gdop, const int maxMeshIterations, LinearSolver linearSolver) :
               gdop(gdop), maxMeshIterations(maxMeshIterations), linearSolver(linearSolver){}

std::string getLinearSolverName(LinearSolver solver) {
    switch (solver) {
        // TODO: add pardiso project sparse solver
        // TODO: add SPRAL solver (bsd licensed)
        case LinearSolver::MUMPS: return "MUMPS";
        case LinearSolver::MA27: return "MA27";
        case LinearSolver::MA57: return "MA57";
        case LinearSolver::MA77: return "MA77";
        case LinearSolver::MA86: return "MA86";
        case LinearSolver::MA97: return "MA97";
        case LinearSolver::PARDISO: return "pardisomkl";
        default: return "MUMPS";
    }
}

int Solver::solve() {
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    // numeric jacobian and hessian
    // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    // app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");

    // test derivatives
    // app->Options()->SetStringValue("derivative_test", "second-order");
    // app->Options()->SetNumericValue("derivative_test_tol", 1e-2);

    app->Options()->SetNumericValue("tol", tolerance);
    app->Options()->SetNumericValue("acceptable_tol", tolerance * 1e3);
    app->Options()->SetStringValue("mu_strategy", "adaptive");

    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("timing_statistics", "yes");

    app->Options()->SetStringValue("linear_solver", getLinearSolverName(linearSolver));
    app->Options()->SetStringValue("hsllib", "/home/linus/masterarbeit/ThirdParty-HSL/.libs/libcoinhsl.so.2.2.5");

    app->Options()->SetStringValue("output_file", "ipopt.out");

    ApplicationReturnStatus status = app->Initialize();
    status = app->OptimizeTNLP(gdop);

    int iteration = 0;
    while (iteration < maxMeshIterations) {
        const double sigma = 2.5;
        auto intervals = basicStochasticStrategy(sigma);
        if (sz(intervals) == 0)
            return status;

        // interpolate x and u, update the mesh
        refine(intervals);

        // small overhead because of the reinitializing of the GDOP, but thus all variables are initialized correctly
        // and invariant constant during each optimization remain const
        this->gdop = new GDOP(std::move(gdop->problem), gdop->mesh, gdop->rk, InitVars::CALLBACK);
        gdop->x_cb = cbValues;
        status = app->OptimizeTNLP(gdop);
        iteration++;
    }
    return status;
}

std::vector<int> Solver::basicStochasticStrategy(const double sigma) const {
    std::vector<int> intervals = {};
    std::vector<std::vector<double>> absIntSum = {};
    std::vector<double> means;
    std::vector<double> stds;

    for (int u = 0; u < gdop->problem.sizeU; u++) {
        absIntSum.emplace_back();
        for (int i = 0; i < gdop->mesh.intervals; i++) {
            double sum = 0;
            for (int j = -1; j < gdop->rk.steps; j++) {
                if (i != 0 || j != -1) {
                    double u1 = gdop->optimum[u + gdop->offX + i * gdop->offXUBlock + (j + 1) * gdop->offXU];
                    double u2 = gdop->optimum[u + gdop->offX + i * gdop->offXUBlock + j * gdop->offXU];
                    if (u1 > u2) {
                        sum += u1 - u2;
                    } else {
                        sum += u2 - u1;
                    }
                }
            }
            absIntSum[u].push_back(sum);
        }
        std::vector<double> absIntSumCopy = absIntSum[u];
        std::sort(absIntSumCopy.begin(), absIntSumCopy.end());
        std::vector<double> absSum95 = std::vector<double>(absIntSumCopy.begin(),
                                                           absIntSumCopy.begin() + int(0.95 * sz(absIntSumCopy)));
        const double mean = calculateMean(absSum95);
        means.push_back(mean);
        stds.push_back(calculateStdDev(absSum95, mean));
    }

    for (int i = 0; i < gdop->mesh.intervals; i++) {
        bool containsInterval = false;
        for (int u = 0; u < gdop->problem.sizeU; u++) {
            if (absIntSum[u][i] > means[u] + sigma * stds[u] || containsInterval) {
                intervals.push_back(i);
                containsInterval = true;
            }
        }
    }
    return intervals;
}

void Solver::refine(std::vector<int> &intervals) {
    gdop->mesh.update(intervals);
    int varCount = (gdop->problem.sizeX + gdop->problem.sizeU) * gdop->rk.steps * gdop->mesh.intervals +
                   gdop->problem.sizeP;
    cbValues.reserve(varCount);

    // INTERPOLATE / EXTRAPOLATE OTHER VALUES
    // TODO: init cbValues with range as 0 -> later set values
    int index = 0;
    for (int i = 0; i < gdop->mesh.intervals; i++) {
        if (intervals[index] == i) {
            for (int v = 0; v < gdop->offXU; v++) {
                std::vector<double> localVars;
                if (i > 0) {
                    /*
                    if (i == 0 && v < gdop->offX) {
                        localVars.push_back(gdop->problem.x0[v]);
                    }
                    if (i == 0) {
                        int jstart = 0;
                    }
                    else {
                        int jstart = -1;*/

                    for (int j = -1; j < gdop->rk.steps; j++) {
                        localVars.push_back(gdop->optimum[v + i * gdop->offXUBlock + j * gdop->offXU]);
                    }
                    auto const vals = gdop->rk.interpolate(localVars);
                }
            }
            index++;
        }
    }
    cbValues.push_back(0);
}
