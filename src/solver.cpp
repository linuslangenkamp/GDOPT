#include "solver.h"
#include "IpIpoptApplication.hpp"

// TODO: Output solution to csv, refine

Solver::Solver(const SmartPtr<GDOP>& gdop, const int maxMeshIterations, LinearSolver linearSolver) :
               gdop(gdop), maxMeshIterations(maxMeshIterations), linearSolver(linearSolver){}

std::string getLinearSolverName(LinearSolver solver) {
    switch (solver) {
        // TODO: add pardiso project sparse solver
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

int Solver::solve() const {
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

    int iteration = 0;
    bool stopRefinement = false;
    while (iteration <= maxMeshIterations || stopRefinement){
        status = app->OptimizeTNLP(gdop);
        const double sigma = 2.5;
        auto intervals = basicStochasticStrategy(sigma);
        // interpolate x and u
        // rebuild gdop -> run OptimizeTNLP again
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
                    }
                    else {
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