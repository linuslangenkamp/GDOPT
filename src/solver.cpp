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

    int meshIteration = 0;
    while (meshIteration < maxMeshIterations) {
        const double sigma = 2.5;
        auto markedIntervals = basicStochasticStrategy(sigma);
        if (sz(markedIntervals) == 0)
            return status;

        // interpolate x and u, update the mesh
        refine(markedIntervals);

        // small overhead because of the reinitializing of the GDOP, but thus all variables are initialized correctly
        // and invariant constant during each optimization remain const
        gdop = new GDOP(gdop->problem, gdop->mesh, gdop->rk, InitVars::CALLBACK);;
        gdop->x_cb = cbValues;

        // optimize again
        status = app->OptimizeTNLP(gdop);

        meshIteration++;
    }
    return status;
}

std::vector<int> Solver::basicStochasticStrategy(const double sigma) const {
    std::vector<int> markedIntervals = {};
    std::vector<std::vector<double>> absIntSum = {};
    std::vector<double> means;
    std::vector<double> stds;

    for (int u = 0; u < gdop->problem->sizeU; u++) {
        absIntSum.emplace_back();
        for (int i = 0; i < gdop->mesh.intervals; i++) {
            double sum = 0;
            for (int j = -1; j < gdop->rk.steps; j++) {
                if (!((i == 0 && j == -1) || (i == gdop->mesh.intervals - 1 && j == gdop->rk.steps - 1)))  {
                    double u1 = gdop->optimum.at(u + gdop->offX + i * gdop->offXUBlock + (j + 1) * gdop->offXU);
                    double u2 = gdop->optimum.at(u + gdop->offX + i * gdop->offXUBlock + j * gdop->offXU);
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
        for (int u = 0; u < gdop->problem->sizeU; u++) {
            if (absIntSum[u][i] > means[u] + sigma * stds[u] || containsInterval) {
                markedIntervals.push_back(i);
                containsInterval = true;
            }
        }
    }
    return markedIntervals;
}

void Solver::refine(std::vector<int> &markedIntervals) {
    const int oldIntervalLen = gdop->mesh.intervals;
    gdop->mesh.update(markedIntervals); // create newMesh here
    int newOffXUTotal = (gdop->problem->sizeX + gdop->problem->sizeU) * gdop->rk.steps * gdop->mesh.intervals;
    int newNumberVars = newOffXUTotal + gdop->problem->sizeP;
    cbValues.resize(newNumberVars, 0.0);

    // interpolate all values on marked intervals
    int index = 0;
    for (int i = 0; i < oldIntervalLen; i++) {
        if (markedIntervals[index] == i) {
            for (int v = 0; v < gdop->offXU; v++) {         // iterate over every var in {x, u} -> interpolate
                std::vector<double> localVars;
                // i > 0 interval cases
                if (i > 0) {
                    for (int j = -1; j < gdop->rk.steps; j++) {
                        localVars.push_back(gdop->optimum[v + i * gdop->offXUBlock + j * gdop->offXU]);
                    }
                    auto const polyVals = gdop->rk.interpolate(localVars);
                    for (int k = 0; k < sz(polyVals); k++) {
                        cbValues[v + (i + index) * gdop->offXUBlock + k * gdop->offXU] = polyVals[k];
                    }
                }
                else {
                    // 0-th interval cases
                    if (v < gdop->offX) {
                        localVars.push_back(gdop->problem->x0[v]);
                        for (int j = 0; j < gdop->rk.steps; j++) {
                            localVars.push_back(gdop->optimum[v + j * gdop->offXU]);
                            auto const polyVals = gdop->rk.interpolate(localVars);
                            for (int k = 0; k < sz(polyVals); k++) {
                                cbValues[v + i * gdop->offXUBlock + k * gdop->offXU] = polyVals[k];
                            }
                        }
                    }
                    else {
                        // 0-th control -> interpolate with order one less (only not rk.steps points, not rk.steps + 1)
                        for (int j = 0; j < gdop->rk.steps; j++) {
                            localVars.push_back(gdop->optimum[v + j * gdop->offXU]);
                            }
                        auto const polyVals = gdop->rk.interpolateFirstControl(localVars);
                        for (int k = 0; k < sz(polyVals); k++) {
                            cbValues[v + i * gdop->offXUBlock + k * gdop->offXU] = polyVals[k];
                        }
                    }
                }
            }
            index++;       // go to next marked interval
        }
        // not marked interval: copy optimal values
        else {
            for (int v = 0; v < gdop->offXU; v++) {
                for (int k = 0; k < gdop->rk.steps; k++) {
                    cbValues[v + (i + index) * gdop->offXUBlock + k * gdop->offXU] =
                    gdop->optimum.at(v + i * gdop->offXUBlock + k * gdop->offXU);
                }
            }
        }
    }
    for (int p = 0; p < gdop->offP; p++) {
        cbValues[newOffXUTotal + p] = gdop->optimum[gdop->offXUTotal + p];
    }
}
