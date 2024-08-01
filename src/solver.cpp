#include <chrono>
#include "solver.h"
#include "IpIpoptApplication.hpp"

Solver::Solver(const SmartPtr<GDOP>& gdop, const int maxMeshIterations, LinearSolver linearSolver, MeshAlgorithm meshAlgorithm) :
               gdop(gdop), maxMeshIterations(maxMeshIterations), linearSolver(linearSolver), meshAlgorithm(meshAlgorithm){}

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
    initSolvingProcess();
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    // numeric jacobian and hessian
    // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    // app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");

    // test derivatives
    // app->Options()->SetStringValue("derivative_test", "second-order");

    app->Options()->SetNumericValue("tol", tolerance);
    app->Options()->SetNumericValue("acceptable_tol", tolerance * 1e3);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetIntegerValue("max_iter", 100000);

    app->Options()->SetIntegerValue("print_level", 5);
    app->Options()->SetStringValue("timing_statistics", "yes");

    app->Options()->SetStringValue("linear_solver", getLinearSolverName(linearSolver));
    app->Options()->SetStringValue("hsllib", "/home/linus/masterarbeit/ThirdParty-HSL/.libs/libcoinhsl.so.2.2.5");

    app->Options()->SetStringValue("output_file", "ipopt.out");

    ApplicationReturnStatus status = app->Initialize();

    // initial optimization
    initialIntervals = gdop->mesh.intervals;
    status = app->OptimizeTNLP(gdop);
    postOptimization();
    while (meshIteration <= maxMeshIterations) {
        auto markedIntervals = detect();

        if (sz(markedIntervals) == 0) {
            finalizeOptimization();
            return status;
        }

        // interpolate x and u, update the mesh
        refine(markedIntervals);

        // small overhead because of the reinitializing of the GDOP, but thus all variables are initialized correctly
        // and invariant constant during each optimization remain const
        gdop = new GDOP(gdop->problem, gdop->mesh, gdop->rk, InitVars::CALLBACK);

        // set new starting values
        gdop->x_cb = cbValues;

        // optimize again
        status = app->OptimizeTNLP(gdop);
        postOptimization();
    }

    finalizeOptimization();
    return status;
}

std::vector<int> Solver::detect() {
    switch (meshAlgorithm) {
        case MeshAlgorithm::NONE: return {};
        case MeshAlgorithm::BASIC: return basicStrategy(2.5);
        case MeshAlgorithm::L2_BOUNDARY_NORM: return l2BoundaryNorm();
        default: return {};
    }
}

std::vector<int> Solver::l2BoundaryNorm() const {
    std::vector<int> markedIntervals = {};

    for (int i = 0; i < gdop->mesh.intervals; i++) {
        for (int u = 0; u < gdop->problem->sizeU; u++) {
            if (i > 0) {
                // TODO: i = 0 -> interpolate with deg 1 less, clean up integrator, last interval new interval
                // -> mark based on tolerances

                std::vector<double> uCoeffs;
                for (int j = -1; j < gdop->rk.steps; j++) {
                    uCoeffs.push_back(gdop->optimum[u + gdop->offX + i * gdop->offXUBlock + j * gdop->offXU]);
                }

                // values of the (1st, 2nd) diff of the interpolating polynomial at 0, c1, c2, ...
                std::vector<double> p_uDiff = gdop->rk.evalLagrangeDiff(uCoeffs);
                std::vector<double> p_uDiff2 = gdop->rk.evalLagrangeDiff2(uCoeffs);

                // squared values of the (1st, 2nd) diff of the interpolating polynomial at c1, c2, ...
                std::vector<double> sq_p_uDiff;
                std::vector<double> sq_p_uDiff2;
                for (int k = 1; k < sz(p_uDiff); k++) {
                    sq_p_uDiff.push_back(p_uDiff[k] * p_uDiff[k]);
                    sq_p_uDiff2.push_back(p_uDiff2[k] * p_uDiff2[k]);
                }

                // int_0^1 d^{1,2}/dt^{1,2} p_u(t)^2 dt - L2 norm of the (1st, 2nd) diff
                double L2Diff1 = gdop->rk.integrate(sq_p_uDiff);
                double L2Diff2 = gdop->rk.integrate(sq_p_uDiff2);
            }

        }

    }
    return markedIntervals;
}

std::vector<int> Solver::basicStrategy(const double sigma) const {
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
            index++;    // go to next marked interval
        }
        // not marked interval: copy optimal values
        else {
            for (int v = 0; v < gdop->offXU; v++) {
                for (int k = 0; k < gdop->rk.steps; k++) {
                    cbValues[v + (i + index) * gdop->offXUBlock + k * gdop->offXU] =
                    gdop->optimum[v + i * gdop->offXUBlock + k * gdop->offXU];
                }
            }
        }
    }
    for (int p = 0; p < gdop->offP; p++) {
        cbValues[newOffXUTotal + p] = gdop->optimum[gdop->offXUTotal + p];
    }
}

void Solver::setTolerance(double d) {
    tolerance = d;
}

void Solver::setExportOptimumPath(const std::string& exportPath) {
    this->exportOptimumPath = exportPath;
    this->exportOptimum = true;
}

void Solver::initSolvingProcess() {
    solveStartTime = std::chrono::high_resolution_clock::now();
}

void Solver::postOptimization() {
    objectiveHistory.push_back(gdop->objective);
    if (exportOptimum) {
        gdop->exportOptimum(exportOptimumPath + "/" + gdop->problem->name + std::to_string(meshIteration) + ".csv");
    }
    meshIteration++;
}

void Solver::printObjectiveHistory() {
    // TODO: Output with scientific notation
    std::cout << "\nMesh refinement history\n" << std::endl;
    std::cout << std::setw(5) << "iteration" << std::setw(17) << "objective" << std::endl;
    std::cout << "--------------------------------" << std::endl;
    for (int iteration = 0; iteration < sz(objectiveHistory); iteration++) {
        std::cout << std::setw(5) << iteration << std::setw(27) << double2Str(objectiveHistory[iteration]) << std::endl;
    }
}
void Solver::finalizeOptimization() {
    if (maxMeshIterations > 0) {
        printMeshStats();
        printObjectiveHistory();
    }
    auto timeTaken = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - solveStartTime).count();
    std::cout << "\nSolving took: " << timeTaken << " seconds" << std::endl;
}

void Solver::printMeshStats() {
    std::cout << "\nNumber of intervals\n" << std::endl;
    std::cout << "Initial:" << std::setw(7) << initialIntervals << std::endl;
    std::cout << "Inserted:  "<< std::setw(4) << gdop->mesh.intervals - initialIntervals << std::endl;
    std::cout << "Final:" << std::setw(9) << gdop->mesh.intervals << std::endl;
}
