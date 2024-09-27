#include "solver.h"

#include <chrono>

#include "IpIpoptApplication.hpp"
#include "gdop_impl.h"

/* BIG TODOS: TYPE | IMPORTANCE | EFFORT from 1 to 5
 * TODO:

    main todos:
    1 add, construct, test more mesh refinement algorithms!              4, 4
    2 OpenModelica interface                                             4, 5
    3 test framework for huge examples / industry relevant               3, 2
    4 play with setting in ipopt / pivoting etc. -> basically done       2, 3
    5 check long double to double cast in evals?!                        1, 2

    delayed:
    6 tf as free variable                                                2, 4.5
    7 vectorized equation, derivatives with local substitutions          2, 4
    -> define vec(f,g), vec(r), vec(a(p))
    8 better initial guess evolutionary algorithms                       1, 3
    9 detection for nominal, linear, quadratic, const hessian            1, 2
    10 constructing a p / hp-method?                                     3, 5

    abandoned:
    10 saving of local hessian and jacobian structures (contained in 4?) 0, 1
    11 creation of local jacobian structure (contained in 4?)            0, 1

    others:
    12 plotting features for path constraints, lagrange terms            1, 1
    13 splitting const jacobian equality / inequality                    1, 1
*/

struct SolverPrivate {
    SmartPtr<GDOP> gdop;
};

Solver::Solver(GDOP* gdop, const int maxMeshIterations, LinearSolver linearSolver, MeshAlgorithm meshAlgorithm)
    : maxMeshIterations(maxMeshIterations), linearSolver(linearSolver), meshAlgorithm(meshAlgorithm) {
    this->_priv = std::make_unique<SolverPrivate>();
    this->_priv->gdop = gdop;
}

Solver::~Solver() = default;

std::string getLinearSolverName(LinearSolver solver) {
    switch (solver) {
        case LinearSolver::MUMPS:
            return "MUMPS";
        case LinearSolver::MA27:
            return "MA27";
        case LinearSolver::MA57:
            return "MA57";
        case LinearSolver::MA77:
            return "MA77";
        case LinearSolver::MA86:
            return "MA86";
        case LinearSolver::MA97:
            return "MA97";
        case LinearSolver::PARDISO:
            return "pardisomkl";
        default:
            return "MUMPS";
    }
}

int Solver::solve() {
    // init solving, pre optimization
    initSolvingProcess();

    // create IPOPT application, set linear solver, tolerances, etc.
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    setSolverFlags(*app);
    ApplicationReturnStatus status = app->Initialize();

    // initial optimization
    status = app->OptimizeTNLP(_priv->gdop);
    postOptimization();

    while (meshIteration <= maxMeshIterations) {
        // detect intervals that have to be refined
        auto markedIntervals = detect();

        if (sz(markedIntervals) == 0) {
            finalizeOptimization();
            return status;
        }

        // interpolate x and u, update the mesh
        refine(markedIntervals);

        // small overhead because of the reinitializing of the GDOP, but thus all
        // variables are initialized correctly and invariant constants during each
        // optimization remain constant
        _priv->gdop = new GDOP(_priv->gdop->problem, _priv->gdop->mesh, _priv->gdop->rk, InitVars::CALLBACK);

        // set new starting values
        _priv->gdop->x_cb = cbValues;

        // update solver flags
        setSolverFlags(*app);

        // optimize again
        status = app->OptimizeTNLP(_priv->gdop);
        postOptimization();
    }

    finalizeOptimization();
    return status;
}

std::vector<int> Solver::detect() const {
    switch (meshAlgorithm) {
        case MeshAlgorithm::NONE:
            return {};
        case MeshAlgorithm::BASIC:
            return basicStrategy();
        case MeshAlgorithm::L2_BOUNDARY_NORM:
            return L2BoundaryNorm();
        default:
            return {};
    }
}

void Solver::setRefinementParameters() {
    switch (meshAlgorithm) {
        case MeshAlgorithm::NONE:
            return;
        case MeshAlgorithm::BASIC:
            return setBasicStrategy();
        case MeshAlgorithm::L2_BOUNDARY_NORM:
            return setL2BoundaryNorm();
        default:
            return;
    }
}

std::vector<int> Solver::L2BoundaryNorm() const {
    const double cDiff = 1;
    const double cDiff2 = cDiff / 2;
    std::set<int> markerSet;

    // init last derivatives u^(d)_{i-1,m} as 0
    std::vector<std::vector<double>> lastDiffs;
    lastDiffs.reserve(_priv->gdop->problem->sizeU);
    for (int u = 0; u < _priv->gdop->problem->sizeU; u++) {
        lastDiffs.push_back({0, 0});
    }

    // calculate max, min of u^(d)
    std::vector<double> maxU;
    std::vector<double> minU;
    maxU.reserve(_priv->gdop->problem->sizeU);
    minU.reserve(_priv->gdop->problem->sizeU);
    for (int i = 0; i < _priv->gdop->mesh.intervals; i++) {
        for (int j = 0; j < _priv->gdop->rk.steps; j++) {
            for (int u = 0; u < _priv->gdop->problem->sizeU; u++) {
                if (i == 0 && j == 0) {
                    maxU.push_back(_priv->gdop->optimum[u + _priv->gdop->offX]);
                    minU.push_back(_priv->gdop->optimum[u + _priv->gdop->offX]);
                }
                else {
                    if (_priv->gdop->optimum[u + _priv->gdop->offX + i * _priv->gdop->offXUBlock + j * _priv->gdop->offXU] > maxU[u]) {
                        maxU[u] = _priv->gdop->optimum[u + _priv->gdop->offX + i * _priv->gdop->offXUBlock + j * _priv->gdop->offXU];
                    }
                    else if (_priv->gdop->optimum[u + _priv->gdop->offX + i * _priv->gdop->offXUBlock + j * _priv->gdop->offXU] < minU[u]) {
                        minU[u] = _priv->gdop->optimum[u + _priv->gdop->offX + i * _priv->gdop->offXUBlock + j * _priv->gdop->offXU];
                    }
                }
            }
        }
    }
    std::vector<double> rangeU;
    rangeU.reserve(_priv->gdop->problem->sizeU);
    for (int u = 0; u < _priv->gdop->problem->sizeU; u++) {
        rangeU.push_back(maxU[u] - minU[u]);
    }

    std::vector<double> boundsDiff;
    std::vector<double> boundsDiff2;
    for (int u = 0; u < _priv->gdop->problem->sizeU; u++) {
        boundsDiff.push_back(cDiff * rangeU[u] / initialIntervals * pow(10, -L2Level));
        boundsDiff2.push_back(cDiff2 * rangeU[u] / initialIntervals * pow(10, -L2Level));
    }

    for (int i = 0; i < _priv->gdop->mesh.intervals; i++) {
        bool cornerTrigger = false;
        for (int u = 0; u < _priv->gdop->problem->sizeU; u++) {
            std::vector<double> uCoeffs;
            if (i == 0) {
                for (int j = 0; j < _priv->gdop->rk.steps; j++) {
                    uCoeffs.push_back(_priv->gdop->optimum[u + _priv->gdop->offX + i * _priv->gdop->offXUBlock + j * _priv->gdop->offXU]);
                }
                uCoeffs.insert(uCoeffs.begin(), Integrator::evalLagrange(_priv->gdop->rk.c, uCoeffs, 0.0));
            }
            else {
                for (int j = -1; j < _priv->gdop->rk.steps; j++) {
                    uCoeffs.push_back(_priv->gdop->optimum[u + _priv->gdop->offX + i * _priv->gdop->offXUBlock + j * _priv->gdop->offXU]);
                }
            }

            // values of the (1st, 2nd) diff of the interpolating polynomial at 0, c1, c2, ...
            std::vector<double> p_uDiff = _priv->gdop->rk.evalLagrangeDiff(uCoeffs);
            std::vector<double> p_uDiff2 = _priv->gdop->rk.evalLagrangeDiff2(uCoeffs);

            // squared values of the (1st, 2nd) diff of the interpolating polynomial at c1, c2, ...
            std::vector<double> sq_p_uDiff;
            std::vector<double> sq_p_uDiff2;
            for (int k = 1; k < sz(p_uDiff); k++) {
                sq_p_uDiff.push_back(p_uDiff[k] * p_uDiff[k]);
                sq_p_uDiff2.push_back(p_uDiff2[k] * p_uDiff2[k]);
            }

            // (int_0^1 (d^{1,2}/dt^{1,2} p_u(t))^2 dt)^0.5 - L2 norm of the (1st, 2nd) diff
            double L2Diff1 = std::sqrt(_priv->gdop->rk.integrate(sq_p_uDiff));
            double L2Diff2 = std::sqrt(_priv->gdop->rk.integrate(sq_p_uDiff2));
            if (i > 0) {
                // difference in derivatives from polynomial of adjacent intervals
                // must not exceed some eps using p1 (+1) error; basically
                // isclose(.) in numpy bib
                double p1ErrorDiff = std::abs(p_uDiff[0] - lastDiffs[u][0]) / (1 + std::max({std::abs(p_uDiff[0]), std::abs(lastDiffs[u][0])}));
                double p1ErrorDiff2 = std::abs(p_uDiff2[0] - lastDiffs[u][1]) / (1 + std::max({std::abs(p_uDiff2[0]), std::abs(lastDiffs[u][1])}));

                if (p1ErrorDiff > L2CornerTol || p1ErrorDiff2 > L2CornerTol) {
                    cornerTrigger = true;
                }
            }
            lastDiffs[u] = {p_uDiff[sz(_priv->gdop->rk.c)], p_uDiff2[sz(_priv->gdop->rk.c)]};

            // detection which intervals should be bisected
            if (L2Diff1 > boundsDiff[u] || L2Diff2 > boundsDiff2[u] || cornerTrigger) {
                // splitting the interval itself
                markerSet.insert(i);

                // on interval / L2 criterion -> forces adjacent intervals to be split
                // as well
                if (L2Diff1 > boundsDiff[u] || L2Diff2 > boundsDiff2[u]) {
                    if (i >= 1) {
                        markerSet.insert(i - 1);
                    }
                    if (i <= _priv->gdop->mesh.intervals - 2) {
                        markerSet.insert(i + 1);
                    }
                }

                // corner criterion (only exists for i > 0) -> forces left adjacent
                // interval split
                if (cornerTrigger) {
                    markerSet.insert(i - 1);
                }

                break;
            }
        }
    }
    return {markerSet.begin(), markerSet.end()};
}

std::vector<int> Solver::basicStrategy() const {
    std::vector<int> markedIntervals = {};
    std::vector<std::vector<double>> absIntSum = {};
    std::vector<double> means;
    std::vector<double> stds;

    for (int u = 0; u < _priv->gdop->problem->sizeU; u++) {
        absIntSum.emplace_back();
        for (int i = 0; i < _priv->gdop->mesh.intervals; i++) {
            double sum = 0;
            for (int j = -1; j < _priv->gdop->rk.steps; j++) {
                if (!((i == 0 && j == -1) || (i == _priv->gdop->mesh.intervals - 1 && j == _priv->gdop->rk.steps - 1))) {
                    double u1 = _priv->gdop->optimum[u + _priv->gdop->offX + i * _priv->gdop->offXUBlock + (j + 1) * _priv->gdop->offXU];
                    double u2 = _priv->gdop->optimum[u + _priv->gdop->offX + i * _priv->gdop->offXUBlock + j * _priv->gdop->offXU];
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
        std::vector<double> absSum95 = std::vector<double>(absIntSumCopy.begin(), absIntSumCopy.begin() + int(0.95 * sz(absIntSumCopy)));
        const double mean = calculateMean(absSum95);
        means.push_back(mean);
        stds.push_back(calculateStdDev(absSum95, mean));
    }

    for (int i = 0; i < _priv->gdop->mesh.intervals; i++) {
        bool containsInterval = false;
        for (int u = 0; u < _priv->gdop->problem->sizeU; u++) {
            if (absIntSum[u][i] > means[u] + basicStrategySigma * stds[u] || containsInterval) {
                markedIntervals.push_back(i);
                containsInterval = true;
            }
        }
    }
    return markedIntervals;
}

void Solver::refine(std::vector<int>& markedIntervals) {
    switch (refinementMethod) {
        case RefinementMethod::LINEAR_SPLINE:
            refineLinear(markedIntervals);
            break;

        case RefinementMethod::POLYNOMIAL:
            refinePolynomial(markedIntervals);
            break;

        default:
            refineLinear(markedIntervals);
            break;
    }
}

void Solver::refineLinear(std::vector<int>& markedIntervals) {
    // does a linear spline on each marked interval -> new values of control, states
    const int oldIntervalLen = _priv->gdop->mesh.intervals;
    _priv->gdop->mesh.update(markedIntervals);  // create new mesh here
    int newOffXUTotal = (_priv->gdop->problem->sizeX + _priv->gdop->problem->sizeU) * _priv->gdop->rk.steps * _priv->gdop->mesh.intervals;
    int newNumberVars = newOffXUTotal + _priv->gdop->problem->sizeP;
    cbValues.resize(newNumberVars, 0.0);

    // interpolate all values on marked intervals
    int index = 0;
    for (int i = 0; i < oldIntervalLen; i++) {
        if (markedIntervals[index] == i) {
            for (int v = 0; v < _priv->gdop->offXU; v++) {  // iterate over every var in {x, u} -> interpolate
                std::vector<double> localVars = {};
                // i > 0 interval cases
                if (i > 0) {
                    for (int j = -1; j < _priv->gdop->rk.steps; j++) {
                        localVars.push_back(_priv->gdop->optimum[v + i * _priv->gdop->offXUBlock + j * _priv->gdop->offXU]);
                    }
                    auto const splineVals = _priv->gdop->rk.evalLinearSplineNewNodes(localVars);
                    for (int k = 0; k < sz(splineVals); k++) {
                        cbValues[v + (i + index) * _priv->gdop->offXUBlock + k * _priv->gdop->offXU] = splineVals[k];
                    }
                }
                else {
                    // 0-th interval cases
                    if (v < _priv->gdop->offX) {
                        localVars.push_back(_priv->gdop->problem->x0[v]);
                        for (int j = 0; j < _priv->gdop->rk.steps; j++) {
                            localVars.push_back(_priv->gdop->optimum[v + j * _priv->gdop->offXU]);
                        }
                        auto const splineVals = _priv->gdop->rk.evalLinearSplineNewNodes(localVars);
                        for (int k = 0; k < sz(splineVals); k++) {
                            cbValues[v + i * _priv->gdop->offXUBlock + k * _priv->gdop->offXU] = splineVals[k];
                        }
                    }
                    else {
                        // 0-th control -> interpolate to get the value at t=0, linear splines
                        for (int j = 0; j < _priv->gdop->rk.steps; j++) {
                            localVars.push_back(_priv->gdop->optimum[v + j * _priv->gdop->offXU]);
                        }
                        localVars.insert(localVars.begin(), Integrator::evalLagrange(_priv->gdop->rk.c, localVars, 0.0));
                        auto const splineVals = _priv->gdop->rk.evalInterpolationNewNodes(localVars);
                        for (int k = 0; k < sz(splineVals); k++) {
                            cbValues[v + i * _priv->gdop->offXUBlock + k * _priv->gdop->offXU] = splineVals[k];
                        }
                    }
                }
            }
            index++;  // go to next marked interval
        }
        // not marked interval: copy optimal values
        else {
            for (int v = 0; v < _priv->gdop->offXU; v++) {
                for (int k = 0; k < _priv->gdop->rk.steps; k++) {
                    cbValues[v + (i + index) * _priv->gdop->offXUBlock + k * _priv->gdop->offXU] =
                        _priv->gdop->optimum[v + i * _priv->gdop->offXUBlock + k * _priv->gdop->offXU];
                }
            }
        }
    }
    for (int p = 0; p < _priv->gdop->offP; p++) {
        cbValues[newOffXUTotal + p] = _priv->gdop->optimum[_priv->gdop->offXUTotal + p];
    }
}

void Solver::refinePolynomial(std::vector<int>& markedIntervals) {
    const int oldIntervalLen = _priv->gdop->mesh.intervals;
    _priv->gdop->mesh.update(markedIntervals);  // create new mesh here
    int newOffXUTotal = (_priv->gdop->problem->sizeX + _priv->gdop->problem->sizeU) * _priv->gdop->rk.steps * _priv->gdop->mesh.intervals;
    int newNumberVars = newOffXUTotal + _priv->gdop->problem->sizeP;
    cbValues.resize(newNumberVars, 0.0);

    // interpolate all values on marked intervals
    int index = 0;
    for (int i = 0; i < oldIntervalLen; i++) {
        if (markedIntervals[index] == i) {
            for (int v = 0; v < _priv->gdop->offXU; v++) {  // iterate over every var in {x, u} -> interpolate
                std::vector<double> localVars = {};
                // i > 0 interval cases
                if (i > 0) {
                    for (int j = -1; j < _priv->gdop->rk.steps; j++) {
                        localVars.push_back(_priv->gdop->optimum[v + i * _priv->gdop->offXUBlock + j * _priv->gdop->offXU]);
                    }
                    auto const polyVals = _priv->gdop->rk.evalInterpolationNewNodes(localVars);
                    for (int k = 0; k < sz(polyVals); k++) {
                        cbValues[v + (i + index) * _priv->gdop->offXUBlock + k * _priv->gdop->offXU] = polyVals[k];
                    }
                }
                else {
                    // 0-th interval cases
                    if (v < _priv->gdop->offX) {
                        localVars.push_back(_priv->gdop->problem->x0[v]);
                        for (int j = 0; j < _priv->gdop->rk.steps; j++) {
                            localVars.push_back(_priv->gdop->optimum[v + j * _priv->gdop->offXU]);
                        }
                        auto const polyVals = _priv->gdop->rk.evalInterpolationNewNodes(localVars);
                        for (int k = 0; k < sz(polyVals); k++) {
                            cbValues[v + i * _priv->gdop->offXUBlock + k * _priv->gdop->offXU] = polyVals[k];
                        }
                    }
                    else {
                        // 0-th control -> interpolate with order one less (only
                        // not rk.steps points, not rk.steps + 1)
                        for (int j = 0; j < _priv->gdop->rk.steps; j++) {
                            localVars.push_back(_priv->gdop->optimum[v + j * _priv->gdop->offXU]);
                        }
                        auto const polyVals = _priv->gdop->rk.interpolateFirstControl(localVars);
                        for (int k = 0; k < sz(polyVals); k++) {
                            cbValues[v + i * _priv->gdop->offXUBlock + k * _priv->gdop->offXU] = polyVals[k];
                        }
                    }
                }
            }
            index++;  // go to next marked interval
        }
        // not marked interval: copy optimal values
        else {
            for (int v = 0; v < _priv->gdop->offXU; v++) {
                for (int k = 0; k < _priv->gdop->rk.steps; k++) {
                    cbValues[v + (i + index) * _priv->gdop->offXUBlock + k * _priv->gdop->offXU] =
                        _priv->gdop->optimum[v + i * _priv->gdop->offXUBlock + k * _priv->gdop->offXU];
                }
            }
        }
    }
    for (int p = 0; p < _priv->gdop->offP; p++) {
        cbValues[newOffXUTotal + p] = _priv->gdop->optimum[_priv->gdop->offXUTotal + p];
    }
}

void Solver::setTolerance(double tol) {
    tolerance = tol;
}

void Solver::setExportOptimumPath(const std::string& exportPath) {
    exportOptimumPath = exportPath;
    exportOptimum = true;
}

void Solver::setExportHessianPath(const std::string& exportPath) {
    exportHessianPath = exportPath;
    exportHessian = true;
}

void Solver::setExportJacobianPath(const std::string& exportPath) {
    exportJacobianPath = exportPath;
    exportJacobian = true;
}

void Solver::setRefinementMethod(RefinementMethod method) {
    refinementMethod = method;
}

void Solver::initSolvingProcess() {
    setRefinementParameters();
    solveStartTime = std::chrono::high_resolution_clock::now();
    initialIntervals = _priv->gdop->mesh.intervals;
    if (exportHessian) {
        _priv->gdop->exportHessianPath = exportHessianPath;
        _priv->gdop->exportHessian = true;
    }
    if (exportJacobian) {
        _priv->gdop->exportJacobianPath = exportJacobianPath;
        _priv->gdop->exportJacobian = true;
    }
}

void Solver::postOptimization() {
    auto io_start = std::chrono::high_resolution_clock::now();
    objectiveHistory.push_back(_priv->gdop->objective);
    if (exportOptimum) {
        _priv->gdop->exportOptimum(exportOptimumPath + "/" + _priv->gdop->problem->name + std::to_string(meshIteration) + ".csv");
    }
    meshIteration++;
    timedeltaIO += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - io_start);
}

void Solver::printObjectiveHistory() {
    std::cout << "\n----------------------------------------------------------------" << std::endl;
    std::cout << "\nMesh refinement history\n" << std::endl;
    std::cout << std::setw(5) << "iteration" << std::setw(20) << "objective" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    for (int iteration = 0; iteration < sz(objectiveHistory); iteration++) {
        std::cout << std::setw(5) << iteration << std::setw(30) << double2Str(objectiveHistory[iteration]) << std::endl;
    }
}

void Solver::createModelInfo() const {
    // writes a tmp file that contains all relevant information of the model like
    // the last mesh iteration and so on contains the metadata of the model
    std::ofstream outFile("/tmp/modelinfo.txt");
    if (!outFile) {
        std::cerr << "Error opening file for writing: "
                  << "/tmp/modelinfo.txt" << std::endl;
        return;
    }
    outFile << "maxMeshIteration, " << meshIteration - 1 << "\n";
    outFile << "totalTimeInSolver, " << solveTotalTimeTaken.count() << "\n";
    outFile << "actualTimeInSolver, " << solveActualTimeTaken.count() << "\n";
    outFile << "totalTimeInIO, " << timedeltaIO.count() << "\n";
    outFile.close();
}

void Solver::finalizeOptimization() {
    std::cout << "\n----------------------------------------------------------------" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;
    std::cout << "\nOutput for optimization of model: " << _priv->gdop->problem->name << std::endl;

    printMeshStats();
    printObjectiveHistory();

    solveTotalTimeTaken = (std::chrono::high_resolution_clock::now() - solveStartTime);
    solveActualTimeTaken = solveTotalTimeTaken - timedeltaIO;
    std::cout << std::fixed << std::setprecision(5);
    std::cout << "\n----------------------------------------------------------------" << std::endl;
    std::cout << "\nTotal time in Solver (w/o I/O): " << std::setw(8) << solveTotalTimeTaken.count() << " seconds" << std::endl;
    std::cout << "Total time in I/O: " << std::setw(21) << timedeltaIO.count() << " seconds" << std::endl;
    std::cout << "Total time in Solver: " << std::setw(18) << solveActualTimeTaken.count() << " seconds" << std::endl;
    std::cout << std::defaultfloat;

    createModelInfo();
}

void Solver::printMeshStats() const {
    std::cout << "\n----------------------------------------------------------------" << std::endl;
    std::cout << "\nNumber of intervals\n" << std::endl;
    std::cout << "Initial:" << std::setw(7) << initialIntervals << std::endl;
    std::cout << "Inserted:  " << std::setw(4) << _priv->gdop->mesh.intervals - initialIntervals << std::endl;
    std::cout << "Final:" << std::setw(9) << _priv->gdop->mesh.intervals << std::endl;
}

void Solver::setSolverFlags(IpoptApplication& app) {
    // numeric jacobian and hessian
    // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    // app->Options()->SetStringValue("jacobian_approximation", "finite-difference-values");

    // test derivatives
    // app->Options()->SetStringValue("derivative_test", "second-order");

    if (meshIteration == 0) {
        // flags for the initial optimization
        app.Options()->SetNumericValue("bound_push", 1e-2);
        app.Options()->SetNumericValue("bound_frac", 1e-2);
    }
    else {
        // flags for every following refined optimization
        app.Options()->SetNumericValue("bound_push", 1e-4);
        app.Options()->SetNumericValue("bound_frac", 1e-4);
    }

    // setting the standard flags
    setStandardSolverFlags(app);
}

void Solver::setStandardSolverFlags(IpoptApplication& app) {
    // these are flags are always set, no matter if a mesh refinement is currently executed

    // mu-update strategy
    // turns out kkt-error adaptive_mu works really well for the provided examples
    app.Options()->SetStringValue("mu_strategy", "adaptive");
    app.Options()->SetStringValue("adaptive_mu_globalization", "kkt-error");

    // iterations and tolereances
    app.Options()->SetNumericValue("tol", tolerance);
    app.Options()->SetNumericValue("acceptable_tol", tolerance * 1e3);
    app.Options()->SetIntegerValue("max_iter", 15000);

    // ipopt dump
    app.Options()->SetIntegerValue("print_level", 5);
    app.Options()->SetStringValue("timing_statistics", "yes");

    // linear solver
    auto const libHSLPath = getenv("LIB_HSL");
    auto const linSolver = getLinearSolverName(linearSolver);
    if (libHSLPath != nullptr) {
        // HSL found -> set chosen solver
        app.Options()->SetStringValue("hsllib", libHSLPath);
        app.Options()->SetStringValue("linear_solver", linSolver);
    }
    else if ((libHSLPath == nullptr and linSolver != "MUMPS")){ 
        // HSL not found but set -> set MUMPS as fallback
        std::cout << "\nEnvironment variable 'LIB_HSL' not found! Fallback to standard linear solver 'MUMPS'\n" << std::endl;
        app.Options()->SetStringValue("linear_solver", "MUMPS");
    }
    else {
        // set chosen solver
        app.Options()->SetStringValue("linear_solver", linSolver);
    }
 
    // scaling
    if (userScaling) {
        app.Options()->SetStringValue("nlp_scaling_method", "user-scaling");
    }
    else {
        app.Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    }

    // constant derivatives reduce the number of function evals
    if (_priv->gdop->problem->quadraticObjLinearConstraints) {
        app.Options()->SetStringValue("hessian_constant", "yes");
    }
    if (_priv->gdop->problem->linearObjective) {
        app.Options()->SetStringValue("grad_f_constant", "yes");
    }
    if (_priv->gdop->problem->linearConstraints) {
        app.Options()->SetStringValue("jac_c_constant", "yes");
        app.Options()->SetStringValue("jac_d_constant", "yes");
    }
}

void Solver::setMeshParameter(const std::string& field, double value) {
    meshParameters.emplace(field, value);
}

void Solver::setL2BoundaryNorm() {
    if (meshParameters.count("level") > 0)
        L2Level = meshParameters["level"];
    if (meshParameters.count("ctol") > 0)
        L2CornerTol = meshParameters["ctol"];
}

void Solver::setBasicStrategy() {
    if (meshParameters.count("sigma") > 0)
        basicStrategySigma = meshParameters["sigma"];
}
