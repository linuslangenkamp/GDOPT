#include <cassert>
#include <algorithm>
#include "gdop.h"
#include "gdop_impl.h"
#include "util.h"

// TODO?: recheck derivatives, indices

bool GDOP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
    // #vars
    n = numberVars;

    // #eqs
    m = sz(problem->A) + sz(problem->R) + (sz(problem->F) + sz(problem->G)) * rk.steps * mesh.intervals;

    // #nnz in jac
    init_jac(nnz_jac_g);

    // # nnz in hessian + creation of sparsity maps and offsets
    init_h(nnz_h_lag);

    // index starts at 0
    index_style = TNLP::C_STYLE;
    return true;
}

void GDOP::init_jac(Index &nnz_jac_g) {
    // eval jacobian nnz
    nnz_jac_g = 0;

    // nnz jac: dynamics
    int nnzDynBlock = 0;
    int containedIndex = 0; // target idx
    for (const auto& dyn : problem->F) {
        nnzDynBlock += sz(dyn->adj.indX) + sz(dyn->adj.indU) + sz(dyn->adj.indP);
        // check if the k-th component of x, that is always contained in this block equation,
        // is contained in grad f_k(.) as well => reduce block nnz by 1
        auto it = std::find(dyn->adj.indX.begin(), dyn->adj.indX.end(), containedIndex);
        if (it != dyn->adj.indX.end())
            nnzDynBlock -= 1;
        containedIndex++;
    }
    nnz_jac_g += mesh.intervals * rk.steps * nnzDynBlock;                // nnz for f(v_{ij}) eval
    nnz_jac_g += mesh.intervals * rk.steps * rk.steps * (sz(problem->F)); // nnz for sum_k ~a_{jk} * x_{ik}
    nnz_jac_g += (mesh.intervals - 1) * rk.steps * sz(problem->F);        // nnz for sum_k ~a_{jk} * (-x_{i-1,m}), i>=1

    // nnz jac: path constraints
    int nnzPathBlock = 0;
    for (const auto& pathConstr : problem->G) {
        nnzPathBlock += sz(pathConstr->adj.indX) + sz(pathConstr->adj.indU) + sz(pathConstr->adj.indP);
    }
    nnz_jac_g += mesh.intervals * rk.steps * nnzPathBlock; // RHS = nnz for g(v_{ij}) eval

    // nnz jac: final constraints
    for (const auto& finalConstr : problem->R) {
        nnz_jac_g += sz(finalConstr->adj.indX) + sz(finalConstr->adj.indU) + sz(finalConstr->adj.indP);
    }

    // nnz jac: algebraic parameter constraints
    for (const auto& paramConstr : problem->A) {
        nnz_jac_g += sz(paramConstr->adj.indP);
    }
}

void GDOP::init_h(Index &nnz_h_lag) {
    nnz_h_lag = 0;

    // dense local hessian adjacency structure, will be reduced to hashmap later
    std::vector<std::vector<int>> denseS0(offXU), denseS0t(offXU), denseS1(offP), denseS1t(offP), denseS2(offP);

    // init local hessian adjacency
    for (int i = 0; i < offXU; i++) {
        denseS0[i] = std::vector<int>(i + 1, 0);
        denseS0t[i] = std::vector<int>(i + 1, 0);
    }
    for (int i = 0; i < offP; i++) {
        denseS1[i] = std::vector<int>(offXU, 0);
        denseS1t[i] = std::vector<int>(offXU, 0);
        denseS2[i] = std::vector<int>(offP, 0);
    }

    // update lagrange, dynamics, path constraint hessian structure
    if (problem->L)
        updateDenseHessianLFG(*problem->L, denseS0, denseS0t, denseS1, denseS1t, denseS2);

    for (const auto& f : problem->F) {
        updateDenseHessianLFG(*f, denseS0, denseS0t, denseS1, denseS1t, denseS2);
    }

    for (const auto& g : problem->G) {
        updateDenseHessianLFG(*g, denseS0, denseS0t, denseS1, denseS1t, denseS2);
    }

    // update mayer, final constraint hessian structure
    if (problem->M)
        updateDenseHessianMR(*problem->M,denseS0t, denseS1t, denseS2);

    for (const auto& r : problem->R) {
        updateDenseHessianMR(*r, denseS0t, denseS1t, denseS2);
    }

    // update algebraic parameter constraint hessian structure
    for (const auto& a : problem->A) {
        updateDenseHessianA(*a, denseS2);
    }

    // at this point all local hessian structures have been obtained
    nnz_h_lag = (fullSum(denseS0) + fullSum(denseS1)) * (mesh.intervals * rk.steps - 1) +
                 fullSum(denseS0t) + fullSum(denseS1t) + fullSum(denseS2);

    createSparseHessian(denseS0, denseS0t, denseS1, denseS1t, denseS2, nnz_h_lag);
}

void GDOP::createSparseHessian(std::vector<std::vector<int>>& denseS0,
                                std::vector<std::vector<int>>& denseS0t,
                                std::vector<std::vector<int>>& denseS1,
                                std::vector<std::vector<int>>& denseS1t,
                                std::vector<std::vector<int>>& denseS2,
                                Index &nnz_h_lag) {
    /**
     it: equation index in COO format
     S0: shift block 0,0 -> i,j: shift by lengthS0 * (i * rk.steps + j)
     S0t: is exact -> it = correct eq index
     S1: is not exact, it = starting index for p-th parameter:  rowLengthS1Block[p] * (i * rk.steps + j)
     S1t: is exact -> it = correct eq index
     S2: is exact -> it = correct eq index
    **/
    int it = 0; // eq index
    for (int i = 0; i < offXU; i++) {
        for (int j = 0; j <= i; j++) {
            if (denseS0[i][j] == 1) {
                S0.insert({{i, j}, it});        // "it" not exact; shift needed based on
                it++;
            }
        }
    }
    lengthS0 = it;
    it *= (mesh.intervals * rk.steps - 1);

    for (int i = 0; i < offXU; i++) {
        for (int j = 0; j <= i; j++) {
            if (denseS0t[i][j] == 1) {
                S0t.insert({{i, j}, it});      // "it" exact; no shift needed later on
                it++;
            }
        }
    }

    for (int i = 0; i < offP; i++) {
        int nnzS1row = 0;
        for (int j = 0; j < offXU; j++) {
            if (denseS1[i][j] == 1) {
                S1.insert({{i, j}, it});       // "it" not exact; shift needed based on row
                it++;
                nnzS1row++;
            }
        }
        rowLengthS1Block.push_back(nnzS1row);
        it += nnzS1row * (mesh.intervals * rk.steps - 2);   // nnzS1row * (mesh.intervals * rk.steps - 1) - nnzS1row

        for (int j = 0; j < offXU; j++) {
            if (denseS1t[i][j] == 1) {
                S1t.insert({{i, j}, it});   // "it" exact; no shift needed later on
                it++;
            }
        }

        for (int j = 0; j < offP; j++) {
            if (denseS2[i][j] == 1) {
                S2.insert({{i, j}, it});    // "it" exact; no shift needed later on
                it++;
            }
        }
    }
    assert(it == nnz_h_lag);
}

void GDOP::updateDenseHessianLFG(const Expression& expr,
                                    std::vector<std::vector<int>>& denseS0,
                                    std::vector<std::vector<int>>& denseS0t,
                                    std::vector<std::vector<int>>& denseS1,
                                    std::vector<std::vector<int>>& denseS1t,
                                    std::vector<std::vector<int>>& denseS2) const {
    // update of Hessian Blocks for Lagrange-term, dynamic and path constraints
    for (auto const [x1, x2] : expr.adjDiff.indXX) {
        denseS0[x1][x2] = 1;
        denseS0t[x1][x2] = 1;
    }
    for (auto const [u, x] : expr.adjDiff.indUX) {
        denseS0[offX + u][x] = 1;
        denseS0t[offX + u][x] = 1;
    }
    for (auto const [u1, u2] : expr.adjDiff.indUU) {
        denseS0[offX + u1][offX + u2] = 1;
        denseS0t[offX + u1][offX + u2] = 1;
    }
    for (auto const [p, x] : expr.adjDiff.indPX) {
        denseS1[p][x] = 1;
        denseS1t[p][x] = 1;
    }
    for (auto const [p, u] : expr.adjDiff.indPU) {
        denseS1[p][offX + u] = 1;
        denseS1t[p][offX + u] = 1;
    }
    for (auto const [p1, p2] : expr.adjDiff.indPP) {
        denseS2[p1][p2] = 1;
    }
}

void GDOP::updateDenseHessianMR(const Expression& expr,
                                std::vector<std::vector<int>>& denseS0t,
                                std::vector<std::vector<int>>& denseS1t,
                                std::vector<std::vector<int>>& denseS2) const {
    // update of Hessian Blocks for Mayer-term and final constraints
    for (auto const [x1, x2] : expr.adjDiff.indXX) {
        denseS0t[x1][x2] = 1;
    }
    for (auto const [u, x] : expr.adjDiff.indUX) {
        denseS0t[offX + u][x] = 1;
    }
    for (auto const [u1, u2] : expr.adjDiff.indUU) {
        denseS0t[offX + u1][offX + u2] = 1;
    }
    for (auto const [p, x] : expr.adjDiff.indPX) {
        denseS1t[p][x] = 1;
    }
    for (auto const [p, u] : expr.adjDiff.indPU) {
        denseS1t[p][offX + u] = 1;
    }
    for (auto const [p1, p2] : expr.adjDiff.indPP) {
        denseS2[p1][p2] = 1;
    }
}

void GDOP::updateDenseHessianA(const ParamExpression& expr, std::vector<std::vector<int>>& denseS2) {
    // update of Hessian Blocks for parametric constraints
    for (auto const [p1, p2] : expr.adjDiff.indPP) {
        denseS2[p1][p2] = 1;
    }
}

bool GDOP::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) {
    assert(n == numberVars);
    assert(m == sz(problem->A) + sz(problem->R) + (sz(problem->F) + sz(problem->G)) * rk.steps * mesh.intervals);

    // var bounds
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < rk.steps; j++) {

            const int xij = i * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
            const int uij = i * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)

            // state bounds
            for (int d = 0; d < problem->sizeX; d++) {
                x_l[xij + d] = problem->lbX[d];
                x_u[xij + d] = problem->ubX[d];
            }

            // control bounds
            for (int d = 0; d < problem->sizeU; d++) {
                x_l[uij + d] = problem->lbU[d];
                x_u[uij + d] = problem->ubU[d];
            }
        }
    }

    // parameter bounds
    for (int d = 0; d < problem->sizeP; d++) {
        x_l[offXUTotal + d] = problem->lbP[d];
        x_u[offXUTotal + d] = problem->ubP[d];
    }

    // equation bounds
    int eq = 0;
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < rk.steps; j++) {
            // dynamic bound must equal 0
            for (int d = 0; d < sz(problem->F); d++) {
                g_l[eq] = 0.0;
                g_u[eq] = 0.0;
                eq++;
            }

            // path constraint bounds
            for (auto const& pathConstr : problem->G) {
                g_l[eq] = pathConstr->lb;
                g_u[eq] = pathConstr->ub;
                eq++;
            }
        }
    }

    // final constraint bounds
    for (auto const& finalConstr : problem->R) {
        g_l[eq] = finalConstr->lb;
        g_u[eq] = finalConstr->ub;
        eq++;
    }

    // parameter constraint bounds
    for (auto const& paramConstr : problem->A) {
        g_l[eq] = paramConstr->lb;
        g_u[eq] = paramConstr->ub;
        eq++;
    }
    assert(m == eq);
    return true;
}

bool GDOP::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index m,
                              bool init_lambda, Number *lambda) {
    // TODO: implement different strategies:
    //  * SOLVE(.) -> given u(t), p guesses and x0: simulate the dynamic with given scheme, newton iteration, linear system
    //  * genetic / bionic algorithm
    assert(n == numberVars);
    if (init_x) {

        std::vector<std::vector<double>> initialStates;
        double tijOld;

        switch (initVars) {
            case InitVars::CONST:
                // states will be constant globally
                for (int i = 0; i < mesh.intervals; i++) {
                    for (int j = 0; j < rk.steps; j++) {
                        const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
                        const int xij = i * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
                        const int uij = i * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)

                        for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                            x[xij + dimX] = problem->x0[dimX];
                        }

                        std::vector<double> uGuess = problem->uInitialGuess(tij);
                        for (int dimU = 0; dimU < problem->sizeU; dimU++) {
                            x[uij + dimU] = uGuess[dimU];
                        }
                    }
                }
                for (int dimP = 0; dimP < problem-> sizeP; dimP++) {
                    x[offXUTotal + dimP] = problem->pInitialGuess[dimP];
                }

                break;

            case InitVars::SOLVE:
                // currently reading in a solution from the frontend (solved by scipy with Radau5 or BDF)

                initialStates = readInitialValues(problem->initialStatesPath);

                for (int i = 0; i < mesh.intervals; i++) {
                    for (int j = 0; j < rk.steps; j++) {
                        const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
                        const int xij = i * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
                        const int uij = i * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)

                        for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                            x[xij + dimX] = initialStates[rk.steps * i + j][dimX + 1];
                        }

                        std::vector<double> uGuess = problem->uInitialGuess(tij);
                        for (int dimU = 0; dimU < problem->sizeU; dimU++) {
                            x[uij + dimU] = uGuess[dimU];
                        }
                    }
                }

                for (int dimP = 0; dimP < problem-> sizeP; dimP++) {
                    x[offXUTotal + dimP] = problem->pInitialGuess[dimP];
                }

                break;

            case InitVars::SOLVE_EXPLICIT_EULER:
                // explicit Euler for initial guess of states (on each sub-interval [t_ij, t_{i, j+1}])
                for (int dimP = 0; dimP < problem-> sizeP; dimP++) {
                    x[offXUTotal + dimP] = problem->pInitialGuess[dimP];
                }

                tijOld = 0;

                for (int i = 0; i < mesh.intervals; i++) {
                    for (int j = 0; j < rk.steps; j++) {
                        const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
                        const double dt = tij - tijOld;
                        const int xij = i * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
                        const int uij = i * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)

                        if (i == 0 && j == 0) {
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                x[xij + dimX] = problem->x0[dimX];
                            }
                        }
                        else {
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                x[xij + dimX] = x[xij - offXU + dimX] + dt * problem->F[dimX]->eval(&x[xij - offXU], &x[uij - offXU], &x[offXUTotal], tijOld);
                            }

                        }

                        const std::vector<double> uGuess = problem->uInitialGuess(tij);
                        for (int dimU = 0; dimU < problem->sizeU; dimU++) {
                            x[uij + dimU] = uGuess[dimU];
                        }

                        tijOld = tij;
                    }
                }

                break;

            case InitVars::SOLVE_EXPLICIT:
                // classic Runge-Kutta scheme for initial guess of states (on each sub-interval [t_ij, t_{i, j+1}])
                for (int dimP = 0; dimP < problem-> sizeP; dimP++) {
                    x[offXUTotal + dimP] = problem->pInitialGuess[dimP];
                }

                tijOld = 0;

                for (int i = 0; i < mesh.intervals; i++) {
                    for (int j = 0; j < rk.steps; j++) {
                        const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
                        const double dt = tij - tijOld;
                        const int xij = i * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
                        const int uij = i * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)

                        if (i == 0 && j == 0) {
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                x[xij + dimX] = problem->x0[dimX];
                            }
                        }
                        else {
                            double k1X[problem->sizeX];
                            double k2X[problem->sizeX];
                            double k3X[problem->sizeX];
                            double k4X[problem->sizeX];

                            // k1
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                k1X[dimX] = problem->F[dimX]->eval(&x[xij - offXU], &x[uij - offXU], &x[offXUTotal], tijOld);
                            }

                            const double t2 = tijOld + 0.5 * dt;
                            std::vector<double> uGuess2 = problem->uInitialGuess(t2);

                            std::vector<double> xTemp(&x[xij - offXU], &x[xij - offXU] + problem->sizeX);

                            // x_i + 0.5 * dt * k1
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                xTemp[dimX] += 0.5 * dt * k1X[dimX];
                            }

                            // k2
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                k2X[dimX] = problem->F[dimX]->eval(xTemp.data(), uGuess2.data(), &x[offXUTotal], t2);
                            }

                            std::vector<double> uGuess3 = problem->uInitialGuess(t2);

                            // x_i + 0.5 * dt * k2
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                xTemp[dimX] = x[xij - offXU + dimX] + 0.5 * dt * k2X[dimX];
                            }

                            // k3
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                k3X[dimX] = problem->F[dimX]->eval(xTemp.data(), uGuess3.data(), &x[offXUTotal], t2);
                            }

                            std::vector<double> uGuess4 = problem->uInitialGuess(tij);

                            // x_i + dt * k3
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                xTemp[dimX] = x[xij - offXU + dimX] + dt * k3X[dimX];
                            }

                            // k4
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                k4X[dimX] = problem->F[dimX]->eval(xTemp.data(), uGuess4.data(), &x[offXUTotal], tij);
                            }

                            // x_{i+1} = x_i + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
                            for (int dimX = 0; dimX < problem->sizeX; dimX++) {
                                x[xij + dimX] = x[xij - offXU + dimX] + dt / 6.0 * (k1X[dimX] + 2.0 * (k2X[dimX] + k3X[dimX]) + k4X[dimX]);
                            }

                        }

                        const std::vector<double> uGuess = problem->uInitialGuess(tij);
                        for (int dimU = 0; dimU < problem->sizeU; dimU++) {
                            x[uij + dimU] = uGuess[dimU];
                        }

                        tijOld = tij;
                    }
                }

                break;

            case InitVars::CALLBACK:
                // strict callback case, optimizer will set initVars to Callback, if a mesh refinement is executed
                for (int i = 0; i < n; i++) {
                    x[i] = x_cb[i];
                }

                break;
        }
    }
    return true;
}

bool GDOP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value) {
    double MAY = 0;
    double LAG = 0;

    // mayer term evaluation
    if (problem->M) {
        MAY = problem->M->eval(&x[offXUTotal-offXU], &x[offXUTotal-offU], &x[offXUTotal], mesh.tf);
    }

    // lagrange term evaluation
    if (problem->L) {
        for (int i = 0; i < mesh.intervals; i++){
            for (int j = 0; j < rk.steps; j++){
                const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
                const int xij = i * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
                const int uij = i * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)
                LAG += mesh.deltaT[i] * rk.b[j] * problem->L->eval(&x[xij], &x[uij], &x[offXUTotal], tij);
            }
        }
    }
    obj_value = MAY + LAG;
    return true;
}

bool GDOP::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) {
    // initialize the dense gradient as 0 for all elements
    std::fill(grad_f, grad_f + numberVars, 0);

    // mayer term derivative
    if (problem->M) {

        const int xnm = offXUTotal - offXU; // index of 1st x var at collocation point (n,m)
        const int unm = offXUTotal - offU;  // index of 1st u var at collocation point (n,m)
        const auto diffMAY = problem->M->evalDiff(&x[xnm], &x[unm], &x[offXUTotal], mesh.tf);

        for (int v = 0; v < sz(diffMAY[0]); v++) {
            grad_f[xnm + problem->M->adj.indX[v]] += diffMAY[0][v];
        }
        for (int v = 0; v < sz(diffMAY[1]); v++) {
            grad_f[unm + problem->M->adj.indU[v]] += diffMAY[1][v];
        }
        for (int v = 0; v < sz(diffMAY[2]); v++) {
            grad_f[offXUTotal + problem->M->adj.indP[v]] += diffMAY[2][v];
        }
    }

    // lagrange term derivative
    if (problem->L) {
        for(int i = 0; i < mesh.intervals; i++){
            for (int j = 0; j < rk.steps; j++){

                const int xij = i * offXUBlock + j * offXU;        // index of 1st x var at collocation point (i,j)
                const int uij = i * offXUBlock + j * offXU + offX; // index of 1st u var at collocation point (i,j)
                const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
                const double quadCoeff = mesh.deltaT[i] * rk.b[j];
                const auto diffLAG = problem->L->evalDiff(&x[xij], &x[uij], &x[offXUTotal], tij);

                for (int v = 0; v < sz(diffLAG[0]); v++) {
                    grad_f[xij + problem->L->adj.indX[v]] += quadCoeff * diffLAG[0][v];
                }
                for (int v = 0; v < sz(diffLAG[1]); v++) {
                    grad_f[uij + problem->L->adj.indU[v]] += quadCoeff * diffLAG[1][v];
                }
                for (int v = 0; v < sz(diffLAG[2]); v++) {
                    grad_f[offXUTotal + problem->L->adj.indP[v]] += quadCoeff * diffLAG[2][v];
                }
            }
        }
    }
    return true;
}

bool GDOP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) {
    int eq = 0;
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < rk.steps; j++) {

            const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
            const int xij = i * offXUBlock + j * offXU;         // first index (dim=0) of x_{ij}
            const int uij = i * offXUBlock + j * offXU + offX;  // first index (dim=0) of u_{ij}
            const int xi1_m = i * offXUBlock - offXU;           // first index (dim=0) of x_{i-1,m}; m=rk.steps

            // dynamic constraints
            // 0 = sum_k ~a_{jk} * (x_{ik} - x_{i-1,m)) - del t_i * f(x_{ij}, u_{ij}, p, t_{ij}), i=0 -> x_{i-1,m) = x0
            for (int d = 0; d < sz(problem->F); d++) {
                double rkLinearComb = 0;
                for (int k = 0; k < rk.steps; k++) {
                    const int xik = i * offXUBlock + k * offXU;  // first index (dim=0) of x_{ik}
                    if (i == 0) {
                        rkLinearComb += rk.Ainv[j][k] * (x[xik + d] - problem->x0[d]);
                    }
                    else {
                        rkLinearComb += rk.Ainv[j][k] * (x[xik + d] - x[xi1_m + d]);
                    }
                }
                g[eq] = rkLinearComb - mesh.deltaT[i] * problem->F[d]->eval(&x[xij], &x[uij], &x[offXUTotal], tij);
                eq++;
            }

            // path constraints: LB_g <= g(x_{ij}, u_{ij}, p, t_{ij}) <= UB_g
            for (auto const& pathConstr : problem->G) {
                g[eq] = pathConstr->eval(&x[xij], &x[uij], &x[offXUTotal], tij);
                eq++;
            }
        }
    }

    // final constraints: LB_r <= r(x_{nm}, u_{nm}, p, t_{nm}) <= UB_r
    const int xnm = offXUTotal - offXU; // first index (dim=0) of x_{n,m}
    const int unm = offXUTotal - offU;  // first index (dim=0) of u_{n,m}
    for (auto const& finalConstr : problem->R) {
        g[eq] = finalConstr->eval(&x[xnm], &x[unm], &x[offXUTotal], mesh.tf);
        eq++;
    }

    // parameter constraints: LB_a <= a(p) <= UB_a
    for (auto const& paramConstr : problem->A) {
        g[eq] = paramConstr->eval(&x[offXUTotal]);
        eq++;
    }
    assert(m == eq);
    return true;
}

int GDOP::init_jac_sparsity(Index *iRow, Index *jCol) {
    // Just iterates over all blocks atm, no use of block struct
    int it = 0;
    int eq = 0;
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < rk.steps; j++) {

            const int xij = i * offXUBlock + j * offXU;         // first index (dim=0) of x_{ij}
            const int uij = i * offXUBlock + j * offXU + offX;  // first index (dim=0) of u_{ij}
            const int xi1_m = i * offXUBlock - offXU;           // first index (dim=0) of x_{i-1,m}; m=rk.steps

            for (int d = 0; d < sz(problem->F); d++) {

                // sum_k ~a_{jk} * (-x_{i-1,m}), i>=1
                if (i > 0) {
                    iRow[it] = eq;
                    jCol[it] = xi1_m + d;
                    it++;
                }

                // sum_k ~a_{jk} * (x_{i,k})
                for (int k = 0; k < rk.steps; k++) {
                    const int xik = i * offXUBlock + k * offXU;
                    iRow[it] = eq;
                    jCol[it] = xik + d;
                    it++;
                }

                // f(v_{ij})_x: care handle duplicates if d is contained in indX!!
                for (int v : problem->F[d]->adj.indX) {
                    if (v != d) {
                        iRow[it] = eq;
                        jCol[it] = xij + v;
                        it++;
                    }
                }

                // f(v_{ij})_u
                for (int v : problem->F[d]->adj.indU) {
                    iRow[it] = eq;
                    jCol[it] = uij + v;
                    it++;
                }

                // f(v_{ij})_p
                for (int v : problem->F[d]->adj.indP) {
                    iRow[it] = eq;
                    jCol[it] = offXUTotal + v;
                    it++;
                }
                eq++;
            }

            for (const auto & constrG : problem->G) {

                // g(v_{ij})_x
                for (int v : constrG->adj.indX) {
                    iRow[it] = eq;
                    jCol[it] = xij + v;
                    it++;
                }

                // g(v_{ij})_u
                for (int v : constrG->adj.indU) {
                    iRow[it] = eq;
                    jCol[it] = uij + v;
                    it++;
                }

                // g(v_{ij})_p
                for (int v : constrG->adj.indP) {
                    iRow[it] = eq;
                    jCol[it] = offXUTotal + v;
                    it++;
                }
                eq++;
            }
        }
    }

    const int xnm = offXUTotal - offXU;
    const int unm = offXUTotal - offU;

    for (const auto & constrR : problem->R) {

        // r(v_{nm})_x
        for (int v : constrR->adj.indX) {
            iRow[it] = eq;
            jCol[it] = xnm + v;
            it++;
        }

        // r(v_{nm})_u
        for (int v : constrR->adj.indU) {
            iRow[it] = eq;
            jCol[it] = unm + v;
            it++;
        }

        // r(v_{nm})_p
        for (int v : constrR->adj.indP) {
            iRow[it] = eq;
            jCol[it] = offXUTotal + v;
            it++;
        }
        eq++;
    }

    for (const auto & constrA : problem->A) {

        // a_p
        for (int v : constrA->adj.indP) {
            iRow[it] = eq;
            jCol[it] = offXUTotal + v;
            it++;
        }
        eq++;
    }
    return eq;
}

void GDOP::get_jac_values(const Number *x, Number *values) {
    int containedIdx = -1;
    int it = 0;
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < rk.steps; j++) {
            const int xij = i * offXUBlock + j * offXU;         // first index (dim=0) of x_{ij}
            const int uij = i * offXUBlock + j * offXU + offX;  // first index (dim=0) of u_{ij}
            const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];

            for (int d = 0; d < sz(problem->F); d++) {

                // sum_k ~a_{jk} * (-x_{i-1,m}), i>=1
                if (i > 0) {
                    values[it] = -rk.invRowSum[j];
                    it++;
                }

                // sum_k ~a_{jk} * (x_{i,k})
                for (int k = 0; k < rk.steps; k++) {
                    values[it] = rk.Ainv[j][k];
                    if (k == j) {
                        containedIdx = it;
                    }
                    it++;
                }

                // eval grad f(v_{ij})
                auto const diffF = problem->F[d]->evalDiff(&x[xij], &x[uij], &x[offXUTotal], tij);

                // f(v_{ij})_x: care handle duplicates if d is contained in indX!!
                for (int v = 0; v < sz(problem->F[d]->adj.indX); v++) {
                    auto idx = problem->F[d]->adj.indX[v];
                    if (idx != d) {
                        values[it] = -mesh.deltaT[i] * diffF[0][v];
                        it++;
                    }
                    else {
                        values[containedIdx] -= mesh.deltaT[i] * diffF[0][v];
                    }
                }

                // f(v_{ij})_u
                for (int v = 0; v < sz(problem->F[d]->adj.indU); v++) {
                    values[it] = -mesh.deltaT[i] * diffF[1][v];
                    it++;
                }

                // f(v_{ij})_p
                for (int v = 0; v < sz(problem->F[d]->adj.indP); v++) {
                    values[it] = -mesh.deltaT[i] * diffF[2][v];
                    it++;
                }
            }

            for (const auto & constrG : problem->G) {

                // eval grad g(v_{ij})
                auto const diffG = constrG->evalDiff(&x[xij], &x[uij], &x[offXUTotal], tij);

                // g(v_{ij})_x
                for (int v = 0; v < sz(constrG->adj.indX); v++) {
                    values[it] = diffG[0][v];
                    it++;
                }

                // g(v_{ij})_u
                for (int v = 0; v < sz(constrG->adj.indU); v++) {
                    values[it] = diffG[1][v];
                    it++;
                }

                // g(v_{ij})_p
                for (int v = 0; v < sz(constrG->adj.indP); v++) {
                    values[it] = diffG[2][v];
                    it++;
                }
            }
        }
    }

    const int xnm = offXUTotal - offXU;
    const int unm = offXUTotal - offU;

    for (const auto & constrR : problem->R) {

        // eval grad r(v_{nm})
        auto const diffR = constrR->evalDiff(&x[xnm], &x[unm], &x[offXUTotal], mesh.tf);

        // r(v_{nm})_x
        for (int v = 0; v < sz(constrR->adj.indX); v++) {
            values[it] = diffR[0][v];
            it++;
        }

        // r(v_{nm})_u
        for (int v = 0; v < sz(constrR->adj.indU); v++) {
            values[it] = diffR[1][v];
            it++;
        }

        // r(v_{nm})_p
        for (int v = 0; v < sz(constrR->adj.indP); v++) {
            values[it] = diffR[2][v];
            it++;
        }
    }

    for (const auto & constrA : problem->A) {

        // eval grad a(p)
        auto const diffA = constrA->evalDiff(&x[offXUTotal]);

        // a_p
        for (int v = 0; v < sz(constrA->adj.indP); v++) {
            values[it] = diffA[v];
            it++;
        }
    }
}

bool GDOP::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol,
                      Number *values) {
    assert(n == numberVars);
    assert(m == sz(problem->A) + sz(problem->R) + (sz(problem->F) + sz(problem->G)) * rk.steps * mesh.intervals);
    if (values == nullptr) {
        int eq = init_jac_sparsity(iRow, jCol);
        if (exportJacobian) {
            exportSparsity(iRow, jCol, nele_jac, {m, n}, exportJacobianPath + "/" + problem->name + "_jacobian.csv");
        }
        assert(eq == m);
    }
    else
    {
        std::fill(values, values + nele_jac, 0);
        get_jac_values(x, values);
    }
    return true;
}

void GDOP::init_h_sparsity(Index *iRow, Index *jCol) {

    // S0 block forall i, j except very last interval (n, m)
    for (const auto &[vars, it]: S0) {
        auto const [v1, v2] = vars;
        for (int i = 0; i < mesh.intervals - 1; i++) {
            for (int j = 0; j < rk.steps; j++) {
                const int idxij = it + lengthS0 * (i * rk.steps + j);
                const int v1ij = i * offXUBlock + j * offXU + v1;
                const int v2ij = i * offXUBlock + j * offXU + v2;
                iRow[idxij] = v1ij;
                jCol[idxij] = v2ij;
            }
        }
         for (int j = 0; j < rk.steps - 1; j++) {
            const int idxij = it + lengthS0 * ((mesh.intervals - 1) * rk.steps + j);
            const int v1ij = (mesh.intervals - 1) * offXUBlock + j * offXU + v1;
            const int v2ij = (mesh.intervals - 1) * offXUBlock + j * offXU + v2;
            iRow[idxij] = v1ij;
            jCol[idxij] = v2ij;
        }
    }

    // S0t for very last interval, note that "it" is the correct array index by construction
    for (const auto &[vars, it]: S0t) {
        auto const [v1, v2] = vars;
        for (int j = 0; j < rk.steps; j++) {
            const int v1ij = (mesh.intervals - 1) * offXUBlock + j * offXU + v1;
            const int v2ij = (mesh.intervals - 1) * offXUBlock + j * offXU + v2;
            iRow[it] = v1ij;
            jCol[it] = v2ij;
        }
    }

    // S1 block forall i, j except very last interval (n, m)
    for (const auto &[vars, it]: S1) {
        auto const [p, v] = vars;
        for (int i = 0; i < mesh.intervals - 1; i++) {
            for (int j = 0; j < rk.steps; j++) {
                const int idxij = it + rowLengthS1Block[p] * (i * rk.steps + j);
                const int pshifted = offXUTotal + p;
                const int vij = i * offXUBlock + j * offXU + v;
                iRow[idxij] = pshifted;
                jCol[idxij] = vij;
            }
        }
        for (int j = 0; j < rk.steps - 1; j++) {
            const int idxij = it + rowLengthS1Block[p] * ((mesh.intervals - 1) * rk.steps + j);
            const int pshifted = offXUTotal + p;
            const int vij = (mesh.intervals - 1) * offXUBlock + j * offXU + v;
            iRow[idxij] = pshifted;
            jCol[idxij] = vij;
        }
    }

    // S1t for very last interval, note that "it" is the correct array index by construction
    for (const auto &[vars, it]: S1t) {
        auto const [p, v] = vars;
        for (int j = 0; j < rk.steps; j++) {
            const int pshifted = offXUTotal + p;
            const int vij = (mesh.intervals - 1) * offXUBlock + j * offXU + v;
            iRow[it] = pshifted;
            jCol[it] = vij;
        }
    }

    // S2
    for (const auto &[vars, it]: S2) {
        auto const [p1, p2] = vars;
        const int p1shifted = offXUTotal + p1;
        const int p2shifted = offXUTotal + p2;
        iRow[it] = p1shifted;
        jCol[it] = p2shifted;
    }
}

void GDOP::evalHessianS0_S1(Number* values, const Number *x, Expression& expr, const double factor, const int xij,
        const int uij, const double tij, const int i, const int j) {
    auto const diff2Expr = expr.evalDiff2(&x[xij], &x[uij], &x[offXUTotal], tij);
    for (int k = 0; k < sz(expr.adjDiff.indXX); k++) {
        auto const vars = expr.adjDiff.indXX[k];
        auto const idx = lengthS0 * (i * rk.steps + j) + S0[vars];
        values[idx] += factor * diff2Expr[0][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indUX); k++) {
        auto const [uvar, xvar] = expr.adjDiff.indUX[k];
        auto const idx = lengthS0 * (i * rk.steps + j) + S0[{offX + uvar, xvar}];
        values[idx] += factor * diff2Expr[1][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indUU); k++) {
        auto const [uvar1, uvar2] = expr.adjDiff.indUU[k];
        auto const idx = lengthS0 * (i * rk.steps + j) + S0[{offX + uvar1, offX + uvar2}];
        values[idx] += factor * diff2Expr[2][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indPX); k++) {
        auto const [pvar, xvar] = expr.adjDiff.indPX[k];
        auto const it = S1[{pvar, xvar}];
        auto const idx = it + rowLengthS1Block[pvar] * (i * rk.steps + j);
        values[idx] += factor * diff2Expr[3][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indPU); k++) {
        auto const [pvar, uvar] = expr.adjDiff.indPU[k];
        auto const it = S1[{pvar, uvar + offX}];
        auto const idx = it + rowLengthS1Block[pvar] * (i * rk.steps + j);
        values[idx] += factor * diff2Expr[4][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indPP); k++) {
        auto const vars = expr.adjDiff.indPP[k];
        auto const it = S2[vars];
        values[it] += factor * diff2Expr[5][k];
    }
}

void GDOP::evalHessianS0t_S1t(Number* values, const Number* x, Expression& expr, const double factor, const int xij,
        const int uij, const double tij) {
    auto const diff2Expr = expr.evalDiff2(&x[xij], &x[uij], &x[offXUTotal], tij);
    for (int k = 0; k < sz(expr.adjDiff.indXX); k++) {
        auto const vars = expr.adjDiff.indXX[k];
        values[S0t[vars]] += factor * diff2Expr[0][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indUX); k++) {
        auto const [uvar, xvar] = expr.adjDiff.indUX[k];
        values[S0t[{offX + uvar, xvar}]] += factor * diff2Expr[1][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indUU); k++) {
        auto const [uvar1, uvar2] = expr.adjDiff.indUU[k];
        values[S0t[{offX + uvar1, offX + uvar2}]] += factor * diff2Expr[2][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indPX); k++) {
        auto const vars = expr.adjDiff.indPX[k];
        values[S1t[vars]] += factor * diff2Expr[3][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indPU); k++) {
        auto const [pvar, uvar] = expr.adjDiff.indPU[k];
        values[S1t[{pvar, uvar + offX}]] += factor * diff2Expr[4][k];
    }
    for (int k = 0; k < sz(expr.adjDiff.indPP); k++) {
        auto const vars = expr.adjDiff.indPP[k];
        values[S2[vars]] += factor * diff2Expr[5][k];
    }
}

void GDOP::evalHessianS2(Number* values, const Number *x, ParamExpression& expr, double factor) {
    auto const diff2Expr = expr.evalDiff2(&x[offXUTotal]);
    for (int k = 0; k < sz(expr.adjDiff.indPP); k++) {
        auto const vars = expr.adjDiff.indPP[k];
        values[S2[vars]] += factor * diff2Expr[k];
    }
}

int GDOP::get_h_values(const Number *x, Number *values, Number obj_factor, const Number *lambda) {
    int eq = 0;
    for (int i = 0; i < mesh.intervals - 1; i++) {

        for (int j = 0; j < rk.steps; j++) {
            const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
            const int xij = i * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
            const int uij = i * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)

            // eval hessian lagrangian
            if (problem->L) {
                const double lFactor = obj_factor * mesh.deltaT[i] * rk.b[j];
                evalHessianS0_S1(values, x, *problem->L, lFactor, xij, uij, tij, i, j);
            }

            // eval hessian dynamics
            for (auto const& f : problem->F) {
                const double fFactor = -lambda[eq] * mesh.deltaT[i];
                evalHessianS0_S1(values, x, *f, fFactor, xij, uij, tij, i, j);
                eq++;
            }

            // eval hessian path constraints
            for (auto const& g : problem->G) {
                const double gFactor = lambda[eq];
                evalHessianS0_S1(values, x, *g, gFactor, xij, uij, tij, i, j);
                eq++;
            }
        }
    }

    // same for the last rk.steps-1 blocks excluding the very last block -> S0t handling
    for (int j = 0; j < rk.steps - 1; j++) {
        const double tij = mesh.grid[mesh.intervals - 1] + rk.c[j] * mesh.deltaT[mesh.intervals - 1];
        const int xij = (mesh.intervals - 1) * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
        const int uij = (mesh.intervals - 1) * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)

        // eval hessian lagrangian
        if (problem->L) {
            const double lagrFactor = obj_factor * mesh.deltaT[(mesh.intervals - 1)] * rk.b[j];
            evalHessianS0_S1(values, x, *problem->L, lagrFactor, xij, uij, tij, (mesh.intervals - 1), j);
        }

        // eval hessian dynamics
        for (auto const& f : problem->F) {
            const double fFactor = lambda[eq] * (-mesh.deltaT[(mesh.intervals - 1)]);
            evalHessianS0_S1(values, x, *f, fFactor, xij, uij, tij, mesh.intervals - 1, j);
            eq++;
        }

        // eval hessian path constraints
        for (auto const& g : problem->G) {
            const double gFactor = lambda[eq];
            evalHessianS0_S1(values, x, *g, gFactor, xij, uij, tij, mesh.intervals - 1, j);
            eq++;
        }
    }

    // S0t, S1t: contains M, L, F, G, R
    const int xnm = (mesh.intervals - 1) * offXUBlock + (rk.steps - 1) * offXU;         // index of 1st x var at collocation point (n,m)
    const int unm = (mesh.intervals - 1) * offXUBlock + (rk.steps - 1) * offXU + offX;  // index of 1st u var at collocation point (n,m)

    // eval hessian mayer
    if (problem->M) {
        const double mayFactor = obj_factor;
        evalHessianS0t_S1t(values, x, *problem->M, mayFactor, xnm, unm, mesh.tf);
    }

    // eval hessian lagrangian
    if (problem->L) {
        const double lagrFactor = obj_factor * mesh.deltaT[(mesh.intervals - 1)] * rk.b[rk.steps - 1];
        evalHessianS0t_S1t(values, x, *problem->L, lagrFactor, xnm, unm, mesh.tf);
    }

    // eval hessian dynamics
    for (auto const& f : problem->F) {
        const double fFactor = lambda[eq] * (-mesh.deltaT[(mesh.intervals - 1)]);
        evalHessianS0t_S1t(values, x, *f, fFactor, xnm, unm, mesh.tf);
        eq++;
    }

    // eval hessian path constraints
    for (auto const& g : problem->G) {
        const double gFactor = lambda[eq];
        evalHessianS0t_S1t(values, x, *g, gFactor, xnm, unm, mesh.tf);
        eq++;
    }

    // eval hessian final constraints
    for (auto const& r : problem->R) {
        const double rFactor = lambda[eq];
        evalHessianS0t_S1t(values, x, *r, rFactor, xnm, unm, mesh.tf);
        eq++;
    }

    // eval hessian algebraic parameter constraints
    for (auto const& a : problem->A) {
        const double aFactor = lambda[eq];
        evalHessianS2(values, x, *a, aFactor);
        eq++;
    }
    return eq;
}

bool GDOP::eval_h(Index n, const Number *x, bool new_x, Number obj_factor, Index m, const Number *lambda,
                  bool new_lambda, Index nele_hess, Index *iRow, Index *jCol, Number *values) {
    if (values == nullptr) {
        init_h_sparsity(iRow, jCol);
        if (exportHessian) {
            exportSparsity(iRow, jCol, nele_hess, {n, n}, exportHessianPath + "/" + problem->name + "_hessian.csv");
        }
    }
    else
    {
        std::fill(values, values + nele_hess, 0);
        int eq = get_h_values(x, values, obj_factor, lambda);
        assert(eq == m); // same as before eq index over lambda must be equal to #constrs
    }
    return true;
}

void GDOP::finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L, const Number *z_U,
                             Index m, const Number *g, const Number *lambda, Number obj_value, const IpoptData *ip_data,
                             IpoptCalculatedQuantities *ip_cq) {
    optimum.assign(x, x + n);
    objective = obj_value;
}

void GDOP::exportOptimum(const std::string& filename) const {
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }
    std::string vars;
    for (int vx = 0; vx < problem->sizeX; vx++) {
        vars += ",x[" + std::to_string(vx) + "]";
    }
    for (int vu = 0; vu < problem->sizeU; vu++) {
        vars += ",u[" + std::to_string(vu) + "]";
    }
    for (int vp = 0; vp < problem->sizeP; vp++) {
        vars += ",p[" + std::to_string(vp) + "]";
    }
    auto header = "time" + vars;
    outFile << header << "\n";

    // time = 0 -> interpolate the control backwards

    std::string values = "0";
    for (int vx = 0; vx < problem->sizeX; vx++) {
        values += "," + double2Str(problem->x0[vx]);
    }
    for (int vu = 0; vu < problem->sizeU; vu++) {
        std::vector<double> uValues = {};
        uValues.reserve(rk.steps);
        for (int j = 0; j < rk.steps; j++) {
            uValues.push_back(optimum[vu + offX + offXU * j]);
        }
        values += "," + double2Str(Integrator::evalLagrange(rk.c, uValues, 0.0));
    }
    for (int vp = 0; vp < problem->sizeP; vp++) {
        values += "," + double2Str(optimum[vp + offXUTotal]);
    }
    outFile << values << "\n";

    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < rk.steps; j++) {
            values = double2Str(mesh.grid[i] + rk.c[j] * mesh.deltaT[i]);
            for (int vx = 0; vx < problem->sizeX; vx++) {
                values += "," + double2Str(optimum[vx + offXU * j + offXUBlock * i]);
            }
            for (int vu = 0; vu < problem->sizeU; vu++) {
                values += "," + double2Str(optimum[vu + offX + offXU * j + offXUBlock * i]);
            }
            for (int vp = 0; vp < problem->sizeP; vp++) {
                values += "," + double2Str(optimum[vp + offXUTotal]);
            }
            outFile << values << "\n";
        }
    }
}

GDOP::GDOP(const std::shared_ptr<const Problem> & problem, Mesh &mesh, Integrator &rk, InitVars initVars) :
                problem(problem), mesh(mesh), rk(rk), initVars(initVars) {}


GDOP *create_gdop(const std::shared_ptr<const Problem> & problem, Mesh &mesh, Integrator &rk, InitVars initVars) {
    return new GDOP(problem, mesh, rk, initVars);
}