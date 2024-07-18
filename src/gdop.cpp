#include <cassert>
#include <algorithm>
#include <iostream>
#include "gdop.h"
#include "util.h"


// TODO: Check entire indices!
bool GDOP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
    // #vars
    n = numberVars;

    // #eqs
    m = sz(problem.A) + sz(problem.R) + (sz(problem.F) + sz(problem.G)) * rk.steps * mesh.intervals;

    // eval jacobian nnz
    nnz_jac_g = 0;

    // nnz jac: dynamics
    int nnzDynBlock = 0;
    int containedIndex = 0; // target idx
    for (const auto& dyn : problem.F) {
        nnzDynBlock += sz(dyn->adj.indX) + sz(dyn->adj.indU) + sz(dyn->adj.indP);
        // intersection of state indices may not be empty!
        // check if the k-th component of x, that is always contained in this block equation,
        // is contained in grad f_k(.) as well => reduce block nnz by 1
        auto it = std::find(dyn->adj.indX.begin(), dyn->adj.indX.end(), containedIndex);
        if (it != dyn->adj.indX.end())
            nnzDynBlock -= 1;
        containedIndex++;
    }
    nnz_jac_g += mesh.intervals * rk.steps * nnzDynBlock;                // nnz for f(v_{ij}) eval
    nnz_jac_g += mesh.intervals * rk.steps * rk.steps * (sz(problem.F)); // nnz for sum_k ~a_{jk} * x_{ik}
    nnz_jac_g += (mesh.intervals - 1) * rk.steps * sz(problem.F);        // nnz for sum_k ~a_{jk} * (-x_{i-1,m}), i>=1

    // nnz jac: path constraints
    int nnzPathBlock = 0;
    for (const auto& pathConstr : problem.G) {
        nnzPathBlock += sz(pathConstr->adj.indX) + sz(pathConstr->adj.indU) + sz(pathConstr->adj.indP);
    }
    nnz_jac_g += mesh.intervals * rk.steps * nnzPathBlock; // RHS = nnz for g(v_{ij}) eval

    // nnz jac: final constraints
    for (const auto& finalConstr : problem.R) {
        nnz_jac_g += sz(finalConstr->adj.indX) + sz(finalConstr->adj.indU) + sz(finalConstr->adj.indP);
    }

    // nnz jac: algebraic parameter constraints
    for (const auto& paramConstr : problem.A) {
        nnz_jac_g += sz(paramConstr->adj.indP);
    }

    // TODO: implement hessian sparsity
    // eval hessian nnz
    nnz_h_lag = 0;

    index_style = TNLP::C_STYLE; // index starts at 0
    return true;
}

bool GDOP::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) {
    assert(n == numberVars);
    assert(m == sz(problem.A) + sz(problem.R) + (sz(problem.F) + sz(problem.G)) * rk.steps * mesh.intervals);

    // var bounds
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < rk.steps; j++) {

            const int xij = i * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
            const int uij = i * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)

            // state bounds
            for (int d = 0; d < problem.sizeX; d++) {
                x_l[xij + d] = problem.lbX[d];
                x_u[xij + d] = problem.ubX[d];
            }

            // control bounds
            for (int d = 0; d < problem.sizeU; d++) {
                x_l[uij + d] = problem.lbU[d];
                x_u[uij + d] = problem.ubU[d];
            }
        }
    }

    // parameter bounds
    for (int d = 0; d < problem.sizeP; d++) {
        x_l[offXUTotal + d] = problem.lbP[d];
        x_u[offXUTotal + d] = problem.ubP[d];
    }

    // equation bounds
    int eq = 0;
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < rk.steps; j++) {
            // dynamic bound must equal 0
            for (int d = 0; d < problem.F.size(); d++) {
                g_l[eq] = 0.0;
                g_u[eq] = 0.0;
                eq++;
            }

            // path constraint bounds
            for (auto const& pathConstr : problem.G) {
                g_l[eq] = pathConstr->lb;
                g_u[eq] = pathConstr->ub;
                eq++;
            }
        }
    }

    // final constraint bounds
    for (auto const& finalConstr : problem.R) {
        g_l[eq] = finalConstr->lb;
        g_u[eq] = finalConstr->ub;
        eq++;
    }

    // parameter constraint bounds
    for (auto const& paramConstr : problem.A) {
        g_l[eq] = paramConstr->lb;
        g_u[eq] = paramConstr->ub;
        eq++;
    }
    assert(m == eq);
    return true;
}

bool GDOP::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index m,
                              bool init_lambda, Number *lambda) {
    // TODO: implement different strategies: SOLVE(), WARM_START
    assert(n == numberVars);
    switch (initVars) {
        case InitVars::CONST:
            for (int i = 0; i < n; i++) {
                x[i] = 1;
            }
    }
    return true;
}

bool GDOP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value) {
    double MAY = 0;
    double LAG = 0;

    // mayer term evaluation
    if (problem.M) {
        MAY = problem.M->eval(&x[offXUTotal-offXU], &x[offXUTotal-offU], &x[offXUTotal], mesh.tf);
    }

    // lagrange term evaluation
    if (problem.L) {
        for (int i = 0; i < mesh.intervals; i++){
            for (int j = 0; j < rk.steps; j++){
                const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
                const int xij = i * offXUBlock + j * offXU;         // index of 1st x var at collocation point (i,j)
                const int uij = i * offXUBlock + j * offXU + offX;  // index of 1st u var at collocation point (i,j)
                LAG += mesh.deltaT[i] * rk.b[j] * problem.L->eval(&x[xij], &x[uij], &x[offXUTotal], tij);
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
    if (problem.M) {

        const int xnm = offXUTotal - offXU; // index of 1st x var at collocation point (n,m)
        const int unm = offXUTotal - offU;  // index of 1st u var at collocation point (n,m)
        const auto diffMAY = problem.M->evalDiff(&x[xnm], &x[unm], &x[offXUTotal], mesh.tf);

        for (int v = 0; v < diffMAY[0].size(); v++) {
            grad_f[xnm + problem.M->adj.indX[v]] += diffMAY[0][v];
        }
        for (int v = 0; v < diffMAY[1].size(); v++) {
            grad_f[unm + problem.M->adj.indU[v]] += diffMAY[1][v];
        }
        for (int v = 0; v < diffMAY[2].size(); v++) {
            grad_f[offXUTotal + problem.M->adj.indP[v]] += diffMAY[2][v];
        }
    }

    // lagrange term derivative
    if (problem.L) {
        for(int i = 0; i < mesh.intervals; i++){
            for (int j = 0; j < rk.steps; j++){

                const int xij = i * offXUBlock + j * offXU;        // index of 1st x var at collocation point (i,j)
                const int uij = i * offXUBlock + j * offXU + offX; // index of 1st u var at collocation point (i,j)
                const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];
                const double quadCoeff = mesh.deltaT[i] * rk.b[j];
                const auto diffLAG = problem.L->evalDiff(&x[xij], &x[uij], &x[offXUTotal], tij);

                for (int v = 0; v < diffLAG[0].size(); v++) {
                    grad_f[xij + problem.L->adj.indX[v]] += quadCoeff * diffLAG[0][v];
                }
                for (int v = 0; v < diffLAG[1].size(); v++) {
                    grad_f[uij + problem.L->adj.indU[v]] += quadCoeff * diffLAG[1][v];
                }
                for (int v = 0; v < diffLAG[2].size(); v++) {
                    grad_f[offXUTotal + problem.L->adj.indP[v]] += quadCoeff * diffLAG[2][v];
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
            if (i == 0){
                for (int d = 0; d < problem.F.size(); d++) {
                    double rkLinearComb = 0;
                    for (int k = 0; k < rk.steps; k++) {
                        const int xik = i * offXUBlock + k * offXU;  // first index (dim=0) of x_{ik}
                        rkLinearComb += rk.Ainv[j][k] * (x[xik + d] - problem.x0[d]);
                    }
                    g[eq] = rkLinearComb - mesh.deltaT[i] * problem.F[d]->eval(&x[xij], &x[uij], &x[offXUTotal], tij);
                    eq++;
                }
            }
            else {
                for (int d = 0; d < problem.F.size(); d++) {
                    double rkLinearComb = 0;
                    for (int k = 0; k < rk.steps; k++) {
                        const int xik = i * offXUBlock + k * offXU;  // first index (dim=0) of x_{ik}
                        rkLinearComb += rk.Ainv[j][k] * (x[xik + d] - x[xi1_m + d]);
                    }
                    g[eq] = rkLinearComb - mesh.deltaT[i] * problem.F[d]->eval(&x[xij], &x[uij], &x[offXUTotal], tij);
                    eq++;
                }
            }

            // path constraints: LB_g <= g(x_{ij}, u_{ij}, p, t_{ij}) <= UB_g
            for (auto const& pathConstr : problem.G) {
                g[eq] = pathConstr->eval(&x[xij], &x[uij], &x[offXUTotal], tij);
                eq++;
            }
        }
    }

    // final constraints: LB_r <= r(x_{nm}, u_{nm}, p, t_{nm}) <= UB_r
    const int xnm = offXUTotal - offXU; // first index (dim=0) of x_{n,m}
    const int unm = offXUTotal - offU;  // first index (dim=0) of u_{n,m}
    for (auto const& finalConstr : problem.R) {
        g[eq] = finalConstr->eval(&x[xnm], &x[unm], &x[offXUTotal], mesh.tf);
        eq++;
    }

    // parameter constraints: LB_a <= a(p) <= UB_a
    for (auto const& paramConstr : problem.A) {
        g[eq] = paramConstr->eval(&x[offXUTotal]);
        eq++;
    }
    assert(m == eq);
    return true;
}

bool GDOP::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol,
                      Number *values) {

    assert(n == numberVars);
    assert(m == sz(problem.A) + sz(problem.R) + (sz(problem.F) + sz(problem.G)) * rk.steps * mesh.intervals);
    if (values == nullptr) {
        int it = 0;
        int eq = 0;
        for (int i = 0; i < mesh.intervals; i++) {
            for (int j = 0; j < rk.steps; j++) {
                // TODO: maybe use block structure of jacobian (might be slight speed up)
                const int xij = i * offXUBlock + j * offXU;         // first index (dim=0) of x_{ij}
                const int uij = i * offXUBlock + j * offXU + offX;  // first index (dim=0) of u_{ij}
                const int xi1_m = i * offXUBlock - offXU;           // first index (dim=0) of x_{i-1,m}; m=rk.steps

                for (int d = 0; d < problem.F.size(); d++) {

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
                    for (int v : problem.F[d]->adj.indX) {
                        if (v != d) {
                            iRow[it] = eq;
                            jCol[it] = xij + v;
                            it++;
                        }
                    }

                    // f(v_{ij})_u
                    for (int v : problem.F[d]->adj.indU) {
                        iRow[it] = eq;
                        jCol[it] = uij + v;
                        it++;
                    }

                    // f(v_{ij})_p
                    for (int v : problem.F[d]->adj.indP) {
                        iRow[it] = eq;
                        jCol[it] = offXUTotal + v;
                        it++;
                    }
                    eq++;
                }

                for (const auto & constrG : problem.G) {

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

        for (const auto & constrR : problem.R) {

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

        for (const auto & constrA : problem.A) {

            // a_p
            for (int v : constrA->adj.indP) {
                iRow[it] = eq;
                jCol[it] = offXUTotal + v;
                it++;
            }
            eq++;
        }
        assert(m == eq);
    }
    else
    {
        int containedIdx = -1;
        int it = 0;
        for (int i = 0; i < mesh.intervals; i++) {
            for (int j = 0; j < rk.steps; j++) {
                const int xij = i * offXUBlock + j * offXU;         // first index (dim=0) of x_{ij}
                const int uij = i * offXUBlock + j * offXU + offX;  // first index (dim=0) of u_{ij}
                const double tij = mesh.grid[i] + rk.c[j] * mesh.deltaT[i];

                for (int d = 0; d < problem.F.size(); d++) {

                    // sum_k ~a_{jk} * (-x_{i-1,m}), i>=1
                    if (i > 0) {
                        values[it] = -rk.invRowSum[j];
                        it++;
                    }

                    // sum_k ~a_{jk} * (x_{i,k})
                    for (int k = 0; k < rk.steps; k++) {
                        values[it] = rk.Ainv[j][k];
                        if (k == j)
                            containedIdx = it;
                        it++;
                    }


                    // eval grad f(v_{ij})
                    auto const diffF = problem.F[d]->evalDiff(&x[xij], &x[uij], &x[offXUTotal], tij);

                    // f(v_{ij})_x: care handle duplicates if d is contained in indX!!
                    for (int v : problem.F[d]->adj.indX) {
                        if (v != d) {
                            values[it] = -mesh.deltaT[i] * diffF[0][v];
                            it++;
                        }
                        else {
                            values[containedIdx] += -mesh.deltaT[i] * diffF[0][v];
                        }
                    }

                    // f(v_{ij})_u
                    for (int v : problem.F[d]->adj.indU) {
                        values[it] = -mesh.deltaT[i] * diffF[1][v];
                        it++;
                    }

                    // f(v_{ij})_p
                    for (int v : problem.F[d]->adj.indP) {
                        values[it] = mesh.deltaT[i] * diffF[2][v];
                        it++;
                    }
                }

                for (const auto & constrG : problem.G) {

                    // eval grad g(v_{ij})
                    auto const diffG = constrG->evalDiff(&x[xij], &x[uij], &x[offXUTotal], tij);

                    // g(v_{ij})_x
                    for (int v : constrG->adj.indX) {
                        values[it] = diffG[0][v];
                        it++;
                    }

                    // g(v_{ij})_u
                    for (int v : constrG->adj.indU) {
                        values[it] = diffG[1][v];
                        it++;
                    }

                    // g(v_{ij})_p
                    for (int v : constrG->adj.indP) {
                        values[it] = diffG[2][v];
                        it++;
                    }
                }
            }
        }

        const int xnm = offXUTotal - offXU;
        const int unm = offXUTotal - offU;

        for (const auto & constrR : problem.R) {

            // eval grad r(v_{nm})
            auto const diffR = constrR->evalDiff(&x[xnm], &x[unm], &x[offXUTotal], mesh.tf);

            // r(v_{nm})_x
            for (int v : constrR->adj.indX) {
                values[it] = diffR[0][v];
                it++;
            }

            // r(v_{nm})_u
            for (int v : constrR->adj.indU) {
                values[it] = diffR[1][v];
                it++;
            }

            // r(v_{nm})_p
            for (int v : constrR->adj.indP) {
                values[it] = diffR[2][v];
                it++;
            }
        }

        for (const auto & constrA : problem.A) {

            // eval grad a(p)
            auto const diffA = constrA->evalDiff(&x[offXUTotal]);

            // a_p
            for (int v : constrA->adj.indP) {
                values[it] = diffA[v];
                it++;
            }
        }
    }
    return true;
}

bool GDOP::eval_h(Index n, const Number *x, bool new_x, Number obj_factor, Index m, const Number *lambda,
                  bool new_lambda, Index nele_hess, Index *iRow, Index *jCol, Number *values) {
    // TODO: hessian eval
    return true;
}

void GDOP::finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L, const Number *z_U,
                             Index m, const Number *g, const Number *lambda, Number obj_value, const IpoptData *ip_data,
                             IpoptCalculatedQuantities *ip_cq) {
    // TODO: output solution
}

GDOP::GDOP(Problem problem, Mesh &mesh, Integrator &rk, InitVars initVars) : problem(std::move(problem)), mesh(mesh),
                                                                             rk(rk), initVars(initVars) {
}
