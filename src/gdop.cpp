#include <cassert>
#include <algorithm>
#include "gdop.h"
#include "util.h"

bool GDOP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
    n = numberVars;
    m = sz(problem.A) + sz(problem.R) + (sz(problem.F) + sz(problem.G)) * rk.steps * mesh.intervals;

    // eval jacobian nnz
    nnz_jac_g = 0;

    // nnz jac - dynamics
    int nnzDynBlock = 0;
    int containedIndex = 0; // target idx
    for (const auto& dynConstr : problem.F) {
        nnzDynBlock += sz(dynConstr->adj.indX) + sz(dynConstr->adj.indU) + sz(dynConstr->adj.indP);
        // intersection of state indices may not be empty!
        // check if the k-th component of x, that is always contained in this block equation,
        // is contained in grad f_k(.) as well => reduce block nnz by 1
        auto it = std::find(dynConstr->adj.indX.begin(), dynConstr->adj.indX.end(), containedIndex);
        if (it != dynConstr->adj.indX.end())
            nnzDynBlock -= 1;
        containedIndex++;
    }
    nnz_jac_g += mesh.intervals * rk.steps * nnzDynBlock; // nnz for f(v_{i,j}) eval
    nnz_jac_g += mesh.intervals * rk.steps * rk.steps * (sz(problem.F)); // nnz for sum_k ~a_{j,k} * x_{i,k}
    nnz_jac_g += (mesh.intervals - 1) * rk.steps * sz(problem.F); // nnz for sum_k ~a_{j,k} * (-x_{i-1,m}), i = 1,2,...

    // nnz jac - path
    int nnzPathBlock = 0;
    for (const auto& pathConstr : problem.G) {
        nnzPathBlock += sz(pathConstr->adj.indX) + sz(pathConstr->adj.indU) + sz(pathConstr->adj.indP);
    }
    nnz_jac_g += mesh.intervals * rk.steps * nnzPathBlock; // RHS = nnz for g(v_{i,j}) eval

    // nnz jac - final constraints
    for (const auto& finalConstr : problem.R) {
        nnz_jac_g += sz(finalConstr->adj.indX) + sz(finalConstr->adj.indU) + sz(finalConstr->adj.indP);
    }

    // nnz jac - algebraic parameter constraints
    for (const auto& paramConstr : problem.A) {
        nnz_jac_g += sz(paramConstr->adj.indP);
    }

    // TODO: implement hessian sparsity
    // eval hessian nnz
    nnz_h_lag = -1;

    index_style = TNLP::C_STYLE; // index starts at 0
    return true;
}

bool GDOP::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) {
    // TODO: implement bounds init
    return true;
}

bool GDOP::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index m,
                              bool init_lambda, Number *lambda) {
    // TODO: implement x init
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
        for(Index i = 0; i < mesh.intervals; i++){
            for (Index j = 0; j < rk.steps; j++){
                const int offset = i * offXUBlock + j * offXU;  // index of 1st x var at collocation point (i,j)
                LAG += mesh.deltaT[i] * rk.b[j] * problem.L->eval(&x[offset],
                                                                  &x[offset + offU],
                                                                  &x[offXUTotal],
                                                                  mesh.grid[i] + rk.c[j] * mesh.deltaT[i]);
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
        const auto diffMAY = problem.M->evalDiff(&x[offXUTotal - offXU], &x[offXUTotal - offU], &x[offXUTotal], mesh.tf);
        for (size_t k = 0; k < diffMAY[0].size(); k++) {
            grad_f[offXUTotal - offXU + problem.M->adj.indX[k]] += diffMAY[0][k];
        }
        for (size_t k = 0; k < diffMAY[1].size(); k++) {
            grad_f[offXUTotal - offU + problem.M->adj.indU[k]] += diffMAY[1][k];
        }
        for (size_t k = 0; k < diffMAY[2].size(); k++) {
            grad_f[offXUTotal + problem.M->adj.indP[k]] += diffMAY[2][k];
        }
    }

    // lagrange term derivative
    if (problem.L) {
        for(Index i = 0; i < mesh.intervals; i++){
            const double deltaT_i = mesh.deltaT[i];
            const double grid_i = mesh.grid[i];
            for (Index j = 0; j < rk.steps; j++){
                const int offset = i * offXUBlock + j * offXU; // index of 1st x var at collocation point (i,j)
                const double quadCoeff = mesh.deltaT[i] * rk.b[j];
                const auto diffLAG = problem.L->evalDiff(&x[offset], &x[offset + offU], &x[offXUTotal], grid_i + rk.c[j] * deltaT_i);

                for (size_t k = 0; k < diffLAG[0].size(); k++) {
                    grad_f[offset + problem.L->adj.indX[k]] += quadCoeff * diffLAG[0][k];
                }
                for (size_t k = 0; k < diffLAG[1].size(); k++) {
                    grad_f[offset + offU + problem.L->adj.indU[k]] += quadCoeff * diffLAG[1][k];
                }
                for (size_t k = 0; k < diffLAG[2].size(); k++) {
                    grad_f[offXUTotal + problem.L->adj.indP[k]] += quadCoeff * diffLAG[2][k];
                }
            }
        }
    }
    return true;
}

bool GDOP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) {
    // TODO: constr eval
    return true;
}

bool GDOP::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol,
                      Number *values) {
    // TODO: constr evalDiff
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

GDOP::GDOP(Problem problem, Mesh &mesh, Integrator &rk) : problem(std::move(problem)), mesh(mesh), rk(rk) {
}
