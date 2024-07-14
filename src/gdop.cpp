//
// Created by Linus on 14.07.2024.
//

#include "gdop.h"
#include <cassert>


bool GDOP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) {
    return true;
}

bool GDOP::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) {
    return true;
}

bool GDOP::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index m,
                              bool init_lambda, Number *lambda) {
    return true;
}

bool GDOP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value) {
    double MAY = 0;
    double LAG = 0;
    if (this->problem.M) {
        MAY = this->problem.M->eval(&x[offXUTotal-offXU], &x[offXUTotal-offU], &x[offXUTotal], mesh.tf);
    }
    if (this->problem.L) {
        for( Index i = 0; i < mesh.intervals; i++ ){
            for ( Index j = 0; j < rk.steps; j++ ){
                LAG += mesh.deltaT[i] * rk.b[j] * problem.L->eval(&x[i * offXUBlock + j * offXU],
                                                                  &x[i * offXUBlock + j * offXU + offU],
                                                                  &x[offXUTotal],
                                                                  mesh.grid[i] + rk.c[j] * mesh.deltaT[i]);
            }
        }
    }
    obj_value = MAY + LAG;
    return true;
}

bool GDOP::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) {
    return true;
}

bool GDOP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) {
    return true;
}

bool GDOP::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol,
                      Number *values) {
    return true;
}

bool GDOP::eval_h(Index n, const Number *x, bool new_x, Number obj_factor, Index m, const Number *lambda,
                  bool new_lambda, Index nele_hess, Index *iRow, Index *jCol, Number *values) {
    return true;
}

void GDOP::finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L, const Number *z_U,
                             Index m, const Number *g, const Number *lambda, Number obj_value, const IpoptData *ip_data,
                             IpoptCalculatedQuantities *ip_cq) {
}

GDOP::GDOP(Problem problem, Mesh &mesh, Integrator &rk) : problem(std::move(problem)), mesh(mesh), rk(rk) {
}
