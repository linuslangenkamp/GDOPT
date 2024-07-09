//
// Created by Linus on 08.07.2024.
//

#include "test_nlp.h"
#include <cassert>

bool TestNLP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, TNLP::IndexStyleEnum &index_style) {
    n = 3;
    m = 2;
    nnz_jac_g = 4;
    nnz_h_lag = 3;
    index_style = TNLP::C_STYLE;
    return true;
}

bool TestNLP::get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) {
    assert(n == 3);
    assert(m == 2);
    for( Index i = 0; i < n; i++ )
    {
        x_l[i] = 0.0;
        x_u[i] = 10.0;
    }
    g_l[0] = 0.0;
    g_u[0] = 5.0;
    g_l[1] = 1.0;
    g_u[1] = 12.0;
    return true;
}

bool TestNLP::get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index m,
                                 bool init_lambda, Number *lambda) {
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);
    x[0] = 2.0;
    x[1] = 2.0;
    x[2] = 2.0;
    return true;
}

bool TestNLP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value) {
    assert(n == 3);
    obj_value = x[0] * x[1] + x[2] * x[2];
    return true;
}

bool TestNLP::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) {
    assert(n == 3);
    grad_f[0] = x[1];
    grad_f[1] = x[0];
    grad_f[2] = 2 * x[2];;
    return true;
}

bool TestNLP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) {
    assert(n == 3);
    assert(m == 2);
    g[0] = x[0] + x[1];
    g[1] = x[0] * x[0] + x[2] * x[2];
    return true;
}

bool TestNLP::eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol,
                         Number *values) {
    assert(n == 3);
    assert(m == 2);
    if( values == NULL )
    {
        iRow[0] = 0;
        jCol[0] = 0;
        iRow[1] = 0;
        jCol[1] = 1;
        iRow[2] = 1;
        jCol[2] = 0;
        iRow[3] = 1;
        jCol[3] = 2;
    }
    else
    {
        values[0] = 1;
        values[1] = 1;
        values[2] = 2 * x[0];
        values[3] = 2 * x[2];
    }
    return true;
}

bool TestNLP::eval_h(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
                     const Ipopt::Number *lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index *iRow,
                     Ipopt::Index *jCol, Ipopt::Number *values) {
    if( values == NULL )
    {
        iRow[0] = 1;
        jCol[0] = 0;
        iRow[1] = 2;
        jCol[1] = 0;
        iRow[2] = 2;
        jCol[2] = 2;
        // later with assert some idx == nele_hess -> give structure of hessian
    }
    else
    {
        values[0] = obj_factor;
        values[1] = lambda[1];
        values[2] = 2 * obj_factor;
    }
    return true;
}

void
TestNLP::finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L, const Number *z_U, Index m,
                           const Number *g, const Number *lambda, Number obj_value, const IpoptData *ip_data,
                           IpoptCalculatedQuantities *ip_cq) {

}
