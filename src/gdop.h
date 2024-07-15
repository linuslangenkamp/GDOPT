//
// Created by Linus on 14.07.2024.
//

#ifndef IPOPT_DO_GDOP_H
#define IPOPT_DO_GDOP_H

#include <IpTNLP.hpp>
#include "problem.h"
#include "mesh.h"
#include "integrator.h"


using namespace Ipopt;

class GDOP : public TNLP {
public:
    GDOP(Problem problem, Mesh& mesh, Integrator& rk);

    Problem problem;
    Mesh mesh;
    Integrator rk;

    const int offX = problem.sizeX;
    const int offU = problem.sizeU;
    const int offXU = problem.sizeX + problem.sizeU; // number of vars for one collocation knod
    const int offXUBlock = (problem.sizeX + problem.sizeU) * rk.steps;  // number of vars per interval
    const int offXUTotal = (problem.sizeX + problem.sizeU) * rk.steps * mesh.intervals; // first const parameter variable
    const int numberVars = (problem.sizeX + problem.sizeU) * rk.steps * mesh.intervals + problem.sizeP;

    bool get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) override;

    bool get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) override;

    bool get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index m,
                            bool init_lambda, Number *lambda) override;

    bool eval_f(Index n, const Number *x, bool new_x, Number &obj_value) override;

    bool eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) override;

    bool eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) override;

    bool eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol,
                    Number *values) override;


    bool eval_h(Ipopt::Index n, const Ipopt::Number *x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m,
                const Ipopt::Number *lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index *iRow,
                Ipopt::Index *jCol, Ipopt::Number *values) override;

    void finalize_solution(SolverReturn status, Index n, const Number *x, const Number *z_L, const Number *z_U, Index m,
                           const Number *g, const Number *lambda, Number obj_value, const IpoptData *ip_data,
                           IpoptCalculatedQuantities *ip_cq) override;
};


#endif //IPOPT_DO_GDOP_H
