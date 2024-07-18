#ifndef IPOPT_DO_PROBLEM_H
#define IPOPT_DO_PROBLEM_H

#include <vector>
#include <memory>
#include "expression.h"
#include "constraint.h"

class Problem {
public:
    Problem(int sizeX, int sizeU, int sizeP,
            std::vector<double> x0, std::vector<double> lbX, std::vector<double> ubX,
            std::vector<double> lbU, std::vector<double> ubU,
            std::vector<double> lbP, std::vector<double> ubP,
            std::unique_ptr<Expression> M, std::unique_ptr<Expression> L,
            std::vector<std::unique_ptr<Expression>> F,
            std::vector<std::unique_ptr<Constraint>> G,
            std::vector<std::unique_ptr<Constraint>> R,
            std::vector<std::unique_ptr<paramConstraint>> A);

    const int sizeX;
    const int sizeU;
    const int sizeP;

    std::vector<double> x0;                          // starting value for states
    std::vector<double> lbX;                         // global lower bound on state vars
    std::vector<double> ubX;                         // global upper bound on state vars
    std::vector<double> lbU;                         // global lower bound on control vars
    std::vector<double> ubU;                         // global upper bound on control vars
    std::vector<double> lbP;                         // global lower bound on parameters
    std::vector<double> ubP;                         // global upper bound on parameters

    std::unique_ptr<Expression> M;                   // mayer term
    std::unique_ptr<Expression> L;                   // lagrange term
    std::vector<std::unique_ptr<Expression>> F;      // state dynamics, RHS of ODE
    std::vector<std::unique_ptr<Constraint>> G;      // algebraic path constraints for states, control, parameters and time
    std::vector<std::unique_ptr<Constraint>> R;      // algebraic final constraints
    std::vector<std::unique_ptr<paramConstraint>> A; // algebraic constraints for parameters only:
};


#endif //IPOPT_DO_PROBLEM_H
