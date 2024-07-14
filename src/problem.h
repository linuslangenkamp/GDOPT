//
// Created by Linus on 10.07.2024.
//

#ifndef IPOPT_DO_PROBLEM_H
#define IPOPT_DO_PROBLEM_H

#include <vector>
#include <memory>
#include "expression.h"
#include "constraint.h"

class Problem {
public:
    Problem(int sizeX, int sizeU, int sizeP,
            std::vector<double> startX, std::vector<double> lbX, std::vector<double> ubX,
            std::vector<double> lbU, std::vector<double> ubU,
            std::vector<double> lbP, std::vector<double> ubP,
            std::unique_ptr<Expression> M, std::unique_ptr<Expression> L,
            std::vector<std::unique_ptr<Constraint>> F,
            std::vector<std::unique_ptr<Constraint>> G,
            std::vector<std::unique_ptr<Constraint>> R,
            std::vector<std::unique_ptr<Constraint>> A);

    const int sizeX;
    const int sizeU;
    const int sizeP;

    std::vector<double> startX;
    std::vector<double> lbX;
    std::vector<double> ubX;
    std::vector<double> lbU;
    std::vector<double> ubU;
    std::vector<double> lbP;
    std::vector<double> ubP;

    std::unique_ptr<Expression> M;                  // Mayer term
    std::unique_ptr<Expression> L;                  // Lagrange term
    std::vector<std::unique_ptr<Constraint>> F;     // Differential constraints - state dynamics
    std::vector<std::unique_ptr<Constraint>> G;     // Algebraic path constraints for states, control, parameters and time
    std::vector<std::unique_ptr<Constraint>> R;     // Algebraic final constraints
    std::vector<std::unique_ptr<Constraint>> A;     // Algebraic constraints for parameters only
};


#endif //IPOPT_DO_PROBLEM_H
