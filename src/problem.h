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
    Problem(int sizeX, int sizeU, int sizeP, int dimF, int dimG, int dimR, int dimA,
            std::unique_ptr<Expression> M, std::unique_ptr<Expression> L,
            std::vector<std::unique_ptr<Constraint>> F,
            std::vector<std::unique_ptr<Constraint>> G,
            std::vector<std::unique_ptr<Constraint>> R,
            std::vector<std::unique_ptr<Constraint>> A)
            : sizeX(sizeX), sizeU(sizeU), sizeP(sizeP), dimF(dimF), dimG(dimG), dimR(dimR), dimA(dimA),
              M(std::move(M)), L(std::move(L)), F(std::move(F)), G(std::move(G)), R(std::move(R)), A(std::move(A)) {}

    const int sizeX;
    const int sizeU;
    const int sizeP;
    const int dimF;
    const int dimG;
    const int dimR;
    const int dimA;
    std::unique_ptr<Expression> M;                  // Mayer term
    std::unique_ptr<Expression> L;                  // Lagrange term
    std::vector<std::unique_ptr<Constraint>> F;     // Differential constraints - state dynamics
    std::vector<std::unique_ptr<Constraint>> G;     // Algebraic path constraints for states, control, parameters and time
    std::vector<std::unique_ptr<Constraint>> R;     // Algebraic final constraints
    std::vector<std::unique_ptr<Constraint>> A;     // Algebraic constraints for parameters only

};


#endif //IPOPT_DO_PROBLEM_H
