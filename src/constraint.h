//
// Created by Linus on 10.07.2024.
//

#ifndef IPOPT_DO_CONSTRAINT_H
#define IPOPT_DO_CONSTRAINT_H

#include "expression.h"

class Constraint : public Expression {
public:
    Constraint(Adjacency adj, double lb = 0.0, double ub = 0.0) : Expression(std::move(adj)), lb{lb}, ub{ub} {}

    virtual ~Constraint() = default;

    const double lb;
    const double ub;
};

#endif //IPOPT_DO_CONSTRAINT_H
