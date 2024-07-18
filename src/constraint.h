#ifndef IPOPT_DO_CONSTRAINT_H
#define IPOPT_DO_CONSTRAINT_H

#include "expression.h"

class Constraint : public Expression {
public:
    explicit Constraint(Adjacency adj, double lb = 0.0, double ub = 0.0) : Expression(std::move(adj)), lb{lb}, ub{ub} {}

    virtual ~Constraint() = default;

    const double lb;
    const double ub;
};

class paramConstraint : public paramExpression {
public:
    explicit paramConstraint(paramAdjacency adj, double lb = 0.0, double ub = 0.0) : paramExpression(std::move(adj)), lb{lb}, ub{ub} {}

    virtual ~paramConstraint() = default;

    const double lb;
    const double ub;
};

#endif //IPOPT_DO_CONSTRAINT_H
