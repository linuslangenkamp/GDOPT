#ifndef IPOPT_DO_CONSTRAINT_H
#define IPOPT_DO_CONSTRAINT_H

#include "expression.h"

class Constraint : public Expression {
public:
    explicit Constraint(Adjacency adj, AdjacencyDiff adjDiff, double lb = 0.0, double ub = 0.0)
        : Expression(std::move(adj), std::move(adjDiff)), lb{lb}, ub{ub} {
    }

    virtual ~Constraint() = default;

    const double lb;
    const double ub;
};

class ParamConstraint : public ParamExpression {
public:
    explicit ParamConstraint(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb = 0.0, double ub = 0.0)
        : ParamExpression(std::move(adj), std::move(adjDiff)), lb{lb}, ub{ub} {
    }

    virtual ~ParamConstraint() = default;

    const double lb;
    const double ub;
};

#endif  // IPOPT_DO_CONSTRAINT_H
