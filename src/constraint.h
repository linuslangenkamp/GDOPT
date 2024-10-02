/**
 * GDOPT - General Dynamic Optimizer
 * Copyright (C) 2024  Linus Langenkamp
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/

#ifndef GDOPT_CONSTRAINT_H
#define GDOPT_CONSTRAINT_H

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

#endif  // GDOPT_CONSTRAINT_H
