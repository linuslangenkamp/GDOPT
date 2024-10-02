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

#ifndef GDOPT_PROBLEM_H
#define GDOPT_PROBLEM_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "constraint.h"
#include "expression.h"

class Problem {
public:
    Problem(int sizeX, int sizeU, int sizeP, std::vector<double> x0, std::vector<double> lbX, std::vector<double> ubX,
            std::function<std::vector<double>(double)> uInitialGuess, std::vector<double> lbU, std::vector<double> ubU, std::vector<double> pInitialGuess,
            std::vector<double> lbP, std::vector<double> ubP, std::unique_ptr<Expression> M, std::unique_ptr<Expression> L,
            std::vector<std::unique_ptr<Expression>> F, std::vector<std::unique_ptr<Constraint>> G, std::vector<std::unique_ptr<Constraint>> R,
            std::vector<std::unique_ptr<ParamConstraint>> A, std::string name);

    const int sizeX;
    const int sizeU;
    const int sizeP;

    std::vector<double> x0;                                    // starting value for states
    std::vector<double> lbX;                                   // global lower bound on state vars
    std::vector<double> ubX;                                   // global upper bound on state vars
    std::function<std::vector<double>(double)> uInitialGuess;  // starting values for control vars - not forced, only initial values for optimization
    std::vector<double> lbU;                                   // global lower bound on control vars
    std::vector<double> ubU;                                   // global upper bound on control vars
    std::vector<double> pInitialGuess;                         // starting values for parameters - not forced, only initial values for optimization
    std::vector<double> lbP;                                   // global lower bound on parameters
    std::vector<double> ubP;                                   // global upper bound on parameters

    std::vector<double> nominalsX = {};  // nominal values for the states
    std::vector<double> nominalsU = {};  // nominal values for the inputs
    std::vector<double> nominalsP = {};  // nominal values for the parameters

    std::unique_ptr<Expression> M;                    // mayer term
    std::unique_ptr<Expression> L;                    // lagrange term
    std::vector<std::unique_ptr<Expression>> F;       // state dynamics, RHS of ODE
    std::vector<std::unique_ptr<Constraint>> G;       // algebraic path constraints for states, control, parameters and time
    std::vector<std::unique_ptr<Constraint>> R;       // algebraic final constraints
    std::vector<std::unique_ptr<ParamConstraint>> A;  // algebraic constraints for parameters only:

    double nominalObjective = 1;         // nominal value for the objective
    std::vector<double> nominalsF = {};  // nominal values for the state dynamics
    std::vector<double> nominalsG = {};  // nominal values for the path constraints
    std::vector<double> nominalsR = {};  // nominal values for the final constraints
    std::vector<double> nominalsA = {};  // nominal values for the parametric constraints

    std::string name;
    std::string initialStatesPath;

    bool linearObjective = false;                // true if M and L are linear
    bool linearConstraints = false;              // true if f, g, r and a are all linear
    bool quadraticObjLinearConstraints = false;  // true if f, g, r and a are all linear and M and L are at most quadratic
};

#endif  // GDOPT_PROBLEM_H
