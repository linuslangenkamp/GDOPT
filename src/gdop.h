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

#ifndef GDOP_H
#define GDOP_H

#include <memory>

#include "integrator.h"
#include "mesh.h"
#include "problem.h"

enum class InitVars {
    CONST,                 // states will be constant globally
    SOLVE,                 // solving the ode with the given RadauIIA scheme
    SOLVE_EXPLICIT,        // solving the ode with the classic Runge-Kutta method
    SOLVE_EXPLICIT_EULER,  // solving the ode with the explicit Euler method
    CALLBACK               // callback case for recursive calls (only called from Solver)
};

class GDOP;

GDOP* create_gdop(const std::shared_ptr<const Problem>& problem, Mesh& mesh, Integrator& rk, InitVars initVars);

#endif  // GDOP_H
