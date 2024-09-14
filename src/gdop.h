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
