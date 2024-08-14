//
// Created by linus on 14.08.24.
//

#ifndef GDOP_H
#define GDOP_H

#include <memory>
#include "problem.h"
#include "mesh.h"
#include "integrator.h"

enum class InitVars {
    CONST,
    CALLBACK
};

class GDOP;

GDOP *create_gdop(const std::shared_ptr<const Problem> & problem, Mesh &mesh, Integrator &rk, InitVars initVars);

#endif //GDOP_H
