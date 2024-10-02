#ifndef GDOPT_CONFIG_H
#define GDOPT_CONFIG_H

#include <optional>

#include "solver.h"

extern InitVars INIT_VARS;
extern IntegratorSteps RADAU_INTEGRATOR;
extern RefinementMethod REFINEMENT_METHOD;
extern LinearSolver LINEAR_SOLVER;
extern MeshAlgorithm MESH_ALGORITHM;

extern int INTERVALS;
extern double FINAL_TIME;
extern int MESH_ITERATIONS;
extern double TOLERANCE;
extern bool USER_SCALING;

extern std::string EXPORT_OPTIMUM_PATH;
extern std::string EXPORT_HESSIAN_PATH;
extern std::string EXPORT_JACOBIAN_PATH;
extern std::string INITIAL_STATES_PATH;

extern std::optional<double> SIGMA;
extern std::optional<double> LEVEL;
extern std::optional<double> C_TOL;

std::unordered_map<std::string, std::string> readConfig(const std::string& filename);
void setGlobalStdConfiguration(const std::unordered_map<std::string, std::string>& configMap);

#endif  // GDOPT_CONFIG_H
