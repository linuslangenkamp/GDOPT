#ifndef GDOPT_CONFIG_H
#define GDOPT_CONFIG_H

#include <optional>

#include "solver.h"

extern InitVars INIT_VARS;
extern IntegratorSteps RADAU_INTEGRATOR;
extern RefinementMethod REFINEMENT_METHOD;
extern LinearSolver LINEAR_SOLVER;
extern MeshAlgorithm MESH_ALGORITHM;

extern double TOLERANCE;
extern double FINAL_TIME;
extern int INTERVALS;
extern int MAX_ITERATIONS;
extern int MESH_ITERATIONS;
extern bool USER_SCALING;

extern bool KKT_ERROR_MU_GLOBALIZATION;

extern bool LINEAR_OBJECTIVE;
extern bool LINEAR_CONSTRAINTS;
extern bool QUADRATIC_OBJECTIVE_LINEAR_CONSTRAINTS;

extern std::string EXPORT_OPTIMUM_PATH;
extern std::string EXPORT_HESSIAN_PATH;
extern std::string EXPORT_JACOBIAN_PATH;
extern std::string INITIAL_STATES_PATH;

extern std::optional<int> IPOPT_PRINT_LEVEL;
extern std::optional<double> SIGMA;
extern std::optional<double> LEVEL;
extern std::optional<double> C_TOL;

std::unordered_map<std::string, std::string> readConfig(const std::string& filename);
void setGlobalStandardConfiguration(const std::unordered_map<std::string, std::string>& configMap);

#endif  // GDOPT_CONFIG_H
