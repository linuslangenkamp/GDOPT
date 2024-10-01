#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>

#include "solver.h"

// declaration of global variables

// always set
InitVars INIT_VARS;
IntegratorSteps RADAU_INTEGRATOR;
RefinementMethod REFINEMENT_METHOD;
LinearSolver LINEAR_SOLVER;
MeshAlgorithm MESH_ALGORITHM;

int INTERVALS;
double FINAL_TIME;
int MESH_ITERATIONS;
double TOLERANCE;
bool USER_SCALING;

// pseudo optional, if null -> empty string ""
std::string EXPORT_OPTIMUM_PATH;
std::string EXPORT_HESSIAN_PATH;
std::string EXPORT_JACOBIAN_PATH;
std::string INITIAL_STATES_PATH;

// actual optionals
std::optional<double> SIGMA;
std::optional<double> LEVEL;
std::optional<double> C_TOL;

const std::unordered_map<std::string, InitVars> initVarsMap = {{"CONST", InitVars::CONST},
                                                               {"SOLVE", InitVars::SOLVE},
                                                               {"SOLVE_EXPLICIT", InitVars::SOLVE_EXPLICIT},
                                                               {"SOLVE_EXPLICIT_EULER", InitVars::SOLVE_EXPLICIT_EULER},
                                                               {"CALLBACK", InitVars::CALLBACK}};

const std::unordered_map<std::string, LinearSolver> linearSolverMap = {
    {"MUMPS", LinearSolver::MUMPS}, {"MA27", LinearSolver::MA27}, {"MA57", LinearSolver::MA57},      {"MA77", LinearSolver::MA77},
    {"MA86", LinearSolver::MA86},   {"MA97", LinearSolver::MA97}, {"PARDISO", LinearSolver::PARDISO}};

const std::unordered_map<std::string, RefinementMethod> refinementMethodMap = {{"POLYNOMIAL", RefinementMethod::POLYNOMIAL},
                                                                               {"LINEAR_SPLINE", RefinementMethod::LINEAR_SPLINE}};

const std::unordered_map<std::string, MeshAlgorithm> meshAlgorithmMap = {
    {"NONE", MeshAlgorithm::NONE}, {"BASIC", MeshAlgorithm::BASIC}, {"L2_BOUNDARY_NORM", MeshAlgorithm::L2_BOUNDARY_NORM}};

InitVars stringToInitVars(const std::string &str) {
    auto it = initVarsMap.find(str);
    if (it != initVarsMap.end()) {
        return it->second;
    }
    else {
        throw std::invalid_argument("Invalid InitVars configuration argument: " + str);
    }
}

LinearSolver stringToLinearSolver(const std::string &str) {
    auto it = linearSolverMap.find(str);
    if (it != linearSolverMap.end()) {
        return it->second;
    }
    else {
        throw std::invalid_argument("Invalid LinearSolver configuration argument: " + str);
    }
}

RefinementMethod stringToRefinementMethod(const std::string &str) {
    auto it = refinementMethodMap.find(str);
    if (it != refinementMethodMap.end()) {
        return it->second;
    }
    else {
        throw std::invalid_argument("Invalid RefinementMethod configuration argument: " + str);
    }
}

MeshAlgorithm stringToMeshAlgorithm(const std::string &str) {
    auto it = meshAlgorithmMap.find(str);
    if (it != meshAlgorithmMap.end()) {
        return it->second;
    }
    else {
        throw std::invalid_argument("Invalid MeshAlgorithm configuration argument: " + str);
    }
}

std::string trim(const std::string &str) {
    size_t first = str.find_first_not_of(" \t\"");
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(" \t\"");
    return str.substr(first, last - first + 1);
}

std::unordered_map<std::string, std::string> readConfig(const std::string &filename) {
    std::unordered_map<std::string, std::string> configMap;
    std::ifstream infile(filename);
    std::string line;

    while (std::getline(infile, line)) {
        line = trim(line);

        // skip
        if (line.empty() || line[0] == '#' || line.front() == '[' && line.back() == ']')
            continue;

        std::istringstream iss(line);
        std::string key, value;

        if (std::getline(iss, key, ' ') && std::getline(iss, value)) {
            key = trim(key);
            value = trim(value);
            configMap[key] = value;
        }
    }

    // printout model parameters to stdout
    std::cout << "Model parameters:\n";
    for (auto [arg, val] : configMap) {
        std::cout << arg << " -> " << val << "\n";
    }

    return configMap;
}

void setGlobalStdConfiguration(const std::unordered_map<std::string, std::string> &configMap) {
    // sets the global variables from the given configuration (excluding runtime parameters)

    // basic types
    INTERVALS = std::stoi(configMap.at("INTERVALS"));
    MESH_ITERATIONS = std::stoi(configMap.at("MESH_ITERATIONS"));
    FINAL_TIME = std::stod(configMap.at("FINAL_TIME"));
    TOLERANCE = std::stod(configMap.at("TOLERANCE"));
    USER_SCALING = configMap.at("USER_SCALING") == "true";

    // enums
    INIT_VARS = stringToInitVars(configMap.at("INIT_VARS"));
    LINEAR_SOLVER = stringToLinearSolver(configMap.at("LINEAR_SOLVER"));
    REFINEMENT_METHOD = stringToRefinementMethod(configMap.at("REFINEMENT_METHOD"));
    MESH_ALGORITHM = stringToMeshAlgorithm(configMap.at("MESH_ALGORITHM"));
    RADAU_INTEGRATOR = (IntegratorSteps)std::stoi(configMap.at("RADAU_INTEGRATOR"));

    // optional string arguments
    if ((configMap.find("EXPORT_OPTIMUM_PATH") != configMap.end())) {
        EXPORT_OPTIMUM_PATH = configMap.at("EXPORT_OPTIMUM_PATH");
    }
    if ((configMap.find("EXPORT_HESSIAN_PATH") != configMap.end())) {
        EXPORT_HESSIAN_PATH = configMap.at("EXPORT_HESSIAN_PATH");
    }
    if ((configMap.find("EXPORT_JACOBIAN_PATH") != configMap.end())) {
        EXPORT_JACOBIAN_PATH = configMap.at("EXPORT_JACOBIAN_PATH");
    }
    if ((configMap.find("INITIAL_STATES_PATH") != configMap.end())) {
        INITIAL_STATES_PATH = configMap.at("INITIAL_STATES_PATH");
    }

    // optional double mesh flags
    if ((configMap.find("SIGMA") != configMap.end())) {
        SIGMA = std::stod(configMap.at("SIGMA"));
    }
    if ((configMap.find("LEVEL") != configMap.end())) {
        LEVEL = std::stod(configMap.at("LEVEL"));
    }
    if ((configMap.find("C_TOL") != configMap.end())) {
        C_TOL = std::stod(configMap.at("C_TOL"));
    }
}
