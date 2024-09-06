
// CODEGEN FOR MODEL "analyticHypersensitive"

// includes
#define _USE_MATH_DEFINES
#include "analyticHypersensitiveGeneratedParams.h"
#include <cmath>
#include <string>
#include "constants.h"
#include <problem.h>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"


// runtime parameters and global constants


// lagrange term
class LagrangeanalyticHypersensitive : public Expression {
public:
	static std::unique_ptr<LagrangeanalyticHypersensitive> create() {
		Adjacency adj{{0}, {0}, {}};
		AdjacencyDiff adjDiff{{{0, 0}}, {}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<LagrangeanalyticHypersensitive>(new LagrangeanalyticHypersensitive(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.5*(pow(u[0], 2) + pow(x[0], 2));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{1.0*x[0]}, {1.0*u[0]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{1.0}, {}, {1.0}, {}, {}, {}};
	}
private:
	LagrangeanalyticHypersensitive(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0analyticHypersensitive : public Expression {
public:
	static std::unique_ptr<F0analyticHypersensitive> create() {
		Adjacency adj{{0}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F0analyticHypersensitive>(new F0analyticHypersensitive(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return u[0] - x[0];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	F0analyticHypersensitive(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// final constraints
class R0analyticHypersensitive : public Constraint {
public:
	static std::unique_ptr<R0analyticHypersensitive> create() {
		Adjacency adj{{0}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<R0analyticHypersensitive>(new R0analyticHypersensitive(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 1.0 - x[0];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	R0analyticHypersensitive(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


std::vector<double> uInitialGuess(double t) {
	 return {0};
};

Problem createProblem_analyticHypersensitive() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0analyticHypersensitive::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0analyticHypersensitive::create());
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            1, 1, 0,  // #vars
            {1.5},  // x0
            {MINUS_INFINITY},  // lb x
            {PLUS_INFINITY},  // ub x
            &uInitialGuess,  // u0 initial guesses for optimization
            {MINUS_INFINITY},  // lb u
            {PLUS_INFINITY},  // ub u
            {},  // p0 initial guesses for optimization
            {},  // lb p
            {},  // ub p
            {},
            LagrangeanalyticHypersensitive::create(),
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "analyticHypersensitive");
            
    #ifdef INITIAL_STATES_PATH
    problem.initialStatesPath = INITIAL_STATES_PATH "/initialValues.csv";
    #endif
    
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_analyticHypersensitive());
    InitVars initVars = INIT_VARS;
    Integrator rk = Integrator::radauIIA(RADAU_INTEGRATOR);
    Mesh mesh = Mesh::createEquidistantMesh(INTERVALS, FINAL_TIME);
    LinearSolver linearSolver = LINEAR_SOLVER;
    MeshAlgorithm meshAlgorithm = MESH_ALGORITHM;
    int meshIterations = MESH_ITERATIONS;

    Solver solver = Solver(create_gdop(problem, mesh, rk, initVars), meshIterations, linearSolver, meshAlgorithm);

    // set solver flags
    #ifdef EXPORT_OPTIMUM_PATH
    solver.setExportOptimumPath(EXPORT_OPTIMUM_PATH);
    #endif
    
    #ifdef EXPORT_HESSIAN_PATH
    solver.setExportHessianPath(EXPORT_HESSIAN_PATH);
    #endif
    
    #ifdef EXPORT_JACOBIAN_PATH
    solver.setExportJacobianPath(EXPORT_JACOBIAN_PATH);
    #endif
    
    #ifdef TOLERANCE
    solver.setTolerance(TOLERANCE);
    #endif
    
    // set solver mesh parameters
    #ifdef LEVEL
    solver.setMeshParameter("level", LEVEL);
    #endif
    
    #ifdef C_TOL
    solver.setMeshParameter("ctol", C_TOL);
    #endif
    
    #ifdef SIGMA
    solver.setMeshParameter("sigma", SIGMA);
    #endif
    
    // optimize
    int status = solver.solve();
    return status;
}        
        