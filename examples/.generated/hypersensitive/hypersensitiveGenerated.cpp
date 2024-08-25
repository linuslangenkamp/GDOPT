
// CODEGEN FOR MODEL "hypersensitive"

// includes
#define _USE_MATH_DEFINES
#include "hypersensitiveGeneratedParams.h"
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
class Lagrangehypersensitive : public Expression {
public:
	static std::unique_ptr<Lagrangehypersensitive> create() {
		Adjacency adj{{0}, {0}, {}};
		AdjacencyDiff adjDiff{{{0, 0}}, {}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<Lagrangehypersensitive>(new Lagrangehypersensitive(std::move(adj), std::move(adjDiff)));
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
	Lagrangehypersensitive(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0hypersensitive : public Expression {
public:
	static std::unique_ptr<F0hypersensitive> create() {
		Adjacency adj{{0}, {0}, {}};
		AdjacencyDiff adjDiff{{{0, 0}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F0hypersensitive>(new F0hypersensitive(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return u[0] - pow(x[0], 3);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-3*pow(x[0], 2)}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-6*x[0]}, {}, {}, {}, {}, {}};
	}
private:
	F0hypersensitive(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// final constraints
class R0hypersensitive : public Constraint {
public:
	static std::unique_ptr<R0hypersensitive> create() {
		Adjacency adj{{0}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<R0hypersensitive>(new R0hypersensitive(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 1.5 - x[0];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	R0hypersensitive(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


Problem createProblem_hypersensitive() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0hypersensitive::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0hypersensitive::create());
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            1, 1, 0,  // #vars
            {1},  // x0
            {MINUS_INFINITY},  // lb x
            {PLUS_INFINITY},  // ub x
            {0},  // u0 initial guesses for optimization
            {MINUS_INFINITY},  // lb u
            {PLUS_INFINITY},  // ub u
            {},  // p0 initial guesses for optimization
            {},  // lb p
            {},  // ub p
            {},
            Lagrangehypersensitive::create(),
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "hypersensitive");
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_hypersensitive());
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
        