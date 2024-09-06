
// CODEGEN FOR MODEL "benchmark"

// includes
#define _USE_MATH_DEFINES
#include "benchmarkGeneratedParams.h"
#include <cmath>
#include <string>
#include "constants.h"
#include <problem.h>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"


// runtime parameters and global constants


// mayer term
class Mayerbenchmark : public Expression {
public:
	static std::unique_ptr<Mayerbenchmark> create() {
		Adjacency adj{{3}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<Mayerbenchmark>(new Mayerbenchmark(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return x[3];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	Mayerbenchmark(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0benchmark : public Expression {
public:
	static std::unique_ptr<F0benchmark> create() {
		Adjacency adj{{1}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F0benchmark>(new F0benchmark(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return x[1];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	F0benchmark(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F1benchmark : public Expression {
public:
	static std::unique_ptr<F1benchmark> create() {
		Adjacency adj{{2}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 2}}, {}, {}, {}, {}};
		return std::unique_ptr<F1benchmark>(new F1benchmark(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -8 + 16*t - u[0]*x[2];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-u[0]}, {-x[2]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {-1}, {}, {}, {}, {}};
	}
private:
	F1benchmark(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F2benchmark : public Expression {
public:
	static std::unique_ptr<F2benchmark> create() {
		Adjacency adj{{}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F2benchmark>(new F2benchmark(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return u[0];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	F2benchmark(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F3benchmark : public Expression {
public:
	static std::unique_ptr<F3benchmark> create() {
		Adjacency adj{{0, 1, 2}, {0}, {}};
		AdjacencyDiff adjDiff{{{0, 0}, {1, 1}, {2, 1}, {2, 2}}, {{0, 1}, {0, 2}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F3benchmark>(new F3benchmark(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return pow(x[0], 2) + pow(x[1], 2) + 0.0005*pow(-8 + 16*t + x[1] - 0.1*pow(u[0], 2)*x[2], 2);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], 2);
        const double x1 = -8 + 16*t + x[1] - 0.1*x0*x[2];
		return {std::vector<double>{2*x[0], 0.001*x1 + 2*x[1], -0.0001*x0*x1}, {-0.0002*x1*u[0]*x[2]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = -0.0002*u[0]*x[2];
        const double x1 = pow(u[0], 2);
        const double x2 = 0.0002*(-8 + 16*t + x[1] - 0.1*x1*x[2]);
        const double x3 = 2e-05*pow(u[0], 3)*x[2] - x2*u[0];
		return {std::vector<double>{2, 2.001, -0.0001*x1, 1e-05*pow(u[0], 4)}, {x0, x3}, {4e-05*x1*pow(x[2], 2) - x2*x[2]}, {}, {}, {}};
	}
private:
	F3benchmark(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


std::vector<double> uInitialGuess(double t) {
	 return {7.5};
};

Problem createProblem_benchmark() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0benchmark::create());
    F.push_back(F1benchmark::create());
    F.push_back(F2benchmark::create());
    F.push_back(F3benchmark::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            4, 1, 0,  // #vars
            {0, -1, -sqrt(5), 0},  // x0
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},  // lb x
            {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},  // ub x
            &uInitialGuess,  // u0 initial guesses for optimization
            {-4},  // lb u
            {10},  // ub u
            {},  // p0 initial guesses for optimization
            {},  // lb p
            {},  // ub p
            Mayerbenchmark::create(),
            {},
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "benchmark");
            
    #ifdef INITIAL_STATES_PATH
    problem.initialStatesPath = INITIAL_STATES_PATH "/initialValues.csv";
    #endif
    
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_benchmark());
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
        