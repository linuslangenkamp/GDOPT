
// CODEGEN FOR MODEL "vanDerPol"

// includes
#define _USE_MATH_DEFINES
#include "vanDerPolGeneratedParams.h"
#include <cmath>
#include <string>
#include "constants.h"
#include <problem.h>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"


// runtime parameters and global constants
const double RP = RP_VALUE;


// lagrange term
class LagrangevanDerPol : public Expression {
public:
	static std::unique_ptr<LagrangevanDerPol> create() {
		Adjacency adj{{0, 1}, {0}, {}};
		AdjacencyDiff adjDiff{{{0, 0}, {1, 1}}, {}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<LagrangevanDerPol>(new LagrangevanDerPol(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return pow(u[0], 2)*RP_VALUE + pow(x[0], 2) + pow(x[1], 2);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{2*x[0], 2*x[1]}, {2*u[0]*RP_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{2, 2}, {}, {2*RP_VALUE}, {}, {}, {}};
	}
private:
	LagrangevanDerPol(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0vanDerPol : public Expression {
public:
	static std::unique_ptr<F0vanDerPol> create() {
		Adjacency adj{{0, 1}, {0}, {}};
		AdjacencyDiff adjDiff{{{1, 0}, {1, 1}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F0vanDerPol>(new F0vanDerPol(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return u[0] - x[1] + x[0]*(1 - pow(x[1], 2));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{1 - pow(x[1], 2), -1 - 2*x[0]*x[1]}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-2*x[1], -2*x[0]}, {}, {}, {}, {}, {}};
	}
private:
	F0vanDerPol(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F1vanDerPol : public Expression {
public:
	static std::unique_ptr<F1vanDerPol> create() {
		Adjacency adj{{0}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F1vanDerPol>(new F1vanDerPol(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return x[0];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	F1vanDerPol(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


Problem createProblem_vanDerPol() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0vanDerPol::create());
    F.push_back(F1vanDerPol::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            2, 1, 0,  // #vars
            {0, 1},  // x0
            {MINUS_INFINITY, MINUS_INFINITY},  // lb x
            {PLUS_INFINITY, PLUS_INFINITY},  // ub x
            {0},  // u0 initial guesses for optimization
            {MINUS_INFINITY},  // lb u
            {0.8},  // ub u
            {},  // p0 initial guesses for optimization
            {},  // lb p
            {},  // ub p
            {},
            LagrangevanDerPol::create(),
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "vanDerPol");
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_vanDerPol());
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
        