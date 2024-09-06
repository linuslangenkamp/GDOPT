
// CODEGEN FOR MODEL "invertedPendulum"

// includes
#define _USE_MATH_DEFINES
#include "invertedPendulumGeneratedParams.h"
#include <cmath>
#include <string>
#include "constants.h"
#include <problem.h>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"


// runtime parameters and global constants
const double Ms = MS_VALUE;
const double Mp = MP_VALUE;
const double R = R_VALUE;


// lagrange term
class LagrangeinvertedPendulum : public Expression {
public:
	static std::unique_ptr<LagrangeinvertedPendulum> create() {
		Adjacency adj{{2}, {}, {}};
		AdjacencyDiff adjDiff{{{2, 2}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<LagrangeinvertedPendulum>(new LagrangeinvertedPendulum(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return pow(sin((1.0/2.0)*x[2]), 2);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = (1.0/2.0)*x[2];
		return {std::vector<double>{sin(x0)*cos(x0)}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = (1.0/2.0)*x[2];
		return {std::vector<double>{(-1.0/2.0)*pow(sin(x0), 2) + (1.0/2.0)*pow(cos(x0), 2)}, {}, {}, {}, {}, {}};
	}
private:
	LagrangeinvertedPendulum(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0invertedPendulum : public Expression {
public:
	static std::unique_ptr<F0invertedPendulum> create() {
		Adjacency adj{{1}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F0invertedPendulum>(new F0invertedPendulum(std::move(adj), std::move(adjDiff)));
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
	F0invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F1invertedPendulum : public Expression {
public:
	static std::unique_ptr<F1invertedPendulum> create() {
		Adjacency adj{{2, 3}, {0}, {}};
		AdjacencyDiff adjDiff{{{2, 2}, {3, 2}, {3, 3}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F1invertedPendulum>(new F1invertedPendulum(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = sin(x[2]);
        const double x1 = x0*MP_VALUE;
		return u[0] + (-9.81*x1*cos(x[2]) - x1*pow(x[3], 2)*R_VALUE)/(MS_VALUE + pow(x0, 2)*MP_VALUE);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = sin(x[2]);
        const double x1 = pow(x0, 2)*MP_VALUE;
        const double x2 = MS_VALUE + x1;
        const double x3 = 9.81*MP_VALUE;
        const double x4 = cos(x[2]);
        const double x5 = x0*x4;
        const double x6 = R_VALUE*MP_VALUE;
        const double x7 = x6*pow(x[3], 2);
        const double x8 = pow(x2, -1);
		return {std::vector<double>{x8*(9.81*x1 - x4*x7 - pow(x4, 2)*x3) - 2*(-x0*x7 - x3*x5)*x5*MP_VALUE/pow(x2, 2), -2*x0*x6*x8*x[3]}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = sin(x[2]);
        const double x1 = cos(x[2]);
        const double x2 = x0*x1;
        const double x3 = x2*MP_VALUE;
        const double x4 = R_VALUE*MP_VALUE;
        const double x5 = x4*pow(x[3], 2);
        const double x6 = x0*x5;
        const double x7 = pow(x0, 2);
        const double x8 = x7*MP_VALUE;
        const double x9 = MS_VALUE + x8;
        const double x10 = pow(x9, -1);
        const double x11 = pow(x1, 2);
        const double x12 = 9.81*MP_VALUE;
        const double x13 = -x6 - x2*x12;
        const double x14 = x13*x11;
        const double x15 = pow(x9, -2);
        const double x16 = 2*x15;
        const double x17 = x7*pow(MP_VALUE, 2);
        const double x18 = 4*x15;
        const double x19 = x1*x[3];
        const double x20 = 2*x4*x10;
		return {std::vector<double>{x10*(39.24*x3 + x6) - x14*x16*MP_VALUE - x3*x18*(9.81*x8 - x1*x5 - x12*x11) + x8*x13*x16 + 8*x14*x17/pow(x9, 3), -x20*x19 + x19*x18*x17*R_VALUE, -x0*x20}, {}, {}, {}, {}, {}};
	}
private:
	F1invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F2invertedPendulum : public Expression {
public:
	static std::unique_ptr<F2invertedPendulum> create() {
		Adjacency adj{{3}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F2invertedPendulum>(new F2invertedPendulum(std::move(adj), std::move(adjDiff)));
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
	F2invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F3invertedPendulum : public Expression {
public:
	static std::unique_ptr<F3invertedPendulum> create() {
		Adjacency adj{{2, 3}, {0}, {}};
		AdjacencyDiff adjDiff{{{2, 2}, {3, 2}, {3, 3}}, {{0, 2}}, {}, {}, {}, {}};
		return std::unique_ptr<F3invertedPendulum>(new F3invertedPendulum(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = cos(x[2]);
        const double x1 = sin(x[2]);
        const double x2 = 9.81*x1;
		return (x2 - x0*(u[0] + (-x0*x2*MP_VALUE - x1*pow(x[3], 2)*R_VALUE*MP_VALUE)/(MS_VALUE + pow(x1, 2)*MP_VALUE)))/R_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = cos(x[2]);
        const double x1 = 9.81*x0;
        const double x2 = sin(x[2]);
        const double x3 = pow(x2, 2)*MP_VALUE;
        const double x4 = MS_VALUE + x3;
        const double x5 = pow(x4, -1);
        const double x6 = x2*MP_VALUE;
        const double x7 = pow(x[3], 2)*R_VALUE;
        const double x8 = -x1*x6 - x6*x7;
        const double x9 = 2*x0*x6;
        const double x10 = pow(R_VALUE, -1);
		return {std::vector<double>{x10*(x1 - x0*(x5*(9.81*x3 - 9.81*pow(x0, 2)*MP_VALUE - x0*x7*MP_VALUE) - x8*x9/pow(x4, 2)) + x2*(u[0] + x5*x8)), x5*x9*x[3]}, {-x0*x10}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = cos(x[2]);
        const double x1 = sin(x[2]);
        const double x2 = x0*MP_VALUE;
        const double x3 = x2*x1;
        const double x4 = pow(x[3], 2)*R_VALUE;
        const double x5 = x1*x4*MP_VALUE;
        const double x6 = pow(x1, 2);
        const double x7 = x6*MP_VALUE;
        const double x8 = MS_VALUE + x7;
        const double x9 = pow(x8, -1);
        const double x10 = pow(x0, 2);
        const double x11 = x10*MP_VALUE;
        const double x12 = 2*x11;
        const double x13 = 9.81*x1;
        const double x14 = -x5 - x2*x13;
        const double x15 = pow(x8, -2);
        const double x16 = x15*x14;
        const double x17 = 2*x7;
        const double x18 = x6*x10*pow(MP_VALUE, 2);
        const double x19 = -9.81*x11 + 9.81*x7 - x2*x4;
        const double x20 = 4*x15;
        const double x21 = 2*x1;
        const double x22 = x2*x21;
        const double x23 = pow(R_VALUE, -1);
        const double x24 = x1*x23;
        const double x25 = x9*x[3];
		return {std::vector<double>{x23*(-x13 + x0*(u[0] + x9*x14) - x0*(-x12*x16 + x17*x16 + x9*(39.24*x3 + x5) - x3*x20*x19 + 8*x14*x18/pow(x8, 3)) + x21*(-x22*x16 + x9*x19)), x25*x12 - x25*x17 - x20*x18*x[3], x9*x22}, {x24}, {}, {}, {}, {}};
	}
private:
	F3invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


std::vector<double> uInitialGuess(double t) {
	 return {0};
};

Problem createProblem_invertedPendulum() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0invertedPendulum::create());
    F.push_back(F1invertedPendulum::create());
    F.push_back(F2invertedPendulum::create());
    F.push_back(F3invertedPendulum::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            4, 1, 0,  // #vars
            {0, 0, -0.001 + acos(-1), 0},  // x0
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},  // lb x
            {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},  // ub x
            &uInitialGuess,  // u0 initial guesses for optimization
            {-2.5},  // lb u
            {2.5},  // ub u
            {},  // p0 initial guesses for optimization
            {},  // lb p
            {},  // ub p
            {},
            LagrangeinvertedPendulum::create(),
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "invertedPendulum");
            
    #ifdef INITIAL_STATES_PATH
    problem.initialStatesPath = INITIAL_STATES_PATH "/initialValues.csv";
    #endif
    
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_invertedPendulum());
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
        