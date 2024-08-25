
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
        const double x1 = pow(x0, 2);
        const double x2 = x1*MP_VALUE;
        const double x3 = MS_VALUE + x2;
        const double x4 = 9.81*MP_VALUE;
        const double x5 = cos(x[2]);
        const double x6 = x0*x5;
        const double x7 = pow(x[3], 2);
        const double x8 = x7*R_VALUE;
        const double x9 = x0*x8;
        const double x10 = (-x4*x6 - x9*MP_VALUE)/pow(x3, 2);
        const double x11 = pow(x3, -1);
        const double x12 = x0*x11*MP_VALUE;
		return {std::vector<double>{x11*(9.81*x2 - x4*pow(x5, 2) - x5*x8*MP_VALUE) - 2*x6*x10*MP_VALUE, -2*x12*x[3]*R_VALUE}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = sin(x[2]);
        const double x1 = cos(x[2]);
        const double x2 = x0*x1;
        const double x3 = x2*MP_VALUE;
        const double x4 = pow(x[3], 2);
        const double x5 = x0*x4;
        const double x6 = x5*R_VALUE;
        const double x7 = x6*MP_VALUE;
        const double x8 = pow(x0, 2);
        const double x9 = x8*MP_VALUE;
        const double x10 = MS_VALUE + x9;
        const double x11 = pow(x10, -1);
        const double x12 = pow(x1, 2);
        const double x13 = 9.81*x2;
        const double x14 = -x7 - x13*MP_VALUE;
        const double x15 = 2*x14;
        const double x16 = pow(x10, -2);
        const double x17 = x16*MP_VALUE;
        const double x18 = x15*x16;
        const double x19 = pow(x10, -3);
        const double x20 = x14*x19;
        const double x21 = x8*pow(MP_VALUE, 2);
        const double x22 = x1*x4;
        const double x23 = x22*R_VALUE;
        const double x24 = 9.81*x12;
        const double x25 = 9.81*x8;
        const double x26 = x16*(-x23*MP_VALUE - x24*MP_VALUE + x25*MP_VALUE);
        const double x27 = 4*x3;
        const double x28 = -x26 + x20*x27;
        const double x29 = pow(x0, 3);
        const double x30 = 4*x1;
        const double x31 = x16*(-x13 - x6);
        const double x32 = 2*x31;
        const double x33 = x11*(-x23 - x24 + x25) - x2*x18 - x3*x32 - x8*x26 + x30*x20*x29*MP_VALUE;
        const double x34 = x11*MP_VALUE;
        const double x35 = x21*x16;
        const double x36 = -x34*x22 + 2*x35*x22;
        const double x37 = 2*x[3];
        const double x38 = x34*R_VALUE;
        const double x39 = x0*x37;
        const double x40 = x39*R_VALUE;
        const double x41 = x40*x17;
        const double x42 = x29*x17;
        const double x43 = -x40*x11 + x42*x37*R_VALUE;
        const double x44 = -x34*x39;
        const double x45 = x15*x19;
        const double x46 = -x31 + x8*x45;
        const double x47 = x5*x17;
        const double x48 = x4*x42 - x5*x11;
		return {std::vector<double>{x11*(39.24*x3 + x7) - x26*x27 + x9*x18 - x15*x12*x17 + 8*x20*x21*x12, -x1*x38*x37 + x30*x35*x[3]*R_VALUE, -2*x0*x38}, {}, {}, {}, {}, {}};
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
        const double x3 = pow(x2, 2);
        const double x4 = x3*MP_VALUE;
        const double x5 = MS_VALUE + x4;
        const double x6 = pow(x5, -1);
        const double x7 = x2*x1;
        const double x8 = pow(x[3], 2);
        const double x9 = x8*R_VALUE;
        const double x10 = x2*x9;
        const double x11 = -x10*MP_VALUE - x7*MP_VALUE;
        const double x12 = u[0] + x6*x11;
        const double x13 = x11/pow(x5, 2);
        const double x14 = x0*MP_VALUE;
        const double x15 = 2*x2*x14;
        const double x16 = pow(R_VALUE, -1);
        const double x17 = x0*x16;
		return {std::vector<double>{x16*(x1 - x0*(-x15*x13 + x6*(9.81*x4 - 9.81*pow(x0, 2)*MP_VALUE - x9*x14)) + x2*x12), x6*x15*x[3]}, {-x17}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = cos(x[2]);
        const double x1 = sin(x[2]);
        const double x2 = x0*MP_VALUE;
        const double x3 = x2*x1;
        const double x4 = pow(x[3], 2);
        const double x5 = x4*R_VALUE;
        const double x6 = x1*x5;
        const double x7 = x6*MP_VALUE;
        const double x8 = pow(x1, 2);
        const double x9 = x8*MP_VALUE;
        const double x10 = MS_VALUE + x9;
        const double x11 = pow(x10, -1);
        const double x12 = pow(x0, 2);
        const double x13 = x12*MP_VALUE;
        const double x14 = pow(x10, -2);
        const double x15 = 9.81*x1;
        const double x16 = x0*x15;
        const double x17 = -x7 - x16*MP_VALUE;
        const double x18 = x14*x17;
        const double x19 = 2*x18;
        const double x20 = x17/pow(x10, 3);
        const double x21 = pow(MP_VALUE, 2);
        const double x22 = x21*x12;
        const double x23 = x8*x22;
        const double x24 = x0*x5;
        const double x25 = 9.81*x12;
        const double x26 = 9.81*x8;
        const double x27 = -x24*MP_VALUE - x25*MP_VALUE + x26*MP_VALUE;
        const double x28 = x27*x14;
        const double x29 = u[0] + x11*x17;
        const double x30 = x0*x29;
        const double x31 = x1*x18;
        const double x32 = 2*x31;
        const double x33 = -x2*x32 + x27*x11;
        const double x34 = 2*x1;
        const double x35 = pow(R_VALUE, -1);
        const double x36 = x1*x35;
        const double x37 = 4*x20;
        const double x38 = pow(x1, 3);
        const double x39 = x2*x38;
        const double x40 = -x16 - x6;
        const double x41 = x40*x14;
        const double x42 = x2*x34;
        const double x43 = x0*(-x0*x32 + x11*(-x24 - x25 + x26) + x37*x39 - x41*x42 - x8*x28);
        const double x44 = x40*x11 - x8*x18;
        const double x45 = pow(R_VALUE, -2);
        const double x46 = -x45*(9.81*x0 - x0*x33 + x1*x29);
        const double x47 = x4*x11;
        const double x48 = x4*x14;
        const double x49 = 2*x8;
        const double x50 = x48*x49;
        const double x51 = x9*x47;
        const double x52 = x14*x[3];
        const double x53 = x11*x[3];
        const double x54 = 2*x53;
        const double x55 = -x52*x42;
        const double x56 = x53*x34;
        const double x57 = 2*x52;
        const double x58 = x0*x45;
        const double x59 = x0*x35;
        const double x60 = 2*x20;
        const double x61 = x49*x20;
        const double x62 = -x58*x18 - x2*x48*x36;
        const double x63 = x38*MP_VALUE;
        const double x64 = x1*x47;
        const double x65 = x63*x48;
        const double x66 = x58*x44;
		return {std::vector<double>{x35*(-x15 + x30 - x0*(x11*(39.24*x3 + x7) - x13*x19 + 8*x23*x20 - 4*x3*x28 + x9*x19) + x34*x33), -4*x52*x23 + x54*x13 - x9*x54, x42*x11}, {x36}, {}, {}, {}, {}};
	}
private:
	F3invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
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
            {0},  // u0 initial guesses for optimization
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
        