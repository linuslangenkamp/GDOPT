
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
		return {std::vector<double>{-1.0/2.0*pow(sin(x0), 2) + (1.0/2.0)*pow(cos(x0), 2)}, {}, {}, {}, {}, {}};
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
        const double x1 = MP_VALUE*pow(x0, 2) + MS_VALUE;
		return (-MP_VALUE*R_VALUE*x0*pow(x[3], 2) - 4.9050000000000002*MP_VALUE*sin(2*x[2]) + u[0]*x1)/x1;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = sin(x[2]);
        const double x1 = pow(x[3], 2);
        const double x2 = R_VALUE*x1;
        const double x3 = x0*x2;
        const double x4 = 2*x[2];
        const double x5 = 4.9050000000000002*sin(x4);
        const double x6 = pow(x0, 2);
        const double x7 = MP_VALUE*x6 + MS_VALUE;
        const double x8 = (-MP_VALUE*x3 - MP_VALUE*x5 + u[0]*x7)/pow(x7, 2);
        const double x9 = MP_VALUE*cos(x[2]);
        const double x10 = 2*x0*x9;
        const double x11 = 1.0/x7;
        const double x12 = MP_VALUE*x0*x11;
		return {std::vector<double>{-x10*x8 + x11*(-9.8100000000000005*MP_VALUE*cos(x4) + u[0]*x10 - x2*x9), -2*R_VALUE*x12*x[3]}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = sin(x[2]);
        const double x1 = pow(x0, 2);
        const double x2 = MP_VALUE*x1;
        const double x3 = MS_VALUE + x2;
        const double x4 = 1.0/x3;
        const double x5 = pow(x[3], 2);
        const double x6 = R_VALUE*x5;
        const double x7 = x0*x6;
        const double x8 = MP_VALUE*x7;
        const double x9 = 2*x[2];
        const double x10 = sin(x9);
        const double x11 = cos(x[2]);
        const double x12 = pow(x11, 2);
        const double x13 = 2*u[0];
        const double x14 = MP_VALUE*x13;
        const double x15 = 4.9050000000000002*x10;
        const double x16 = -MP_VALUE*x15 + u[0]*x3 - x8;
        const double x17 = 2*x16;
        const double x18 = pow(x3, -2);
        const double x19 = MP_VALUE*x18;
        const double x20 = x17*x18;
        const double x21 = x0*x11;
        const double x22 = 9.8100000000000005*cos(x9);
        const double x23 = x11*x6;
        const double x24 = x18*(-MP_VALUE*x22 - MP_VALUE*x23 + x14*x21);
        const double x25 = MP_VALUE*x21;
        const double x26 = 4*x25;
        const double x27 = pow(MP_VALUE, 2);
        const double x28 = pow(x3, -3);
        const double x29 = x16*x28;
        const double x30 = x13*x18;
        const double x31 = -x24 - x25*x30 + x26*x29;
        const double x32 = u[0]*x1;
        const double x33 = x18*(-x15 + x32 - x7);
        const double x34 = pow(x0, 3);
        const double x35 = 4*x11;
        const double x36 = MP_VALUE*x29*x34*x35 - x1*x24 - x20*x21 - 2*x25*x33 + x4*(x13*x21 - x22 - x23);
        const double x37 = MP_VALUE*x4;
        const double x38 = x11*x5;
        const double x39 = 2*x1;
        const double x40 = x18*x27;
        const double x41 = -x37*x38 + x38*x39*x40;
        const double x42 = R_VALUE*x[3];
        const double x43 = 2*x42;
        const double x44 = 2*x0;
        const double x45 = x37*x44;
        const double x46 = x42*x44;
        const double x47 = x19*x46;
        const double x48 = x19*x34;
        const double x49 = -x4*x46 + x43*x48;
        const double x50 = -x45*x[3];
        const double x51 = x17*x28;
        const double x52 = x1*x51 - x18*x32 - x33;
        const double x53 = x0*x5;
        const double x54 = x19*x53;
        const double x55 = -x4*x53 + x48*x5;
		return {std::vector<double>{8*x1*x12*x27*x29 - x12*x17*x19 + x2*x20 - x24*x26 + x4*(19.620000000000001*MP_VALUE*x10 + x12*x14 - x13*x2 + x8), x1*x35*x40*x42 - x11*x37*x43, -R_VALUE*x45}, {}, {}, {}, {}, {}};
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
        const double x0 = sin(x[2]);
        const double x1 = MP_VALUE*pow(x0, 2) + MS_VALUE;
		return (9.8100000000000005*x0*x1 - (-MP_VALUE*R_VALUE*x0*pow(x[3], 2) - 4.9050000000000002*MP_VALUE*sin(2*x[2]) + u[0]*x1)*cos(x[2]))/(R_VALUE*x1);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = cos(x[2]);
        const double x1 = 1.0/R_VALUE;
        const double x2 = x0*x1;
        const double x3 = sin(x[2]);
        const double x4 = pow(x3, 2);
        const double x5 = MP_VALUE*x4;
        const double x6 = MS_VALUE + x5;
        const double x7 = 9.8100000000000005*x3;
        const double x8 = pow(x[3], 2);
        const double x9 = R_VALUE*x8;
        const double x10 = x3*x9;
        const double x11 = 2*x[2];
        const double x12 = 4.9050000000000002*sin(x11);
        const double x13 = -MP_VALUE*x10 - MP_VALUE*x12 + u[0]*x6;
        const double x14 = -x0*x13 + x6*x7;
        const double x15 = x14/pow(x6, 2);
        const double x16 = MP_VALUE*x3;
        const double x17 = 2*x16;
        const double x18 = u[0]*x0;
        const double x19 = MP_VALUE*x0;
        const double x20 = 1.0/x6;
        const double x21 = x1*x20;
        const double x22 = x1*x15;
		return {std::vector<double>{-x15*x17*x2 + x21*(19.620000000000001*x0*x5 + 9.8100000000000005*x0*x6 - x0*(-9.8100000000000005*MP_VALUE*cos(x11) + x17*x18 - x19*x9) + x13*x3), 2*x19*x20*x3*x[3]}, {-x2}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = sin(x[2]);
        const double x1 = pow(x0, 3);
        const double x2 = 19.620000000000001*MP_VALUE;
        const double x3 = cos(x[2]);
        const double x4 = pow(x[3], 2);
        const double x5 = R_VALUE*x4;
        const double x6 = x0*x5;
        const double x7 = MP_VALUE*x6;
        const double x8 = 2*x[2];
        const double x9 = sin(x8);
        const double x10 = pow(x3, 2);
        const double x11 = MP_VALUE*x10;
        const double x12 = 2*x11;
        const double x13 = pow(x0, 2);
        const double x14 = u[0]*x13;
        const double x15 = 2*MP_VALUE;
        const double x16 = MP_VALUE*x13;
        const double x17 = MS_VALUE + x16;
        const double x18 = x0*x17;
        const double x19 = 9.8100000000000005*x18;
        const double x20 = 4.9050000000000002*x9;
        const double x21 = -MP_VALUE*x20 + u[0]*x17 - x7;
        const double x22 = x21*x3;
        const double x23 = u[0]*x0;
        const double x24 = x23*x3;
        const double x25 = 9.8100000000000005*cos(x8);
        const double x26 = x3*x5;
        const double x27 = -MP_VALUE*x25 - MP_VALUE*x26 + x15*x24;
        const double x28 = 2*x0;
        const double x29 = 1.0/R_VALUE;
        const double x30 = 1.0/x17;
        const double x31 = x29*x30;
        const double x32 = x19 - x22;
        const double x33 = pow(x17, -2);
        const double x34 = x32*x33;
        const double x35 = x29*x34;
        const double x36 = 2*x16;
        const double x37 = 9.8100000000000005*x3;
        const double x38 = x0*x21 + 19.620000000000001*x16*x3 + x17*x37 - x27*x3;
        const double x39 = x0*x29;
        const double x40 = MP_VALUE*x3;
        const double x41 = x33*x40;
        const double x42 = pow(x17, -3);
        const double x43 = x29*x42;
        const double x44 = pow(MP_VALUE, 2)*x10;
        const double x45 = x13*x44;
        const double x46 = x29*x33;
        const double x47 = x38*x46;
        const double x48 = -u[0]*x3 + 9.8100000000000005*x0;
        const double x49 = x15*x3;
        const double x50 = x33*x39*x49;
        const double x51 = 4*x32;
        const double x52 = x31*(x23 + x37) + x39*x40*x42*x51 - x47 - x48*x50;
        const double x53 = x14 - x20 - x6;
        const double x54 = 9.8100000000000005*x1 - x3*x53;
        const double x55 = x1*x40;
        const double x56 = -x13*x47 - 2*x3*x34*x39 + x31*(x0*x53 + 29.43*x13*x3 - x3*(2*x24 - x25 - x26)) + x43*x51*x55 - x50*x54;
        const double x57 = x11*x4;
        const double x58 = x16*x4;
        const double x59 = pow(R_VALUE, -2);
        const double x60 = x30*x59;
        const double x61 = 2*x13;
        const double x62 = x4*x46;
        const double x63 = x3*x59;
        const double x64 = MP_VALUE*x28*x63;
        const double x65 = x34*x64 - x38*x60 - x44*x61*x62;
        const double x66 = x33*x[3];
        const double x67 = x30*x[3];
        const double x68 = x3*x30;
        const double x69 = x28*x68;
        const double x70 = -x28*x40*x66;
        const double x71 = -x1*x49*x66 + x69*x[3];
        const double x72 = x46*x48;
        const double x73 = 2*x32;
        const double x74 = x43*x73;
        const double x75 = x46*x54;
        const double x76 = -x13*x72 + x13*x74 - x75;
        const double x77 = x39*x4;
        const double x78 = x34*x59;
        const double x79 = -x41*x77 - x48*x60 + x78;
        const double x80 = x13*x78 - x54*x60 - x55*x62 + x68*x77;
		return {std::vector<double>{-x12*x35 + x31*(58.859999999999999*x0*x11 - x1*x2 - x19 + x22 + x27*x28 - x3*(u[0]*x12 - x14*x15 + x2*x9 + x7)) + 8*x32*x43*x45 + x35*x36 - 4*x38*x39*x41, x12*x67 - x36*x67 - 4*x45*x66, MP_VALUE*x69}, {x39}, {}, {}, {}, {}};
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
            {0, 0, -0.001 + M_PI, 0},  // x0
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},  // lb x
            {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},  // ub x
            {-2.5},  // lb u
            {2.5},  // ub u
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
        