
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


// runtime parameters
const double Parameter_Ms = PARAMETER_MS_VALUE;
const double Parameter_Mp = PARAMETER_MP_VALUE;
const double Parameter_R = PARAMETER_R_VALUE;


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
        const double s1 = (1.0/2.0)*x[2];
		return {std::vector<double>{sin(s1)*cos(s1)}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s0 = (1.0/2.0)*x[2];
		return {std::vector<double>{-1.0/2.0*pow(sin(s0), 2) + (1.0/2.0)*pow(cos(s0), 2)}, {}, {}, {}, {}, {}};
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
        const double s0 = sin(x[2]);
        const double s1 = Parameter_Mp*pow(s0, 2) + Parameter_Ms;
		return (u[0]*s1 - Parameter_Mp*Parameter_R*pow(x[3], 2)*s0 - 4.9050000000000002*Parameter_Mp*sin(2*x[2]))/s1;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s18 = sin(x[2]);
        const double s19 = Parameter_Mp*pow(s18, 2) + Parameter_Ms;
        const double s20 = 1.0/s19;
        const double s21 = 2*x[2];
        const double s22 = cos(x[2]);
        const double s23 = Parameter_Mp*s22;
        const double s24 = Parameter_R*pow(x[3], 2);
        const double s25 = Parameter_Mp*s18;
		return {std::vector<double>{-2*s18*s23*(u[0]*s19 - 4.9050000000000002*Parameter_Mp*sin(s21) - s24*s25)/pow(s19, 2) + s20*(2*u[0]*Parameter_Mp*s18*s22 - 9.8100000000000005*Parameter_Mp*cos(s21) - s23*s24), -2*Parameter_R*x[3]*s20*s25}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s2 = 2*x[2];
        const double s3 = sin(s2);
        const double s4 = sin(x[2]);
        const double s5 = pow(s4, 2);
        const double s6 = 2*u[0];
        const double s7 = cos(x[2]);
        const double s8 = pow(s7, 2);
        const double s9 = Parameter_R*pow(x[3], 2);
        const double s10 = s4*s9;
        const double s11 = s4*s7;
        const double s12 = Parameter_Mp*s5;
        const double s13 = Parameter_Ms + s12;
        const double s14 = 1.0/s13;
        const double s15 = Parameter_Mp*s14;
        const double s16 = 2*s14;
        const double s17 = 2*Parameter_R*s15;
		return {std::vector<double>{s15*(s10 + 4*s11*s15*(-s11*s6 + s7*s9 + 9.8100000000000005*cos(s2)) - s16*(-u[0]*s13 + Parameter_Mp*s10 + 4.9050000000000002*Parameter_Mp*s3)*(4*s12*s14*s8 + s5 - s8) + 19.620000000000001*s3 - s5*s6 + s6*s8), x[3]*s17*s7*(s12*s16 - 1), -s17*s4}, {}, {}, {}, {}, {}};
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
        const double s0 = sin(x[2]);
        const double s1 = Parameter_Mp*pow(s0, 2) + Parameter_Ms;
		return (9.8100000000000005*s0*s1 - (u[0]*s1 - Parameter_Mp*Parameter_R*pow(x[3], 2)*s0 - 4.9050000000000002*Parameter_Mp*sin(2*x[2]))*cos(x[2]))/(Parameter_R*s1);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s26 = sin(x[2]);
        const double s27 = Parameter_Mp*pow(s26, 2);
        const double s28 = Parameter_Ms + s27;
        const double s29 = 9.8100000000000005*s28;
        const double s30 = cos(x[2]);
        const double s31 = 2*x[2];
        const double s32 = Parameter_Mp*s26;
        const double s33 = Parameter_R*pow(x[3], 2);
        const double s34 = u[0]*s28 - 4.9050000000000002*Parameter_Mp*sin(s31) - s32*s33;
        const double s35 = 1.0/Parameter_R;
        const double s36 = s30*s35;
        const double s37 = 2*s32;
        const double s38 = 1.0/s28;
		return {std::vector<double>{s35*s38*(s26*s34 + 19.620000000000001*s27*s30 + s29*s30 - s30*(2*u[0]*Parameter_Mp*s26*s30 - Parameter_Mp*s30*s33 - 9.8100000000000005*Parameter_Mp*cos(s31))) - s36*s37*(s26*s29 - s30*s34)/pow(s28, 2), x[3]*s30*s37*s38}, {-s36}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s2 = 1.0/Parameter_R;
        const double s3 = sin(x[2]);
        const double s4 = pow(s3, 2);
        const double s5 = Parameter_Mp*s4;
        const double s6 = Parameter_Ms + s5;
        const double s7 = 1.0/s6;
        const double s8 = 2*x[2];
        const double s9 = sin(s8);
        const double s10 = 2*u[0];
        const double s11 = cos(x[2]);
        const double s12 = pow(s11, 2);
        const double s13 = Parameter_R*pow(x[3], 2);
        const double s14 = s13*s3;
        const double s15 = Parameter_Mp*s11;
        const double s16 = -s10*s11*s3 + s11*s13 + 9.8100000000000005*cos(s8);
        const double s17 = 4*s7;
        const double s18 = s12*s5;
        const double s19 = -s12 + s4;
        const double s20 = 9.8100000000000005*s6;
        const double s21 = -u[0]*s6 + Parameter_Mp*s14 + 4.9050000000000002*Parameter_Mp*s9;
        const double s22 = s11*s21 + s20*s3;
        const double s23 = s15*s3;
        const double s24 = s2*s3;
        const double s25 = 2*s7;
		return {std::vector<double>{s2*s7*(58.859999999999999*Parameter_Mp*s12*s3 - 2*Parameter_Mp*s16*s3 + 2*Parameter_Mp*s22*s7*(s17*s18 + s19) - 19.620000000000001*Parameter_Mp*pow(s3, 3) - s15*(s10*s12 - s10*s4 + s14 + 19.620000000000001*s9) - s17*s23*(s11*s20 + 19.620000000000001*s11*s5 + s15*s16 - s21*s3) - s22), -Parameter_Mp*x[3]*s25*(s18*s25 + s19), s23*s25}, {s24}, {}, {}, {}, {}};
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
            {0, 0, 3.1405926535897932, 0},  // x0
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
        