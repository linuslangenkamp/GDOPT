
// CODEGEN FOR MODEL "oilShalePyrolysis"

// includes
#define _USE_MATH_DEFINES
#include "oilShalePyrolysisGeneratedParams.h"
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
class MayeroilShalePyrolysis : public Expression {
public:
	static std::unique_ptr<MayeroilShalePyrolysis> create() {
		Adjacency adj{{1}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<MayeroilShalePyrolysis>(new MayeroilShalePyrolysis(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -x[1];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	MayeroilShalePyrolysis(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0oilShalePyrolysis : public Expression {
public:
	static std::unique_ptr<F0oilShalePyrolysis> create() {
		Adjacency adj{{0, 1}, {0}, {}};
		AdjacencyDiff adjDiff{{{1, 0}}, {{0, 0}, {0, 1}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F0oilShalePyrolysis>(new F0oilShalePyrolysis(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], -1);
		return -exp(8.86 - 10215.3784219002*x0)*x[0] - x[0]*x[1]*(exp(18.75 - 14190.8212560386*x0) + exp(20.7 - 15599.8389694042*x0) + exp(23.67 - 17008.8566827697*x0));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], -1);
        const double x1 = exp(8.86 - 10215.3784219002*x0);
        const double x2 = exp(20.7 - 15599.8389694042*x0);
        const double x3 = exp(18.75 - 14190.8212560386*x0);
        const double x4 = exp(23.67 - 17008.8566827697*x0);
        const double x5 = x2 + x3 + x4;
        const double x6 = pow(u[0], -2);
		return {std::vector<double>{-x1 - x5*x[1], -x5*x[0]}, {-10215.3784219002*x1*x6*x[0] - x[0]*x[1]*(15599.8389694042*x2*x6 + 14190.8212560386*x3*x6 + 17008.8566827697*x4*x6)}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], -1);
        const double x1 = exp(23.67 - 17008.8566827697*x0);
        const double x2 = pow(u[0], -2);
        const double x3 = exp(18.75 - 14190.8212560386*x0);
        const double x4 = exp(20.7 - 15599.8389694042*x0);
        const double x5 = 17008.8566827697*x2*x1 + 14190.8212560386*x2*x3 + 15599.8389694042*x2*x4;
        const double x6 = exp(8.86 - 10215.3784219002*x0);
        const double x7 = -10215.3784219002*x2*x6 - x5*x[1];
        const double x8 = -x5*x[0];
        const double x9 = pow(u[0], -4);
        const double x10 = pow(u[0], -3);
        const double x11 = x6*x[0];
		return {std::vector<double>{-(x1 + x3 + x4)}, {x7, x8}, {20430.7568438003*x11*x10 - 104353956.302623*x9*x11 - x[0]*x[1]*(-34017.7133655395*x1*x10 + 289301205.655*x1*x9 - 28381.6425120773*x3*x10 + 201379407.920838*x3*x9 - 31199.6779388084*x4*x10 + 243354975.871341*x4*x9)}, {}, {}, {}};
	}
private:
	F0oilShalePyrolysis(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F1oilShalePyrolysis : public Expression {
public:
	static std::unique_ptr<F1oilShalePyrolysis> create() {
		Adjacency adj{{0, 1}, {0}, {}};
		AdjacencyDiff adjDiff{{{1, 0}}, {{0, 0}, {0, 1}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F1oilShalePyrolysis>(new F1oilShalePyrolysis(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], -1);
		return exp(8.86 - 10215.3784219002*x0)*x[0] - exp(24.25 - 18820.4508856683*x0)*x[1] + exp(23.67 - 17008.8566827697*x0)*x[0]*x[1];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], -1);
        const double x1 = exp(8.86 - 10215.3784219002*x0);
        const double x2 = exp(23.67 - 17008.8566827697*x0);
        const double x3 = x2*x[1];
        const double x4 = exp(24.25 - 18820.4508856683*x0);
        const double x5 = pow(u[0], -2);
        const double x6 = x5*x[0];
		return {std::vector<double>{x1 + x3, -x4 + x2*x[0]}, {10215.3784219002*x1*x6 + 17008.8566827697*x3*x6 - 18820.4508856683*x4*x5*x[1]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], -1);
        const double x1 = exp(23.67 - 17008.8566827697*x0);
        const double x2 = pow(u[0], -2);
        const double x3 = 17008.8566827697*x2;
        const double x4 = exp(8.86 - 10215.3784219002*x0);
        const double x5 = 10215.3784219002*x2*x4 + x1*x3*x[1];
        const double x6 = x1*x[0];
        const double x7 = exp(24.25 - 18820.4508856683*x0);
        const double x8 = -18820.4508856683*x2*x7 + x3*x6;
        const double x9 = pow(u[0], -3);
        const double x10 = x6*x[1];
        const double x11 = pow(u[0], -4);
        const double x12 = x7*x[1];
        const double x13 = x4*x[0];
		return {std::vector<double>{x1}, {x5, x8}, {289301205.655*x11*x10 - 354209371.539852*x12*x11 + 104353956.302623*x13*x11 - 34017.7133655395*x9*x10 + 37640.9017713366*x9*x12 - 20430.7568438003*x9*x13}, {}, {}, {}};
	}
private:
	F1oilShalePyrolysis(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F2oilShalePyrolysis : public Expression {
public:
	static std::unique_ptr<F2oilShalePyrolysis> create() {
		Adjacency adj{{0, 1}, {0}, {}};
		AdjacencyDiff adjDiff{{{1, 0}}, {{0, 0}, {0, 1}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F2oilShalePyrolysis>(new F2oilShalePyrolysis(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], -1);
		return exp(24.25 - 18820.4508856683*x0)*x[1] + exp(18.75 - 14190.8212560386*x0)*x[0]*x[1];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], -1);
        const double x1 = exp(18.75 - 14190.8212560386*x0);
        const double x2 = x1*x[1];
        const double x3 = exp(24.25 - 18820.4508856683*x0);
        const double x4 = pow(u[0], -2);
		return {std::vector<double>{x2, x3 + x1*x[0]}, {14190.8212560386*x2*x4*x[0] + 18820.4508856683*x4*x3*x[1]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = pow(u[0], -1);
        const double x1 = exp(18.75 - 14190.8212560386*x0);
        const double x2 = pow(u[0], -2);
        const double x3 = 14190.8212560386*x2;
        const double x4 = x1*x3*x[1];
        const double x5 = x1*x[0];
        const double x6 = exp(24.25 - 18820.4508856683*x0);
        const double x7 = 18820.4508856683*x2*x6 + x3*x5;
        const double x8 = pow(u[0], -4);
        const double x9 = x5*x[1];
        const double x10 = pow(u[0], -3);
        const double x11 = x6*x[1];
		return {std::vector<double>{x1}, {x4, x7}, {-37640.9017713366*x11*x10 + 354209371.539852*x8*x11 + 201379407.920838*x8*x9 - 28381.6425120773*x9*x10}, {}, {}, {}};
	}
private:
	F2oilShalePyrolysis(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F3oilShalePyrolysis : public Expression {
public:
	static std::unique_ptr<F3oilShalePyrolysis> create() {
		Adjacency adj{{0, 1}, {0}, {}};
		AdjacencyDiff adjDiff{{{1, 0}}, {{0, 0}, {0, 1}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F3oilShalePyrolysis>(new F3oilShalePyrolysis(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return exp(20.7 - 15599.8389694042*pow(u[0], -1))*x[0]*x[1];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double x0 = exp(20.7 - 15599.8389694042*pow(u[0], -1));
        const double x1 = x0*x[0];
		return {std::vector<double>{x0*x[1], x1}, {15599.8389694042*x1*x[1]/pow(u[0], 2)}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = exp(20.7 - 15599.8389694042*pow(u[0], -1));
        const double x1 = 15599.8389694042/pow(u[0], 2);
        const double x2 = x0*x1*x[1];
        const double x3 = x0*x[0];
        const double x4 = x1*x3;
        const double x5 = x3*x[1];
		return {std::vector<double>{x0}, {x2, x4}, {243354975.871341*x5/pow(u[0], 4) - 31199.6779388084*x5/pow(u[0], 3)}, {}, {}, {}};
	}
private:
	F3oilShalePyrolysis(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


std::vector<double> uInitialGuess(double t) {
	 return {700};
};

Problem createProblem_oilShalePyrolysis() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0oilShalePyrolysis::create());
    F.push_back(F1oilShalePyrolysis::create());
    F.push_back(F2oilShalePyrolysis::create());
    F.push_back(F3oilShalePyrolysis::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            4, 1, 0,  // #vars
            {1, 0, 0, 0},  // x0
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},  // lb x
            {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},  // ub x
            &uInitialGuess,  // u0 initial guesses for optimization
            {698.15},  // lb u
            {748.15},  // ub u
            {},  // p0 initial guesses for optimization
            {},  // lb p
            {},  // ub p
            MayeroilShalePyrolysis::create(),
            {},
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "oilShalePyrolysis");
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_oilShalePyrolysis());
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
        