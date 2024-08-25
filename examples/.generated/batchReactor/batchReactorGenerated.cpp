
// CODEGEN FOR MODEL "batchReactor"

// includes
#define _USE_MATH_DEFINES
#include "batchReactorGeneratedParams.h"
#include <cmath>
#include <string>
#include "constants.h"
#include <problem.h>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"


// runtime parameters and global constants
const double eff = EFF_VALUE;


// mayer term
class MayerbatchReactor : public Expression {
public:
	static std::unique_ptr<MayerbatchReactor> create() {
		Adjacency adj{{1}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<MayerbatchReactor>(new MayerbatchReactor(std::move(adj), std::move(adjDiff)));
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
	MayerbatchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0batchReactor : public Expression {
public:
	static std::unique_ptr<F0batchReactor> create() {
		Adjacency adj{{0}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 0}}, {{0, 0}}, {}, {}, {}};
		return std::unique_ptr<F0batchReactor>(new F0batchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -x[0]*(u[0]*EFF_VALUE + (1.0/2.0)*pow(u[0], 2));
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-(u[0]*EFF_VALUE + (1.0/2.0)*pow(u[0], 2))}, {-x[0]*(EFF_VALUE + u[0])}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double x0 = -(EFF_VALUE + u[0]);
        const double x1 = -u[0];
        const double x2 = -x[0];
		return {std::vector<double>{}, {x0}, {x2}, {}, {}, {}};
	}
private:
	F0batchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F1batchReactor : public Expression {
public:
	static std::unique_ptr<F1batchReactor> create() {
		Adjacency adj{{0}, {0}, {}};
		AdjacencyDiff adjDiff{{}, {{0, 0}}, {}, {}, {}, {}};
		return std::unique_ptr<F1batchReactor>(new F1batchReactor(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return u[0]*x[0]*EFF_VALUE;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{u[0]*EFF_VALUE}, {x[0]*EFF_VALUE}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {EFF_VALUE}, {}, {}, {}, {}};
	}
private:
	F1batchReactor(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


Problem createProblem_batchReactor() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0batchReactor::create());
    F.push_back(F1batchReactor::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            2, 1, 0,  // #vars
            {1, 0},  // x0
            {0, 0},  // lb x
            {1, 1},  // ub x
            {0},  // lb u
            {5},  // ub u
            {},  // lb p
            {},  // ub p
            MayerbatchReactor::create(),
            {},
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "batchReactor");
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_batchReactor());
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
        