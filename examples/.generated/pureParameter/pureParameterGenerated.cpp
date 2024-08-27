
// CODEGEN FOR MODEL "pureParameter"

// includes
#define _USE_MATH_DEFINES
#include "pureParameterGeneratedParams.h"
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
class MayerpureParameter : public Expression {
public:
	static std::unique_ptr<MayerpureParameter> create() {
		Adjacency adj{{}, {}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<MayerpureParameter>(new MayerpureParameter(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -(0.152873968104274*p[0] + 1.7168252587759*p[10] + 1.88173719133442*p[11] + 1.63174344016537*p[12] + 0.526654031353459*p[13] + 0.252286564167189*p[14] + 0.444153979792149*p[1] + 0.964924808424995*p[2] + 0.648268388616711*p[3] + 0.285426759974794*p[4] + 1.96804465856293*p[5] + 1.01193586952027*p[6] + 0.144581537101869*p[7] + 1.50505417146067*p[8] + 0.0238087490909578*p[9]);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {-0.152873968104274, -0.444153979792149, -0.964924808424995, -0.648268388616711, -0.285426759974794, -1.96804465856293, -1.01193586952027, -0.144581537101869, -1.50505417146067, -0.0238087490909578, -1.7168252587759, -1.88173719133442, -1.63174344016537, -0.526654031353459, -0.252286564167189}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	MayerpureParameter(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// dynamic constraints
class F0pureParameter : public Expression {
public:
	static std::unique_ptr<F0pureParameter> create() {
		Adjacency adj{{}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F0pureParameter>(new F0pureParameter(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	F0pureParameter(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


// parametric constraints
class A0pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A0pureParameter> create() {
		ParamAdjacency adj{{0}};
		ParamAdjacencyDiff adjDiff{{{0, 0}}};
		return std::unique_ptr<A0pureParameter>(new A0pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[0] + pow(p[0], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[0])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A0pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A1pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A1pureParameter> create() {
		ParamAdjacency adj{{1}};
		ParamAdjacencyDiff adjDiff{{{1, 1}}};
		return std::unique_ptr<A1pureParameter>(new A1pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[1] + pow(p[1], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[1])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A1pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A2pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A2pureParameter> create() {
		ParamAdjacency adj{{2}};
		ParamAdjacencyDiff adjDiff{{{2, 2}}};
		return std::unique_ptr<A2pureParameter>(new A2pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[2] + pow(p[2], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[2])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A2pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A3pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A3pureParameter> create() {
		ParamAdjacency adj{{3}};
		ParamAdjacencyDiff adjDiff{{{3, 3}}};
		return std::unique_ptr<A3pureParameter>(new A3pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[3] + pow(p[3], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[3])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A3pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A4pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A4pureParameter> create() {
		ParamAdjacency adj{{4}};
		ParamAdjacencyDiff adjDiff{{{4, 4}}};
		return std::unique_ptr<A4pureParameter>(new A4pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[4] + pow(p[4], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[4])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A4pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A5pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A5pureParameter> create() {
		ParamAdjacency adj{{5}};
		ParamAdjacencyDiff adjDiff{{{5, 5}}};
		return std::unique_ptr<A5pureParameter>(new A5pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[5] + pow(p[5], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[5])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A5pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A6pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A6pureParameter> create() {
		ParamAdjacency adj{{6}};
		ParamAdjacencyDiff adjDiff{{{6, 6}}};
		return std::unique_ptr<A6pureParameter>(new A6pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[6] + pow(p[6], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[6])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A6pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A7pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A7pureParameter> create() {
		ParamAdjacency adj{{7}};
		ParamAdjacencyDiff adjDiff{{{7, 7}}};
		return std::unique_ptr<A7pureParameter>(new A7pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[7] + pow(p[7], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[7])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A7pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A8pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A8pureParameter> create() {
		ParamAdjacency adj{{8}};
		ParamAdjacencyDiff adjDiff{{{8, 8}}};
		return std::unique_ptr<A8pureParameter>(new A8pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[8] + pow(p[8], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[8])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A8pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A9pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A9pureParameter> create() {
		ParamAdjacency adj{{9}};
		ParamAdjacencyDiff adjDiff{{{9, 9}}};
		return std::unique_ptr<A9pureParameter>(new A9pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[9] + pow(p[9], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[9])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A9pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A10pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A10pureParameter> create() {
		ParamAdjacency adj{{10}};
		ParamAdjacencyDiff adjDiff{{{10, 10}}};
		return std::unique_ptr<A10pureParameter>(new A10pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[10] + pow(p[10], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[10])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A10pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A11pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A11pureParameter> create() {
		ParamAdjacency adj{{11}};
		ParamAdjacencyDiff adjDiff{{{11, 11}}};
		return std::unique_ptr<A11pureParameter>(new A11pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[11] + pow(p[11], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[11])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A11pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A12pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A12pureParameter> create() {
		ParamAdjacency adj{{12}};
		ParamAdjacencyDiff adjDiff{{{12, 12}}};
		return std::unique_ptr<A12pureParameter>(new A12pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[12] + pow(p[12], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[12])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A12pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A13pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A13pureParameter> create() {
		ParamAdjacency adj{{13}};
		ParamAdjacencyDiff adjDiff{{{13, 13}}};
		return std::unique_ptr<A13pureParameter>(new A13pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[13] + pow(p[13], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[13])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A13pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A14pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A14pureParameter> create() {
		ParamAdjacency adj{{14}};
		ParamAdjacencyDiff adjDiff{{{14, 14}}};
		return std::unique_ptr<A14pureParameter>(new A14pureParameter(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double* p) override {
		return -(-p[14] + pow(p[14], 2));
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-(-1 + 2*p[14])};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{-2};
	}
private:
	A14pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class A15pureParameter : public ParamConstraint {
public:
	static std::unique_ptr<A15pureParameter> create() {
		ParamAdjacency adj{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}};
		ParamAdjacencyDiff adjDiff{{}};
		return std::unique_ptr<A15pureParameter>(new A15pureParameter(std::move(adj), std::move(adjDiff), MINUS_INFINITY, 5));
	}

	double eval(const double* p) override {
		return 0.989614145245674*p[0] + 1.83381414008353*p[10] + 1.9969247970382*p[11] + 0.272265247538793*p[12] + 0.832453337096074*p[13] + 1.54098370456447*p[14] + 1.79396354140524*p[1] + 0.657159420149244*p[2] + 0.440374400300925*p[3] + 1.8349598556623*p[4] + 0.231431022250533*p[5] + 0.0386835389577318*p[6] + 1.79708257968256*p[7] + 1.00215126546455*p[8] + 1.25414419413018*p[9];
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{0.989614145245674, 1.79396354140524, 0.657159420149244, 0.440374400300925, 1.8349598556623, 0.231431022250533, 0.0386835389577318, 1.79708257968256, 1.00215126546455, 1.25414419413018, 1.83381414008353, 1.9969247970382, 0.272265247538793, 0.832453337096074, 1.54098370456447};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{};
	}
private:
	A15pureParameter(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


Problem createProblem_pureParameter() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0pureParameter::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    A.push_back(A0pureParameter::create());
    A.push_back(A1pureParameter::create());
    A.push_back(A2pureParameter::create());
    A.push_back(A3pureParameter::create());
    A.push_back(A4pureParameter::create());
    A.push_back(A5pureParameter::create());
    A.push_back(A6pureParameter::create());
    A.push_back(A7pureParameter::create());
    A.push_back(A8pureParameter::create());
    A.push_back(A9pureParameter::create());
    A.push_back(A10pureParameter::create());
    A.push_back(A11pureParameter::create());
    A.push_back(A12pureParameter::create());
    A.push_back(A13pureParameter::create());
    A.push_back(A14pureParameter::create());
    A.push_back(A15pureParameter::create());

    Problem problem(
            1, 0, 15,  // #vars
            {0},  // x0
            {MINUS_INFINITY},  // lb x
            {PLUS_INFINITY},  // ub x
            {},  // u0 initial guesses for optimization
            {},  // lb u
            {},  // ub u
            {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5},  // p0 initial guesses for optimization
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  // lb p
            {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},  // ub p
            MayerpureParameter::create(),
            {},
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "pureParameter");
    return problem;
};

int main() {
    auto problem = std::make_shared<const Problem>(createProblem_pureParameter());
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
        