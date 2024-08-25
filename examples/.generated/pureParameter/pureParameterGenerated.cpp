
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
		return -(1.49859411224726*p[0] + 1.60048027388215*p[10] + 1.44191387106179*p[11] + 1.74363772891805*p[12] + 1.60599325211726*p[13] + 0.239668663322606*p[14] + 1.13615952323541*p[1] + 0.116810012013373*p[2] + 0.899662676217498*p[3] + 0.908114831692897*p[4] + 1.37417836465884*p[5] + 1.20454058585417*p[6] + 0.555418042584864*p[7] + 1.90677403041815*p[8] + 1.9791205060385*p[9]);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {-1.49859411224726, -1.13615952323541, -0.116810012013373, -0.899662676217498, -0.908114831692897, -1.37417836465884, -1.20454058585417, -0.555418042584864, -1.90677403041815, -1.9791205060385, -1.60048027388215, -1.44191387106179, -1.74363772891805, -1.60599325211726, -0.239668663322606}};
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
		return p[0]*(-1 + p[0]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[0]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[1]*(-1 + p[1]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[1]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[2]*(-1 + p[2]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[2]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[3]*(-1 + p[3]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[3]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[4]*(-1 + p[4]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[4]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[5]*(-1 + p[5]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[5]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[6]*(-1 + p[6]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[6]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[7]*(-1 + p[7]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[7]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[8]*(-1 + p[8]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[8]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[9]*(-1 + p[9]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[9]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[10]*(-1 + p[10]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[10]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[11]*(-1 + p[11]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[11]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[12]*(-1 + p[12]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[12]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[13]*(-1 + p[13]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[13]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return p[14]*(-1 + p[14]);
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{-1 + 2*p[14]};
	}

	std::vector<double> evalDiff2(const double* p) override {
		return std::vector<double>{2};
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
		return 0.727693606498591*p[0] + 0.224646561664913*p[10] + 1.09107248153227*p[11] + 0.0023843031091233*p[12] + 0.999940542009322*p[13] + 0.191483926010037*p[14] + 0.0843435964646715*p[1] + 0.321462017870338*p[2] + 1.20118337865561*p[3] + 0.11695612101188*p[4] + 0.0286320009805114*p[5] + 1.36077417000029*p[6] + 1.1800241898134*p[7] + 1.58524236389105*p[8] + 0.18995163968575*p[9];
	}

	std::vector<double> evalDiff(const double* p) override {
		return std::vector<double>{0.727693606498591, 0.0843435964646715, 0.321462017870338, 1.20118337865561, 0.11695612101188, 0.0286320009805114, 1.36077417000029, 1.1800241898134, 1.58524236389105, 0.18995163968575, 0.224646561664913, 1.09107248153227, 0.0023843031091233, 0.999940542009322, 0.191483926010037};
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
        