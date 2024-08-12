
// CODEGEN FOR MODEL "pureParameter"

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "pureParameterGenerated.h"
#include "constants.h"


class MayerpureParameter : public Expression {
public:
	static std::unique_ptr<MayerpureParameter> create() {
		Adjacency adj{{}, {}, {0, 1}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<MayerpureParameter>(new MayerpureParameter(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -3*p[0] - 2*p[1];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {-3, -2}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	MayerpureParameter(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


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


class G0pureParameter : public Constraint {
public:
	static std::unique_ptr<G0pureParameter> create() {
		Adjacency adj{{}, {}, {0, 1}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {{0, 0}, {1, 1}}};
		return std::unique_ptr<G0pureParameter>(new G0pureParameter(std::move(adj), std::move(adjDiff), 1, 1));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return pow(p[0], 2) + pow(p[1], 2);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {2*p[0], 2*p[1]}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {2, 2}};
	}
private:
	G0pureParameter(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


Problem createProblem_pureParameter() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0pureParameter::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    G.push_back(G0pureParameter::create());
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            1, 0, 2,  // #vars
            {0},  // x0
            {MINUS_INFINITY},  // lb x
            {PLUS_INFINITY},  // ub x
            {},  // lb u
            {},  // ub u
            {-1, -1},  // lb p
            {1, 1},  // ub p
            MayerpureParameter::create(),
            {},
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "pureParameter");
    return problem;
};
