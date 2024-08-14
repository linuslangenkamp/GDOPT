
// CODEGEN FOR MODEL "simpleParameter"

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "simpleParameterGenerated.h"
#include "constants.h"



class MayersimpleParameter : public Expression {
public:
	static std::unique_ptr<MayersimpleParameter> create() {
		Adjacency adj{{}, {}, {0, 1}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {{1, 0}}};
		return std::unique_ptr<MayersimpleParameter>(new MayersimpleParameter(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return p[0]*p[1];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {p[1], p[0]}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {1}};
	}
private:
	MayersimpleParameter(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F0simpleParameter : public Expression {
public:
	static std::unique_ptr<F0simpleParameter> create() {
		Adjacency adj{{}, {}, {0}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F0simpleParameter>(new F0simpleParameter(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 2*p[0]*t;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {2*t}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	F0simpleParameter(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class G0simpleParameter : public Constraint {
public:
	static std::unique_ptr<G0simpleParameter> create() {
		Adjacency adj{{0}, {0}, {1}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<G0simpleParameter>(new G0simpleParameter(std::move(adj), std::move(adjDiff), 0.2, 0.25));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return u[0] + p[1] - x[0];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-1}, {1}, {1}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	G0simpleParameter(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


Problem createProblem_simpleParameter() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0simpleParameter::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    G.push_back(G0simpleParameter::create());
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            1, 1, 2,  // #vars
            {0},  // x0
            {MINUS_INFINITY},  // lb x
            {PLUS_INFINITY},  // ub x
            {-1},  // lb u
            {1},  // ub u
            {MINUS_INFINITY, MINUS_INFINITY},  // lb p
            {PLUS_INFINITY, PLUS_INFINITY},  // ub p
            MayersimpleParameter::create(),
            {},
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "simpleParameter");
    return problem;
};
