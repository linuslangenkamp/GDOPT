
// CODEGEN FOR MODEL "satellite"

#include <cmath>
#include <string>
#include "satellite.h"
#include "constants.h"


const double I1 = 1000000;
const double I2 = 833333;
const double I3 = 916667;
const double T1S = 550;
const double T2S = 50;
const double T3S = 550;


class Mayersatellite : public Expression {
public:
	static std::unique_ptr<Mayersatellite> create() {
		Adjacency adj{{0, 1, 2, 3, 4, 5, 6}, {}, {}};
		AdjacencyDiff adjDiff{{{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<Mayersatellite>(new Mayersatellite(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return pow(x[4], 2) + pow(x[5], 2) + pow(x[6], 2) + pow(x[0] - 0.70106000000000002, 2) + pow(x[1] - 0.092299999999999993, 2) + pow(x[2] - 0.56098000000000003, 2) + pow(x[3] - 0.43047000000000002, 2);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{2*x[0] - 1.40212, 2*x[1] - 0.18459999999999999, 2*x[2] - 1.1219600000000001, 2*x[3] - 0.86094000000000004, 2*x[4], 2*x[5], 2*x[6]}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{2, 2, 2, 2, 2, 2, 2}, {}, {}, {}, {}, {}};
	}
private:
	Mayersatellite(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class Lagrangesatellite : public Expression {
public:
	static std::unique_ptr<Lagrangesatellite> create() {
		Adjacency adj{{}, {0, 1, 2}, {}};
		AdjacencyDiff adjDiff{{}, {}, {{0, 0}, {1, 1}, {2, 2}}, {}, {}, {}};
		return std::unique_ptr<Lagrangesatellite>(new Lagrangesatellite(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.5*pow(u[0], 2) + 0.5*pow(u[1], 2) + 0.5*pow(u[2], 2);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {1.0*u[0], 1.0*u[1], 1.0*u[2]}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {1.0, 1.0, 1.0}, {}, {}, {}};
	}
private:
	Lagrangesatellite(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F0satellite : public Expression {
public:
	static std::unique_ptr<F0satellite> create() {
		Adjacency adj{{1, 2, 3, 4, 5, 6}, {}, {}};
		AdjacencyDiff adjDiff{{{4, 3}, {5, 2}, {6, 1}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F0satellite>(new F0satellite(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.5*x[1]*x[6] - 0.5*x[2]*x[5] + 0.5*x[3]*x[4];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{0.5*x[6], -0.5*x[5], 0.5*x[4], 0.5*x[3], -0.5*x[2], 0.5*x[1]}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{0.5, -0.5, 0.5}, {}, {}, {}, {}, {}};
	}
private:
	F0satellite(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F1satellite : public Expression {
public:
	static std::unique_ptr<F1satellite> create() {
		Adjacency adj{{0, 2, 3, 4, 5, 6}, {}, {}};
		AdjacencyDiff adjDiff{{{4, 2}, {5, 3}, {6, 0}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F1satellite>(new F1satellite(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -0.5*x[0]*x[6] + 0.5*x[2]*x[4] + 0.5*x[3]*x[5];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.5*x[6], 0.5*x[4], 0.5*x[5], 0.5*x[2], 0.5*x[3], -0.5*x[0]}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{0.5, 0.5, -0.5}, {}, {}, {}, {}, {}};
	}
private:
	F1satellite(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F2satellite : public Expression {
public:
	static std::unique_ptr<F2satellite> create() {
		Adjacency adj{{0, 1, 3, 4, 5, 6}, {}, {}};
		AdjacencyDiff adjDiff{{{4, 1}, {5, 0}, {6, 3}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F2satellite>(new F2satellite(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return 0.5*x[0]*x[5] - 0.5*x[1]*x[4] + 0.5*x[3]*x[6];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{0.5*x[5], -0.5*x[4], 0.5*x[6], -0.5*x[1], 0.5*x[0], 0.5*x[3]}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.5, 0.5, 0.5}, {}, {}, {}, {}, {}};
	}
private:
	F2satellite(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F3satellite : public Expression {
public:
	static std::unique_ptr<F3satellite> create() {
		Adjacency adj{{0, 1, 2, 4, 5, 6}, {}, {}};
		AdjacencyDiff adjDiff{{{4, 0}, {5, 1}, {6, 2}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F3satellite>(new F3satellite(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return -0.5*x[0]*x[4] - 0.5*x[1]*x[5] - 0.5*x[2]*x[6];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.5*x[4], -0.5*x[5], -0.5*x[6], -0.5*x[0], -0.5*x[1], -0.5*x[2]}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-0.5, -0.5, -0.5}, {}, {}, {}, {}, {}};
	}
private:
	F3satellite(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F4satellite : public Expression {
public:
	static std::unique_ptr<F4satellite> create() {
		Adjacency adj{{5, 6}, {0}, {}};
		AdjacencyDiff adjDiff{{{6, 5}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F4satellite>(new F4satellite(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return (T1S*u[0] + x[5]*x[6]*(I2 - I3))/I1;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{x[6]*(I2 - I3)/I1, x[5]*(I2 - I3)/I1}, {T1S/I1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{(I2 - I3)/I1}, {}, {}, {}, {}, {}};
	}
private:
	F4satellite(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F5satellite : public Expression {
public:
	static std::unique_ptr<F5satellite> create() {
		Adjacency adj{{4, 6}, {1}, {}};
		AdjacencyDiff adjDiff{{{6, 4}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F5satellite>(new F5satellite(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return (T2S*u[1] + x[4]*x[6]*(-I1 + I3))/I2;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{x[6]*(-I1 + I3)/I2, x[4]*(-I1 + I3)/I2}, {T2S/I2}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{-(I1 - I3)/I2}, {}, {}, {}, {}, {}};
	}
private:
	F5satellite(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class F6satellite : public Expression {
public:
	static std::unique_ptr<F6satellite> create() {
		Adjacency adj{{4, 5}, {2}, {}};
		AdjacencyDiff adjDiff{{{5, 4}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<F6satellite>(new F6satellite(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return (T3S*u[2] + x[4]*x[5]*(I1 - I2))/I3;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{x[5]*(I1 - I2)/I3, x[4]*(I1 - I2)/I3}, {T3S/I3}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{(I1 - I2)/I3}, {}, {}, {}, {}, {}};
	}
private:
	F6satellite(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


Problem createProblem_satellite() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0satellite::create());
    F.push_back(F1satellite::create());
    F.push_back(F2satellite::create());
    F.push_back(F3satellite::create());
    F.push_back(F4satellite::create());
    F.push_back(F5satellite::create());
    F.push_back(F6satellite::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            7, 3, 0,  // #vars
            {0, 0, 0, 1, 0.01, 0.005, 0.001},  // x0
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},  // lb x
            {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},  // ub x
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},  // lb u
            {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},  // ub u
            {},  // lb p
            {},  // ub p
            Mayersatellite::create(),
            Lagrangesatellite::create(),
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "satellite");
    return problem;
};
