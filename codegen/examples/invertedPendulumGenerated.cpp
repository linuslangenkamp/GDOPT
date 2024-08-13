
// CODEGEN FOR MODEL "invertedPendulum"

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "invertedPendulumGenerated.h"
#include "constants.h"


class LagrangeinvertedPendulum : public Expression {
public:
	static std::unique_ptr<LagrangeinvertedPendulum> create() {
		Adjacency adj{{2}, {}, {}};
		AdjacencyDiff adjDiff{{{2, 2}}, {}, {}, {}, {}, {}};
		return std::unique_ptr<LagrangeinvertedPendulum>(new LagrangeinvertedPendulum(std::move(adj), std::move(adjDiff)));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return pow(x[2], 2);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{2*x[2]}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{2}, {}, {}, {}, {}, {}};
	}
private:
	LagrangeinvertedPendulum(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


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
        const double s1 = pow(s0, 2) + 100;
		return (u[0]*s1 + pow(x[3], 2)*s0 + 4.9050000000000002*sin(2*x[2]))/s1;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s16 = sin(x[2]);
        const double s17 = pow(s16, 2) + 100;
        const double s18 = 1.0/s17;
        const double s19 = 2*x[2];
        const double s20 = cos(x[2]);
        const double s21 = pow(x[3], 2);
        const double s22 = 2*s16;
        const double s23 = s20*s22;
		return {std::vector<double>{s18*(u[0]*s23 + s20*s21 + 9.8100000000000005*cos(s19)) - s23*(u[0]*s17 + s16*s21 + 4.9050000000000002*sin(s19))/pow(s17, 2), x[3]*s18*s22}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s2 = sin(x[2]);
        const double s3 = pow(s2, 2);
        const double s4 = s3 + 100;
        const double s5 = 1.0/s4;
        const double s6 = 2*x[2];
        const double s7 = sin(s6);
        const double s8 = pow(x[3], 2);
        const double s9 = s2*s8;
        const double s10 = 2*u[0];
        const double s11 = cos(x[2]);
        const double s12 = pow(s11, 2);
        const double s13 = s11*s2;
        const double s14 = 4*s5;
        const double s15 = 2*s5;
		return {std::vector<double>{s5*(2*u[0]*s12 - s10*s3 - s13*s14*(s10*s13 + s11*s8 + 9.8100000000000005*cos(s6)) + 2*s5*(u[0]*s4 + 4.9050000000000002*s7 + s9)*(s12*s14*s3 - s12 + s3) - 19.620000000000001*s7 - s9), x[3]*s11*s15*(-s15*s3 + 1), s15*s2}, {}, {}, {}, {}, {}};
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
		return (-100.25*u[0]*cos(x[2]) + 0.24999999999999997*u[0]*cos(3*x[2]) - 0.49999999999999994*pow(x[3], 2)*sin(2*x[2]) - 990.80999999999995*s0)/(pow(s0, 2) + 100);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s22 = sin(x[2]);
        const double s23 = pow(s22, 2) + 100;
        const double s24 = 1.0/s23;
        const double s25 = cos(x[2]);
        const double s26 = 3*x[2];
        const double s27 = pow(x[3], 2);
        const double s28 = 2*x[2];
        const double s29 = 100.25*s25;
        const double s30 = cos(s26);
        const double s31 = sin(s28);
		return {std::vector<double>{-2*s22*s25*(-u[0]*s29 + 0.24999999999999997*u[0]*s30 - 990.80999999999995*s22 - 0.49999999999999994*s27*s31)/pow(s23, 2) + s24*(100.25*u[0]*s22 - 0.74999999999999989*u[0]*sin(s26) - 990.80999999999995*s25 - 0.99999999999999989*s27*cos(s28)), -0.99999999999999989*x[3]*s24*s31}, {s24*(-s29 + 0.24999999999999997*s30)}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s1 = sin(x[2]);
        const double s2 = pow(s1, 2);
        const double s3 = 1.0/(s2 + 100);
        const double s4 = pow(x[3], 2);
        const double s5 = 2*x[2];
        const double s6 = sin(s5);
        const double s7 = s4*s6;
        const double s8 = 3*x[2];
        const double s9 = cos(s8);
        const double s10 = cos(x[2]);
        const double s11 = pow(s10, 2);
        const double s12 = 4*s3;
        const double s13 = 0.24999999999999997*s9;
        const double s14 = 100.25*s10;
        const double s15 = u[0]*s14 + 990.80999999999995*s1;
        const double s16 = 2*s3;
        const double s17 = 100.25*s1;
        const double s18 = 0.74999999999999989*sin(s8);
        const double s19 = cos(s5);
        const double s20 = s1*s10;
        const double s21 = s3*(s16*s20*(-s13 + s14) + s17 - s18);
		return {std::vector<double>{s3*(-2.2499999999999996*u[0]*s9 + s12*s20*(-u[0]*s17 + u[0]*s18 + 990.80999999999995*s10 + 0.99999999999999989*s19*s4) + s15 - s16*(-u[0]*s13 + s15 + 0.49999999999999994*s7)*(s11*s12*s2 - s11 + s2) + 1.9999999999999998*s7), 1.9999999999999998*x[3]*s3*(s1*s10*s3*s6 - s19), -0.99999999999999989*s3*s6}, {s21}, {}, {}, {}, {}};
	}
private:
	F3invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {}
};


class R0invertedPendulum : public Constraint {
public:
	static std::unique_ptr<R0invertedPendulum> create() {
		Adjacency adj{{0}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<R0invertedPendulum>(new R0invertedPendulum(std::move(adj), std::move(adjDiff), 5, 5));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return x[0];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	R0invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class R1invertedPendulum : public Constraint {
public:
	static std::unique_ptr<R1invertedPendulum> create() {
		Adjacency adj{{1}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<R1invertedPendulum>(new R1invertedPendulum(std::move(adj), std::move(adjDiff), 0, 0));
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
	R1invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class R2invertedPendulum : public Constraint {
public:
	static std::unique_ptr<R2invertedPendulum> create() {
		Adjacency adj{{2}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<R2invertedPendulum>(new R2invertedPendulum(std::move(adj), std::move(adjDiff), 0, 0));
	}

	double eval(const double *x, const double *u, const double *p, double t) override {
		return x[2];
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{1}, {}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}
private:
	R2invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


class R3invertedPendulum : public Constraint {
public:
	static std::unique_ptr<R3invertedPendulum> create() {
		Adjacency adj{{3}, {}, {}};
		AdjacencyDiff adjDiff{{}, {}, {}, {}, {}, {}};
		return std::unique_ptr<R3invertedPendulum>(new R3invertedPendulum(std::move(adj), std::move(adjDiff), 0, 0));
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
	R3invertedPendulum(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {}
};


Problem createProblem_invertedPendulum() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0invertedPendulum::create());
    F.push_back(F1invertedPendulum::create());
    F.push_back(F2invertedPendulum::create());
    F.push_back(F3invertedPendulum::create());
    
    std::vector<std::unique_ptr<Constraint>> G;
    
    
    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0invertedPendulum::create());
    R.push_back(R1invertedPendulum::create());
    R.push_back(R2invertedPendulum::create());
    R.push_back(R3invertedPendulum::create());
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    

    Problem problem(
            4, 1, 0,  // #vars
            {0, 0, 1, 0},  // x0
            {0, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},  // lb x
            {5, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},  // ub x
            {-0.5},  // lb u
            {0.5},  // ub u
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
