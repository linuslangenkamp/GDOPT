
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
        const double s1 = 0.5*pow(s0, 2) + 1.0;
		return (u[0]*s1 - 0.5*pow(x[3], 2)*s0 - 2.4525000000000001*sin(2*x[2]))/s1;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s16 = sin(x[2]);
        const double s17 = 0.5*pow(s16, 2) + 1.0;
        const double s18 = 1.0/s17;
        const double s19 = 2*x[2];
        const double s20 = cos(x[2]);
        const double s21 = 0.5*pow(x[3], 2);
        const double s22 = 1.0*s16;
		return {std::vector<double>{s18*(1.0*u[0]*s16*s20 - s20*s21 - 4.9050000000000002*cos(s19)) - s20*s22*(u[0]*s17 - s16*s21 - 2.4525000000000001*sin(s19))/pow(s17, 2), -x[3]*s18*s22}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s2 = sin(x[2]);
        const double s3 = pow(s2, 2);
        const double s4 = 0.5*s3 + 1.0;
        const double s5 = 1.0/s4;
        const double s6 = 2*x[2];
        const double s7 = sin(s6);
        const double s8 = 1.0*s3;
        const double s9 = cos(x[2]);
        const double s10 = pow(s9, 2);
        const double s11 = 1.0*s10;
        const double s12 = 0.5*pow(x[3], 2);
        const double s13 = s12*s2;
        const double s14 = 1.0*s2;
        const double s15 = s5*s9;
		return {std::vector<double>{s5*(u[0]*s11 - u[0]*s8 + s13 + 2.0*s15*s2*(-u[0]*s14*s9 + s12*s9 + 4.9050000000000002*cos(s6)) - s5*(-u[0]*s4 + s13 + 2.4525000000000001*s7)*(2.0*s10*s3*s5 - s11 + s8) + 9.8100000000000005*s7), 1.0*x[3]*s15*(s3*s5 - 1), -s14*s5}, {}, {}, {}, {}, {}};
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
		return (-1.125*u[0]*cos(x[2]) + 0.125*u[0]*cos(3*x[2]) + 0.25*pow(x[3], 2)*sin(2*x[2]) + 14.715*s0)/(0.5*pow(s0, 2) + 1.0);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s22 = sin(x[2]);
        const double s23 = 0.5*pow(s22, 2) + 1.0;
        const double s24 = 1.0/s23;
        const double s25 = cos(x[2]);
        const double s26 = 3*x[2];
        const double s27 = pow(x[3], 2);
        const double s28 = 2*x[2];
        const double s29 = 1.125*s25;
        const double s30 = cos(s26);
        const double s31 = sin(s28);
		return {std::vector<double>{-1.0*s22*s25*(-u[0]*s29 + 0.125*u[0]*s30 + 14.715*s22 + 0.25*s27*s31)/pow(s23, 2) + s24*(1.125*u[0]*s22 - 0.375*u[0]*sin(s26) + 14.715*s25 + 0.5*s27*cos(s28)), 0.5*x[3]*s24*s31}, {s24*(-s29 + 0.125*s30)}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s1 = sin(x[2]);
        const double s2 = pow(s1, 2);
        const double s3 = 1.0/(0.5*s2 + 1.0);
        const double s4 = pow(x[3], 2);
        const double s5 = 2*x[2];
        const double s6 = sin(s5);
        const double s7 = s4*s6;
        const double s8 = 3*x[2];
        const double s9 = cos(s8);
        const double s10 = 0.125*s9;
        const double s11 = cos(x[2]);
        const double s12 = 1.125*s11;
        const double s13 = -u[0]*s12 + 14.715*s1;
        const double s14 = pow(s11, 2);
        const double s15 = 2.0*s3;
        const double s16 = 1.125*s1;
        const double s17 = 0.375*sin(s8);
        const double s18 = cos(s5);
        const double s19 = s1*s11;
        const double s20 = s3*(s16 - s17 + 1.0*s19*s3*(-s10 + s12));
        const double s21 = 0.5*s3*s6;
		return {std::vector<double>{s3*(-1.125*u[0]*s9 - s13 - s15*s19*(u[0]*s16 - u[0]*s17 + 14.715*s11 + 0.5*s18*s4) + s3*(u[0]*s10 + s13 + 0.25*s7)*(s14*s15*s2 - 1.0*s14 + 1.0*s2) - 1.0*s7), x[3]*s3*(1.0*s18 - s19*s21), s21}, {s20}, {}, {}, {}, {}};
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
