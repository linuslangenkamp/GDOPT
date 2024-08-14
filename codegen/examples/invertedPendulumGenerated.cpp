
// CODEGEN FOR MODEL "invertedPendulum"

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "invertedPendulumGenerated.h"
#include "constants.h"

const double mS = 1;
const double mP = 0.5;
const double r = 1;


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
        const double s1 = mP*pow(s0, 2) + mS;
		return (u[0]*s1 - mP*r*pow(x[3], 2)*s0 - 4.9050000000000002*mP*sin(2*x[2]))/s1;
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s18 = sin(x[2]);
        const double s19 = mP*pow(s18, 2) + mS;
        const double s20 = 1.0/s19;
        const double s21 = 2*x[2];
        const double s22 = cos(x[2]);
        const double s23 = mP*s22;
        const double s24 = r*pow(x[3], 2);
        const double s25 = mP*s18;
		return {std::vector<double>{-2*s18*s23*(u[0]*s19 - 4.9050000000000002*mP*sin(s21) - s24*s25)/pow(s19, 2) + s20*(2*u[0]*mP*s18*s22 - 9.8100000000000005*mP*cos(s21) - s23*s24), -2*r*x[3]*s20*s25}, {1}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s2 = 2*x[2];
        const double s3 = sin(s2);
        const double s4 = sin(x[2]);
        const double s5 = pow(s4, 2);
        const double s6 = 2*u[0];
        const double s7 = cos(x[2]);
        const double s8 = pow(s7, 2);
        const double s9 = r*pow(x[3], 2);
        const double s10 = s4*s9;
        const double s11 = s4*s7;
        const double s12 = mP*s5;
        const double s13 = mS + s12;
        const double s14 = 1.0/s13;
        const double s15 = mP*s14;
        const double s16 = 2*s14;
        const double s17 = 2*r*s15;
		return {std::vector<double>{s15*(s10 + 4*s11*s15*(-s11*s6 + s7*s9 + 9.8100000000000005*cos(s2)) - s16*(-u[0]*s13 + mP*s10 + 4.9050000000000002*mP*s3)*(4*s12*s14*s8 + s5 - s8) + 19.620000000000001*s3 - s5*s6 + s6*s8), x[3]*s17*s7*(s12*s16 - 1), -s17*s4}, {}, {}, {}, {}, {}};
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
        const double s1 = mP*pow(s0, 2) + mS;
		return (9.8100000000000005*s0*s1 - (u[0]*s1 - mP*r*pow(x[3], 2)*s0 - 4.9050000000000002*mP*sin(2*x[2]))*cos(x[2]))/(r*s1);
	}

	std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        const double s26 = sin(x[2]);
        const double s27 = mP*pow(s26, 2);
        const double s28 = mS + s27;
        const double s29 = 9.8100000000000005*s28;
        const double s30 = cos(x[2]);
        const double s31 = 2*x[2];
        const double s32 = mP*s26;
        const double s33 = r*pow(x[3], 2);
        const double s34 = u[0]*s28 - 4.9050000000000002*mP*sin(s31) - s32*s33;
        const double s35 = 1.0/r;
        const double s36 = s30*s35;
        const double s37 = 2*s32;
        const double s38 = 1.0/s28;
		return {std::vector<double>{s35*s38*(s26*s34 + 19.620000000000001*s27*s30 + s29*s30 - s30*(2*u[0]*mP*s26*s30 - mP*s30*s33 - 9.8100000000000005*mP*cos(s31))) - s36*s37*(s26*s29 - s30*s34)/pow(s28, 2), x[3]*s30*s37*s38}, {-s36}, {}};
	}

	std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        const double s2 = 1.0/r;
        const double s3 = sin(x[2]);
        const double s4 = pow(s3, 2);
        const double s5 = mP*s4;
        const double s6 = mS + s5;
        const double s7 = 1.0/s6;
        const double s8 = 2*x[2];
        const double s9 = sin(s8);
        const double s10 = 2*u[0];
        const double s11 = cos(x[2]);
        const double s12 = pow(s11, 2);
        const double s13 = r*pow(x[3], 2);
        const double s14 = s13*s3;
        const double s15 = mP*s11;
        const double s16 = -s10*s11*s3 + s11*s13 + 9.8100000000000005*cos(s8);
        const double s17 = 4*s7;
        const double s18 = s12*s5;
        const double s19 = -s12 + s4;
        const double s20 = 9.8100000000000005*s6;
        const double s21 = -u[0]*s6 + mP*s14 + 4.9050000000000002*mP*s9;
        const double s22 = s11*s21 + s20*s3;
        const double s23 = s15*s3;
        const double s24 = s2*s3;
        const double s25 = 2*s7;
		return {std::vector<double>{s2*s7*(58.859999999999999*mP*s12*s3 - 2*mP*s16*s3 + 2*mP*s22*s7*(s17*s18 + s19) - 19.620000000000001*mP*pow(s3, 3) - s15*(s10*s12 - s10*s4 + s14 + 19.620000000000001*s9) - s17*s23*(s11*s20 + 19.620000000000001*s11*s5 + s15*s16 - s21*s3) - s22), -mP*x[3]*s25*(s18*s25 + s19), s23*s25}, {s24}, {}, {}, {}, {}};
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
