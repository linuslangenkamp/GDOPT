#include <string>
#include "expression.h"
#include "constraint.h"
#include "problem.h"
#include "constants.h"

/**
 * Hypersensitive OCP with an analytic solution from
 * "Symplectic algorithms with mesh refinement for a hypersensitive optimal control problem"
 * Peng, et. al., 2014
 * x*(t) = (1.5 * exp((t-2*t_f)*2^0.5) - exp((t-t_f)*2^0.5) + exp((-t-t_f)*2^0.5) - 1.5 * exp(-t*2^0.5))  / (exp(-t_f * 2 * 2^0.5) - 1)
 * u*(t) = x*(t) + d/dt x*(t)
 * f(x*) = int_0^tf 1/2 (x*(t)^2 + u*(t)^2) dt
 * Analytic solutions:
 * tf = 25,    f(x*) = 1.673097038856277579275453
 * tf = 100,   f(x*) = 1.673097038856279454302744
 * tf = 1000,  f(x*) = 1.673097038856279454302744
 * tf = 10000, f(x*) = 1.673097038856279454302744
**/

class Lagrange : public Expression {
public:
    static std::unique_ptr<Lagrange> create() {
        Adjacency adj{
                {0},
                {0},
                {}
        };
        AdjacencyDiff adjDiff{
                {{0, 0}},
                {},
                {{0, 0}},
                {},
                {},
                {}
        };
        return std::unique_ptr<Lagrange>(new Lagrange(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return 0.5 * (x[0]*x[0] + u[0]*u[0]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{x[0]}, {u[0]}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{1}, {}, {1}, {}, {}, {}};
    }
private:
    Lagrange(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F0 : public Expression {
public:
    static std::unique_ptr<F0> create() {
        Adjacency adj{
                {0},
                {0},
                {}
        };
        AdjacencyDiff adjDiff{
                {},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<F0>(new F0(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -x[0] + u[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-1}, {1}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }

private:
    F0(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class R0 : public Constraint {
public:
    static std::unique_ptr<R0> create() {
        Adjacency adj{
                {0},
                {},
                {}
        };
        AdjacencyDiff adjDiff{
                {},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<R0>(new R0(std::move(adj), std::move(adjDiff), 0, 0));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return 1. - x[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-1}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }
private:
    R0(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

Problem createProblem_analyticHypersensitive() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0::create());

    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0::create());

    Problem problem(
            1, 1, 0,
            {1.5},
            {MINUS_INFINITY}, {PLUS_INFINITY},
            {MINUS_INFINITY}, {PLUS_INFINITY},
            {}, {},
            {},
            Lagrange::create(),
            std::move(F),
            {},
            std::move(R),
            {},
            "analyticHypersensitive");
    return problem;
};
