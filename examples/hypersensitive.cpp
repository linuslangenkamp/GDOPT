#include <string>
#include "expression.h"
#include "constraint.h"
#include "problem.h"
#include "constants.h"

class LagrangeHyp : public Expression {
public:
    static std::unique_ptr<LagrangeHyp> create() {
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
        return std::unique_ptr<LagrangeHyp>(new LagrangeHyp(std::move(adj), std::move(adjDiff)));
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
    LagrangeHyp(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F0Hyp : public Expression {
public:
    static std::unique_ptr<F0Hyp> create() {
        Adjacency adj{
                {0},
                {0},
                {}
        };
        AdjacencyDiff adjDiff{
                            {{0, 0}},
                            {},
                            {},
                            {},
                            {},
                            {}
        };
        return std::unique_ptr<F0Hyp>(new F0Hyp(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -x[0] * x[0] * x[0] + u[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-3*x[0]*x[0]}, {1}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-6*x[0]}, {}, {}, {}, {}, {}};
    }

private:
    F0Hyp(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class R0Hyp : public Constraint {
public:
    static std::unique_ptr<R0Hyp> create() {
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
        return std::unique_ptr<R0Hyp>(new R0Hyp(std::move(adj), std::move(adjDiff), 0, 0));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return 1.5 - x[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-1}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }
private:
    R0Hyp(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

Problem createProblem_hypersensitive() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0Hyp::create());

    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0Hyp::create());

    Problem problem(
            1, 1, 0,
            {1},
            {MINUS_INFINITY}, {PLUS_INFINITY},
            {MINUS_INFINITY}, {PLUS_INFINITY},
            {}, {},
            {},
            LagrangeHyp::create(),
            std::move(F),
            {},
            std::move(R),
            {},
            "Hypersensitive");
    return problem;
};
