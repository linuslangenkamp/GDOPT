#include <string>
#include "trivialBangBang.h"
#include "expression.h"
#include "problem.h"
#include "constants.h"


class MayerBB : public Expression {
public:
    static std::unique_ptr<MayerBB> create() {
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
        return std::unique_ptr<MayerBB>(new MayerBB(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -x[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-1.}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }
private:
    MayerBB(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F0BB : public Expression {
public:
    static std::unique_ptr<F0BB> create() {
        Adjacency adj{
                {1},
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
        return std::unique_ptr<F0BB>(new F0BB(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[1];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{1.0}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }

private:
    F0BB(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F1BB : public Expression {
public:
    static std::unique_ptr<F1BB> create() {
        Adjacency adj{
                {},
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
        return std::unique_ptr<F1BB>(new F1BB(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return u[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {1.0}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }

private:
    F1BB(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};


class G0BB : public Constraint {
public:
    static std::unique_ptr<G0BB> create() {
        Adjacency adj{
                {1},
                {0},
                {}
        };
        AdjacencyDiff adjDiff{
                {},
                {{0, 1}},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<G0BB>(new G0BB(std::move(adj), std::move(adjDiff), -30.0, 30.0));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[1] * u[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{u[0]}, {x[1]}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {1.0}, {}, {}, {}, {}};
    }
private:
    G0BB(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

class R0BB : public Constraint {
public:
    static std::unique_ptr<R0BB> create() {
        Adjacency adj{
                {1},
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
        return std::unique_ptr<R0BB>(new R0BB(std::move(adj), std::move(adjDiff), 0, 0));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[1];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{1.0}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }
private:
    R0BB(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

Problem createProblem_trivalBangBang() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0BB::create());
    F.push_back(F1BB::create());

    std::vector<std::unique_ptr<Constraint>> G;
    G.push_back(G0BB::create());

    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0BB::create());

    Problem problem(
            2, 1, 0,
            {0, 0},
            {MINUS_INFINITY, MINUS_INFINITY}, {PLUS_INFINITY, PLUS_INFINITY},
            {-10}, {10},
            {}, {},
            MayerBB::create(),
            {},
            std::move(F),
            std::move(G),
            std::move(R),
            {},
            "trivialBangBang");
    return problem;
};
