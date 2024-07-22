#include "expression.h"
#include "problem.h"
#include "constants.h"

constexpr double OMEGA_0 = 1;

class LagrangeP : public Expression {
public:
    static std::unique_ptr<LagrangeP> create() {
        Adjacency adj{
                {0, 1},
                {},
                {}
        };
        AdjacencyDiff adjDiff{
                        {{0, 0}, {1, 1}},
                        {},
                        {},
                        {},
                        {},
                        {}
        };
        return std::unique_ptr<LagrangeP>(new LagrangeP(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return 0.5 * (x[0]*x[0] + OMEGA_0*OMEGA_0*x[1]*x[1]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{x[0], OMEGA_0*OMEGA_0*x[1]}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{1, OMEGA_0*OMEGA_0}, {}, {}, {}, {}, {}};
    }
private:
    LagrangeP(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F0P : public Expression {
public:
    static std::unique_ptr<F0P> create() {
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
        return std::unique_ptr<F0P>(new F0P(std::move(adj), std::move(adjDiff)));
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
    F0P(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F1P : public Expression {
public:
    static std::unique_ptr<F1P> create() {
        Adjacency adj{
                    {0, 1},
                    {},
                    {0}
        };
        AdjacencyDiff adjDiff{
                                {},
                                {},
                                {},
                                {{0, 1}},
                                {},
                                {}
        };
        return std::unique_ptr<F1P>(new F1P(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -2*p[0]*x[1] - OMEGA_0*OMEGA_0*x[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-OMEGA_0*OMEGA_0, -2*p[0]}, {}, {-2*x[1]}};
    }
    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {-2}, {}, {}};
    }

private:
    F1P(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

Problem createProblem_parameterSweep() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0P::create());
    F.push_back(F1P::create());

    Problem problem(
            2, 0, 1,
            {1, 1},
            {MINUS_INFINITY, MINUS_INFINITY}, {PLUS_INFINITY, PLUS_INFINITY},
            {}, {},
            {0.0}, {0.7},
            {},
            LagrangeP::create(),
            std::move(F),
            {},
            {},
            {});
    return problem;
};
