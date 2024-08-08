#include <string>
#include <cmath>
#include "expression.h"
#include "constraint.h"
#include "problem.h"
#include "constants.h"

/*
 * see:  Tamimi, J.: Development of the Efficient Algorithms for Model Predictive Control of Fast Systems.
 * PhD Thesis, Technische Universität Ilmenau, VDI Verlag, 2011 Chapter: 7.4
 * also OpenModelica DynamicOptimization testsuite
 * tf = 100 with optimal objective f(x*) = 0.463944283849
*/
const double I1 = 1000000.;
const double I2 = 833333.;
const double I3 = 916667.;
const double T1S = 550.;
const double T2S = 50.;
const double T3S = 550.;

class LagrangeSat : public Expression {
public:
    static std::unique_ptr<LagrangeSat> create() {
        Adjacency adj{
                {},
                {0, 1, 2},
                {}
        };
        AdjacencyDiff adjDiff{
                {},
                {},
                {{0, 0}, {1, 1}, {2, 2}},
                {},
                {},
                {}
        };
        return std::unique_ptr<LagrangeSat>(new LagrangeSat(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return 0.5 * (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {u[0], u[1], u[2]}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {1, 1, 1}, {}, {}, {}};
    }
private:
    LagrangeSat(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class MayerSat : public Expression {
public:
    static std::unique_ptr<MayerSat> create() {
        Adjacency adj{
                {0, 1, 2, 3, 4, 5, 6},
                {},
                {}
        };
        AdjacencyDiff adjDiff{
                {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<MayerSat>(new MayerSat(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return (x[0] - 0.70106)*(x[0] - 0.70106) + (x[1] - 0.0923)*(x[1] - 0.0923) + (x[2] - 0.56098)*(x[2] - 0.56098) +
                (x[3] - 0.43047)*(x[3] - 0.43047) + x[4] * x[4] + x[5] * x[5] + x[6] * x[6];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{2 * (x[0] - 0.70106), 2 * (x[1] - 0.0923), 2 * (x[2] - 0.56098), 2 * (x[3] - 0.43047),
        2 * x[4], 2 * x[5], 2 * x[6]}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{2, 2, 2, 2, 2, 2, 2}, {}, {}, {}, {}, {}};
    }
private:
    MayerSat(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F0Sat : public Expression {
public:
    static std::unique_ptr<F0Sat> create() {
        Adjacency adj{
                {1, 2, 3, 4, 5, 6},
                {},
                {}
        };
        AdjacencyDiff adjDiff{
                {{4, 3}, {5, 2}, {6, 1}},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<F0Sat>(new F0Sat(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return 0.5 * (x[4] * x[3] - x[5] * x[2] + x[6] * x[1]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{0.5 * x[6], -0.5 * x[5], 0.5 * x[4], 0.5 * x[3], -0.5 * x[2], 0.5 * x[1]}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{0.5, -0.5, 0.5}, {}, {}, {}, {}, {}};
    }

private:
    F0Sat(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F1Sat : public Expression {
public:
    static std::unique_ptr<F1Sat> create() {
        Adjacency adj{
                {0, 2, 3, 4, 5, 6},
                {},
                {}
        };
        AdjacencyDiff adjDiff{
                {{4, 2}, {5, 3}, {6, 0}},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<F1Sat>(new F1Sat(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return 0.5 * (x[4] * x[2] + x[5] * x[3] - x[6] * x[0]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-0.5 * x[6], 0.5 * x[4], 0.5 * x[5], 0.5 * x[2], 0.5 * x[3], -0.5 * x[0]}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{0.5, 0.5, -0.5}, {}, {}, {}, {}, {}};
    }

private:
    F1Sat(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F2Sat : public Expression {
public:
    static std::unique_ptr<F2Sat> create() {
        Adjacency adj{
                {0, 1, 3, 4, 5, 6},
                {},
                {}
        };
        AdjacencyDiff adjDiff{
                {{4, 1}, {5, 0}, {6, 3}},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<F2Sat>(new F2Sat(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return 0.5 * (-x[4] * x[1] + x[5] * x[0] + x[6] * x[3]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{0.5 * x[5], -0.5 * x[4], 0.5 * x[6], -0.5 * x[1], 0.5 * x[0], 0.5 * x[3]}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-0.5, 0.5, 0.5}, {}, {}, {}, {}, {}};
    }

private:
    F2Sat(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F3Sat : public Expression {
public:
    static std::unique_ptr<F3Sat> create() {
        Adjacency adj{
                {0, 1, 2, 4, 5, 6},
                {},
                {}
        };
        AdjacencyDiff adjDiff{
                {{4, 0}, {5, 1}, {6, 2}},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<F3Sat>(new F3Sat(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -0.5 * (x[4] * x[0] + x[5] * x[1] + x[6] * x[2]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-0.5 * x[4], -0.5 * x[5], -0.5 * x[6], -0.5 * x[0], -0.5 * x[1], -0.5 * x[2]}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-0.5, -0.5, -0.5}, {}, {}, {}, {}, {}};
    }

private:
    F3Sat(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F4Sat : public Expression {
public:
    static std::unique_ptr<F4Sat> create() {
        Adjacency adj{
                {5, 6},
                {0},
                {}
        };
        AdjacencyDiff adjDiff{
                {{6, 5}},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<F4Sat>(new F4Sat(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return ((I2 - I3) * x[5] * x[6] + T1S * u[0]) / I1 ;
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{(I2 - I3) / I1 * x[6], (I2 - I3) / I1 * x[5]}, {T1S / I1}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{(I2 - I3) / I1}, {}, {}, {}, {}, {}};
    }

private:
    F4Sat(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F5Sat : public Expression {
public:
    static std::unique_ptr<F5Sat> create() {
        Adjacency adj{
                {4, 6},
                {1},
                {}
        };
        AdjacencyDiff adjDiff{
                {{6, 4}},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<F5Sat>(new F5Sat(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return ((I3 - I1) * x[6] * x[4] + T2S * u[1]) / I2;
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{(I3 - I1) / I2 * x[6], (I3 - I1) / I2 * x[4]}, {T2S / I2}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{(I3 - I1) / I2}, {}, {}, {}, {}, {}};
    }

private:
    F5Sat(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F6Sat : public Expression {
public:
    static std::unique_ptr<F6Sat> create() {
        Adjacency adj{
                {4, 5},
                {2},
                {}
        };
        AdjacencyDiff adjDiff{
                {{5, 4}},
                {},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<F6Sat>(new F6Sat(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return ((I1 - I2) * x[4] * x[5] + T3S * u[2]) / I3;
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{(I1 - I2) / I3 * x[5], (I1 - I2) / I3 * x[4]}, {T3S / I3}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{(I1 - I2) / I3}, {}, {}, {}, {}, {}};
    }

private:
    F6Sat(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

Problem createProblem_satellite() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0Sat::create());
    F.push_back(F1Sat::create());
    F.push_back(F2Sat::create());
    F.push_back(F3Sat::create());
    F.push_back(F4Sat::create());
    F.push_back(F5Sat::create());
    F.push_back(F6Sat::create());

    Problem problem(
            7, 3, 0,
            {0, 0, 0, 1, 0.01, 0.005, 0.001},
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY},
            {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},
            {0, 0, 0}, {1, 1, 1},
            {}, {},
            MayerSat::create(),
            LagrangeSat::create(),
            std::move(F),
            {},
            {},
            {},
            "satellite");
    return problem;
};
