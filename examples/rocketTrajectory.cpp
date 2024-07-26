#include <string>
#include <cmath>
#include "expression.h"
#include "constraint.h"
#include "problem.h"
#include "constants.h"

const double earthX = 0;
const double earthY = 0;
const double moonX = 15;
const double moonY = 0;
const double c1 = 0.25;
const double c2 = 0.15;


class LagrangeRT : public Expression {
public:
    static std::unique_ptr<LagrangeRT> create() {
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
        return std::unique_ptr<LagrangeRT>(new LagrangeRT(std::move(adj), std::move(adjDiff)));
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
    LagrangeRT(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F0RT : public Expression {
public:
    static std::unique_ptr<F0RT> create() {
        Adjacency adj{
                {2},
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
        return std::unique_ptr<F0RT>(new F0RT(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[2];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{1.0}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }

private:
    F0RT(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F1RT : public Expression {
public:
    static std::unique_ptr<F1RT> create() {
        Adjacency adj{
                {3},
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
        return std::unique_ptr<F1RT>(new F1RT(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[3];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{1.0}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }

private:
    F1RT(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};


class F2RT : public Expression {
public:
    static std::unique_ptr<F2RT> create() {
        Adjacency adj{
                {0, 1},
                {0, 1},
                {}
        };
        AdjacencyDiff adjDiff{
                {{0,0}, {1, 0}, {1, 1}},
                {},
                {{1, 1}, {1, 0}},
                {},
                {},
                {}
        };
        return std::unique_ptr<F2RT>(new F2RT(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return c1*(earthX - x[0])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -1.5) + c2*(moonX - x[0])*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -1.5) + u[0]*cos(u[1]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{3.0*c1*pow(earthX - x[0], 2)*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -2.5) - c1*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -1.5) + 3.0*c2*pow(moonX - x[0], 2)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -2.5) - c2*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -1.5),
                                    3.0*(c1*(earthX - x[0])*(earthY - x[1])*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 2.5) + c2*(moonX - x[0])*(moonY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 2.5))*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -2.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -2.5)},
                                    {cos(u[1]), -u[0]*sin(u[1])}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-9.0*c1*earthX*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -2.5) + 9.0*c1*x[0]*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -2.5) + 15.0*c1*pow(earthX - x[0], 3)*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -3.5) - 9.0*c2*moonX*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -2.5) + 9.0*c2*x[0]*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -2.5) + 15.0*c2*pow(moonX - x[0], 3)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -3.5),
                                    pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -6.0)*(15.0*c1*pow(earthX - x[0], 2)*(earthY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 2.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 6.0) - 3.0*c1*(earthY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 3.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 6.0) + 15.0*c2*pow(moonX - x[0], 2)*(moonY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 2.5) - 3.0*c2*(moonY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 3.5)),
                                    pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -6.0)*(15.0*c1*(earthX - x[0])*pow(earthY - x[1], 2)*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 2.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 6.0) - 3.0*c1*(earthX - x[0])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 3.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 6.0) + 15.0*c2*(moonX - x[0])*pow(moonY - x[1], 2)*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 2.5) - 3.0*c2*(moonX - x[0])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 3.5))},
                                    {}, {-u[0]*cos(u[1]), -sin(u[1])}, {}, {}, {}};
    }

private:
    F2RT(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F3RT : public Expression {
public:
    static std::unique_ptr<F3RT> create() {
        Adjacency adj{
                {0, 1},
                {0, 1},
                {}
        };
        AdjacencyDiff adjDiff{
                {{0,0}, {1, 0}, {1, 1}},
                {},
                {{1, 1}, {1, 0}},
                {},
                {},
                {}
        };
        return std::unique_ptr<F3RT>(new F3RT(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return c1*(earthY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -1.5) + c2*(moonY - x[1])*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -1.5) + u[0]*sin(u[1]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{3.0*(c1*(earthX - x[0])*(earthY - x[1])*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 2.5) + c2*(moonX - x[0])*(moonY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 2.5))*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -2.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -2.5),
                                    3.0*c1*pow(earthY - x[1], 2)*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -2.5) - c1*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -1.5) + 3.0*c2*pow(moonY - x[1], 2)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -2.5) - c2*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -1.5)},
                                    {sin(u[1]), u[0]*cos(u[1])}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -6.0)*(15.0*c1*pow(earthX - x[0], 2)*(earthY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 2.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 6.0) - 3.0*c1*(earthY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 3.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 6.0) + 15.0*c2*pow(moonX - x[0], 2)*(moonY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 2.5) - 3.0*c2*(moonY - x[1])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 3.5)),
                                    pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -6.0)*(15.0*c1*(earthX - x[0])*pow(earthY - x[1], 2)*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 2.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 6.0) - 3.0*c1*(earthX - x[0])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 3.5)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 6.0) + 15.0*c2*(moonX - x[0])*pow(moonY - x[1], 2)*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 2.5) - 3.0*c2*(moonX - x[0])*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), 6.0)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), 3.5)),
                                    -9.0*c1*earthY*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -2.5) + 9.0*c1*x[1]*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -2.5) + 15.0*c1*pow(earthY - x[1], 3)*pow(pow(earthX - x[0], 2) + pow(earthY - x[1], 2), -3.5) - 9.0*c2*moonY*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -2.5) + 9.0*c2*x[1]*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -2.5) + 15.0*c2*pow(moonY - x[1], 3)*pow(pow(moonX - x[0], 2) + pow(moonY - x[1], 2), -3.5)},
                                    {}, {-u[0]*sin(u[1]), cos(u[1])}, {}, {}, {}};
    }

private:
    F3RT(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class G0RT : public Constraint {
public:
    static std::unique_ptr<G0RT> create() {
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
        return std::unique_ptr<G0RT>(new G0RT(std::move(adj), std::move(adjDiff), 1, PLUS_INFINITY));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return (x[0]*x[0] + x[1]*x[1]) / 1.3;
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{2.0/1.3 * x[0], 2.0/1.3 * x[1]}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{2.0/1.3, 2.0/1.3}, {}, {}, {}, {}, {}};
    }
private:
    G0RT(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

class G1RT : public Constraint {
public:
    static std::unique_ptr<G1RT> create() {
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
        return std::unique_ptr<G1RT>(new G1RT(std::move(adj), std::move(adjDiff), 0.25, PLUS_INFINITY));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return ((x[0]-moonX)*(x[0]-moonX) + (x[1]-moonY)*(x[1]-moonY));
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{2.0*(x[0]-moonX), 2.0*(x[1]-moonY)}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{2.0, 2.0}, {}, {}, {}, {}, {}};
    }
private:
    G1RT(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

class R0RT : public Constraint {
public:
    static std::unique_ptr<R0RT> create() {
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
        return std::unique_ptr<R0RT>(new R0RT(std::move(adj), std::move(adjDiff), 15 + 0.7 * pow(2, 0.5), 15 + 0.7 * pow(2, 0.5)));
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
    R0RT(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

class R1RT : public Constraint {
public:
    static std::unique_ptr<R1RT> create() {
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
        return std::unique_ptr<R1RT>(new R1RT(std::move(adj), std::move(adjDiff), 0.7 * pow(2, 0.5), 0.7 * pow(2, 0.5)));
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
    R1RT(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

class R2RT : public Constraint {
public:
    static std::unique_ptr<R2RT> create() {
        Adjacency adj{
                {2},
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
        return std::unique_ptr<R2RT>(new R2RT(std::move(adj), std::move(adjDiff), 0.5 / pow(2, 0.5), 0.5 / pow(2, 0.5)));
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
    R2RT(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

class R3RT : public Constraint {
public:
    static std::unique_ptr<R3RT> create() {
        Adjacency adj{
                {3},
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
        return std::unique_ptr<R3RT>(new R3RT(std::move(adj), std::move(adjDiff), -0.5 / pow(2, 0.5), -0.5 / pow(2, 0.5)));
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
    R3RT(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {
    }
};

Problem createProblem_rocketTrajectory() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0RT::create());
    F.push_back(F1RT::create());
    F.push_back(F2RT::create());
    F.push_back(F3RT::create());

    std::vector<std::unique_ptr<Constraint>> G;
    G.push_back(G0RT::create());
    G.push_back(G1RT::create());

    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0RT::create());
    R.push_back(R1RT::create());
    R.push_back(R2RT::create());
    R.push_back(R3RT::create());

    Problem problem(
            4, 2, 0,
            {-1., 0.5, 0, 0},
            {MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY, MINUS_INFINITY}, {PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY, PLUS_INFINITY},
            {0, -PI}, {3, PI},
            {}, {},
            {},
            LagrangeRT::create(),
            std::move(F),
            std::move(G),
            std::move(R),
            {},
            "RocketTrajectory");
    return problem;
};
