#include <string>
#include "batchReactor.h"
#include "expression.h"

/**
 * Batch Reactor from Parallel Multiple-Shooting and Collocation Optimization with OpenModelica,
 * Bachmann, Ochel, et. al., 2012

 optimal objectives: f(x*) = -0.57354505750936147 with n = 1e5, m = 3
                             -0.57354505720920057 with n = 5e5, m = 3
                             -0.57354505683399926 with n = 1e6, m = 3, time = 140.658s, MA57
                             -0.57354505633373065 with n = 1e6, m = 5, time = 327.169s, MA57
                             -0.57354505723421423 with n = 2e5, m = 7, time = 127.079s, MA57
                             -0.57354505670893241 with n = 5e5, m = 7, time = 223.389s, MA57

 model batchReactor
  Real x1(start=1, fixed=true, min=0, max=1);
  Real x2(start=0, fixed=true, min=0, max=1);
  Real may = -x2 annotation(isMayer=true);
  input Real u(min=0, max=5);
equation
  der(x1) = -(u + u^2/2) * x1;
  der(x2) = u * x1;
  annotation(experiment(StartTime=0, StopTime=1, Tolerance=1e-14),
__OpenModelica_simulationFlags(solver="optimization", optimizerNP="3"),
__OpenModelica_commandLineOptions="+g=Optimica");
end BatchReactor;
**/

class MayerBR : public Expression {
public:
    static std::unique_ptr<MayerBR> create() {
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
        return std::unique_ptr<MayerBR>(new MayerBR(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -x[1];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-1.}, {}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {}, {}, {}, {}, {}};
    }
private:
    MayerBR(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F0BR : public Expression {
public:
    static std::unique_ptr<F0BR> create() {
        Adjacency adj{
                {0},
                {0},
                {}
        };
        AdjacencyDiff adjDiff{
                {},
                {{0, 0}},
                {{0, 0}},
                {},
                {},
                {}
        };
        return std::unique_ptr<F0BR>(new F0BR(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -(u[0] + u[0]*u[0] / 2) * x[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-(u[0] + u[0]*u[0] / 2)}, {-(1 + u[0]) * x[0]}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {-(1 + u[0])}, {-x[0]}, {}, {}, {}};
    }

private:
    F0BR(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

class F1BR : public Expression {
public:
    static std::unique_ptr<F1BR> create() {
        Adjacency adj{
                {0},
                {0},
                {}
        };
        AdjacencyDiff adjDiff{
                {},
                {{0, 0}},
                {},
                {},
                {},
                {}
        };
        return std::unique_ptr<F1BR>(new F1BR(std::move(adj), std::move(adjDiff)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return u[0] * x[0];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{u[0]}, {x[0]}, {}};
    }

    std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{}, {1.}, {}, {}, {}, {}};
    }

private:
    F1BR(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {
    }
};

Problem createProblem_batchReactor() {

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0BR::create());
    F.push_back(F1BR::create());

    Problem problem(
            2, 1, 0,
            {1, 0},
            {0, 0}, {1, 1},
            {0}, {5},
            {}, {},
            MayerBR::create(),
            {},
            std::move(F),
            {},
            {},
            {},
            "BatchReactor");
    return problem;
};
