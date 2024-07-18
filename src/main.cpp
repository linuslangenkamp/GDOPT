#include <iostream>
#include <memory>
#include <IpIpoptApplication.hpp>
#include "test_nlp.h"
#include "integrator.h"
#include "expression.h"
#include "constraint.h"
#include "mesh.h"
#include "problem.h"
#include "constants.h"
#include "gdop.h"

using namespace Ipopt;

class Mayer : public Expression {
public:
    static std::unique_ptr<Mayer> create() {
        Adjacency adj{
                {0, 1},
                {},
                {}
        };
        return std::unique_ptr<Mayer>(new Mayer(std::move(adj)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[0] * (0.23 + x[1] * x[1]);
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{0.23 + x[1] * x[1], 2 * x[1] * x[0]}, {}, {}};
    }
private:
    Mayer(Adjacency adj) : Expression(std::move(adj)) {
    }
};

class Lagrange : public Expression {
public:
    static std::unique_ptr<Lagrange> create() {
        Adjacency adj{
                {0, 1},
                {},
                {}
        };
        return std::unique_ptr<Lagrange>(new Lagrange(std::move(adj)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[0]*x[0] - x[1] - 0.35;
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{2*x[0], -1}, {}, {}};
    }
private:
    Lagrange(Adjacency adj) : Expression(std::move(adj)) {
    }
};

class F0 : public Expression {
public:
    static std::unique_ptr<F0> create() {
        Adjacency adj{
                {0, 1},
                {0, 1},
                {}
        };
        return std::unique_ptr<F0>(new F0(std::move(adj)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[0] + x[1]*x[1] - u[0]*u[0] - u[1];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{1., 2*x[1]}, {-2*u[0], -1}, {}};
    }
private:
    F0(Adjacency adj) : Expression(std::move(adj)) {
    }
};

class F1 : public Expression {
public:
    static std::unique_ptr<F1> create() {
        Adjacency adj{
                {1},
                {0, 1},
                {}
        };
        return std::unique_ptr<F1>(new F1(std::move(adj)));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -x[1] - u[0] - u[1];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-1}, {-1, -1}, {}};
    }
private:
    F1(Adjacency adj) : Expression(std::move(adj)) {
    }
};

class G0 : public Constraint {
public:
    static std::unique_ptr<G0> create() {
        Adjacency adj{
                {0},
                {1},
                {}
        };
        return std::unique_ptr<G0>(new G0(std::move(adj), MIN_DOUBLE, 3.5));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[0]*x[0] + u[1]*u[1];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{2*x[0]}, {2*u[1]}, {}};
    }
private:
    G0(Adjacency adj, double lb, double ub) : Constraint(std::move(adj), lb, ub) {
    }
};

class R0 : public Constraint {
public:
    static std::unique_ptr<R0> create() {
        Adjacency adj{
                {0, 1},
                {},
                {}
        };
        return std::unique_ptr<R0>(new R0(std::move(adj), 1, MAX_DOUBLE));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -x[0]*x[0] + x[1]*x[1];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-2*x[0], 2*x[1]}, {}, {}};
    }
private:
    R0(Adjacency adj, double lb, double ub) : Constraint(std::move(adj), lb, ub) {
    }
};

int main() {
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps3);
    Mesh mesh = Mesh::createEquidistantMesh(50, 2);

    std::vector<std::unique_ptr<Expression>> F;
    F.push_back(F0::create());
    F.push_back(F1::create());

    std::vector<std::unique_ptr<Constraint>> G;
    G.push_back(G0::create());

    std::vector<std::unique_ptr<Constraint>> R;
    R.push_back(R0::create());

    Problem problem(
            2, 2, 0,
            {0, 1}, {MIN_DOUBLE, MIN_DOUBLE}, {MAX_DOUBLE, MAX_DOUBLE},
            {0, 0}, {2, 2},
            {}, {},
            Mayer::create(),
            Lagrange::create(),
            std::move(F),
            std::move(G),
            std::move(R),
            {}
    );
    GDOP gdop(std::move(problem), mesh, rk);
    int n, m, nJac, nHes;
    GDOP::IndexStyleEnum index_style;
    gdop.get_nlp_info(n, m, nJac, nHes, index_style);

    const std::vector<double> x(600, 2.3);
    std::vector<double> grad_f(600, -2.3);
    double obj_value = 0.0;

    gdop.eval_f(0, x.data(), true, obj_value);
    gdop.eval_grad_f(0, x.data(), true, grad_f.data());
    std::vector<double> g(451, 0);
    gdop.eval_g(0, x.data(), true, 451, g.data());

    SmartPtr<TestNLP> testNLP = new TestNLP;

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    app->Options()->SetNumericValue("tol", 1e-14);
    //app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("output_file", "ipopt.out");

    ApplicationReturnStatus status;
    status = app->Initialize();
    if( status != Solve_Succeeded )
    {
        std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
        return (int) status;
    }

    status = app->OptimizeTNLP(testNLP);

    if( status == Solve_Succeeded )
    {
        std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
    }
    else
    {
        std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
    }

    return (int) status;

    return 0;
}
