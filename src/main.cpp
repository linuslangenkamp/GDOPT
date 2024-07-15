#include <iostream>
#include <IpIpoptApplication.hpp>
#include <memory>
#include "test_nlp.h"
#include "integrator.h"
#include "expression.h"
#include "constraint.h"
#include "mesh.h"
#include "problem.h"
#include "constants.h"
#include "gdop.h"

using namespace Ipopt;

class F0 : public Constraint {
public:
    static std::unique_ptr<F0> create() {
        Adjacency adj{
                {0, 1},
                {0, 1},
                {}
        };
        return std::unique_ptr<F0>(new F0(std::move(adj), 0., 0.));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return x[0] + x[1]*x[1] - u[0]*u[0] - u[1];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{1., 2*x[1]}, {-2*u[0], -1}, {}};
    }
private:
    F0(Adjacency adj, double lb, double ub) : Constraint(std::move(adj), lb, ub) {
    }
};

class F1 : public Constraint {
public:
    static std::unique_ptr<F1> create() {
        Adjacency adj{
                {1},
                {0, 1},
                {}
        };
        return std::unique_ptr<F1>(new F1(std::move(adj), 0., 0.));
    }

    double eval(const double *x, const double *u, const double *p, double t) override {
        return -x[1] - u[0] - u[1];
    }

    std::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {
        return {std::vector<double>{-1}, {-1, -1}, {}};
    }
private:
    F1(Adjacency adj, double lb, double ub) : Constraint(std::move(adj), lb, ub) {
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


int main() {
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps3);
    Mesh mesh = Mesh::createEquidistantMesh(50, 2);

    std::vector<std::unique_ptr<Constraint>> F;
    F.push_back(F0::create());
    F.push_back(F1::create());

    Problem problem(
            2, 2, 0,
            {0, 1}, {MIN_DOUBLE, MIN_DOUBLE}, {MAX_DOUBLE, MAX_DOUBLE},
            {0, 0}, {2, 2},
            {}, {},
            Mayer::create(),
            Lagrange::create(),
            std::move(F),
            {},
            {},
            {}
    );
    GDOP gdop(std::move(problem), mesh, rk);

    const std::vector<double> x(600, 1.5);
    std::vector<double> grad_f(600, -2.3);
    double obj_value = 0.0;

    gdop.eval_f(0, x.data(), true, obj_value);
    gdop.eval_grad_f(0, x.data(), true, grad_f.data());
    /*
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

    return (int) status;*/

    return 0;
}
