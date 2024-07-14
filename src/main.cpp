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
    static F0 create() {
        Adjacency adj{
                {0, 1},
                {0, 1},
                {}
        };
        return F0(std::move(adj), 0., 0.);
    }

    double eval(double *x, double *u, double *p, double t) override {
        return x[0] + x[1]*x[1] - u[0]*u[0] - u[1];
    }

    std::array<std::vector<double>, 3> evalDiff(double *x, double *u, double *p, double t) override {
        return {std::vector<double>{1., 2*x[1]}, {-2*u[0], -1}, {}};
    }
private:
    F0(Adjacency adj, double lb, double ub) : Constraint(std::move(adj), lb, ub) {
    }
};

class F1 : public Constraint {
public:
    static F1 create() {
        Adjacency adj{
                {1},
                {0, 1},
                {}
        };
        return F1(std::move(adj), 0., 0.);
    }

    double eval(double *x, double *u, double *p, double t) override {
        return -x[1] - u[0] - u[1];
    }

    std::array<std::vector<double>, 3> evalDiff(double *x, double *u, double *p, double t) override {
        return {std::vector<double>{-1}, {-1, -1}, {}};
    }
private:
    F1(Adjacency adj, double lb, double ub) : Constraint(std::move(adj), lb, ub) {
    }
};

class Lagrange : public Expression {
public:
    static Lagrange create() {
        Adjacency adj{
                {0, 1},
                {},
                {}
        };
        return Lagrange(std::move(adj));
    }

    double eval(double *x, double *u, double *p, double t) override {
        return x[0] - x[1];
    }

    std::array<std::vector<double>, 3> evalDiff(double *x, double *u, double *p, double t) override {
        return {std::vector<double>{1, -1}, {}, {}};
    }
private:
    Lagrange(Adjacency adj) : Expression(std::move(adj)) {
    }
};

int main() {
    Integrator rk = Integrator::radauIIA(IntegratorSteps::Steps3);
    Mesh mesh = Mesh::createEquidistantMesh(50, 2);

    auto f0 = F0::create();
    auto f1 = F1::create();

    std::vector<std::unique_ptr<Constraint>> F;
    F.push_back(std::make_unique<F0>(std::move(f0)));
    F.push_back(std::make_unique<F1>(std::move(f1)));

    auto lagrange = Lagrange::create();
    Problem problem(
            2, 2, 0,
            {0, 1}, {MIN_DOUBLE, MIN_DOUBLE}, {MAX_DOUBLE, MAX_DOUBLE},
            {0, 0}, {2, 2},
            {}, {},
            nullptr,
            std::make_unique<Lagrange>(std::move(lagrange)),
            std::move(F),
            {},
            {},
            {}
    );
    GDOP gdop(std::move(problem), mesh, rk);

    std::vector<double> vars = {1, 4, 5, 2};
    auto out = gdop.problem.F[1]->eval(&vars[0], &vars[2], &vars[4], 2);

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

    return (int) status;
     */
    return 0;
}
