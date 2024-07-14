#include <iostream>
#include <IpIpoptApplication.hpp>
#include <memory>
#include "test_nlp.h"
#include "integrator.h"
#include "expression.h"
#include "constraint.h"


using namespace Ipopt;

class TestConstraint : public Constraint {
public:
    static TestConstraint create() {
        Adjacency adj{
                {0, 2},
                {0},
                {}
        };
        return TestConstraint(std::move(adj), 0., 10.);
    }

    double eval(double *x, double *u, double *p, double t) override {
        return x[0] + x[2] * u[0];
    }

    std::array<std::vector<double>, 3> evalDiff(double *x, double *u, double *p, double t) override {
        return {std::vector<double>{1., u[0]}, {x[2]}, {}};
    }
private:
    TestConstraint(Adjacency adj, double lb, double ub) : Constraint(std::move(adj), lb, ub) {
    }
};

int main() {
    Integrator RK = Integrator::radauIIA(IntegratorSteps::Steps3);
    auto testConstraint = TestConstraint::create();
    std::vector<double> vars = {1, 4, 5, 2, 5, 2};
    auto out = testConstraint.evalDiff(&vars[0], &vars[3], &vars[6], 2);

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
