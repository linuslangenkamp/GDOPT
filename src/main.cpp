#include <iostream>
#include <IpIpoptApplication.hpp>
#include "test_nlp.h"
#include "integrator.h"
#include "problem.h"
#include "expression.h"
#include "constraint.h"

using namespace Ipopt;

class TestConstraint : public Constraint {
public:
    static TestConstraint create() {
        Adjacency adj{
                {0, 2},
                {},
                {}
        };

        return TestConstraint(std::move(adj), 0., 10.);
    }


    double eval(const std::vector<double> &x, const std::vector<double> &u, const std::vector<double> &p,
         double t) override {
        return x[0] + x[2];
    }

    std::array<std::vector<double>, 3> evalDiff(const std::vector<double> &x, const std::vector<double> &u,
                                                const std::vector<double> &p, double t) override {
        return {std::vector<double>{1., 1.}, {}, {}};
    }
private:
    TestConstraint(Adjacency adj, double lb, double ub) : Constraint(std::move(adj), lb, ub) {

    }

};

int main() {
    Integrator RK = Integrator::radauIIA(IntegratorSteps::Steps3);

    auto testConstraint = TestConstraint::create();

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
