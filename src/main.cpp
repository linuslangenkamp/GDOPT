#include <iostream>
#include <IpIpoptApplication.hpp>
#include "test_nlp.h"
#include "integrator.h"

using namespace Ipopt;

int main() {
    Integrator RK = Integrator::radauIIA(IntegratorSteps::Steps3);
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
}
