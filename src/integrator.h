#ifndef IPOPT_DO_INTEGRATOR_H
#define IPOPT_DO_INTEGRATOR_H

#include <vector>

enum class IntegratorSteps {
    Steps1,
    Steps2,
    Steps3,
    Steps4,
    Steps5,
    Steps6,
    Steps7
};

class Integrator {
public:
    static Integrator radauIIA(IntegratorSteps steps);

    const std::vector<double> c;
    const std::vector<std::vector<double>> A;
    const std::vector<std::vector<double>> Ainv;
    const std::vector<double> invRowSum;
    const std::vector<double> b;
    const int steps;
    const std::vector<std::vector<double>> firstLagrangeBasis;
    const std::vector<std::vector<double>> lagrangeBasis;

    std::vector<std::vector<double>> firstBasisPolynomial();       // inner basis poly, needed for 1st interval of u
    std::vector<std::vector<double>> basisPolynomial();         // standard collocation polynomial
    std::vector<std::vector<double>> basisPolynomialDiff();
    std::vector<std::vector<double>> basisPolynomialDiff2();
    std::vector<double> interpolateFirstControl(std::vector<double> &uValues);
    std::vector<double> interpolate(std::vector<double>&);
private:
    Integrator(const std::vector<double>& c,
               const std::vector<std::vector<double>>& A,
               const std::vector<std::vector<double>>& Ainv,
               const std::vector<double>& invRowSum,
               int steps);
};

#endif //IPOPT_DO_INTEGRATOR_H
