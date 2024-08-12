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
    Steps7,
    Steps8,
    Steps9,
    Steps10,
    Steps11,
    Steps12,
    Steps13,
    Steps14,
    Steps15,
    Steps16,
    Steps17,
    Steps18,
    Steps19,
    Steps20,
    Steps21,
    Steps22,
    Steps23,
    Steps24,
    Steps25,
    Steps26,
    Steps27,
    Steps28,
    Steps29,
    Steps30,
    Steps31,
    Steps32,
    Steps33,
    Steps34,
    Steps35,
    Steps36
};

class Integrator {
public:
    static Integrator radauIIA(IntegratorSteps steps);

    const std::vector<double> c;
    const std::vector<double> c0; // c including 0 point, i.e. [0, c1, c2, ..., cm]
    const std::vector<std::vector<double>> A;
    const std::vector<std::vector<double>> Ainv;
    const std::vector<double> invRowSum;
    const std::vector<double> b;
    const int steps;

    double integrate(std::vector<double> &);

    // eval lagrange poly based on coefficients, values at some x; O(n^2)
    static double evalLagrange(std::vector<double>, std::vector<double>&, double);

    // interpolation stuff, lagrange basis at c_j/2, c_j/2 + 1/2 for j=0...m
    std::vector<std::vector<double>> interpolationFirstBasisPolynomial();       // inner basis poly, needed for 1st interval of u
    std::vector<std::vector<double>> interpolationBasisPolynomial();         // standard collocation polynomial
    std::vector<double> interpolateFirstControl(std::vector<double> &uValues);
    std::vector<double> interpolate(std::vector<double>&);
    const std::vector<std::vector<double>> interpolationFirstLagrangeBasis;
    const std::vector<std::vector<double>> interpolationLagrangeBasis;

    // all basis coefficients at all c_j for p_u, p_u', p_u''
    std::vector<std::vector<double>> basisPolynomialDiff();
    std::vector<std::vector<double>> basisPolynomialDiff2();
    std::vector<double> evalLagrangeDiff(std::vector<double>&);
    std::vector<double> evalLagrangeDiff2(std::vector<double>&);
    const std::vector<std::vector<double>> lagrangeBasisDiff;
    const std::vector<std::vector<double>> lagrangeBasisDiff2;
private:
    Integrator(const std::vector<double>& c,
               const std::vector<std::vector<double>>& A,
               const std::vector<std::vector<double>>& Ainv,
               const std::vector<double>& invRowSum,
               int steps);
};

#endif //IPOPT_DO_INTEGRATOR_H
