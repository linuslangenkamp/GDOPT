//
// Created by Linus on 09.07.2024.
//
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
private:
    Integrator(const std::vector<double>& c,
               const std::vector<std::vector<double>>& A,
               const std::vector<std::vector<double>>& Ainv,
               const std::vector<double>& invRowSum,
               int steps);
};

#endif //IPOPT_DO_INTEGRATOR_H
