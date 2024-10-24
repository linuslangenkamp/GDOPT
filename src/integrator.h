/**
 * GDOPT - General Dynamic Optimizer
 * Copyright (C) 2024  Linus Langenkamp
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **/

#ifndef GDOPT_INTEGRATOR_H
#define GDOPT_INTEGRATOR_H

#include <vector>

enum class IntegratorSteps {
    Steps1 = 1,
    Steps2 = 2,
    Steps3 = 3,
    Steps4 = 4,
    Steps5 = 5,
    Steps6 = 6,
    Steps7 = 7,
    Steps8 = 8,
    Steps9 = 9,
    Steps10 = 10,
    Steps11 = 11,
    Steps12 = 12,
    Steps13 = 13,
    Steps14 = 14,
    Steps15 = 15,
    Steps16 = 16,
    Steps17 = 17,
    Steps18 = 18,
    Steps19 = 19,
    Steps20 = 20,
    Steps21 = 21,
    Steps22 = 22,
    Steps23 = 23,
    Steps24 = 24,
    Steps25 = 25,
    Steps26 = 26,
    Steps27 = 27,
    Steps28 = 28,
    Steps29 = 29,
    Steps30 = 30,
    Steps31 = 31,
    Steps32 = 32,
    Steps33 = 33,
    Steps34 = 34,
    Steps35 = 35,
    Steps36 = 36,
    Steps37 = 37,
    Steps38 = 38,
    Steps39 = 39,
    Steps40 = 40,
    Steps41 = 41,
    Steps42 = 42,
    Steps43 = 43,
    Steps44 = 44,
    Steps45 = 45,
    Steps46 = 46,
    Steps47 = 47,
    Steps48 = 48,
    Steps49 = 49,
    Steps50 = 50,
    Steps51 = 51,
    Steps52 = 52,
    Steps53 = 53,
    Steps54 = 54,
    Steps55 = 55,
    Steps56 = 56,
    Steps57 = 57,
    Steps58 = 58,
    Steps59 = 59,
    Steps60 = 60,
    Steps61 = 61,
    Steps62 = 62,
    Steps63 = 63,
    Steps64 = 64,
    Steps65 = 65,
    Steps66 = 66,
    Steps67 = 67,
    Steps68 = 68,
    Steps69 = 69,
    Steps70 = 70,
    Steps71 = 71,
    Steps72 = 72,
    Steps73 = 73,
    Steps74 = 74,
    Steps75 = 75,
    Steps76 = 76,
    Steps77 = 77,
    Steps78 = 78
};

class Integrator {
public:
    static Integrator radauIIA(IntegratorSteps steps);
    static Integrator testIntegrator(const double a);
    static Integrator testCollocation2Step(const double c1);
    const std::vector<double> c;
    const std::vector<double> c0;  // c including 0 point, i.e. [0, c1, c2, ..., cm]
    const std::vector<std::vector<double>> A;
    const std::vector<std::vector<double>> Ainv;
    const std::vector<double> invRowSum;
    const std::vector<double> b;
    const std::vector<double> cBisection;  // [1/2 * j + 1/2 * c_i] for j=0, 1 and i = 1, ..., m
    const int steps;

    double integrate(std::vector<double>&);

    // eval lagrange poly based on coefficients, values at some x; O(n^2)
    static double evalLagrange(std::vector<double>, std::vector<double>&, double);

    // interpolation stuff, lagrange basis at c_j/2, c_j/2 + 1/2 for j=0...m
    std::vector<std::vector<double>> interpolationFirstBasisPolynomial();  // inner basis poly, needed for 1st interval of u
    std::vector<std::vector<double>> interpolationBasisPolynomial();       // standard collocation polynomial
    std::vector<double> interpolateFirstControl(std::vector<double>& uValues);
    std::vector<double> evalInterpolationNewNodes(std::vector<double>& values);
    const std::vector<std::vector<double>> interpolationFirstLagrangeBasis;
    const std::vector<std::vector<double>> interpolationLagrangeBasis;

    std::vector<double> evalLinearSplineNewNodes(std::vector<double>& values);

    // all basis coefficients at all c_j for p_u, p_u', p_u''
    std::vector<double> evalLagrangeDiff(std::vector<double>&);
    std::vector<double> evalLagrangeDiff2(std::vector<double>&);

    std::vector<std::vector<double>> basisPolynomialDiff();
    std::vector<std::vector<double>> basisPolynomialDiff2();
    const std::vector<std::vector<double>> lagrangeBasisDiff;
    const std::vector<std::vector<double>> lagrangeBasisDiff2;

private:
    Integrator(const std::vector<double>& c, const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& Ainv,
               const std::vector<double>& invRowSum, int steps);
};

#endif  // GDOPT_INTEGRATOR_H
