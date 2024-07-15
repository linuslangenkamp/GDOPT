//
// Created by Linus on 10.07.2024.
//

#ifndef IPOPT_DO_EXPRESSION_H
#define IPOPT_DO_EXPRESSION_H

#include <utility>
#include <vector>
#include <array>

struct Adjacency {
    const std::vector<int> indX;    // Indices X
    const std::vector<int> indU;    // Indices U
    const std::vector<int> indP;    // Indices P
};

class Expression {
public:
    explicit Expression(Adjacency adj) : adj{std::move(adj)} {}

    virtual double eval(const double* x, const double* u, const double* p, double t) = 0;

    // returns {evalDiff(indicesX), evalDiff(indicesU), evalDiff(indicesP)} - same sorting as adj.indices
    virtual std::array<std::vector<double>, 3> evalDiff(const double* x, const double* u, const double* p, double t) = 0;

    const Adjacency adj;
};

#endif //IPOPT_DO_EXPRESSION_H
