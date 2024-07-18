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

    // returns {evalDiff(indX), evalDiff(indU), evalDiff(indP)} - same sorting as adj.indices
    virtual std::array<std::vector<double>, 3> evalDiff(const double* x, const double* u, const double* p, double t) = 0;

    const Adjacency adj;
};


// similar version of adj, expr and constraints for purely parametric expression
struct paramAdjacency {
    const std::vector<int> indP;    // Indices P
};

class paramExpression {
public:
    explicit paramExpression(paramAdjacency adj) : adj{std::move(adj)} {}

    virtual double eval(const double* p) = 0;

    // returns evalDiff(indicesP) - same sorting as adj.indP
    virtual std::vector<double> evalDiff(const double* p) = 0;

    const paramAdjacency adj;
};


#endif //IPOPT_DO_EXPRESSION_H
