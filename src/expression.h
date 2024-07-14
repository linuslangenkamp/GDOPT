//
// Created by Linus on 10.07.2024.
//

#ifndef IPOPT_DO_EXPRESSION_H
#define IPOPT_DO_EXPRESSION_H

#include <utility>
#include <vector>
#include <array>

struct Adjacency {
    const std::vector<int> indicesX;
    const std::vector<int> indicesU;
    const std::vector<int> indicesP;
};

class Expression {
public:
    explicit Expression(Adjacency adj) : adj{std::move(adj)} {}

    virtual double eval(double *x, double *u, double *p, double t) = 0;

    virtual std::array<std::vector<double>, 3> evalDiff(double *x, double *u, double *p, double t) = 0;
protected:
    const Adjacency adj;
};

#endif //IPOPT_DO_EXPRESSION_H
