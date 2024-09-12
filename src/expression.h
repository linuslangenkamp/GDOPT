#ifndef IPOPT_DO_EXPRESSION_H
#define IPOPT_DO_EXPRESSION_H

#include <array>
#include <utility>
#include <vector>

struct Adjacency {
    const std::vector<int> indX; // Indices X
    const std::vector<int> indU; // Indices U
    const std::vector<int> indP; // Indices P
};

// firstly diff w.r.t. 1st index, after that w.r.t. 2nd index
// only need lower triangular hessian. thus 1st idx >= 2nd idx for diff xx, diff uu, diff pp
struct AdjacencyDiff {
    const std::vector<std::tuple<int, int>> indXX; // Indices X, X :: x1 >= x2 must hold
    const std::vector<std::tuple<int, int>> indUX; // Indices U, X
    const std::vector<std::tuple<int, int>> indUU; // Indices U, U :: u1 >= u2 must hold
    const std::vector<std::tuple<int, int>> indPX; // Indices P, X
    const std::vector<std::tuple<int, int>> indPU; // Indices P, U
    const std::vector<std::tuple<int, int>> indPP; // Indices P, P :: p1 >= p2 must hold
};

class Expression {
  public:
    explicit Expression(Adjacency adj, AdjacencyDiff adjDiff) : adj{std::move(adj)}, adjDiff{std::move(adjDiff)} {
    }

    virtual double eval(const double* x, const double* u, const double* p, double t) = 0;

    // returns {evalDiff(indX), evalDiff(indU), evalDiff(indP)} - same sorting as adj!!
    virtual std::array<std::vector<double>, 3> evalDiff(const double* x, const double* u, const double* p, double t) = 0;

    // return {evalDiff2(indXX), evalDiff2(indUX), evalDiff2(idxUU), evalDiff2(indPX), evalDiff2(indPU), evalDiff2(idxPP)}
    //  - same sorting as adj!!
    virtual std::array<std::vector<double>, 6> evalDiff2(const double* x, const double* u, const double* p, double t) = 0;

    const Adjacency adj;
    const AdjacencyDiff adjDiff;
};

// similar version of adj, expr and constraints for purely parametric expression
struct ParamAdjacency {
    const std::vector<int> indP; // Indices P
};

struct ParamAdjacencyDiff {
    const std::vector<std::tuple<int, int>> indPP; // Indices P, P
};

class ParamExpression {
  public:
    explicit ParamExpression(ParamAdjacency adj, ParamAdjacencyDiff adjDiff) : adj{std::move(adj)}, adjDiff{std::move(adjDiff)} {
    }

    virtual double eval(const double* p) = 0;

    // returns evalDiff(indicesP) - same sorting as adj!!
    virtual std::vector<double> evalDiff(const double* p) = 0;

    // return evalDiff2(indicesPP) - same sorting as adj!!
    virtual std::vector<double> evalDiff2(const double* p) = 0;

    const ParamAdjacency adj;
    const ParamAdjacencyDiff adjDiff;
};

#endif // IPOPT_DO_EXPRESSION_H
