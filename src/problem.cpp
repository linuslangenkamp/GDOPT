#include <string>
#include "problem.h"

Problem::Problem(int sizeX, int sizeU, int sizeP,
                 std::vector<double> x0, std::vector<double> lbX, std::vector<double> ubX,
                 std::vector<double> lbU, std::vector<double> ubU,
                 std::vector<double> lbP, std::vector<double> ubP,
                 std::unique_ptr<Expression> M, std::unique_ptr<Expression> L,
                 std::vector<std::unique_ptr<Expression>> F,
                 std::vector<std::unique_ptr<Constraint>> G,
                 std::vector<std::unique_ptr<Constraint>> R,
                 std::vector<std::unique_ptr<ParamConstraint>> A,
                 std::string name)
        : sizeX(sizeX), sizeU(sizeU), sizeP(sizeP),
          x0(std::move(x0)), lbX(std::move(lbX)), ubX(std::move(ubX)),
          lbU(std::move(lbU)), ubU(std::move(ubU)),
          lbP(std::move(lbP)), ubP(std::move(ubP)),
          M(std::move(M)), L(std::move(L)), F(std::move(F)), G(std::move(G)), R(std::move(R)), A(std::move(A)),
          name(std::move(name)) {
}
