#include <string>
#include "problem.h"

Problem::Problem(int sizeX, int sizeU, int sizeP,
                 std::vector<double> x0, std::vector<double> lbX, std::vector<double> ubX,
                 std::function<std::vector<double>(double)> uInitialGuess, std::vector<double> lbU, std::vector<double> ubU,
                 std::vector<double> pInitialGuess, std::vector<double> lbP, std::vector<double> ubP,
                 std::unique_ptr<Expression> M, std::unique_ptr<Expression> L,
                 std::vector<std::unique_ptr<Expression>> F,
                 std::vector<std::unique_ptr<Constraint>> G,
                 std::vector<std::unique_ptr<Constraint>> R,
                 std::vector<std::unique_ptr<ParamConstraint>> A,
                 std::string name)
        : sizeX(sizeX), sizeU(sizeU), sizeP(sizeP),
          x0(std::move(x0)), lbX(std::move(lbX)), ubX(std::move(ubX)),
          uInitialGuess(std::move(uInitialGuess)), lbU(std::move(lbU)), ubU(std::move(ubU)),
          pInitialGuess(std::move(pInitialGuess)), lbP(std::move(lbP)), ubP(std::move(ubP)),
          M(std::move(M)), L(std::move(L)), F(std::move(F)), G(std::move(G)), R(std::move(R)), A(std::move(A)),
          name(std::move(name)) {
}
