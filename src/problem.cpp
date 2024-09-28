/**
 * GDOPT - General Dynamic Optimization Problem Optimizer
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

#include "problem.h"

#include <string>

Problem::Problem(int sizeX, int sizeU, int sizeP, std::vector<double> x0, std::vector<double> lbX, std::vector<double> ubX,
                 std::function<std::vector<double>(double)> uInitialGuess, std::vector<double> lbU, std::vector<double> ubU, std::vector<double> pInitialGuess,
                 std::vector<double> lbP, std::vector<double> ubP, std::unique_ptr<Expression> M, std::unique_ptr<Expression> L,
                 std::vector<std::unique_ptr<Expression>> F, std::vector<std::unique_ptr<Constraint>> G, std::vector<std::unique_ptr<Constraint>> R,
                 std::vector<std::unique_ptr<ParamConstraint>> A, std::string name)
    : sizeX(sizeX),
      sizeU(sizeU),
      sizeP(sizeP),
      x0(std::move(x0)),
      lbX(std::move(lbX)),
      ubX(std::move(ubX)),
      uInitialGuess(std::move(uInitialGuess)),
      lbU(std::move(lbU)),
      ubU(std::move(ubU)),
      pInitialGuess(std::move(pInitialGuess)),
      lbP(std::move(lbP)),
      ubP(std::move(ubP)),
      M(std::move(M)),
      L(std::move(L)),
      F(std::move(F)),
      G(std::move(G)),
      R(std::move(R)),
      A(std::move(A)),
      name(std::move(name)) {
}
