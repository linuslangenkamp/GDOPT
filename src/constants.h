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

#ifndef IPOPT_DO_CONSTANTS_H
#define IPOPT_DO_CONSTANTS_H

#include <limits>

const double PLUS_INFINITY = std::numeric_limits<double>::infinity();
const double MINUS_INFINITY = -std::numeric_limits<double>::infinity();

#endif  // IPOPT_DO_CONSTANTS_H
