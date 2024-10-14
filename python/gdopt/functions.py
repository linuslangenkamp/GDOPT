###############################################################################
# GDOPT - General Dynamic Optimizer
# Copyright (C) 2024  Linus Langenkamp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

from .expressions import *


# custom piecewise function
def piecewise(*args):

    # *args = (val1, condition1), (val2, condition2), ...

    return Piecewise(*args, (0, True))


# standard guess functions
def guessConstant(const):

    # u(t) = const, can use guess=const directly as well

    return const


def guessLinear(u0, uf):

    # u0 = u(0), uf = u(tf)

    return u0 + TIME_SYMBOL * (uf - u0) / FINAL_TIME_SYMBOL


def guessQuadratic(u0, um, uf):

    # u0 = u(0), um = u(tf/2), uf = u(tf)

    return (
        FINAL_TIME_SYMBOL**2 * u0
        + FINAL_TIME_SYMBOL * TIME_SYMBOL * (-3 * u0 - uf + 4 * um)
        + 2 * TIME_SYMBOL**2 * (u0 + uf - 2 * um)
    )


def guessExponential(u0, uf):

    # u0 = u(0), uf = u(tf)

    return u0 * (uf / u0) ** (TIME_SYMBOL / FINAL_TIME_SYMBOL)


def guessPiecewise(*args):

    # *args = (val1, condition1), (val2, condition2), ...

    return Piecewise(*args, (0, True))


def Abs(expr):

    # custom Abs as piecewise function because symengine.Abs() doesn't work in codegen?!

    return Piecewise((-expr, expr < 0), (expr, True))


def backendReturnCode(code):
    code = int(code)
    IPOPT_RETURN_CODES = {
        0: "Optimal_Solution_Found",
        1: "Solved_To_Acceptable_Level",
        2: "Infeasible_Problem_Detected",
        3: "Search_Direction_Becomes_Too_Small",
        4: "Diverging_Iterates",
        5: "User_Requested_Stop",
        6: "Feasible_Point_Found",
        7: "Maximum_Iterations_Exceeded",
        8: "Restoration_Failed",
        9: "Error_In_Step_Computation",
        10: "Maximum_CpuTime_Exceeded",
        11: "Not_Enough_Degrees_Of_Freedom",
        12: "Invalid_Problem_Definition",
        13: "Invalid_Option",
        14: "Invalid_Number_Detected",
        -1: "Internal_Error",
    }
    if code in IPOPT_RETURN_CODES:
        return IPOPT_RETURN_CODES[code]
    else:
        return f"Unknown return code: {code}"
