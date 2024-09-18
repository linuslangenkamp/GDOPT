from optimization.expressions import *

# standard guess functions
def guessConstant(const):

    # u(t) = const, can use guess=const directly as well

    return const

def guessLinear(u0, uf):
    
    # u0 = u(0), uf = u(tf)

    return u0 + TIME_SYMBOL * (uf - u0) / FINAL_TIME_SYMBOL

def guessQuadratic(u0, um, uf):

    # u0 = u(0), um = u(tf/2), uf = u(tf)

    return FINAL_TIME_SYMBOL**2 * u0 + FINAL_TIME_SYMBOL * TIME_SYMBOL * (-3 * u0 - uf + 4 * um) + 2 * TIME_SYMBOL**2 * (u0 + uf - 2 * um)

def guessExponential(u0, uf):

    # u0 = u(0), uf = u(tf)

    return u0 * (uf / u0) ** (TIME_SYMBOL / FINAL_TIME_SYMBOL)

def guessPiecewise(*args):

    # *args = (val1, condition1), (val2, condition2), ...

    return Piecewise(*args, (0, True))
