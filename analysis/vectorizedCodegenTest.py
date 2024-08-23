from sympy import *


DerivativeMap = {}


class Function():

    def __init__(self, expr, adjPure, adjComp, symbol=None):
        self.expr = expr
        self.adjPure = adjPure
        self.adjComp = adjComp
        if symbol is not None:
            self.symbol = symbols(symbol)
        else:
            self.symbol = symbol

    def diff(self, x):
        der = 0
        for g in self.adjComp:
            if (g, x) not in DerivativeMap:
                DerivativeMap[(g, x)] = g.diff(x)
            dfdg = diff(self.expr, g.symbol)
            der += dfdg * DerivativeMap[(g, x)]
        if x in self.adjPure:
            der += diff(self.expr, x)
        return der


x1, x2 = symbols("x1 x2")
g1 = Function(x1**2 + x2**2, {x1, x2}, set(), symbol="g1")
g2 = Function(g1.symbol**2, set(), {g1}, symbol="g2")
f = Function(sin(g2.symbol) + x1, {x1}, {g1, g2})

print(f.diff(x1))
print(DerivativeMap)