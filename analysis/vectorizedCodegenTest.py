from mpmath import extend
from sympy import *


DerivativeMap = {}
AllSymbols = []

class Function:

    def __init__(self, expr, adjPure, adjComp, symbol=None):
        self.expr = expr
        self.adjPure = adjPure
        self.adjComp = adjComp
        self.adjSymbols = adjPure.union(v.symbol for v in adjComp)
        if symbol is not None:
            self.symbol = symbols(symbol)
            AllSymbols.append(self.symbol)
        else:
            self.symbol = symbol

    def diff(self, x):
        der = 0

        for g in self.adjComp:
            if (g, x) not in DerivativeMap:
                # TODO: actually not g.adjPure, g.adjComp, but only the vars that are contained in dg/dx
                DerivativeMap[(g, x)] = Function(g.diff(x), g.adjPure, g.adjComp, symbol=f"d{g.symbol}d{x.name}")
            dfdg = diff(self.expr, g.symbol)
            der += dfdg * DerivativeMap[(g, x)].symbol

        if x in self.adjPure:
            der += diff(self.expr, x)

        DerivativeMap[(self, x)] = Function(der, self.adjPure, self.adjComp, symbol=f"d{self.symbol}d{x.name}")
        return der


# -> dg2dx1*cos(g2 + g3) + dg3dx1*cos(g2 + g3) + 1 do this with maps that capture every subexpression
"""
proposed approach:

* for all given functions f(x1, ..., xn) in F:
    * find all subexpression g_i in f
    * add f~(x1, ..., xn, g1, ..., gm) = f.subst(x,g) to F~

* for all given functions g in G:
    * calculate the derivative of g w.r.t to the adjacent variables (x*, g*) for all i 
      and add them to a map (g, var): der(g, var)

* for all given functions f~ in F~:
    * calculate the derivative of f w.r.t to the adjacent variables (x*, g*)
    
* get cse of the set of all derivatives

* return
"""

x1, x2 = symbols("x1 x2")
variables = [x1, x2]

g1 = Function(x1**2 + x2**2, {x1, x2}, set(), symbol="g1")
g2 = Function(g1.symbol**2, set(), {g1}, symbol="g2")
g3 = Function(x1**4, {x1}, set(), symbol="g3")
f1 = Function(sin(g2.symbol + g3.symbol) + x1, {x1}, {g2, g3})
f2 = Function(sin(g2.symbol + g3.symbol), set(), {g2, g3})

functionsG = [g1, g2, g3]
for g in functionsG:
    for v in g.adjSymbols:
        g.diff(v)

functionsF = [f1, f2]
derivatives = []

for f in functionsF:
    for v in variables: # TODO: here only adj vars -> need adj info from 1st cse call where g1... are generated
        derivatives.append(f.diff(v))

partialExpression = numbered_symbols(prefix='s')
print(cse(derivatives, partialExpression))