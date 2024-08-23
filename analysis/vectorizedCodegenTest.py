from sympy import symbols, diff, exp, sin

# Dictionary to store first and second derivatives
DerivativeMap = {}
Derivative2Map = {}

class Function:
    def __init__(self, expr, dependencies, direct_symbols=None, symbol=None):
        """
        expr: The symbolic expression of the function.
        dependencies: A list of other Function objects this function depends on.
        direct_symbols: A set of symbols (e.g., x1, x2) that this function directly depends on.
        symbol: The symbolic representation of the function itself.
        """
        self.expr = expr
        self.dependencies = dependencies  # List of g's that f depends on
        self.direct_symbols = direct_symbols if direct_symbols else set()  # Direct dependencies on x1, x2, etc.
        self.symbol = symbols(symbol) if symbol else None

    def diff(self, x):
        if (self, x) in DerivativeMap:
            return DerivativeMap[(self, x)]

        der = 0

        # Chain rule for dependencies
        for g in self.dependencies:
            dgdx = g.diff(x).expr if (g, x) not in DerivativeMap else DerivativeMap[(g, x)].expr
            dfdg = diff(self.expr, g.symbol)
            der += dfdg * dgdx

        # Direct derivative with respect to x
        if x in self.direct_symbols or x in self.expr.free_symbols:
            der += diff(self.expr, x)

        DerivativeMap[(self, x)] = Function(der, self.dependencies, self.direct_symbols)
        return DerivativeMap[(self, x)]

    def diff2(self, x, y):
        if (self, x, y) in Derivative2Map:
            return Derivative2Map[(self, x, y)]

        der2 = 0

        # Chain rule for second derivatives
        for g in self.dependencies:
            df_dg = self.diff(g.symbol).expr
            dg_dx = g.diff(x).expr
            dg_dy = g.diff(y).expr

            df2_dg2 = diff(df_dg, g.symbol)
            term1 = df2_dg2 * dg_dx * dg_dy

            dg2_dxdy = diff(g.expr, x, y)
            term2 = df_dg * dg2_dxdy

            der2 += term1 + term2

        # Direct second derivative contribution
        if x in self.direct_symbols or x in self.expr.free_symbols:
            der2 += diff(self.diff(x).expr, y)

        Derivative2Map[(self, x, y)] = Function(der2, self.dependencies, self.direct_symbols)
        return Derivative2Map[(self, x, y)]

# Define symbols
x1, x2 = symbols("x1 x2")

# Define functions
g1 = Function(x1**2 + x2**2, [], direct_symbols={x1, x2}, symbol="g1")
g2 = Function(g1.symbol**2 + x1, [g1], direct_symbols={x1}, symbol="g2")
g3 = Function(x1**4, [], direct_symbols={x1}, symbol="g3")
f1 = Function(sin(g2.symbol + g3.symbol) + g2.symbol * x1, [g2, g3], direct_symbols={x1}, symbol="f1")
f2 = Function(exp(g2.symbol + g3.symbol) + x2**2, [g2, g3], direct_symbols={x2}, symbol="f2")

# Compute first derivatives for g's
for g in [g1, g2, g3]:
    for v in [x1, x2]:
        print(f"∂{g.symbol}/∂{v} =", g.diff(v).expr)

# Compute first derivatives for f's
for f in [f1, f2]:
    for v in [x1, x2]:
        print(f"∂{f.symbol}/∂{v} =", f.diff(v).expr)

# Compute second derivatives (Hessian) for g's
for g in [g1, g2, g3]:
    for v1 in [x1, x2]:
        for v2 in [x1, x2]:
            print(f"∂²{g.symbol}/∂{v1}∂{v2} =", g.diff2(v1, v2).expr)

# Compute second derivatives (Hessian) for f's
for f in [f1, f2]:
    for v1 in [x1, x2]:
        for v2 in [x1, x2]:
            print(f"∂²{f.symbol}/∂{v1}∂{v2} =", f.diff2(v1, v2).expr)
