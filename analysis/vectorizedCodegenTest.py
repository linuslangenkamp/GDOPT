#### really experimental
from sympy import *

# Dictionary to store first and second derivatives
DerivativeMap = {}
Derivative2Map = {}

class Function:
    def __init__(self, expr, dependencies, direct_symbols=None, symbol=None):
        self.expr = expr
        self.dependencies = dependencies  # List of g's that f depends on
        self.direct_symbols = direct_symbols if direct_symbols else set()
        if isinstance(symbol, Symbol):
            self.symbol = symbol
        else:
            self.symbol = symbols(symbol) if symbol else None

    def diff(self, x):
        if (self, x) in DerivativeMap:
            return DerivativeMap[(self, x)]

        der = 0

        # Chain rule for dependencies
        for g in self.dependencies:
            dgdx = g.diff(x)
            dgdx_expr = dgdx.expr if isinstance(dgdx, Function) else dgdx  # Handle constants
            print(dgdx_expr)

            dfdg = diff(self.expr, g.symbol) if isinstance(g, Function) else diff(self.expr, g)
            der += dfdg * dgdx_expr

        # Direct derivative with respect to x
        if x in self.direct_symbols or x in self.expr.free_symbols:
            der += diff(self.expr, x)

        # Handle constants directly
        if isinstance(der, (Float, int, Rational, Symbol)):
            DerivativeMap[(self, x)] = der
            return der

        # Otherwise, wrap it in a Function object
        DerivativeMap[(self, x)] = Function(der, self.dependencies, self.direct_symbols)
        return DerivativeMap[(self, x)]

    def diff2(self, x, y):
        if (self, x, y) in Derivative2Map:
            return Derivative2Map[(self, x, y)]

        der2 = 0

        # Chain rule for second derivatives
        for g in self.dependencies:
            df_dg = self.diff(g.symbol).expr if isinstance(self.diff(g.symbol), Function) else self.diff(g.symbol)
            dg_dx = g.diff(x).expr if isinstance(g.diff(x), Function) else g.diff(x)
            dg_dy = g.diff(y).expr if isinstance(g.diff(y), Function) else g.diff(y)

            df2_dg2 = diff(df_dg, g.symbol) if isinstance(g, Function) else diff(df_dg, g)
            term1 = df2_dg2 * dg_dx * dg_dy

            if isinstance(x, Function) and isinstance(y, Function):
                dg2_dxdy = diff(g.expr, x.symbol, y.symbol)
            elif isinstance(x, Function):
                dg2_dxdy = diff(df_dg, x.symbol, y)
            elif isinstance(y, Function):
                dg2_dxdy = diff(df_dg, x, y.symbol)
            else:
                dg2_dxdy = diff(df_dg, x, y)

            term2 = df_dg * dg2_dxdy

            der2 += term1 + term2

        # Direct second derivative contribution
        if x in self.direct_symbols or x in self.expr.free_symbols:
            direct_der2 = diff(self.diff(x).expr if isinstance(self.diff(x), Function) else self.diff(x), y.symbol if isinstance(y, Function) else y)
            der2 += direct_der2

        # Handle constants directly
        if isinstance(der2, (Float, int, Rational, Symbol)):
            Derivative2Map[(self, x, y)] = der2
            return der2

        # Otherwise, wrap it in a Function object
        Derivative2Map[(self, x, y)] = Function(der2, self.dependencies, self.direct_symbols)
        return Derivative2Map[(self, x, y)]

# Define symbols
x1, x2, x3, x4, x0, u0 = symbols("x1 x2 x3 x4 x0 u0")
Vars = [x1, x2, x3, x4, x0, u0]
# Define functions
g1 = Function(x1**2 + x2**2, [], direct_symbols={x1, x2}, symbol="g1")
g2 = Function(g1.symbol**2 + x1, [g1], direct_symbols={x1}, symbol="g2")
g3 = Function(x1**4, [], direct_symbols={x1}, symbol="g3")
f1 = Function(sin(g2.symbol + g3.symbol) + g2.symbol * x1, [g2, g3], direct_symbols={x1}, symbol="f1")
f2 = Function(exp(g2.symbol + g3.symbol) + x2**2, [g2, g3], direct_symbols={x2}, symbol="f2")

expr = -0.247230109968751*x3**2 + 5.05613423123327e-5*(16973.4306615008*x2*(x2/x1)**0.283877349159248*(1 - 0.791745582846156*(1/x2)**0.214714714714715)*(4305.15564366859*u0*x0/(pi*(0.0825*u0*x0/pi + 1.44204839461401*x0*x0/pi)) + 930.653957571921)*sqrt((1/x0)**0.785285285285285 - 0.889800866961904*(1/x2)**0.892642642642643)/(sqrt((x2/x1)**0.283877349159248*(4305.15564366859*u0*x0/(pi*(0.0825*u0*x0/pi + 1.44204839461401*x0*x1/pi)) + 930.653957571921)/(u0*x0/(pi*(0.0825*u0*x0/pi + 1.44204839461401*x0*x1/pi)) + 0.328427455094938)**0.283877349159248)*(u0*x0/(pi*(0.0825*u0*x0/pi + 1.44204839461401*x0*x1/pi)) + 0.328427455094938)**0.283877349159248) - 602154.981787439*(1.21364877581381*x1**0.283877349159248 - 1)*sqrt(-x1**2/(0.381107883006746*x3**2 + 1)**7.04529616724739 + 0.255587584313281))/x3
partialExpression = numbered_symbols(prefix='s')
subst, result = cse(expr, partialExpression)

G = []
for sub in subst:
    dependencies = []
    direct_symbols = set()
    for e in sub[1].free_symbols:
        if "s" in e.name:
            whichS = int(e.name[1:])
            dependencies.append(G[whichS])
            print(whichS)
            Vars.append(e)
        else:
            direct_symbols.add(e)
    G.append(Function(sub[1], dependencies, direct_symbols, symbol=sub[0]))

F = []
dependencies = []
direct_symbols = set()
for e in result[0].free_symbols:
    if "s" in e.name:
        dependencies.append(e)
    else:
        direct_symbols.add(e)

F.append(Function(result[0], dependencies, direct_symbols, symbol="expr"))
print(F[0].direct_symbols)
print(subst)

# Compute first derivatives for g's
for g in G:
    for v in g.dependencies + list(g.direct_symbols):
        print(f"∂{g.symbol}/∂{v} =", g.diff(v))

# Compute first derivatives for f's
for f in F:
    print(f.expr)
    for v in f.dependencies + list(f.direct_symbols):
        print(f"∂{f.symbol}/∂{v} =", f.diff(v).expr)

# Compute second derivatives (Hessian) for g's
for g in G:
    for v1 in g.dependencies + list(g.direct_symbols):
        for v2 in g.dependencies + list(g.direct_symbols):
            derivative = g.diff2(v1, v2)
            if isinstance(derivative, Function):
                print(f"∂²{g.symbol}/∂{v1}∂{v2} =",derivative.expr )

# Compute second derivatives (Hessian) for f's
for f in F:
    for v1 in f.dependencies + list(f.direct_symbols):
        for v2 in f.dependencies + list(f.direct_symbols):
            print(f"∂²{f.symbol}/∂{v1}∂{v2} =", f.diff2(v1, v2))
