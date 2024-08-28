from sympy import *
from sympy import Integer

x, y, t = symbols("x y t")

def calcL2(f, var):
    return integrate(diff(f, var)**2, (var, 0, 1))

# hier annahme, dass das polynom auf dem interval [0, 1] == [t_i, t_i + h] sich nicht verändert
# bei sprüngen wird das integral eh nie verschwinden. Überlege Sprungfunktion auf kleiner werdenden interval mit e.g.
# kubischer interpolation, steigung wird nur noch höher

f = -469.993*x**2+531.287*x-143.155

# [t_i + j/2**k * h, t_i + (j+1)/2**k * h] -> [0, 1] ohne einschränkung
# ableitung des transformierten polynom p: d/dt p(t)
# berechne (d/dt p(t)) ** 2
# berechne int_0^1 (d/dt p(t)) ** 2 dt

# hier einfach t_i = 0, h = 1

for k in range(5):
    step = 2**k
    for j in range(step):
        inner = t / Integer(step) + Integer(j) / Integer(step)
        fbar = f.subs(x, inner)
        L2 = calcL2(fbar, t)
        print(f"L2({inner}) = {float(L2)}") #  sum L2 * step = L2 iteration 0
        # jede iteration halbiert mindestens den L2 score, falls das polynom sich nicht mehr ändert!!
        # if the solution is similar enough the convergence is extremely fast
        # if not its not, assume jump or p(x) on [t_i, t_i + h] and in next iteration its p(x) rescaled to [t_i, t_i + h/2]
        # and p(t_i + h/2) const on [t_i + h/2, t_i + h]
    print()

