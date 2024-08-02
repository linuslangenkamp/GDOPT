### testing detection criteria for mesh refinement strategies

import pandas as pd
import matplotlib.pyplot as plt
import sympy as sy
from mpmath import *

path = "/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData/batchReactorRefinement"
model = "BatchReactor"
it = 0
specifCol = 'u0'
interval = [0, 1]
intervals = 100
steps = 3

x = sy.symbols('x')
function = x**(steps-1) * (x-1)**steps
der = sy.expand(sy.diff(function, x, steps-1))
pDer = sy.Poly(der, x)
coefficients = [int(pDer.coeff_monomial(x**(steps-k))) for k in range(steps+1)]
c = [[mp.mpf(1)] * steps, polyroots(coefficients)][1]
c = [float(elem) for elem in c]
c.insert(0, 0)
df = pd.read_csv(path + "/" + model + str(it) + ".csv" , sep=",")
print(df.head())
u0 =  list(df[specifCol])
time = list(df["time"])
u0.insert(0, u0[0])
time.insert(0, 0)
timeBasePoints = [time[i] for i in range(intervals * steps) if i%steps == 0]
t = sy.symbols('t')

def basis(j):
	expr = 1
	for k in range(len(c)):
		if k != j:
			expr *= (t - c[k]) / (c[j] - c[k])
	return expr
	
def p(i):
	poly = 0
	for k in range(len(c)):
		poly += u0[steps * i + k] * basis(k)
	return poly

def relativeRescale(iterable):
	m = max(iterable)
	return [elem / m for elem in iterable]

plt.figure(figsize=(10, 6))
zdashL2 = []
zdashAbs = []
z2dashL2 = []
z2dashAbs = []
# .evalf(subs={t: t0}
p1, p1Diff2 = 0, 0
differencesDiff, differencesDiff2 = [], []

for i in range(intervals):
	#zdashAbs.append(sy.integrate(abs(sy.diff(p(i), t)), (t, 0, 1)))
	zdashL2.append(abs(sy.integrate((sy.diff(p(i), t))**2, (t, 0, 1)))**0.5)
	#z2dashAbs.append(sy.integrate(abs(sy.diff(p(i), t, t)), (t, 0, 1)))
	z2dashL2.append(sy.integrate((sy.diff(p(i), t, t))**2, (t, 0, 1))**0.5)
	oldP1 = p1
	p1 = sy.diff(p(i), t).evalf(subs={t: 1})
	differencesDiff.append(abs(sy.diff(p(i), t).evalf(subs={t: 0}) - oldP1))
	
	oldP1Diff2 = p1Diff2
	p1Diff2 = sy.diff(p(i), t, t).evalf(subs={t: 1})
	differencesDiff2.append(abs(sy.diff(p(i), t, t).evalf(subs={t: 0}) - oldP1Diff2))

print(timeBasePoints)
#plt.plot(timeBasePoints, zdashAbs, label="absolute")
plt.plot(timeBasePoints, (zdashL2), label="L2")
plt.plot(timeBasePoints, (z2dashL2), label="L2'")
plt.plot(timeBasePoints, (differencesDiff), label="Δp'")
plt.plot(timeBasePoints, (differencesDiff2), label="Δp''")
plt.xlabel('Time base points')
plt.ylabel("Z Value norm")
plt.xlim(interval[0], interval[1])
plt.legend()
plt.grid(True)
plt.show()
