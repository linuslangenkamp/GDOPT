### testing detection criteria for mesh refinement strategies

import pandas as pd
import matplotlib.pyplot as plt
import sympy as sy


path = "/mnt/c/Users/Linus/Desktop/Studium/Master/Masterarbeit/VariableData/trivialBangBang"
model = "trivialBangBang"
it = 0
specifCol = 'u0'
interval = [0, 1]
intervals = 37
steps = 3
c = [0, 0.155051025721, 0.644948974278, 1]

	
df = pd.read_csv(path + "/" + model + str(it) + ".csv" , sep=",")
print(df.head())
u0 =  list(df[specifCol])
time = list(df["time"])
# hack 0th element to be = to 1st
u0.insert(0, u0[0])
time.insert(0, 0)
timeBasePoints = [time[i] for i in range(intervals * steps) if i%3 == 0]
print(timeBasePoints)
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


plt.figure(figsize=(10, 6))
zdashL2 = []
zdashAbs = []
z2dashL2 = []
z2dashAbs = []
# .evalf(subs={t: t0}
for i in range(intervals):
	zdashAbs.append(sy.integrate(abs(sy.diff(p(i), t)), (t, 0, 1)))
	zdashL2.append(abs(sy.integrate((sy.diff(p(i), t))**2, (t, 0, 1)))**0.5)
	#z2dashAbs.append(sy.integrate(abs(sy.diff(p(i), t, t)), (t, 0, 1)))
	#z2dashL2.append(sy.integrate((sy.diff(p(i), t, t))**2, (t, 0, 1))**0.5)
	
print(zdashL2)
plt.plot(timeBasePoints, zdashAbs, label="absolute")
plt.plot(timeBasePoints, zdashL2, label="L2")
plt.xlabel('Time base points')
plt.ylabel("Z Value norm")
plt.xlim(interval[0], interval[1])
plt.legend()
plt.grid(True)
plt.show()
