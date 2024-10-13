from gdopt import *

m = Model("lorenzAttractor")

rho = m.addRuntimeParameter(default=28, symbol="rho")
sigma = m.addRuntimeParameter(default=10, symbol="sigma")
beta = m.addRuntimeParameter(default=8/3, symbol="beta")

x = m.addState(start=1)
y = m.addState(start=-1)
z = m.addState(start=0)

u = m.addControl(lb=-1, ub=1)

m.addDynamic(x, u + sigma * (y - x))
m.addDynamic(y, x * (rho - z) - y - u / 2)
m.addDynamic(z, x * y - beta * z - u / 2)

m.addMayer(x)

m.generate()

m.optimize(tf=15, steps=75, rksteps=9)

m.plot(dots=Dots.BASE)
