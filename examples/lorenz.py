from gdopt import *

m = Model("lorenzAttractor")

rho = m.addRuntimeParameter(default=28, symbol="rho")
sigma = m.addRuntimeParameter(default=10, symbol="sigma")
beta = m.addRuntimeParameter(default=8/3, symbol="beta")

x = m.addState(start=1, nominal=15)
y = m.addState(start=-1, nominal=15)
z = m.addState(start=0, nominal=15)

u = m.addControl(lb=-1, ub=1, nominal=0.5)

m.addDynamic(x, u + sigma * (y - x), nominal=15)
m.addDynamic(y, x * (rho - z) - y - u / 2, nominal=15)
m.addDynamic(z, x * y - beta * z - u / 2, nominal=15)

m.addMayer(x, nominal=15)

m.generate()

m.optimize(tf=15, steps=200, rksteps=5,
    flags={"tolerance": 1e-14, "linearSolver": LinearSolver.MA57},
    meshFlags={"algorithm": MeshAlgorithm.L2_BOUNDARY_NORM, "iterations": 5},)

m.plot(dots=Dots.ALL)
m.plotInputsAndRefinement()