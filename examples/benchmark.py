from optimization import *

model = Model("benchmark")

x1 = model.addState(start=0)
x2 = model.addState(start=-1)
x3 = model.addState(start=-sqrt(5))
x4 = model.addState(start=0)

u = model.addInput(lb=-4, ub=10, start=10)

model.addDynamic(x1, x2)
model.addDynamic(x2, -x3 * u + 16 * t - 8)
model.addDynamic(x3, u)
model.addDynamic(x4, x1**2 + x2**2 + 0.0005 * (x2 + 16*t - 8 - 0.1 * x3 * u**2)**2)

model.addMayer(x4, Objective.MINIMIZE)

model.generate()

model.optimize(tf=1, steps=5000, rksteps=1,
               flags={"outputPath": "/tmp",
                      "tolerance": 1e-12,
                      "linearSolver": LinearSolver.MA57,
                      "exportJacobianPath": "/tmp"},
               meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM,
                          "meshIterations": 0})

model.plot()
