from optimization import *


model = Model("hypersensitive")

x = model.addX(start=1.)
u = model.addU()

model.addF(x, -x**3 + u)
model.addR(1.5 - x, eq=0)
model.addL(0.5 * (x**2 + u**2))

model.generate()

model.optimize(tf=10000, steps=200, rksteps=7,
               flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MA57},
               meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM,
                          "meshIterations": 10,
                          "meshLevel": 0})

model.plot(interval=[9980, 10000], dots=True)