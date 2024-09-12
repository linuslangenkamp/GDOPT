from optimization import *

model = Model("vanDerPol")

x1 = model.addState(start=0)
x2 = model.addState(start=1)

u = model.addInput(ub=0.8)

rp = model.addRuntimeParameter(default=1, symbol="RP")

model.addDynamic(x1, (1 - x2**2) * x1 - x2 + u)
model.addDynamic(x2, x1)

model.addLagrange(x1**2 + x2**2 + rp * u**2)

model.generate()

model.optimize(
    tf=10,
    steps=20,
    rksteps=3,
    flags={"outputPath": "/tmp", "linearSolver": LinearSolver.MA57},
    meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM, "meshIterations": 10, "meshLevel": -0.5, "meshCTol": 0.5},
)

model.parametricPlot(x1, x2, dots=Dots.ALL)

model.setValue(rp, 0.1)

model.optimize(resimulate=True)
model.plot(dots=Dots.ALL)
