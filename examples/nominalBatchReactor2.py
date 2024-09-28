from gdopt import *

model = Model("nominalBatchReactor2")

x1 = model.addState(symbol="Reactant", start=1e10, nominal=1e10)
x2 = model.addState(symbol="Product", start=0)

u = model.addInput(symbol="u", lb=0, ub=5, guess=1)

model.addDynamic(x1, -(u + u**2 / 2) * x1, nominal=1e10)
model.addDynamic(x2, u * x1 / 1e10)

model.addMayer(x2, Objective.MAXIMIZE)

model.hasLinearObjective()

model.generate()

model.optimize(
    tf=1,
    steps=500,
    rksteps=3,
    flags={
        "outputPath": "/tmp",
        "linearSolver": LinearSolver.MA57,
        "initVars": InitVars.SOLVE,
        "tolerance": 1e-12,
        "exportJacobianPath": "/tmp",
    },
    meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM, "meshIterations": 0},
)

model.plotInputs(dots=Dots.BASE)
