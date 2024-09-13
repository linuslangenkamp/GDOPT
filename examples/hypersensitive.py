from optimization import *

model = Model("hypersensitive")

x = model.addState(start=1)

u = model.addInput()

model.addDynamic(x, -(x**3) + u)

model.addFinal(1.5 - x, eq=0)

model.addLagrange(0.5 * (x**2 + u**2))

model.generate()

model.optimize(
    tf=10000,
    steps=30,
    rksteps=9,
    flags={"outputPath": "/tmp", "linearSolver": LinearSolver.MA57},
    meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM, "meshIterations": 20, "meshLevel": 0},
)

model.plotInputs(interval=[9980, 10000], dots=Dots.ALL)
model.plotMeshRefinement(interval=[9980, 10000])
