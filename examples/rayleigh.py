from gdopt import *

model = Model("rayleigh")

x = model.addState(start=-5)
y = model.addState(start=-5)

u = model.addInput()

model.addDynamic(x, y)
model.addDynamic(y, -x + y * (1.4 - 0.14 * y**2) + 4 * u)

model.addLagrange(x**2 + u**2)

model.generate()

model.optimize(
    tf=4.5,
    steps=50,
    rksteps=3,
    flags={"linearSolver": LinearSolver.MA57},
    meshFlags={
        "algorithm": MeshAlgorithm.L2_BOUNDARY_NORM,
        "refinementMethod": RefinementMethod.LINEAR_SPLINE,
        "iterations": 5,
    },
)

model.plot(dots=Dots.BASE)