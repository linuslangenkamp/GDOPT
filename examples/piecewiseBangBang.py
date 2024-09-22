from optimization import *

model = Model("piecewiseBangBang")

x1 = model.addState(start=0)
x2 = model.addState(start=0)

u = model.addInput(lb=-100, ub=100)

# x1'' = u
model.addDynamic(x1, x2)
model.addDynamic(x2, u)

model.addPath(x2 * u * piecewise((1, t < 0.25), (0, t >= 0.25)), lb=-30, ub=30) # constraint only has to hold for time < 0.25

model.addFinal(x2, eq=0)

model.addMayer(x1, Objective.MAXIMIZE)

model.generate()

# optimizer attributes can be set directly as well
model.meshAlgorithm = MeshAlgorithm.L2_BOUNDARY_NORM
model.meshIterations = 3
model.outputFilePath = "/tmp"

model.optimize(tf=0.5, steps=150, rksteps=3)

model.plot(dots=Dots.BASE)
model.plotMeshRefinement()
