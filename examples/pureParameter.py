from optimization import *


model = Model("pureParameter")

# by creating purely parametric models, it is possible to solve
# standard NLPs

p1 = model.addParameter(symbol="p1", lb=-1, ub=1)
p2 = model.addParameter(symbol="p2", lb=-1, ub=1)

model.addParametric(p1**2 + p2**2, lb=1, ub=1)

model.addMayer(3*p1 + 2*p2, Objective.MAX)

model.generate()

model.optimize(flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MA57,
                      "tolerance": 1e-14},
               meshFlags={})

model.printResultParameters()
