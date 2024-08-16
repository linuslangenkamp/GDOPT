from optimization import *


model = Model("pureParameter")

# by creating purely parametric models, it is possible to solve
# standard NLPs

p1 = model.addP(lb=-1, ub=1)
p2 = model.addP(lb=-1, ub=1)

model.addA(p1**2 + p2**2, lb=1, ub=1)    

model.addM(3*p1 + 2*p2, Objective.MAX)

model.generate()

model.optimize(steps=1, rksteps=1, tf=0,
               flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MA57,
                      "tolerance": 1e-14},
               meshFlags={})

model.printResults()
