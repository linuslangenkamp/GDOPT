from optimization import *
import random

model = Model("pureParameter")

# by creating purely parametric models, it is possible to solve
# standard NLPs: in this case a knapsack problem

Items = 15

P = [model.addParameter(symbol=f"p{i}", lb=0, ub=1, start=0.5) for i in range(Items)]

# force binary decision variables
for p in P:
    model.addParametric(p * (p - 1), eq=0)

Weights = [random.uniform(0, 2) for i in range(Items)]
Values = [random.uniform(0, 2) for i in range(Items)]

model.addParametric(sum(Weights[i] * P[i] for i in range(Items)), ub=5)

model.addMayer(sum(Values[i] * P[i] for i in range(Items)), Objective.MAX)

model.generate()

model.optimize(flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MA57,
                      "tolerance": 1e-14},
               meshFlags={})

model.printResultParameters()
