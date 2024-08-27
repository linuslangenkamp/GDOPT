from optimization import *
import random

model = Model("pureParameter")

# by creating purely parametric models, it is possible to solve
# standard NLPs: in this case a knapsack problem as a NLP

Items = 15

P = [model.addBinaryParameter(symbol=f"p{i}") for i in range(Items)]

# force binary decision variables

Weights = [random.uniform(0, 2) for i in range(Items)]
Values = [random.uniform(0, 2) for i in range(Items)]

model.addParametric(sum(Weights[i] * P[i] for i in range(Items)), ub=5)

model.addMayer(sum(Values[i] * P[i] for i in range(Items)), Objective.MAXIMIZE)

model.generate()

model.optimize(flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MA57,
                      "tolerance": 1e-14},
               meshFlags={})

model.printResultParameters()
