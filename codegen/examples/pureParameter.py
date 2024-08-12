from optimization import *


model = Model("pureParameter")

# by creating purely parametric models, it is possible to solve
# standard NLPs

p1 = model.addP(lb=-1, ub=1)
p2 = model.addP(lb=-1, ub=1)

model.addG(p1**2 + p2**2, lb=1, ub=1)    

model.addM(3*p1 + 2*p2, Objective.MAX)

model.generate()
