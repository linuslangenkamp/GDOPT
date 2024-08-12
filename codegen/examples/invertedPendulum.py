from optimization import *

# variant of https://apmonitor.com/do/index.php/Main/ControlTypes
# choose implicit Euler, since its the most stable

model = Model("invertedPendulum")

m1 = model.addConst(10)
m2 = model.addConst(1)

eps = model.addConst(m2 / (m1 + m2))

y = model.addState(start=-1)
v = model.addState(start=0)
theta = model.addState(start=0)
q = model.addState(start=0)

u = model.addControl(lb=-1, ub=1)

model.addDynamic(y, v)
model.addDynamic(v, -eps*theta + u)

model.addDynamic(theta, q)
model.addDynamic(q, theta - u)

# stationary final state
model.addFinal(y, lb=0, ub=0)
model.addFinal(v, lb=0, ub=0)
model.addFinal(theta, lb=0, ub=0)
model.addFinal(q, lb=0, ub=0)

# minimize action u**2 or angle theta**2
model.addLagrange(theta**2, Objective.MINIMIZE)

model.generate()
