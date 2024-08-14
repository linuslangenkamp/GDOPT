from optimization import *

# choose implicit Euler, since its the most stable, tf > 7s

model = Model("invertedPendulum")

Ms = model.addConst(1)
Mp = model.addConst(0.5)
R = model.addConst(1)
G = model.addConst(-9.81)

s = model.addState(start=0)
v = model.addState(start=0)
phi = model.addState(start=PI-0.001)
omega = model.addState(start=0)

u = model.addControl(lb=-2.5, ub=2.5)

dvdt = u + (sin(-phi)*Mp*R*omega**2 - cos(phi)*sin(-phi)*Mp*G) / (Ms + Mp * sin(-phi)**2) 

model.addDynamic(s, v)
model.addDynamic(v, dvdt)

model.addDynamic(phi, omega)
model.addDynamic(omega, (sin(-phi)*G - cos(phi)*dvdt)/R)

model.addLagrange(sin(phi/2)**2, Objective.MINIMIZE)

model.generate()
