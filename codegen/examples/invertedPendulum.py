from optimization import *

# variant of https://rmc.dlr.de/rm/de/staff/extcms/images/rmc/users/zimm_di/lecture2022/Lecture12.pdf
# choose implicit Euler, since its the most stable, ideal tf>=5s

model = Model("invertedPendulum")
Ms = model.addConst(10)
Mp = model.addConst(1)
R = model.addConst(0.5)
G = model.addConst(-9.81)

s = model.addState(start=0)
v = model.addState(start=0)
phi = model.addState(start=0)
omega = model.addState(start=0)

u = model.addControl(lb=-0.2, ub=0.2)

dvdt = (sin(phi)*Mp*R*omega**2 - cos(phi)*sin(phi)*Mp*G) / (Ms + Mp * sin(phi)**2) + u

model.addDynamic(s, v)
model.addDynamic(v, dvdt)

model.addDynamic(phi, omega)
model.addDynamic(omega, (sin(phi)*G - cos(phi)*dvdt)/R)

# stationary final state
model.addFinal(s, eq=1)
model.addFinal(v, eq=0)
model.addFinal(phi, eq=0)
model.addFinal(omega, eq=0)

# minimize action u**2 or angle phi**2
# Remark: integral phi**2 dt only changes by ~1% whether you choose
# phi**2 or u**2 as the langrange term
model.addLagrange(phi**2, Objective.MINIMIZE)

model.generate()
