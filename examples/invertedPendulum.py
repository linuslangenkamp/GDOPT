from optimization import *

# choose implicit Euler, since it's the most stable

model = Model("invertedPendulum")

Ms = model.addRuntimeParameter(default=1, symbol="Ms")
Mp = model.addRuntimeParameter(default=0.5, symbol="Mp")
R = model.addRuntimeParameter(default=1, symbol="R")
G = -9.81

s = model.addState(start=0)
v = model.addState(start=0)

phi = model.addState(start=pi-0.001)
omega = model.addState(start=0)

u = model.addInput(lb=-2.5, ub=2.5)

dvdt = u + (sin(-phi)*Mp*R*omega**2 - cos(phi)*sin(-phi)*Mp*G) / (Ms + Mp * sin(-phi)**2) 

model.addDynamic(s, v)
model.addDynamic(v, dvdt)

model.addDynamic(phi, omega)
model.addDynamic(omega, (sin(-phi)*G - cos(phi)*dvdt)/R)

model.addLagrange(sin(phi/2)**2, Objective.MINIMIZE)

model.generate()

model.setValue(Ms, 1.5)

model.optimize(steps=5000, rksteps=1, tf=12,
               flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MA57},
               meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM,
                          "meshIterations": 0})

model.plot(meshIteration=0)
