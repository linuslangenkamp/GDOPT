from optimization import *
import random


model = Model("generalizedBatchReactor")

N = 100

x = []
for v in range(N):
    if v == 0:
        x.append(model.addState(symbol=f"obj", start=0))
    else:
        x.append(model.addState(symbol=f"x{v}", start=1/(N-1)))

u = model.addInput(symbol="u", lb=0, ub=5, start=0)

energy = model.addState(symbol="energy", start=0) # total energy consumed by the control

coeffs = []
for origin in range(N):
    oList = []
    for target in range(N):
        if target <= origin or target - origin > 5:
            oList.append(0)
        else:
            oList.append(random.uniform(0, 1))
    coeffs.append(oList)

for origin in range(N):
    for target in range(origin):
        coeffs[origin][origin] -= coeffs[target][origin]

for v in range(N):
    if v == 0:
        model.addDynamic(x[0], sum(u * coeffs[0][k] * x[k] for k in range(N)) - 15 * x[0] * u**2)
    elif 0 < v < N-1:
        model.addDynamic(x[v], sum(u * coeffs[v][k] * x[k] for k in range(N)))
    else:
        model.addDynamic(x[N-1], sum(u * coeffs[N-1][k] * x[k] for k in range(N)) + 15 * x[0] * u**2)

model.addDynamic(energy, u**4)
model.addFinal(energy, ub=0.5)

model.addMayer(x[0], Objective.MAXIMIZE)

model.generate()

model.optimize(tf=1, steps=25, rksteps=3,
               flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MA57},
               meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM,
                          "meshIterations": 0})

model.plotInputs(dots=True)
model.plot(specifCols=["obj"], dots=True)
