from optimization import *
import random

# this model, although it contains not that many non-zeros, will become dense because of fill-in effects
# maybe (actually) use a non-direct linear solver, like a GMRES / MINRES / BICGSTAB?!
# runtime scales ~quadratically with the number of states N

model = Model("generalizedBatchReactor")

N = 50
CHAIN_SIZE = 15

x = []
for v in range(N):
    if v == 0:
        x.append(model.addState(symbol=f"obj", start=0))
    else:
        x.append(model.addState(symbol=f"x{v}", start=1 / (N - 1)))

u = model.addInput(symbol="u", lb=0, ub=5, guess=1.5 - t)

energy = model.addState(symbol="energy", start=0)  # total energy consumed by the control

EXP_E = model.addRuntimeParameter(default=4, symbol="EXPONENT_ENERGY")
DEPL = model.addRuntimeParameter(default=50, symbol="DEPLETION_COEFF")

coeffs = []
for x1 in range(N):
    oList = []
    for x2 in range(N):
        if x2 <= x1 or x2 - x1 > CHAIN_SIZE:
            oList.append(0)
        else:
            oList.append(random.uniform(0, 1))
    coeffs.append(oList)

for x1 in range(N):
    for x2 in range(x1):
        coeffs[x1][x1] -= coeffs[x2][x1]

for v in range(N):
    if v == 0:
        model.addDynamic(x[0], sum(u * coeffs[0][k] * x[k] for k in range(N)) - DEPL * x[0] * u**2)
    elif 0 < v < N - 1:
        model.addDynamic(x[v], sum(u * coeffs[v][k] * x[k] for k in range(N)))
    else:
        model.addDynamic(x[N - 1], sum(u * coeffs[N - 1][k] * x[k] for k in range(N)) + DEPL * x[0] * u**2)

model.addDynamic(energy, u**EXP_E)

model.addFinal(energy, ub=0.5)

model.addMayer(x[0], Objective.MAXIMIZE)

model.generate()

model.optimize(
    tf=1,
    steps=25,
    rksteps=3,
    flags={"outputPath": "/tmp", "linearSolver": LinearSolver.MA57, "exportJacobianPath": "/tmp"},
    meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM, "meshIterations": 5},
)

model.plot(specifCols=["obj", "u", "energy"], dots=Dots.ALL)
model.plotMeshRefinement()
# model.plotSparseMatrix(MatrixType.JACOBIAN)
