from optimization import *

"""
* Batch Reactor from Parallel Multiple-Shooting and Collocation Optimization with OpenModelica,
    * Bachmann, Ochel, et. al., 2012

optimal objectives: f(x*) = -0.57354505750936147 with n = 1e5, m = 3
-0.57354505720920057 with n = 5e5, m = 3
-0.57354505683399926 with n = 1e6, m = 3, time = 140.658s, MA57
-0.57354505633373065 with n = 1e6, m = 5, time = 327.169s, MA57
-0.57354505723421423 with n = 2e5, m = 7, time = 127.079s, MA57
-0.57354505670893241 with n = 5e5, m = 7, time = 223.389s, MA57

model batchReactor
Real x1(start=1, fixed=true, min=0, max=1);
Real x2(start=0, fixed=true, min=0, max=1);
Real may = -x2 annotation(isMayer=true);
input Real u(min=0, max=5);
equation
der(x1) = -(u + u^2/2) * x1;
der(x2) = u * x1;
annotation(experiment(StartTime=0, StopTime=1, Tolerance=1e-14),
           __OpenModelica_simulationFlags(solver="optimization", optimizerNP="3"),
           __OpenModelica_commandLineOptions="+g=Optimica");
end BatchReactor;
"""

model = Model("batchReactor")

x1 = model.addState(symbol="x1", start=1, lb=0, ub=1)
x2 = model.addState(symbol="x2", start=0, lb=0, ub=1)

u = model.addInput(symbol="u", lb=0, ub=5)

model.addMayer(x2, Objective.MAXIMIZE)

model.addDynamic(x1, -(u + u**2 / 2) * x1)
model.addDynamic(x2, u * x1)

model.generate()

model.optimize(tf=1, steps=25, rksteps=3,
               flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MUMPS,
                      "exportJacobianPath": "/tmp"},
               meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM,
                          "meshIterations": 10})

model.plot(dots=True)
model.plotSparseMatrix(MatrixType.JACOBIAN)
model.printResults()
model.plotPointCumulativeCount(meshIteration="all")