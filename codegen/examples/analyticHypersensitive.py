from codegen.optimization import *

"""
* Hypersensitive OCP with an analytic solution from
* "Symplectic algorithms with mesh refinement for a hypersensitive optimal control problem"
* Peng, et. al., 2014
* x*(t) = (1.5 * exp((t-2*t_f)*2^0.5) - exp((t-t_f)*2^0.5) + exp((-t-t_f)*2^0.5) - 1.5 * exp(-t*2^0.5))  / (exp(-t_f * 2 * 2^0.5) - 1)
* u*(t) = x*(t) + d/dt x*(t)
* f(x*) = int_0^tf 1/2 (x*(t)^2 + u*(t)^2) dt
* Analytic solutions:
* tf = 25,    f(x*) = 1.673097038856277579275453
* tf = 100,   f(x*) = 1.673097038856279454302744
* tf = 1000,  f(x*) = 1.673097038856279454302744
* tf = 10000, f(x*) = 1.673097038856279454302744
"""

model = Model("analyticHypersensitive")

x = model.addX(start=1.5)
u = model.addU()

model.addF(x, -x + u)
model.addR(1. - x, eq=0)
model.addL(0.5 * (x**2 + u**2))

model.generate()

model.optimize(tf=10000, steps=200, rksteps=7,
               flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MA57},
               meshFlags={"meshAlgorithm": MeshAlgorithm.L2_BOUNDARY_NORM,
                          "meshIterations": 10,
                          "meshLevel": 0})

model.plot(interval=[9980, 10000], dots=True)