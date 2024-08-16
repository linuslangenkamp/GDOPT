from optimization import *

### I = [t0=0, tf=1]
### 1 interval with Step2 O(h^3) Radau scheme is sufficient for exact sol.
### p1* = 5/8, p2* = -5/8 with u(1)=-1 (lb)
### f(x*) = -25/64

model = Model("simpleParameter")

x1 = model.addX(start=0)

u1 = model.addU(lb=-1, ub=1)

p1 = model.addP()
p2 = model.addP()

# x1' = 2t * p1, x1(0) = 0 => x1 = p1 * t^2
model.addF(x1, 2*t * p1)

# 0.2 <= p2 - p1 * t^2  + u(t) <= 0.25
model.addG(p2 - x1 + u1, lb=0.2, ub=0.25)    

model.addM(p1 * p2)

model.generate()

model.optimize(steps=1, rksteps=2, tf=1,
               flags={"outputPath": "/tmp",
                      "linearSolver": LinearSolver.MUMPS,
                      "tolerance": 1e-14},
               meshFlags={})

model.printResults()
