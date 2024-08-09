from optimization import *


model = Model("satellite")

# add const  : constant = model.addConst(constValue, optional<symbol in code>)
# add var    : variable = model.addVar(Vartyp<State, Input, Parameter>, optional<symbol for debug>, start if state, Optional: lb, Optional: ub)
# Alias: Input = Control = Continous

# add dynamic constr: model.addDynamic(x', f) for ODE: x' = f
# add path constr   : model.addPath(g, Optional: lb, Optional: ub) 
# add final constr  : model.addFinal(r, Optional: lb, Optional: ub) 
# add param constr  : model.addParametric(a, Optional: lb, Optional: ub) 
# Alias: addDynamic = addOde = addF

I1 = model.addConst(1000000)
I2 = model.addConst(833333)
I3 = model.addConst(916667)

T1S = model.addConst(550)
T2S = model.addConst(50)
T3S = model.addConst(550)

M1 = model.addConst(0.70106)
M2 = model.addConst(0.0923)
M3 = model.addConst(0.56098)
M4 = model.addConst(0.43047)

x1 = model.addVar(State(start=0))
x2 = model.addVar(State(start=0))
x3 = model.addVar(State(start=0))
x4 = model.addVar(State(start=1))
x5 = model.addVar(State(start=0.01))
x6 = model.addVar(State(start=0.005))
x7 = model.addVar(State(start=0.001))

u1 = model.addVar(Control())
u2 = model.addVar(Control())
u3 = model.addVar(Control())

model.addDynamic(x1, 0.5 * (x5 * x4 - x6 * x3 + x7 * x2))
model.addDynamic(x2, 0.5 * (x5 * x3 + x6 * x4 - x7 * x1))
model.addDynamic(x3, 0.5 * (-x5 * x2 + x6 * x1 + x7 * x4))
model.addDynamic(x4, -0.5 * (x5 * x1 + x6 * x2 + x7 * x3))

model.addDynamic(x5, ((I2 - I3) * x6 * x7 + T1S * u1) / I1)
model.addDynamic(x6, ((I3 - I1) * x7 * x5 + T2S * u2) / I2)
model.addDynamic(x7, ((I1 - I2) * x5 * x6 + T3S * u3) / I3)

model.addMayer((x1 - M1)**2 + (x2 - M2)**2 + (x3 - M3)**2 + 
                (x4 - M4)**2 + x5**2 + x6**2 + x7**2, Objective.MINIMIZE)
                
model.addLagrange(0.5 * (u1**2 + u2**2 + u3**2), Objective.MINIMIZE)


model.generate("satelliteGenerated")
