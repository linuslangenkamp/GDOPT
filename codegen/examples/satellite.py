from optimization import *


model = Model("satellite")

I1 = model.addConst(Constant("I1", 1000000))
I2 = model.addConst(Constant("I2", 833333))
I3 = model.addConst(Constant("I3", 916667))
T1S = model.addConst(Constant("T1S", 550))
T2S = model.addConst(Constant("T2S", 50))
T3S = model.addConst(Constant("T3S", 550))

x1 = model.addVar(State('e1', start=0))
x2 = model.addVar(State('e2', start=0))
x3 = model.addVar(State('e3', start=0))
x4 = model.addVar(State('e4', start=1))
x5 = model.addVar(State('w1', start=0.01))
x6 = model.addVar(State('w2', start=0.005))
x7 = model.addVar(State('w3', start=0.001))
u1 = model.addVar(Control('T1'))
u2 = model.addVar(Control('T2'))
u3 = model.addVar(Control('T3'))

model.addDynamic(x1, 0.5 * (x5 * x4 - x6 * x3 + x7 * x2))
model.addDynamic(x2, 0.5 * (x5 * x3 + x6 * x4 - x7 * x1))
model.addDynamic(x3, 0.5 * (-x5 * x2 + x6 * x1 + x7 * x4))
model.addDynamic(x4, -0.5 * (x5 * x1 + x6 * x2 + x7 * x3))

model.addDynamic(x5, ((I2 - I3) * x6 * x7 + T1S * u1) / I1)
model.addDynamic(x6, ((I3 - I1) * x7 * x5 + T2S * u2) / I2)
model.addDynamic(x7, ((I1 - I2) * x5 * x6 + T3S * u3) / I3)

model.addMayer((x1 - 0.70106)*(x1 - 0.70106) + (x2 - 0.0923)*(x2 - 0.0923) + (x3 - 0.56098)*(x3 - 0.56098) + 
                (x4 - 0.43047)*(x4 - 0.43047) + x5 * x5 + x6 * x6 + x7 * x7)
model.addLagrange(0.5 * (u1**2 + u2**2 + u3**2))


model.generate("example")
