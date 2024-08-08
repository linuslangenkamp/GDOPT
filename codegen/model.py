from sympy import Symbol
from enum import Enum


class Vartype(Enum):
    STATE = 1
    INPUT = 2
    PARAMETER = 3


class Variable(Symbol):
    id_counter = 0 
    
    def __new__(cls, symbol, lb=-float("inf"), ub=float("inf")):
        obj = Symbol.__new__(cls, symbol)
        obj.lb = lb
        obj.ub = ub
        obj.type = None
        obj.id = cls.id_counter 
        cls.id_counter += 1
        return obj


class State(Variable):
    
    def __new__(cls, symbol, start, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol, lb, ub)
        obj.type = Vartype.STATE
        obj.start = start
        return obj


class Input(Variable):
    
    def __new__(cls, symbol, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol, lb, ub)
        obj.type = Vartype.INPUT
        return obj
        
        
class Parameter(Variable):
    
    def __new__(cls, symbol, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol, lb, ub)
        obj.type = Vartype.PARAMETER
        return obj


class Expression:
    
    def __init__(self, expr):
        self.expr = expr


class DynExpression(Expression):
    
    def __init__(self, diffVar, expr):
        super().__init__(expr)
        self.diffVar = diffVar


class Constraint(Expression):
    
    def __init__(self, expr, lb=-float("inf"), ub=float("inf")):
        super().__init__(expr)
        self.lb = lb
        self.ub = ub


class ParametricConstraint(Expression):
    
    def __init__(self, expr, lb=-float("inf"), ub=float("inf")):
        super().__init__(expr)
        self.lb = lb
        self.ub = ub


class Model:
    
    def __init__(self):
        self.xVars = []
        self.uVars = []
        self.pVars = []
        self.M = None
        self.L = None
        self.F = []
        self.G = []
        self.R = []
        self.A = []
        self.alias = {}
        
    def addVar(self, variable):
        if variable.type == Vartype.STATE:
            self.xVars.append(variable)
        elif variable.type == Vartype.INPUT:
            self.uVars.append(variable)
        elif variable.type == Vartype.PARAMETER:
            self.pVars.append(variable)
            
    def addMayer(self, expr):
        self.M = expr

    def addLagrange(self, expr):
        self.L = expr
    
    def addDynamic(self, diffVar, expr):
        self.F.append(DynExpression(diffVar, expr))
    
    def addPath(self, expr, lb=-float("inf"), ub=float("inf")):
        self.G.append(Constraint(expr, lb=lb, ub=ub))
    
    def addFinal(self, expr, lb=-float("inf"), ub=float("inf")):
        self.R.append(Constraint(expr, lb=lb, ub=ub))
        
    def addParametric(self, expr, lb=-float("inf"), ub=float("inf")):
        self.A.append(ParametricConstraint(expr, lb=lb, ub=ub))

    def codegen(self):
        self.xVars.sort(key=lambda v: v.id, reverse=False)
        self.uVars.sort(key=lambda v: v.id, reverse=False)
        self.pVars.sort(key=lambda v: v.id, reverse=False)
        self.F.sort(key=lambda dyn: dyn.diffVar.id, reverse=False)
		


u1 = Input('u1', lb=0, ub=1)
x1 = State('x1', start=1)
x2 = State('x2', start=1)

p1 = Parameter('p1')

model = Model()
model.addVar(x1)
model.addVar(x2)
model.addVar(u1)
model.addVar(p1)

model.addDynamic(x2, x1 * u1)
model.addDynamic(x1, x2 + u1)

model.addPath(u1 + x2 + x1, lb=0, ub=3)

model.codegen()

