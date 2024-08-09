from sympy import Symbol, diff, sin, ccode
from enum import Enum

class Variable(Symbol):
    id_counter = 0

    def __new__(cls, symbol, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol)
        obj.lb = lb
        obj.ub = ub
        return obj
        
class State(Variable):
    id_counter = 0

    def __new__(cls, symbol, start, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol, lb, ub)
        obj.id = cls.id_counter
        obj.symbol = f'x[{obj.id}]'
        obj.start = start
        cls.id_counter += 1
        return obj

class Input(Variable):
    id_counter = 0

    def __new__(cls, symbol, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol, lb, ub)
        obj.id = cls.id_counter
        obj.symbol = f'u[{obj.id}]'
        cls.id_counter += 1
        return obj

class Parameter(Variable):
    id_counter = 0

    def __new__(cls, symbol, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol, lb, ub)
        obj.id = cls.id_counter
        obj.symbol = f'p[{obj.id}]'
        cls.id_counter += 1
        return obj

class Expression:
    def __init__(self, expr):
        self.expr = expr
    
    def codegen(self, name, variables):
        cEval = ccode(self.expr)

        # Generate second-order derivatives only for lower triangular part
        cEvalDiff2 = [[], [], [], [], [], []]
        adjDiff_indices = [[], [], [], [], [], []]  # indXX, indUX, indUU, indPX, indPU, indPP

        for i, var1 in enumerate(variables):
            for j, var2 in enumerate(variables):
                if type(var1) == State and type(var2) == State and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[0].append(ccode(der))
                        adjDiff_indices[0].append((var1.id, var2.id))  # indXX
                elif type(var1) == Input and type(var2) == State:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[1].append(ccode(der))
                        adjDiff_indices[1].append((var1.id, var2.id))  # indUX
                elif type(var1) == Input and type(var2) == Input and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[2].append(ccode(der))
                        adjDiff_indices[2].append((var1.id, var2.id))  # indUU
                elif type(var1) == Parameter and type(var2) == State:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[3].append(ccode(der))
                        adjDiff_indices[3].append((var1.id, var2.id))  # indPX
                elif type(var1) == Parameter and type(var2) == Input:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[4].append(ccode(der))
                        adjDiff_indices[4].append((var1.id, var2.id))  # indPU
                elif type(var1) == Parameter and type(var2) == Parameter and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[5].append(ccode(der))
                        adjDiff_indices[5].append((var1.id, var2.id))  # indPP

        adj_indices = [[], [], []]  # indX, indU, indP
        cEvalDiff = [[], [], []]
        for i, var in enumerate(variables):
            der = diff(self.expr, var)
            if der != 0:
                if type(var) == State:
                    adj_indices[0].append(var.id)
                    cEvalDiff[0].append(ccode(der))
                elif type(var) == Input:
                    adj_indices[1].append(var.id)
                    cEvalDiff[1].append(ccode(der))
                elif type(var) == Parameter:
                    adj_indices[2].append(var.id)
                    cEvalDiff[2].append(ccode(der))
                    
        adj = "{{{}, {}, {}}}".format(
            "{{{}}}".format(", ".join(map(str, adj_indices[0]))),
            "{{{}}}".format(", ".join(map(str, adj_indices[1]))),
            "{{{}}}".format(", ".join(map(str, adj_indices[2])))
        )
        adjDiff = "{{{}, {}, {}, {}, {}, {}}}".format(
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[0])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[1])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[2])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[3])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[4])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[5]))
        )
        
        out = f"class {name} : public Expression {{\n"
        out += "public:\n"
        out += f"\tstatic std::unique_ptr<{name}> create() {{\n"
        out += f"\t\tAdjacency adj{adj};\n"
        out += f"\t\tAdjacencyDiff adjDiff{adjDiff};\n"
        out += f"\t\treturn std::unique_ptr<{name}>(new {name}(std::move(adj), std::move(adjDiff)));\n"
        out += "\t}\n\n"
        
        out += "\tdouble eval(const double *x, const double *u, const double *p, double t) override {\n"
        out += f"\t\treturn {cEval};\n"
        out += "\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {{\n"
        out += "\t\treturn {{"
        out += "{{{}}}, ".format(", ".join(cEvalDiff[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff[1]))
        out += "{{{}}}".format(", ".join(cEvalDiff[2]))
        out += "}};\n\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {{\n"
        out += "\t\treturn {{"
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[1]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[2]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[3]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[4]))
        out += "{{{}}}".format(", ".join(cEvalDiff2[5]))
        out += "}};\n\t}\n"
        
        out += "private:\n"
        out += f"\t{name}(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {{}}\n"
        out += "};\n"
        return out
        out += "};\n"
        return out

class DynExpression(Expression):
    
    def __init__(self, diffVar, expr):
        super().__init__(expr)
        self.diffVar = diffVar

class Constraint(Expression):
    
    def __init__(self, expr, lb=-float("inf"), ub=float("inf")):
        super().__init__(expr)
        self.lb = lb
        self.ub = ub

    def codegen(self, name, variables):
        cEval = ccode(self.expr)

        # Generate second-order derivatives only for lower triangular part
        cEvalDiff2 = [[], [], [], [], [], []]
        adjDiff_indices = [[], [], [], [], [], []]  # indXX, indUX, indUU, indPX, indPU, indPP

        for i, var1 in enumerate(variables):
            for j, var2 in enumerate(variables):
                if type(var1) == State and type(var2) == State and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[0].append(ccode(der))
                        adjDiff_indices[0].append((var1.id, var2.id))  # indXX
                elif type(var1) == Input and type(var2) == State:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[1].append(ccode(der))
                        adjDiff_indices[1].append((var1.id, var2.id))  # indUX
                elif type(var1) == Input and type(var2) == Input and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[2].append(ccode(der))
                        adjDiff_indices[2].append((var1.id, var2.id))  # indUU
                elif type(var1) == Parameter and type(var2) == State:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[3].append(ccode(der))
                        adjDiff_indices[3].append((var1.id, var2.id))  # indPX
                elif type(var1) == Parameter and type(var2) == Input:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[4].append(ccode(der))
                        adjDiff_indices[4].append((var1.id, var2.id))  # indPU
                elif type(var1) == Parameter and type(var2) == Parameter and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[5].append(ccode(der))
                        adjDiff_indices[5].append((var1.id, var2.id))  # indPP


        adj_indices = [[], [], []]  # indX, indU, indP
        cEvalDiff = [[], [], []]
        for i, var in enumerate(variables):
            der = diff(self.expr, var)
            if der != 0:
                if type(var) == State:
                    adj_indices[0].append(var.id)
                    cEvalDiff[0].append(ccode(der))
                elif type(var) == Input:
                    adj_indices[1].append(var.id)
                    cEvalDiff[1].append(ccode(der))
                elif type(var) == Parameter:
                    adj_indices[2].append(var.id)
                    cEvalDiff[2].append(ccode(der))

        adj = "{{{}, {}, {}}}".format(
            "{{{}}}".format(", ".join(map(str, adj_indices[0]))),
            "{{{}}}".format(", ".join(map(str, adj_indices[1]))),
            "{{{}}}".format(", ".join(map(str, adj_indices[2])))
        )
        adjDiff = "{{{}, {}, {}, {}, {}, {}}}".format(
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[0])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[1])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[2])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[3])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[4])),
            "{{{}}}".format(", ".join(f"{{{i}, {j}}}" for i, j in adjDiff_indices[5]))
        )
        
        lb = self.lb if self.lb != -float('inf') else "MINUS_INFINITY"
        ub = self.ub if self.ub != float('inf') else "PLUS_INFINITY"
        
        out = f"class {name} : public Constraint {{\n"
        out += "public:\n"
        out += f"\tstatic std::unique_ptr<{name}> create() {{\n"
        out += f"\t\tAdjacency adj{adj};\n"
        out += f"\t\tAdjacencyDiff adjDiff{adjDiff};\n"
        out += f"\t\treturn std::unique_ptr<{name}>(new {name}(std::move(adj), std::move(adjDiff), {lb}, {ub}));\n"
        out += "\t}\n\n"
        
        out += "\tdouble eval(const double *x, const double *u, const double *p, double t) override {\n"
        out += f"\t\treturn {cEval};\n"
        out += "\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {{\n"
        out += "\t\treturn {{"
        out += "{{{}}}, ".format(", ".join(cEvalDiff[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff[1]))
        out += "{{{}}}".format(", ".join(cEvalDiff[2]))
        out += "}};\n\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {{\n"
        out += "\t\treturn {{"
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[1]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[2]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[3]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[4]))
        out += "{{{}}}".format(", ".join(cEvalDiff2[5]))
        out += "}};\n\t}\n"
        
        out += "private:\n"
        out += f"\t{name}(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {{}}\n"
        out += "};\n"
        return out

class ParametricConstraint(Expression):
    
    def __init__(self, expr, lb=-float("inf"), ub=float("inf")):
        super().__init__(expr)
        self.lb = lb
        self.ub = ub
        
    def codegen(self, name, variables):
        cEval = ccode(self.expr)

        # Generate second-order derivatives only for lower triangular part
        cEvalDiff2 = []
        adjDiff_indices = []  # indPP

        for i, var1 in enumerate(variables):
            for j, var2 in enumerate(variables):
                if type(var1) == Parameter and type(var2) == Parameter and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2.append(ccode(der))
                        adjDiff_indices.append((var1.id, var2.id))  # indPP


        adj_indices = []  # indX, indU, indP
        cEvalDiff = []
        for i, var in enumerate(variables):
            der = diff(self.expr, var)
            if der != 0:
                if type(var) == Parameter:
                    adj_indices.append(var.id)
                    cEvalDiff.append(ccode(der))
                    
        adj = f"{{{', '.join(map(str, adj_indices))}}}"
        adjDiff = "{{{}}}".format(", ".join("{{{}}}".format(", ".join(map(str, tpl))) for tpl in adjDiff_indices))
        
        
        lb = self.lb if self.lb != -float('inf') else "MINUS_INFINITY"
        ub = self.ub if self.ub != float('inf') else "PLUS_INFINITY"
           
        out = f"class {name} : public ParamConstraint {{\n"
        out += "public:\n"
        out += f"\tstatic std::unique_ptr<{name}> create() {{\n"
        out += f"\t\tParamAdjacency adj{adj};\n"
        out += f"\t\tParamAdjacencyDiff adjDiff{adjDiff};\n"
        out += f"\t\treturn std::unique_ptr<{name}>(new {name}(std::move(adj), std::move(adjDiff), {lb}, {ub}));\n"
        out += "\t}\n\n"
        
        out += "\tdouble eval(const double* p) override {\n"
        out += f"\t\treturn {cEval};\n"
        out += "\t}\n\n"
        
        out += f"\tstd::vector<double> evalDiff(const double* p) override {{\n"
        out += "\t\treturn "
        out += f"{{{', '.join(cEvalDiff)}}}"
        out += ";\n\t}\n\n"
        
        out += f"\tstd::vector<double> evalDiff2(const double* p) override {{\n"
        out += "\t\treturn "
        out += "{{{}}}".format(", ".join(cEvalDiff2))
        out += ";\n\t}\n"
        
        out += "private:\n"
        out += f"\t{name}(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {{}}\n"
        out += "};\n"
        return out
        
class Model:
    
    def __init__(self, name):
        self.xVars = []
        self.uVars = []
        self.pVars = []
        self.M = None
        self.L = None
        self.F = []
        self.G = []
        self.R = []
        self.A = []
        self.name = name
        
    def addVar(self, variable):
        if type(variable) == State:
            self.xVars.append(variable)
        elif type(variable) == Input:
            self.uVars.append(variable)
        elif type(variable) == Parameter:
            self.pVars.append(variable)
        return variable
            
    def addMayer(self, expr):
        self.M = Expression(expr)

    def addLagrange(self, expr):
        self.L = Expression(expr)
    
    def addDynamic(self, diffVar, expr):
        self.F.append(DynExpression(diffVar, expr))
    
    def addPath(self, expr, lb=-float("inf"), ub=float("inf")):
        self.G.append(Constraint(expr, lb=lb, ub=ub))
    
    def addFinal(self, expr, lb=-float("inf"), ub=float("inf")):
        self.R.append(Constraint(expr, lb=lb, ub=ub))
        
    def addParametric(self, expr, lb=-float("inf"), ub=float("inf")):
        if set(self.pVars).issuperset(expr.free_symbols):
            self.A.append(ParametricConstraint(expr, lb=lb, ub=ub))
        else:
            raise("Parametric constraints only allow parametric variables")
    
    def generate(self):
        print(self.xVars)
        if len(self.F) != len(self.xVars):
            raise ValueError("#states != #differential equations") 
        self.F.sort(key=lambda eq: eq.diffVar.id, reverse=False)
        allVars = self.xVars + self.uVars + self.pVars
        print(allVars)

### GLOBAL VAR DEFINITIONS
Continous = Input
Control = Input
t = Symbol("t")
###


model = Model("Test")

x1 = model.addVar(State('x1', start=0.0))
x2 = model.addVar(State('x2', start=0.0))
u1 = model.addVar(Control('u1'))
p1 = model.addVar(Parameter('p1'))
p2 = model.addVar(Parameter('p2'))

variables = model.xVars + model.uVars + model.pVars

model.addDynamic(x2, x1 + x2*u1 + sin(p1)*x1 + p2**2)
model.addDynamic(x1, x2 + x1**2  + t )
model.addPath(x1**2 + x2**2, lb=2)
model.addLagrange(x1 + x2*u1)
model.addParametric(p1**2 + p2**2, lb=0, ub=1)
cpp_code = model.G[0].codegen("G", variables)

cpp_code = model.A[0].codegen("AAAA", variables)

print(cpp_code)

print(model.F)
model.generate()
print(model.F)
model.L.codegen("test", variables)

