from sympy import Symbol, diff, sin, ccode
from enum import Enum

class Vartype(Enum):
    STATE = 1
    INPUT = 2
    PARAMETER = 3

class Variable(Symbol):
    id_counter = 0

    def __new__(cls, symbol, var_type, lb=-float("inf"), ub=float("inf"), start=None):
        obj = super().__new__(cls, symbol)
        obj.lb = lb
        obj.ub = ub
        obj.type = var_type
        obj.start = start
        return obj
        
class State(Variable):
    id_counter = 0

    def __new__(cls, symbol, start, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol, Vartype.STATE, lb, ub, start)
        obj.id = cls.id_counter
        obj.symbol = f'x[{obj.id}]'
        cls.id_counter += 1
        return obj

class Input(Variable):
    id_counter = 0

    def __new__(cls, symbol, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol, Vartype.INPUT, lb, ub)
        obj.id = cls.id_counter
        obj.symbol = f'u[{obj.id}]'
        cls.id_counter += 1
        return obj

class Parameter(Variable):
    id_counter = 0

    def __new__(cls, symbol, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol, Vartype.PARAMETER, lb, ub)
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
                if var1.type == Vartype.STATE and var2.type == Vartype.STATE and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[0].append(ccode(der))
                        adjDiff_indices[0].append((var1.id, var2.id))  # indXX
                elif var1.type == Vartype.INPUT and var2.type == Vartype.STATE:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[1].append(ccode(der))
                        adjDiff_indices[1].append((var1.id, var2.id))  # indUX
                elif var1.type == Vartype.INPUT and var2.type == Vartype.INPUT and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[2].append(ccode(der))
                        adjDiff_indices[2].append((var1.id, var2.id))  # indUU
                elif var1.type == Vartype.PARAMETER and var2.type == Vartype.STATE:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[3].append(ccode(der))
                        adjDiff_indices[3].append((var1.id, var2.id))  # indPX
                elif var1.type == Vartype.PARAMETER and var2.type == Vartype.INPUT:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[4].append(ccode(der))
                        adjDiff_indices[4].append((var1.id, var2.id))  # indPU
                elif var1.type == Vartype.PARAMETER and var2.type == Vartype.PARAMETER and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[5].append(ccode(der))
                        adjDiff_indices[5].append((var1.id, var2.id))  # indPP

        # Create the adjacency list for non-zero first derivatives
        adj_indices = [[], [], []]  # indX, indU, indP
        cEvalDiff = [[], [], []]
        for i, var in enumerate(variables):
            der = diff(self.expr, var)
            if der != 0:
                if var.type == Vartype.STATE:
                    adj_indices[0].append(var.id)
                    cEvalDiff[0].append(ccode(der))
                elif var.type == Vartype.INPUT:
                    adj_indices[1].append(var.id)
                    cEvalDiff[1].append(ccode(der))
                elif var.type == Vartype.PARAMETER:
                    adj_indices[2].append(var.id)
                    cEvalDiff[2].append(ccode(der))

        # Format adjacency information for C++ code
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
                if var1.type == Vartype.STATE and var2.type == Vartype.STATE and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[0].append(ccode(der))
                        adjDiff_indices[0].append((var1.id, var2.id))  # indXX
                elif var1.type == Vartype.INPUT and var2.type == Vartype.STATE:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[1].append(ccode(der))
                        adjDiff_indices[1].append((var1.id, var2.id))  # indUX
                elif var1.type == Vartype.INPUT and var2.type == Vartype.INPUT and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[2].append(ccode(der))
                        adjDiff_indices[2].append((var1.id, var2.id))  # indUU
                elif var1.type == Vartype.PARAMETER and var2.type == Vartype.STATE:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[3].append(ccode(der))
                        adjDiff_indices[3].append((var1.id, var2.id))  # indPX
                elif var1.type == Vartype.PARAMETER and var2.type == Vartype.INPUT:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[4].append(ccode(der))
                        adjDiff_indices[4].append((var1.id, var2.id))  # indPU
                elif var1.type == Vartype.PARAMETER and var2.type == Vartype.PARAMETER and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        cEvalDiff2[5].append(ccode(der))
                        adjDiff_indices[5].append((var1.id, var2.id))  # indPP

        
        # Create the adjacency list for non-zero first derivatives
        adj_indices = [[], [], []]  # indX, indU, indP
        cEvalDiff = [[], [], []]
        for i, var in enumerate(variables):
            der = diff(self.expr, var.symbol)
            if der != 0:
                if var.type == Vartype.STATE:
                    adj_indices[0].append(var.id)
                    cEvalDiff[0].append(ccode(der))
                elif var.type == Vartype.INPUT:
                    adj_indices[1].append(var.id)
                    cEvalDiff[1].append(ccode(der))
                elif var.type == Vartype.PARAMETER:
                    adj_indices[2].append(var.id)
                    cEvalDiff[2].append(ccode(der))

        # Format adjacency information for C++ code
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
        
        out = f"class {name} : public Constraint {{\n"
        out += "public:\n"
        out += f"\tstatic std::unique_ptr<{name}> create() {{\n"
        out += f"\t\tAdjacency adj{adj};\n"
        out += f"\t\tAdjacencyDiff adjDiff{adjDiff};\n"
        out += f"\t\treturn std::unique_ptr<{name}>(new {name}(std::move(adj), std::move(adjDiff), {self.lb}, {self.ub}));\n"
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
        if variable.type == Vartype.STATE:
            self.xVars.append(variable)
        elif variable.type == Vartype.INPUT:
            self.uVars.append(variable)
        elif variable.type == Vartype.PARAMETER:
            self.pVars.append(variable)
            
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
        self.A.append(ParametricConstraint(expr, lb=lb, ub=ub))
    
    def generate(self):
        self.F.sort(key=lambda dyn: dyn.diffVar.id, reverse=False)
# Define symbols for states, inputs, and parameters

# Define variable instances
x1 = State('x1', start=0.0)
x2 = State('x2', start=0.0)
u1 = Input('u1')
p1 = Parameter('p1')
p2 = Parameter('p2')
x3 = State('x3', start=0.0)

# List of variables
variables = [x1, x2, x3, u1, p1, p2]

print(x1.symbol, x1)


model = Model("Test")
model.addDynamic(x2, x1 + x2*u1 + sin(p1)*x1 + p2**2)
model.addDynamic(x1, x2)
model.addPath(x1**2 + x2**2, lb=2, ub=10)
model.addLagrange(x1 + x2*u1 + sin(x3 + x2 - p2))
cpp_code = model.L.codegen("Constr", variables)
print(cpp_code)
print(model.F)
model.generate()
print(model.F)
model.L.codegen("test", variables)

