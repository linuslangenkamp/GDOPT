import os

from sympy import *
import time as timer
from enum import Enum
import pandas as pd
import matplotlib.pyplot as plt 

# TODO: better structure in codegen, add subfolders and so on
# TODO: only diff if first diff != 0
# add vectorized eval of RHS = [f, g]^T, vectorized evalDiff, evalDiff2?
# or with colored jacobian
# TODO: make framework more robust, cleaner
# TODO: Maximize -> show


class InvalidModel(Exception):
    pass


class Objective(Enum):
    MINIMIZE = 1
    MAXIMIZE = 2
    MAX = MAXIMIZE
    MIN = MINIMIZE


class LinearSolver(Enum):
    MUMPS = 1
    MA27 = 2
    MA57 = 3
    MA77 = 4
    MA86 = 5
    MA97 = 6
    PARDISO = 7
    
    
class MeshAlgorithm(Enum):
    NONE = 1
    BASIC = 2
    L2_BOUNDARY_NORM = 3


class Variable(Symbol):
    id_counter = 0
    
    def __new__(cls, symbol, lb=-float("inf"), ub=float("inf")):
        obj = super().__new__(cls, symbol)
        obj.lb = lb
        obj.ub = ub
        return obj


class State(Variable):
    id_counter = 0

    def __new__(cls, start, symbol=None, lb=-float("inf"), ub=float("inf")):
        if symbol == None:
            symbol = f'x[{cls.id_counter}]'
        obj = super().__new__(cls, symbol, lb, ub)
        obj.id = cls.id_counter
        obj.symbol = f'x[{obj.id}]'
        obj.start = start
        cls.id_counter += 1
        return obj


class Input(Variable):
    id_counter = 0

    def __new__(cls, symbol=None, lb=-float("inf"), ub=float("inf")):
        if symbol == None:
            symbol = f'u[{cls.id_counter}]'
        obj = super().__new__(cls, symbol, lb, ub)
        obj.id = cls.id_counter
        obj.symbol = f'u[{obj.id}]'
        cls.id_counter += 1
        return obj


class Parameter(Variable):
    id_counter = 0

    def __new__(cls, symbol=None, lb=-float("inf"), ub=float("inf")):
        if symbol == None:
            symbol = f'p[{cls.id_counter}]'
        obj = super().__new__(cls, symbol, lb, ub)
        obj.id = cls.id_counter
        obj.symbol = f'p[{obj.id}]'
        cls.id_counter += 1
        return obj


class RuntimeParameter(Variable):
    id_counter = 0

    def __new__(cls, default, symbol=None, lb=-float("inf"), ub=float("inf")):
        if symbol == None:
            symbol = f'Parameter{cls.id_counter}'
        else:
            symbol = f'Parameter_{symbol}'
        obj = super().__new__(cls, symbol, lb, ub)
        obj.id = cls.id_counter
        obj.value = default
        obj.symbol = symbol
        cls.id_counter += 1
        return obj

    def setValue(self, value):
        self.value = value

# sorting of elements for adjacency structures
class_order = {
    'State': 0,
    'Input': 1,
    'Parameter': 2
}

def sort_vars(elem):
    return (class_order[elem.__class__.__name__], elem.id)

class Expression:
    def __init__(self, expr):
        self.expr = simplify(expr)
        self.adj = []
        try:
            for sym in expr.free_symbols:
                if type(sym) == State or type(sym) == Input or type(sym) == Parameter:
                    self.adj.append(sym)
            self.adj.sort(key=sort_vars)
        except:
            self.adj = []
    
    def codegen(self, name):
        partialExpression = numbered_symbols(prefix='s')
        
        subst0, exprEval = cse(self.expr, partialExpression)
        cEval = ccode(exprEval[0])
        
        # Generate second-order derivatives only for lower triangular part
        cEvalDiff2 = [[], [], [], [], [], []]
        adjDiff_indices = [[], [], [], [], [], []]  # indXX, indUX, indUU, indPX, indPU, indPP
        allDiffs2 = []
        diffVars2 = []
        
        for i, var1 in enumerate(self.adj):
            for j, var2 in enumerate(self.adj):
                if not ((type(var1) == State and type(var2) == State and i < j) or (type(var1) == Input and type(var2) == Input and i < j) or (type(var1) == Parameter and type(var2) == Parameter and i < j)):
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        allDiffs2.append(der)
                        diffVars2.append((var1, var2))
                
        subst2, substExpr2 = cse(allDiffs2, symbols=partialExpression)


        for i, expr in enumerate(substExpr2):
            var1, var2 = diffVars2[i]
            if type(var1) == State and type(var2) == State:
                cEvalDiff2[0].append(ccode(expr))
                adjDiff_indices[0].append((var1.id, var2.id))  # indXX
            elif type(var1) == Input and type(var2) == State:
                cEvalDiff2[1].append(ccode(expr))
                adjDiff_indices[1].append((var1.id, var2.id))  # indUX
            elif type(var1) == Input and type(var2) == Input:
                cEvalDiff2[2].append(ccode(expr))
                adjDiff_indices[2].append((var1.id, var2.id))  # indUU
            elif type(var1) == Parameter and type(var2) == State:
                cEvalDiff2[3].append(ccode(expr))
                adjDiff_indices[3].append((var1.id, var2.id))  # indPX
            elif type(var1) == Parameter and type(var2) == Input:
                cEvalDiff2[4].append(ccode(expr))
                adjDiff_indices[4].append((var1.id, var2.id))  # indPU
            elif type(var1) == Parameter and type(var2) == Parameter:
                cEvalDiff2[5].append(ccode(expr))
                adjDiff_indices[5].append((var1.id, var2.id))  # indPP

        adj_indices = [[], [], []]  # indX, indU, indP
        cEvalDiff = [[], [], []]
        allDiffs = []
        diffVars = []
        
        for i, var in enumerate(self.adj):
            der = diff(self.expr, var)
            if der != 0:
                allDiffs.append(der)
                diffVars.append(var)
                
        subst1, substExpr = cse(allDiffs, symbols=partialExpression)

        for v, expr in enumerate(substExpr):
            var = diffVars[v]
            if type(var) == State:
                adj_indices[0].append(var.id)
                cEvalDiff[0].append(ccode(expr))
            elif type(var) == Input:
                adj_indices[1].append(var.id)
                cEvalDiff[1].append(ccode(expr))
            elif type(var) == Parameter:
                adj_indices[2].append(var.id)
                cEvalDiff[2].append(ccode(expr))
        
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
        
        for s in subst0:
            out += "        const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
            
        out += f"\t\treturn {cEval};\n"
        out += "\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {{\n"
        
        for s in subst1:
            out += "        const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
            
        out += "\t\treturn {std::vector<double>"
        out += "{{{}}}, ".format(", ".join(cEvalDiff[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff[1]))
        out += "{{{}}}".format(", ".join(cEvalDiff[2]))
        out += "};\n\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {{\n"
        for s in subst2:
            out += "        const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
            
        out += "\t\treturn {std::vector<double>"
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[1]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[2]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[3]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[4]))
        out += "{{{}}}".format(", ".join(cEvalDiff2[5]))
        out += "};\n\t}\n"
        
        out += "private:\n"
        out += f"\t{name}(Adjacency adj, AdjacencyDiff adjDiff) : Expression(std::move(adj), std::move(adjDiff)) {{}}\n"
        out += "};\n\n\n"
        return out


class DynExpression(Expression):
    
    def __init__(self, diffVar, expr):
        super().__init__(expr)
        self.diffVar = diffVar


class Constraint(Expression):
    
    def __init__(self, expr, lb=-float("inf"), ub=float("inf"), eq=None):
        super().__init__(expr)
        if eq != None:
            self.lb = eq
            self.ub = eq
        else:
            self.lb = lb
            self.ub = ub

    def codegen(self, name):
        partialExpression = numbered_symbols(prefix='s')
        
        subst0, exprEval = cse(self.expr, partialExpression)
        cEval = ccode(exprEval[0])
        
        # Generate second-order derivatives only for lower triangular part
        cEvalDiff2 = [[], [], [], [], [], []]
        adjDiff_indices = [[], [], [], [], [], []]  # indXX, indUX, indUU, indPX, indPU, indPP
        allDiffs2 = []
        diffVars2 = []
        
        for i, var1 in enumerate(self.adj):
            for j, var2 in enumerate(self.adj):
                if not ((type(var1) == State and type(var2) == State and i < j) or (type(var1) == Input and type(var2) == Input and i < j) or (type(var1) == Parameter and type(var2) == Parameter and i < j)):
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        allDiffs2.append(der)
                        diffVars2.append((var1, var2))
                
        subst2, substExpr2 = cse(allDiffs2, symbols=partialExpression)


        for i, expr in enumerate(substExpr2):
            var1, var2 = diffVars2[i]
            if type(var1) == State and type(var2) == State:
                cEvalDiff2[0].append(ccode(expr))
                adjDiff_indices[0].append((var1.id, var2.id))  # indXX
            elif type(var1) == Input and type(var2) == State:
                cEvalDiff2[1].append(ccode(expr))
                adjDiff_indices[1].append((var1.id, var2.id))  # indUX
            elif type(var1) == Input and type(var2) == Input:
                cEvalDiff2[2].append(ccode(expr))
                adjDiff_indices[2].append((var1.id, var2.id))  # indUU
            elif type(var1) == Parameter and type(var2) == State:
                cEvalDiff2[3].append(ccode(expr))
                adjDiff_indices[3].append((var1.id, var2.id))  # indPX
            elif type(var1) == Parameter and type(var2) == Input:
                cEvalDiff2[4].append(ccode(expr))
                adjDiff_indices[4].append((var1.id, var2.id))  # indPU
            elif type(var1) == Parameter and type(var2) == Parameter:
                cEvalDiff2[5].append(ccode(expr))
                adjDiff_indices[5].append((var1.id, var2.id))  # indPP

        adj_indices = [[], [], []]  # indX, indU, indP
        cEvalDiff = [[], [], []]
        allDiffs = []
        diffVars = []
        
        for i, var in enumerate(self.adj):
            der = diff(self.expr, var)
            if der != 0:
                allDiffs.append(der)
                diffVars.append(var)
                
        subst1, substExpr = cse(allDiffs, symbols=partialExpression)

        for v, expr in enumerate(substExpr):
            var = diffVars[v]
            if type(var) == State:
                adj_indices[0].append(var.id)
                cEvalDiff[0].append(ccode(expr))
            elif type(var) == Input:
                adj_indices[1].append(var.id)
                cEvalDiff[1].append(ccode(expr))
            elif type(var) == Parameter:
                adj_indices[2].append(var.id)
                cEvalDiff[2].append(ccode(expr))
        
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
        
        for s in subst0:
            out += "        const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
            
        out += f"\t\treturn {cEval};\n"
        out += "\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {{\n"
        
        for s in subst1:
            out += "        const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
            
        out += "\t\treturn {std::vector<double>"
        out += "{{{}}}, ".format(", ".join(cEvalDiff[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff[1]))
        out += "{{{}}}".format(", ".join(cEvalDiff[2]))
        out += "};\n\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {{\n"
        for s in subst2:
            out += "        const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
            
        out += "\t\treturn {std::vector<double>"
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[1]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[2]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[3]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff2[4]))
        out += "{{{}}}".format(", ".join(cEvalDiff2[5]))
        out += "};\n\t}\n"
        
        out += "private:\n"
        out += f"\t{name}(Adjacency adj, AdjacencyDiff adjDiff, double lb, double ub) : Constraint(std::move(adj), std::move(adjDiff), lb, ub) {{}}\n"
        out += "};\n\n\n"
        return out



class ParametricConstraint(Expression):
    
    def __init__(self, expr, lb=-float("inf"), ub=float("inf"), eq=None):
        super().__init__(expr)
        if eq != None:
            self.lb = eq
            self.ub = eq
        else:
            self.lb = lb
            self.ub = ub

        
    def codegen(self, name):
        partialExpression = numbered_symbols(prefix='s')
        
        subst0, substExpr = cse(self.expr, partialExpression)
        cEval = ccode(substExpr[0])
        
        # Generate second-order derivatives only for lower triangular part
        cEvalDiff2 = []
        adjDiff_indices = []  # indPP
        allDiffs2 = []
        
        for i, var1 in enumerate(self.adj):
            for j, var2 in enumerate(self.adj):
                if type(var1) == Parameter and type(var2) == Parameter and i >= j:
                    der = diff(self.expr, var1, var2)
                    if der != 0:
                        allDiffs2.append(der)
                        adjDiff_indices.append((var1.id, var2.id))  # indPP
                        
        subst2, substExpr2 = cse(allDiffs2, partialExpression)
        
        for expr in substExpr2:
            cEvalDiff2.append(ccode(expr))
        
        adj_indices = []  # indX, indU, indP
        cEvalDiff = []
        allDiffs = []
        
        for i, var in enumerate(self.adj):
            if type(var) == Parameter:
                der = diff(self.expr, var)
                if der != 0:
                    adj_indices.append(var.id)
                    allDiffs.append(der)
        
        subst1, substExpr = cse(allDiffs, partialExpression)
        
        for expr in substExpr:
            cEvalDiff.append(ccode(expr))
        
        adj = f"{{{', '.join(map(str, adj_indices))}}}"
        adjDiff = "{{{}}}".format(", ".join("{{{}}}".format(", ".join(map(str, tpl))) for tpl in adjDiff_indices))
        
        
        lb = self.lb if self.lb != -float('inf') else "MINUS_INFINITY"
        ub = self.ub if self.ub != float('inf') else "PLUS_INFINITY"
           
        out = f"class {name} : public ParamConstraint {{\n"
        out += "public:\n"
        out += f"\tstatic std::unique_ptr<{name}> create() {{\n"
        out += f"\t\tParamAdjacency adj{{{adj}}};\n"
        out += f"\t\tParamAdjacencyDiff adjDiff{{{adjDiff}}};\n"
        out += f"\t\treturn std::unique_ptr<{name}>(new {name}(std::move(adj), std::move(adjDiff), {lb}, {ub}));\n"
        out += "\t}\n\n"
        
        out += "\tdouble eval(const double* p) override {\n"
        
        for s in subst0:
            out += "        const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
                
        out += f"\t\treturn {cEval};\n"
        out += "\t}\n\n"
        
        out += f"\tstd::vector<double> evalDiff(const double* p) override {{\n"
        
        for s in subst1:
            out += "        const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
            
        out += "\t\treturn std::vector<double>"
        out += f"{{{', '.join(cEvalDiff)}}}"
        out += ";\n\t}\n\n"
        
        out += f"\tstd::vector<double> evalDiff2(const double* p) override {{\n"
        
        for s in subst2:
            out += "        const double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
            
        out += "\t\treturn std::vector<double>"
        out += "{{{}}}".format(", ".join(cEvalDiff2))
        out += ";\n\t}\n"
        
        out += "private:\n"
        out += f"\t{name}(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {{}}\n"
        out += "};\n\n\n"
        return out


class Model:
    
    def __init__(self, name):
        self.creationTime = timer.process_time()
        self.xVars = []
        self.uVars = []
        self.pVars = []
        self.rpVars = []
        self.M = None
        self.L = None
        self.F = []
        self.G = []
        self.R = []
        self.A = []
        self.name = name
        
        # additional stuff for running the model
        self.tf = None
        self.steps = None
        self.rksteps = None
        self.outputFilePath = None
        self.linearSolver = LinearSolver.MUMPS
        self.meshAlgorithm = MeshAlgorithm.NONE
        self.meshIterations = 0
        self.tolerance = 1e-14
        self.exportHessianPath = None
        self.exportJacobianPath = None
        self.meshLevel = None
        self.meshCTol = None
        self.meshSigma = None
        
        
    def addVar(self, variable):
        
        # adds any variable to the model, must be assigned previously
        
        if type(variable) == State:
            self.xVars.append(variable)
        elif type(variable) == Input:
            self.uVars.append(variable)
        elif type(variable) == Parameter:
            self.pVars.append(variable)
        elif type(variable) == RuntimeParameter:
            self.rpVars.append(variable)
        return variable
    
    def addState(self, start, symbol=None, lb=-float("inf"), ub=float("inf")):
        
        # adds a state x for optimization, must a have starting value
        # specify lb or ub if needed
        
        variable = State(start, symbol=symbol, lb=lb, ub=ub)
        self.xVars.append(variable)
        return variable
    
    def addInput(self, symbol=None, lb=-float("inf"), ub=float("inf")):
        
        # adds an input / control u for optimization
        # specify lb or ub if needed
        
        variable = Input(symbol=symbol, lb=lb, ub=ub)
        self.uVars.append(variable)
        return variable
        
    def addParameter(self, symbol=None, lb=-float("inf"), ub=float("inf")):
        
        # adds a parameter p for optimization
        # specify lb or ub if needed
        
        variable = Parameter(symbol=symbol, lb=lb, ub=ub)
        self.pVars.append(variable)
        return variable
    
    def addRuntimeParameter(self, default, symbol=None):
        
        # adds a runtime parameter: can be changed for optimization
        # will be handled symbolically -> thus quite slow
        # for actual constants simply write pythonic var = value
        
        variable = RuntimeParameter(default, symbol=symbol)
        self.rpVars.append(variable)
        return variable

    def addMayer(self, expr, obj=Objective.MINIMIZE):
        
        # adds the mayer term: min/max expr(tf)
        
        if self.M:
            raise InvalidModel("Mayer term already set")
        else:
            if obj == Objective.MAXIMIZE:
                self.M = Expression(-1 * expr)
            else:
                self.M = Expression(expr)
                
    def addLagrange(self, expr, obj=Objective.MINIMIZE):
        
        # adds the lagrange term: min/max integral_0_tf expr dt
        
        if self.L:
            raise InvalidModel("Lagrange term already set")
        else:
            if obj == Objective.MAXIMIZE:
                self.L = Expression(-1 * expr)
            else:
                self.L = Expression(expr)
    
    def addDynamic(self, diffVar, expr):
        
        # adds a dynamic constraint: diffVar' = expr
        
        for f in self.F:
            if diffVar == f.diffVar:
                fail = f"Equation for diff({diffVar}) has been added already"
                raise InvalidModel(fail)
        self.F.append(DynExpression(diffVar, expr))
    
    def addPath(self, expr, lb=-float("inf"), ub=float("inf"), eq=None):
        
        # adds a path constraint: lb <= g(.(t)) <= ub or g(.(t)) == eq
        
        if eq != None and (lb != -float("inf") or ub != float("inf")):
            raise InvalidModel("Can't set eq and lb or ub.")
        self.G.append(Constraint(expr, lb=lb, ub=ub, eq=eq))
    
    def addFinal(self, expr, lb=-float("inf"), ub=float("inf"), eq=None):
        
        # adds a final constraint: lb <= r(.(tf)) <= ub or r(.(tf)) == eq
        
        if eq != None and (lb != -float("inf") or ub != float("inf")):
            raise InvalidModel("Can't set eq and lb or ub.")
        self.R.append(Constraint(expr, lb=lb, ub=ub, eq=eq))
        
    def addParametric(self, expr, lb=-float("inf"), ub=float("inf"), eq=None):
        
        # adds a parametric constraint: lb <= a(p) <= ub or a(p) == eq
        
        if set(self.pVars).issuperset(expr.free_symbols):
            self.A.append(ParametricConstraint(expr, lb=lb, ub=ub, eq=eq))
        else:
            raise InvalidModel("Parametric constraints only allow parametric variables")
            
    def _addDummy(self):
        
        # adds a dummy state for purely parametric models
        
        x = self.addState(start=0)
        self.addF(x, 0)
    
    def setFinalTime(self, tf):
        self.tf = tf
    
    def setSteps(self, steps):
        self.steps = steps
    
    def setRkSteps(self, rksteps):
        self.rksteps = rksteps
    
    def setOutputPath(self, path):
        self.outputFilePath = path
    
    def setLinearSolver(self, solver):
        self.linearSolver = solver

    def setTolerance(self, tolerance):
        self.tolerance = tolerance
    
    def setExportHessianPath(self, path):
        self.exportHessianPath = path
        
    def setExportJacobianPath(self, path):
        self.exportJacobianPath = path
    
    def setMeshAlgorithm(self, meshAlgorithm):
        self.meshAlgorithm = meshAlgorithm
    
    def setMeshIterations(self, meshIterations):
        self.meshIterations = meshIterations
    
    def setMeshLevel(self, meshLevel):
        self.meshLevel = meshLevel
    
    def setMeshCTol(self, meshCTol):
        self.meshCTol = meshCTol
    
    def setMeshSigma(self, meshSigma):
        self.meshSigma = meshSigma

    def setFlags(self, flags):
        if "outputPath" in flags:
            self.setOutputPath(flags["outputPath"])
        if "linearSolver" in flags:
            self.setLinearSolver(flags["linearSolver"])
        if "tolerance" in flags:
            self.setTolerance(flags["tolerance"])
        if "exportHessianPath" in flags:
            self.setExportHessianPath(flags["exportHessianPath"])
        if "exportJacobianPath" in flags:
            self.setExportJacobianPath(flags["exportJacobianPath"])

    def setMeshFlags(self, meshFlags):
        if "meshAlgorithm" in meshFlags:
            self.setMeshAlgorithm(meshFlags["meshAlgorithm"])
        if "meshIterations" in meshFlags:
            self.setMeshIterations(meshFlags["meshIterations"])
        if "meshLevel" in meshFlags:
            self.setMeshLevel(meshFlags["meshLevel"])
        if "meshCTol" in meshFlags:
            self.setMeshCTol(meshFlags["meshCTol"])
        if "meshSigma" in meshFlags:
            self.setMeshSigma(meshFlags["meshSigma"])

    def generate(self):
        
        # does the entire code generation of the model
        
        if len(self.F) != len(self.xVars):
            raise InvalidModel("#states != #differential equations") 
        elif len(self.xVars) == 0:
            # adding a dummy var: x(0)=0, x'=0 for purely parametric models
            self._addDummy()
            
        # preparations
        
        self.F.sort(key=lambda eq: eq.diffVar.id, reverse=False)
        
        # codegen
        filename = self.name + "Generated"

        print("Starting .cpp codegen ...\n")
        
        OUTPUT = f'''
// CODEGEN FOR MODEL "{self.name}"\n
// includes
#define _USE_MATH_DEFINES
#include "{filename}Params.h"
#include <cmath>
#include <string>
#include "constants.h"
#include <problem.h>
#include "integrator.h"
#include "mesh.h"
#include "gdop.h"
#include "solver.h"
\n\n'''
        
        OUTPUT += "// runtime parameters\n"
        for rp in self.rpVars:
            OUTPUT += f"const double {rp.symbol} = {str(rp.symbol).upper()}_VALUE;\n"
        else:
            OUTPUT += "\n\n"
        
        if self.M:
            OUTPUT += "// mayer term\n"
            OUTPUT += self.M.codegen("Mayer" + self.name)
            print("Mayer: codegen done.\n")

        if self.L:
            OUTPUT += "// lagrange term\n"
            OUTPUT += self.L.codegen("Lagrange" + self.name)
            print("Lagrange: codegen done.\n")
        
        if self.F != []:
            OUTPUT += "// dynamic constraints\n"
        
        for n, f in enumerate(self.F):
            OUTPUT += f.codegen("F" + str(n) + self.name)
            print(f"Dynamic constraint {n}: codegen done.\n")
            
        if self.G != []:
            OUTPUT += "// path constraints\n"
            
        for n, g in enumerate(self.G):
            OUTPUT += g.codegen("G" + str(n) + self.name)
            print(f"Path constraint {n}: codegen done.\n")
        
        if self.R != []:
            OUTPUT += "// final constraints\n"
        
        for n, r in enumerate(self.R):
            OUTPUT += r.codegen("R" + str(n) + self.name)
            print(f"Final constraint {n}: codegen done.\n")
        
        if self.A != []:
            OUTPUT += "// parametric constraints\n"
        
        for n, a in enumerate(self.A):
            OUTPUT += a.codegen("A" + str(n) + self.name)
            print(f"Parametric constraints {n}: codegen done.\n")
            
        pushF = "\n    ".join("F.push_back(" + "F" + str(n) + self.name + "::create());" for n in range(len(self.F)))
        pushG = "\n    ".join("G.push_back(" + "G" + str(n) + self.name + "::create());" for n in range(len(self.G)))
        pushR = "\n    ".join("R.push_back(" + "R" + str(n) + self.name + "::create());" for n in range(len(self.R)))
        pushA = "\n    ".join("A.push_back(" + "A" + str(n) + self.name + "::create());" for n in range(len(self.A)))
        
        OUTPUT += f"""Problem createProblem_{self.name}() {{

    std::vector<std::unique_ptr<Expression>> F;
    {pushF}
    
    std::vector<std::unique_ptr<Constraint>> G;
    {pushG}
    
    std::vector<std::unique_ptr<Constraint>> R;
    {pushR}
    
    std::vector<std::unique_ptr<ParamConstraint>> A;
    {pushA}

    Problem problem(
            {len(self.xVars)}, {len(self.uVars)}, {len(self.pVars)},  // #vars
            {{{', '.join(str(x.start) for x in self.xVars)}}},  // x0
            {{{', '.join(str(x.lb if x.lb != -float('inf') else "MINUS_INFINITY") for x in self.xVars)}}},  // lb x
            {{{', '.join(str(x.ub if x.ub != float('inf') else "PLUS_INFINITY") for x in self.xVars)}}},  // ub x
            {{{', '.join(str(u.lb if u.lb != -float('inf') else "MINUS_INFINITY") for u in self.uVars)}}},  // lb u
            {{{', '.join(str(u.ub if u.ub != float('inf') else "PLUS_INFINITY") for u in self.uVars)}}},  // ub u
            {{{', '.join(str(p.lb if p.lb != -float('inf') else "MINUS_INFINITY") for p in self.pVars)}}},  // lb p
            {{{', '.join(str(p.ub if p.ub != float('inf') else "PLUS_INFINITY") for p in self.pVars)}}},  // ub p
            {"Mayer" + self.name + "::create()" if self.M else "{}"},
            {"Lagrange" + self.name + "::create()" if self.L else "{}"},
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "{self.name}");
    return problem;
}};\n"""
        OUTPUT += f"""
int main() {{
    auto problem = std::make_shared<const Problem>(createProblem_{self.name}());
    InitVars initVars = INIT_VARS;
    Integrator rk = Integrator::radauIIA(RADAU_INTEGRATOR);
    Mesh mesh = Mesh::createEquidistantMesh(INTERVALS, FINAL_TIME);
    LinearSolver linearSolver = LINEAR_SOLVER;
    MeshAlgorithm meshAlgorithm = MESH_ALGORITHM;
    int meshIterations = MESH_ITERATIONS;

    Solver solver = Solver(create_gdop(problem, mesh, rk, initVars), meshIterations, linearSolver, meshAlgorithm);

    // set solver flags
    #ifdef EXPORT_OPTIMUM_PATH
    solver.setExportOptimumPath(EXPORT_OPTIMUM_PATH);
    #endif
    
    #ifdef EXPORT_HESSIAN_PATH
    solver.setExportHessianPath(EXPORT_HESSIAN_PATH);
    #endif
    
    #ifdef EXPORT_JACOBIAN_PATH
    solver.setExportJacobianPath(EXPORT_JACOBIAN_PATH);
    #endif
    
    #ifdef TOLERANCE
    solver.setTolerance(TOLERANCE);
    #endif
    
    // set solver mesh parameters
    #ifdef LEVEL
    solver.setMeshParameter("level", LEVEL);
    #endif
    
    #ifdef C_TOL
    solver.setMeshParameter("ctol", C_TOL);
    #endif
    
    #ifdef SIGMA
    solver.setMeshParameter("sigma", SIGMA);
    #endif
    
    // optimize
    int status = solver.solve();
    return status;
}}        
        """
        print(".cpp: codegen done.\n")

        with open(f'{filename}.cpp', 'w') as file:
            file.write(OUTPUT)

        print(f"Generated model to {filename}Generated.cpp\n")
        print(f"Model creation, derivative calculations, and code generation took {round(timer.process_time() - self.creationTime, 4)} seconds.")
        return 0
        
    def optimize(self, tf=0, steps=1, rksteps=1, flags={}, meshFlags={}):
        
        # generate corresponding main function with flags, mesh, refinement
        # set runtime parameter file from map
        # run the code

        # always with setter to ensure some security
        
        self.setFinalTime(tf)
        self.setSteps(steps)
        self.setRkSteps(rksteps)

        self.setFlags(flags)

        self.setMeshFlags(meshFlags)

        ### main codegen
        filename = self.name + "Generated"
        OUTPUT = "//defines\n\n"
        OUTPUT += "#define INIT_VARS InitVars::CONST\n"
        OUTPUT += f"#define RADAU_INTEGRATOR IntegratorSteps::Steps{self.rksteps}\n"
        OUTPUT += f"#define INTERVALS {self.steps}\n"
        OUTPUT += f"#define FINAL_TIME {self.tf}\n"
        OUTPUT += f"#define LINEAR_SOLVER LinearSolver::{self.linearSolver.name}\n"
        OUTPUT += f"#define MESH_ALGORITHM MeshAlgorithm::{self.meshAlgorithm.name}\n"
        OUTPUT += f"#define MESH_ITERATIONS {self.meshIterations}\n"
        OUTPUT += f"#define TOLERANCE {self.tolerance}\n"

        if self.outputFilePath:
            OUTPUT += f'#define EXPORT_OPTIMUM_PATH "{self.outputFilePath}"\n'
        if self.exportHessianPath:
            OUTPUT += f'#define EXPORT_HESSIAN_PATH "{self.exportHessianPath}"\n'
        if self.exportJacobianPath:
            OUTPUT += f'#define EXPORT_JACOBIAN_PATH "{self.exportJacobianPath}"\n'

        if self.meshSigma:
            OUTPUT += f'#define SIGMA {self.meshSigma}\n'

        if self.meshLevel:
            OUTPUT += f'#define LEVEL {self.meshLevel}\n'

        if self.meshCTol:
            OUTPUT += f'#define C_TOL {self.meshCTol}\n'

        if self.rpVars:
            OUTPUT += "\n// values for runtime parameters\n"

        for rp in self.rpVars:
            OUTPUT += f'#define {str(rp.symbol).upper()}_VALUE {rp.value}\n'

        with open(f'{filename}Params.h', 'w') as file:
            file.write(OUTPUT)

        os.system(f"g++ {filename}.cpp -O3 -I../../src/ -L../../cmake-build-release/src -lipopt_do -o{self.name}") # vorher src lipopt_do

        os.system(f"LD_LIBRARY_PATH=../../cmake-build-release/src/ ./{self.name}")

        return 0
    
    def plot(self, meshIteration=None, interval=None, specifCols=None, dots=False):
        if meshIteration is None:
            meshIteration = self.meshIterations
        if interval is None:
            interval = [0, self.tf]
        df = pd.read_csv(self.outputFilePath + "/" + self.name + str(meshIteration) + ".csv", sep=",")
        plt.rcParams.update({
    'font.serif': ['Times New Roman'],
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'legend.frameon': True,
    'legend.loc': 'best',
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'grid.linewidth': 0.5,
    'figure.figsize': (12, 8),
    'axes.titlepad': 20,
    'axes.labelpad': 10
        })

        if specifCols is None:
            columns_to_plot = df.columns[1:]
        else:
            columns_to_plot = specifCols

        num_plots = len(columns_to_plot)
        fig, axs = plt.subplots(num_plots, 1, figsize=(12, 8 * num_plots), sharex=True)

        if num_plots == 1:
            axs = [axs]

        for idx, column in enumerate(columns_to_plot):
            ax = axs[idx]
            ax.plot(df['time'], df[column], label=column, linewidth=2, linestyle='-', color='steelblue')
            if dots:
                ax.scatter(df['time'], df[column], color='red', s=30, edgecolor='black', alpha=0.8, zorder=5)
            ax.set_xlabel('time')
            ax.set_ylabel(column)
            ax.set_xlim(interval[0], interval[1])
            ax.legend(frameon=True, loc='best')
            ax.grid(True)
            ax.title.set_fontsize(16)

        plt.tight_layout()
        plt.show()

    def getResults(self, meshIteration=0):
        return pd.read_csv(self.outputFilePath + "/" + self.name + str(meshIteration) + ".csv", sep=",")

    def printResults(self, meshIteration=0):
        df = self.getResults(meshIteration)
        print("")
        print(df)

### GLOBAL ALIAS AND GLOBAL VAR DEFINITIONS 

# variables
Continous = Input
Control = Input
Model.addControl = Model.addInput
Model.addContinous = Model.addInput
Model.addX = Model.addState
Model.addU = Model.addInput
Model.addP = Model.addParameter
Model.addRP = Model.addRuntimeParameter

Model.addV = Model.addVar

# objective
Model.addM = Model.addMayer
Model.addL = Model.addLagrange

# constraints
Model.addOde = Model.addDynamic
Model.addF = Model.addDynamic
Model.addG = Model.addPath
Model.addR = Model.addFinal
Model.addA = Model.addParametric

# time symbol
t = Symbol("t")
time = t

# consts
PI = 3.14159265358979323846
pi = PI

"""
SET THIS AS GEANY EXECUTE
export PYTHONPATH=/home/linus/masterarbeit/ipopt_do/codegen && python3 -u "%f"
"""
