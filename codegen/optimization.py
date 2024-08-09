from sympy import *
import time
from enum import Enum

# TODO: only diff w.r.t. to vars that are contained in a given expr
# check adj -> only diff then -> simplify -> check adj -> diff again -> simplify

class InvalidModel(Exception):
    pass


class Objective(Enum):
    MINIMIZE = 1
    MAXIMIZE = 2
    MAX = MAXIMIZE
    MIN = MINIMIZE


class Constant(Symbol):
    id_counter = 0
    
    def __new__(cls, value, symbol=None):
        if symbol == None:
            symbol = f"c{cls.id_counter}"
        obj = super().__new__(cls, symbol)
        obj.symbol = symbol
        obj.value = value
        cls.id_counter += 1
        return obj


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
        obj.symbol = symbol
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


class Expression:
    def __init__(self, expr):
        self.expr = simplify(expr)
    
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
        out += "\t\treturn {std::vector<double>"
        out += "{{{}}}, ".format(", ".join(cEvalDiff[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff[1]))
        out += "{{{}}}".format(", ".join(cEvalDiff[2]))
        out += "};\n\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {{\n"
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
        out += "\t\treturn {std::vector<double>"
        out += "{{{}}}, ".format(", ".join(cEvalDiff[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff[1]))
        out += "{{{}}}".format(", ".join(cEvalDiff[2]))
        out += "};\n\t}\n\n"
        
        out += f"\tstd::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {{\n"
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
        out += "\t\treturn std::vector<double>"
        out += f"{{{', '.join(cEvalDiff)}}}"
        out += ";\n\t}\n\n"
        
        out += f"\tstd::vector<double> evalDiff2(const double* p) override {{\n"
        out += "\t\treturn std::vector<double>"
        out += "{{{}}}".format(", ".join(cEvalDiff2))
        out += ";\n\t}\n"
        
        out += "private:\n"
        out += f"\t{name}(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {{}}\n"
        out += "};\n\n\n"
        return out


class Model:
    
    def __init__(self, name, constants=False):
        self.creationTime = time.process_time()
        self.useConstants = constants
        self.xVars = []
        self.uVars = []
        self.pVars = []
        self.constants = []
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
    
    def addConst(self, constant, symbol=None):
        if self.useConstants:
            c = Constant(constant, symbol)
            self.constants.append(c)
            return c
        else:
            return constant
    
    def addMayer(self, expr, obj=Objective.MINIMIZE):
        if self.M:
            raise InvalidModel("Mayer term already set")
        else:
            if obj == Objective.MAXIMIZE:
                self.M = Expression(-1 * expr)
            else:
                self.M = Expression(expr)
                
    def addLagrange(self, expr, obj=Objective.MINIMIZE):
        if self.L:
            raise InvalidModel("Lagrange term already set")
        else:
            if obj == Objective.MAXIMIZE:
                self.L = Expression(-1 * expr)
            else:
                self.L = Expression(expr)
    
    def addDynamic(self, diffVar, expr):
        for f in self.F:
            if diffVar == f.diffVar:
                fail = f"Equation for diff({diffVar}) has been added already"
                raise InvalidModel(fail)
        self.F.append(DynExpression(diffVar, expr))
    
    def addPath(self, expr, lb=-float("inf"), ub=float("inf")):
        self.G.append(Constraint(expr, lb=lb, ub=ub))
    
    def addFinal(self, expr, lb=-float("inf"), ub=float("inf")):
        self.R.append(Constraint(expr, lb=lb, ub=ub))
        
    def addParametric(self, expr, lb=-float("inf"), ub=float("inf")):
        if set(self.pVars).issuperset(expr.free_symbols):
            self.A.append(ParametricConstraint(expr, lb=lb, ub=ub))
        else:
            raise InvalidModel("Parametric constraints only allow parametric variables")
    
    def generate(self, filename):
        if len(self.F) != len(self.xVars):
            raise InvalidModel("#states != #differential equations") 
        
        # preparations
        
        self.F.sort(key=lambda eq: eq.diffVar.id, reverse=False)
        allVars = self.xVars + self.uVars + self.pVars
        
        # codegen
        
        HEADEROUTPUT = f"""
// CODEGEN FOR MODEL "{self.name}"\n
#ifndef IPOPT_DO_{self.name.upper()}_H
#define IPOPT_DO_{self.name.upper()}_H

#include <problem.h>

Problem createProblem_{self.name}();

#endif //IPOPT_DO_{self.name.upper()}_H
"""
        
        OUTPUT = f'''
// CODEGEN FOR MODEL "{self.name}"\n
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "{self.name}.h"
#include "constants.h"\n\n\n'''
        
        for constant in self.constants:
            OUTPUT += f"const double {constant.symbol} = {constant.value};\n"
        else:
            OUTPUT += "\n\n"
        
        if self.M:
            OUTPUT += self.M.codegen("Mayer" + self.name, allVars)
            
        if self.L:
            OUTPUT += self.L.codegen("Lagrange" + self.name, allVars)
        
        for n, f in enumerate(self.F):
            OUTPUT += f.codegen("F" + str(n) + self.name, allVars)
        
        for n, g in enumerate(self.G):
            OUTPUT += g.codegen("G" + str(n) + self.name, allVars)
        
        for n, r in enumerate(self.R):
            OUTPUT += r.codegen("R" + str(n) + self.name, allVars)
        
        for n, a in enumerate(self.A):
            OUTPUT += a.codegen("A" + str(n) + self.name, allVars)
        
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
        
        with open(f'{filename}.h', 'w') as file:
            file.write(HEADEROUTPUT)
        
        with open(f'{filename}.cpp', 'w') as file:
            file.write(OUTPUT)
        
        print(f"Generated model to {filename}.h and {filename}.cpp")
        print(f"Model creation, derivative calculations, and code generation took {round(time.process_time() - self.creationTime, 4)} seconds.")
        return 0
        
        
### GLOBAL VAR DEFINITIONS AND GLOBAL ALIAS

# variables
Continous = Input
Control = Input
Model.addC = Model.addConst
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

###

"""
SET THIS AS GEANY EXECUTE
export PYTHONPATH=/home/linus/masterarbeit/ipopt_do/codegen && python3 -u "%f"
"""
