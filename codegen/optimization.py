from sympy import *
import time as timer
from enum import Enum

# TODO: only diff if first diff != 0
# add vectorized eval of RHS = [f, g]^T, vectorized evalDiff, evalDiff2?
# or with colored jacobian
# TODO: integrate entire framework and make it more robust
# be careful to 


class InvalidModel(Exception):
    pass


class Objective(Enum):
    MINIMIZE = 1
    MAXIMIZE = 2
    MAX = MAXIMIZE
    MIN = MINIMIZE


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
        obj.default = default
        obj.symbol = symbol
        cls.id_counter += 1
        return obj
        
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
        for sym in expr.free_symbols:
            if type(sym) == State or type(sym) == Input or type(sym) == Parameter:
                self.adj.append(sym)
        try:
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
        
        print("Starting .h codegen ...\n")
        
        HEADEROUTPUT = f"""
// CODEGEN FOR MODEL "{self.name}"\n
#ifndef IPOPT_DO_{self.name.upper()}_H
#define IPOPT_DO_{self.name.upper()}_H

#include <problem.h>

Problem createProblem_{self.name}();

#endif //IPOPT_DO_{self.name.upper()}_H
"""
        
        print(".h: codegen done.\n")
        print("Starting .cpp codegen ...\n")
        
        OUTPUT = f'''
// CODEGEN FOR MODEL "{self.name}"\n
// includes
#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include "{filename}.h"
#include "constants.h"\n\n'''
        
        OUTPUT += "// runtime parameters\n"
        for rp in self.rpVars:
            OUTPUT += f"const double {rp.symbol} = {rp.default};\n"
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

        print(".cpp: codegen done.\n")
        
        with open(f'{filename}.h', 'w') as file:
            file.write(HEADEROUTPUT)
        
        with open(f'{filename}.cpp', 'w') as file:
            file.write(OUTPUT)
        
        print(f"Generated model to {filename}Generated.h and {filename}Generated.cpp\n")
        print(f"Model creation, derivative calculations, and code generation took {round(timer.process_time() - self.creationTime, 4)} seconds.")
        return 0
        
        
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
Model.addConstant = Model.addRuntimeParameter
Model.addC = Model.addRuntimeParameter

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
