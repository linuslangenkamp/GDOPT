from sympy.printing.c import C99CodePrinter

from optimization.variables import *

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
        cEval = toCode(exprEval[0])

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
                cEvalDiff2[0].append(toCode(expr))
                adjDiff_indices[0].append((var1.id, var2.id))  # indXX
            elif type(var1) == Input and type(var2) == State:
                cEvalDiff2[1].append(toCode(expr))
                adjDiff_indices[1].append((var1.id, var2.id))  # indUX
            elif type(var1) == Input and type(var2) == Input:
                cEvalDiff2[2].append(toCode(expr))
                adjDiff_indices[2].append((var1.id, var2.id))  # indUU
            elif type(var1) == Parameter and type(var2) == State:
                cEvalDiff2[3].append(toCode(expr))
                adjDiff_indices[3].append((var1.id, var2.id))  # indPX
            elif type(var1) == Parameter and type(var2) == Input:
                cEvalDiff2[4].append(toCode(expr))
                adjDiff_indices[4].append((var1.id, var2.id))  # indPU
            elif type(var1) == Parameter and type(var2) == Parameter:
                cEvalDiff2[5].append(toCode(expr))
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
                cEvalDiff[0].append(toCode(expr))
            elif type(var) == Input:
                adj_indices[1].append(var.id)
                cEvalDiff[1].append(toCode(expr))
            elif type(var) == Parameter:
                adj_indices[2].append(var.id)
                cEvalDiff[2].append(toCode(expr))

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
            out += "        const double " + toCode(s[0]) + " = " + toCode(s[1]) + ";\n"

        out += f"\t\treturn {cEval};\n"
        out += "\t}\n\n"

        out += f"\tstd::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {{\n"

        for s in subst1:
            out += "        const double " + toCode(s[0]) + " = " + toCode(s[1]) + ";\n"

        out += "\t\treturn {std::vector<double>"
        out += "{{{}}}, ".format(", ".join(cEvalDiff[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff[1]))
        out += "{{{}}}".format(", ".join(cEvalDiff[2]))
        out += "};\n\t}\n\n"

        out += f"\tstd::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {{\n"
        for s in subst2:
            out += "        const double " + toCode(s[0]) + " = " + toCode(s[1]) + ";\n"

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
        cEval = toCode(exprEval[0])

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
                cEvalDiff2[0].append(toCode(expr))
                adjDiff_indices[0].append((var1.id, var2.id))  # indXX
            elif type(var1) == Input and type(var2) == State:
                cEvalDiff2[1].append(toCode(expr))
                adjDiff_indices[1].append((var1.id, var2.id))  # indUX
            elif type(var1) == Input and type(var2) == Input:
                cEvalDiff2[2].append(toCode(expr))
                adjDiff_indices[2].append((var1.id, var2.id))  # indUU
            elif type(var1) == Parameter and type(var2) == State:
                cEvalDiff2[3].append(toCode(expr))
                adjDiff_indices[3].append((var1.id, var2.id))  # indPX
            elif type(var1) == Parameter and type(var2) == Input:
                cEvalDiff2[4].append(toCode(expr))
                adjDiff_indices[4].append((var1.id, var2.id))  # indPU
            elif type(var1) == Parameter and type(var2) == Parameter:
                cEvalDiff2[5].append(toCode(expr))
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
                cEvalDiff[0].append(toCode(expr))
            elif type(var) == Input:
                adj_indices[1].append(var.id)
                cEvalDiff[1].append(toCode(expr))
            elif type(var) == Parameter:
                adj_indices[2].append(var.id)
                cEvalDiff[2].append(toCode(expr))

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
            out += "        const double " + toCode(s[0]) + " = " + toCode(s[1]) + ";\n"

        out += f"\t\treturn {cEval};\n"
        out += "\t}\n\n"

        out += f"\tstd::array<std::vector<double>, 3> evalDiff(const double *x, const double *u, const double *p, double t) override {{\n"

        for s in subst1:
            out += "        const double " + toCode(s[0]) + " = " + toCode(s[1]) + ";\n"

        out += "\t\treturn {std::vector<double>"
        out += "{{{}}}, ".format(", ".join(cEvalDiff[0]))
        out += "{{{}}}, ".format(", ".join(cEvalDiff[1]))
        out += "{{{}}}".format(", ".join(cEvalDiff[2]))
        out += "};\n\t}\n\n"

        out += f"\tstd::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double t) override {{\n"
        for s in subst2:
            out += "        const double " + toCode(s[0]) + " = " + toCode(s[1]) + ";\n"

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
        cEval = toCode(substExpr[0])

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
            cEvalDiff2.append(toCode(expr))

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
            cEvalDiff.append(toCode(expr))

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
            out += "        const double " + toCode(s[0]) + " = " + toCode(s[1]) + ";\n"

        out += f"\t\treturn {cEval};\n"
        out += "\t}\n\n"

        out += f"\tstd::vector<double> evalDiff(const double* p) override {{\n"

        for s in subst1:
            out += "        const double " + toCode(s[0]) + " = " + toCode(s[1]) + ";\n"

        out += "\t\treturn std::vector<double>"
        out += f"{{{', '.join(cEvalDiff)}}}"
        out += ";\n\t}\n\n"

        out += f"\tstd::vector<double> evalDiff2(const double* p) override {{\n"

        for s in subst2:
            out += "        const double " + toCode(s[0]) + " = " + toCode(s[1]) + ";\n"

        out += "\t\treturn std::vector<double>"
        out += "{{{}}}".format(", ".join(cEvalDiff2))
        out += ";\n\t}\n"

        out += "private:\n"
        out += f"\t{name}(ParamAdjacency adj, ParamAdjacencyDiff adjDiff, double lb, double ub) : ParamConstraint(std::move(adj), std::move(adjDiff), lb, ub) {{}}\n"
        out += "};\n\n\n"
        return out


# custom ccode printer for const handling like pi
class CustomCCodePrinter(C99CodePrinter):
    def _print_Pi(self, expr):
        return 'M_PI'

# define global printer
printer = CustomCCodePrinter()

def toCode(expr):
    return printer.doprint(expr)
