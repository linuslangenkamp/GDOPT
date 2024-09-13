import os
from optimization.variables import *
from optimization.expressions import *
from optimization.structures import *
from optimization.radauHandling import *
from scipy.integrate import solve_ivp
import time as timer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# set global precision
pd.set_option("display.precision", 8)


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
        self.addedDummy = False

        # additional stuff for running the model
        self.tf = None
        self.steps = None
        self.rksteps = None
        self.outputFilePath = None
        self.ivpSolver = IVPSolver.Radau
        self.initVars = InitVars.SOLVE
        self.linearSolver = LinearSolver.MUMPS
        self.meshAlgorithm = MeshAlgorithm.NONE
        self.meshIterations = 0
        self.tolerance = 1e-14
        self.exportHessianPath = None
        self.exportJacobianPath = None
        self.initialStatesPath = "/tmp"
        self.meshLevel = None
        self.meshCTol = None
        self.meshSigma = None
        self.quadraticObjective = False
        self.linearObjective = False
        self.linearConstraints = False
        self.autoVariableNominals = False  # TODO: not in use: automatically scale variables and equations that dont have a nominal value
        self.userScaling = False

        # stuff for analyzing the results
        self.resultHistory = {}
        self.modelInfo = {}
        self.xVarNames = []
        self.uVarNames = []
        self.pVarNames = []
        self.rpVarNames = []

    def addState(self, start, symbol=None, lb=-float("inf"), ub=float("inf"), nominal=None):

        # adds a state x for optimization, must a have starting value
        # specify lb or ub if needed

        info = StateStruct(start, symbol=symbol, lb=lb, ub=ub, nominal=nominal)
        variable = Symbol(f"x[{info.id}]")
        varInfo[variable] = info
        self.xVars.append(variable)
        self.checkNominalNone(nominal)
        return variable

    def addInput(self, symbol=None, lb=-float("inf"), ub=float("inf"), guess=0, nominal=None):

        # adds an input / control u for optimization
        # specify lb or ub, start/initialGuess if needed
        # only provide a guess if you are certain that it will benefit

        info = InputStruct(symbol=symbol, lb=lb, ub=ub, initialGuess=guess, nominal=nominal)
        variable = Symbol(f"u[{info.id}]")
        varInfo[variable] = info
        self.uVars.append(variable)
        self.checkNominalNone(nominal)
        return variable

    def addBinaryInput(self, lb=0, ub=1, symbol=None, guess=None, nominal=None):

        # dangerous, might break the code
        # adds an input / control u for optimization
        # specify lb or ub, start/initialGuess if needed
        # to ensure the binary type, adds a parametric equation (u - lb) * (u - ub) == 0 to the model

        if guess is None:
            guess = (lb + ub) / 2
        info = InputStruct(symbol=symbol, lb=lb, ub=ub, initialGuess=guess, nominal=nominal)
        variable = Symbol(f"u[{info.id}]")
        varInfo[variable] = info
        self.uVars.append(variable)
        self.addPath(variable**2 - variable * (lb + ub) + lb * ub, eq=0)  # forces binary type
        self.checkNominalNone(nominal)
        return variable

    def addParameter(self, symbol=None, lb=-float("inf"), ub=float("inf"), guess=0, nominal=None):

        # adds a parameter p for optimization
        # specify lb or ub if needed, start/initialGuess if needed

        info = ParameterStruct(symbol=symbol, lb=lb, ub=ub, initialGuess=guess, nominal=nominal)
        variable = Symbol(f"p[{info.id}]")
        varInfo[variable] = info
        self.pVars.append(variable)
        self.checkNominalNone(nominal)
        return variable

    def addBinaryParameter(self, lb=0, ub=1, symbol=None, guess=None, nominal=None):

        # dangerous, might break the code
        # adds a binary parameter p for optimization
        # specify lb or ub, start/initialGuess if needed
        # to ensure the binary type, adds a parametric equation (p - lb) * (p - ub) == 0 to the model

        if guess is None:
            guess = (lb + ub) / 2
        info = ParameterStruct(symbol=symbol, lb=lb, ub=ub, initialGuess=guess, nominal=nominal)
        variable = Symbol(f"p[{info.id}]")
        varInfo[variable] = info
        self.pVars.append(variable)
        self.addParametric(variable**2 - variable * (lb + ub) + lb * ub, eq=0)  # forces binary type
        self.checkNominalNone(nominal)
        return variable

    def addRuntimeParameter(self, default, symbol):

        # adds a runtime parameter: can be changed for optimization
        # will be handled symbolically -> thus quite slow
        # for actual constants simply write pythonic var = value

        info = RuntimeParameterStruct(default=default, symbol=symbol)
        variable = Symbol(f"{str(info.symbol).upper().strip(' ')}_VALUE")
        varInfo[variable] = info
        self.rpVars.append(variable)
        return variable

    def setValue(self, runtimeParameter, value):

        # sets a runtime parameter value

        varInfo[runtimeParameter].value = value

    def addObjective(self, mayer, lagrange, obj=Objective.MINIMIZE, nominal=None):

        # adds the mayer and lagrange term
        # if mayer and lagrange have a nominal -> sum is used as nominal value

        self.addMayer(mayer, obj, nominal=nominal)
        self.addLagrange(lagrange, obj, nominal=nominal)

    def addMayer(self, expr, obj=Objective.MINIMIZE, nominal=None):

        # adds the mayer term: min/max expr(tf)
        # if mayer and lagrange have a nominal -> sum is used as nominal value

        if self.M:
            raise InvalidModel("Mayer term already set")
        else:
            if obj == Objective.MAXIMIZE:
                self.M = Expression(-1 * expr, nominal=nominal)
                print("Setting Mayer term as -1 * mayer, since maximization is chosen.\n")
            else:
                self.M = Expression(expr, nominal=nominal)
        self.checkNominalNone(nominal)

    def addLagrange(self, expr, obj=Objective.MINIMIZE, nominal=None):

        # adds the lagrange term: min/max integral_0_tf expr dt
        # if mayer and lagrange have a nominal -> sum is used as nominal value

        if self.L:
            raise InvalidModel("Lagrange term already set")
        else:
            if obj == Objective.MAXIMIZE:
                self.L = Expression(-1 * expr, nominal=nominal)
                print("Setting Lagrange term as -1 * lagrange, since maximization is chosen.\n")
            else:
                self.L = Expression(expr, nominal=nominal)
        self.checkNominalNone(nominal)

    def addDynamic(self, diffVar, expr, nominal=None):

        # adds a dynamic constraint: diffVar' = expr

        for f in self.F:
            if diffVar == f.diffVar:
                fail = f"Equation for diff({diffVar}) has been added already"
                raise InvalidModel(fail)
        self.F.append(DynExpression(diffVar, expr, nominal=nominal))
        self.checkNominalNone(nominal)

    def addPath(self, expr, lb=-float("inf"), ub=float("inf"), eq=None, nominal=None):

        # adds a path constraint: lb <= g(.(t)) <= ub or g(.(t)) == eq

        if eq is not None and (lb != -float("inf") or ub != float("inf")):
            raise InvalidModel("Can't set eq and lb or ub.")
        self.G.append(Constraint(expr, lb=lb, ub=ub, eq=eq, nominal=nominal))
        self.checkNominalNone(nominal)

    def addFinal(self, expr, lb=-float("inf"), ub=float("inf"), eq=None, nominal=None):

        # adds a final constraint: lb <= r(.(tf)) <= ub or r(.(tf)) == eq

        if eq is not None and (lb != -float("inf") or ub != float("inf")):
            raise InvalidModel("Can't set eq and lb or ub.")
        self.R.append(Constraint(expr, lb=lb, ub=ub, eq=eq, nominal=nominal))
        self.checkNominalNone(nominal)

    def addParametric(self, expr, lb=-float("inf"), ub=float("inf"), eq=None, nominal=None):

        # adds a parametric constraint: lb <= a(p) <= ub or a(p) == eq

        if set(self.pVars).issuperset(expr.free_symbols):
            self.A.append(ParametricConstraint(expr, lb=lb, ub=ub, eq=eq, nominal=nominal))
        else:
            raise InvalidModel("Parametric constraints only allow parametric variables")
        self.checkNominalNone(nominal)

    def _addDummy(self):

        # adds a dummy state for purely parametric models

        x = self.addState(start=0)
        self.addDynamic(x, 0)
        self.addedDummy = True

    def solveDynamic(self):

        # solves the dynamic system with the given initial u(t), p, x0

        # define the rhs, uCallback function
        rhs = [
            Lambdify([t] + [x for x in self.xVars] + [u for u in self.uVars] + [p for p in self.pVars] + [rp for rp in self.rpVars], self.F[eq].expr)
            for eq in range(len(self.F))
        ]
        uFuncs = [Lambdify([t], varInfo[u].initialGuess) for u in self.uVars]

        # define rhs as an actual function
        ode = lambda T, x: [
            r(
                *(
                    [T]
                    + [elem for elem in x]
                    + [uFuncs[i](T) for i in range(len(self.uVars))]
                    + [varInfo[pVar].initialGuess for pVar in self.pVars]
                    + [varInfo[rpVar].value for rpVar in self.rpVars]
                )
            )
            for r in rhs
        ]

        x0 = [varInfo[x].start for x in self.xVars]
        timeHorizon = [0, self.tf]
        # scipy.solve_ivp
        solution = solve_ivp(ode, timeHorizon, x0, method=self.ivpSolver.name, dense_output=True, rtol=1e-8)

        # get solution at the RadauIIA knots (based on self.rksteps)
        timeVals = generate_radau_knots(timeHorizon, self.steps, self.rksteps)
        stateVals = solution.sol(timeVals)
        return timeVals, stateVals

    def setFinalTime(self, tf: float):
        self.tf = tf

    def setSteps(self, steps: int):
        self.steps = steps

    def setRkSteps(self, rksteps: int):
        self.rksteps = rksteps

    def setOutputPath(self, path: str):
        self.outputFilePath = path

    def setLinearSolver(self, solver: LinearSolver):
        self.linearSolver = solver

    def setInitVars(self, initVars: InitVars):
        self.initVars = initVars

    def setTolerance(self, tolerance: float):
        self.tolerance = tolerance

    def setExportHessianPath(self, path: str):
        self.exportHessianPath = path

    def setExportJacobianPath(self, path: str):
        self.exportJacobianPath = path

    def setInitialStatesPath(self, path: str):
        self.initialStatesPath = path

    def setIVPSolver(self, solver: IVPSolver):

        # set IVP Solver for initial guess of states. (see scipy.solve_ivp)
        # default = IVPSolver.Radau (should be best), others: BDF, LSODA, RK45, DOP853, RK23

        self.ivpSolver = solver

    def setMeshAlgorithm(self, meshAlgorithm: MeshAlgorithm):
        self.meshAlgorithm = meshAlgorithm

    def setMeshIterations(self, meshIterations: int):
        self.meshIterations = meshIterations

    def setMeshLevel(self, meshLevel: float):
        self.meshLevel = meshLevel

    def setMeshCTol(self, meshCTol: float):
        self.meshCTol = meshCTol

    def setMeshSigma(self, meshSigma: float):
        self.meshSigma = meshSigma

    def setUserScaling(self, userScaling: bool):
        self.userScaling = userScaling

    def hasLinearObjective(self):

        # set this if both Mayer and Lagrange are linear

        self.linearObjective = True

    def hasQuadraticObjective(self):

        # set this if both Mayer and Lagrange are quadratic

        self.quadraticObjective = True

    def hasLinearConstraints(self):

        # set this if all constraints, i.e. f, g, r and a are linear

        self.linearConstraints = True

    def setExpressionSimplification(self, simp):

        # turn initial simplification of expression at generation on / off, standard = off, good for large models

        Expression.simplification = simp

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
        if "initialStatesPath" in flags:
            self.setInitialStatesPath(flags["initialStatesPath"])
        if "ivpSolver" in flags:
            self.setIVPSolver(flags["ivpSolver"])
        if "initVars" in flags:
            self.setInitVars(flags["initVars"])
        if "userScaling" in flags:
            self.setUserScaling(flags["userScaling"])

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

    def uInitialGuessCodegen(self):
        out = "std::vector<double> uInitialGuess(double t) {\n"
        out += f"\t return {{{', '.join(str(toCode(varInfo[u].initialGuess)) for u in self.uVars)}}};"
        out += "\n};\n\n"
        return out

    def generateProblemCondition(self):
        lines = [
            "\tproblem.linearObjective = true;" if self.linearObjective else "",
            "\tproblem.linearConstraints = true;" if self.linearConstraints else "",
            "\tproblem.quadraticObjLinearConstraints = true;" if (self.linearObjective or self.quadraticObjective) and self.linearConstraints else "",
        ]
        out = "\n".join(line for line in lines if line)
        if out != "":
            out = "\n" + out + "\n"
        return out

    def objectiveNominal(self):
        if self.M and self.L:
            if self.M.nominal is not None and self.L.nominal is not None:
                return self.M.nominal + self.L.nominal
            elif self.M.nominal is not None:
                return self.M.nominal
            elif self.L.nominal is not None:
                return self.L.nominal
            else:
                return 1
        elif self.M:
            return self.M.nominal if self.M.nominal is not None else 1
        elif self.L:
            return self.L.nominal if self.L.nominal is not None else 1
        else:
            return 1

    def checkNominalNone(self, nominal):
        # forces userScaling to True, if a single nominal is provided in the model -> else its False
        if nominal is not None:
            self.userScaling = True

    def initAnalysis(self):
        with open("/tmp/modelinfo.txt", "r") as file:
            for line in file:
                line = line.strip()
                key, value = line.split(",")
                key = key.strip()
                value = value.strip()
                self.modelInfo[key] = int(value)

    def generate(self):

        # does the entire code generation of the model

        if len(self.F) != len(self.xVars):
            raise InvalidModel("#states != #differential equations")
        elif len(self.xVars) == 0:
            # adding a dummy var: x(0)=0, x'=0 for purely parametric models
            self._addDummy()

        # preparations

        self.F.sort(key=lambda eq: varInfo[eq.diffVar].id, reverse=False)

        # codegen
        filename = self.name + "Generated"

        print("Starting .cpp codegen...\n")

        OUTPUT = f"""
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
\n\n"""

        OUTPUT += "// runtime parameters and global constants\n"
        for rp in self.rpVars:
            info = varInfo[rp]
            OUTPUT += f"const double {info.symbol} = {str(info.symbol).upper().strip(' ')}_VALUE;\n"
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

        if self.F:
            OUTPUT += "// dynamic constraints\n"

        for n, f in enumerate(self.F):
            OUTPUT += f.codegen("F" + str(n) + self.name)
            print(f"Dynamic constraint {n}: codegen done.\n")

        if self.G:
            OUTPUT += "// path constraints\n"

        for n, g in enumerate(self.G):
            OUTPUT += g.codegen("G" + str(n) + self.name)
            print(f"Path constraint {n}: codegen done.\n")

        if self.R:
            OUTPUT += "// final constraints\n"

        for n, r in enumerate(self.R):
            OUTPUT += r.codegen("R" + str(n) + self.name)
            print(f"Final constraint {n}: codegen done.\n")

        if self.A:
            OUTPUT += "// parametric constraints\n"

        for n, a in enumerate(self.A):
            OUTPUT += a.codegen("A" + str(n) + self.name)
            print(f"Parametric constraints {n}: codegen done.\n")

        OUTPUT += self.uInitialGuessCodegen()

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
            {{{', '.join(str(toCode(varInfo[x].start)) for x in self.xVars)}}},  // x0
            {{{', '.join(str(toCode(varInfo[x].lb) if varInfo[x].lb != -float('inf') else "MINUS_INFINITY") for x in self.xVars)}}},  // lb x
            {{{', '.join(str(toCode(varInfo[x].ub) if varInfo[x].ub != float('inf') else "PLUS_INFINITY") for x in self.xVars)}}},  // ub x
            &uInitialGuess,  // u0 initial guesses for optimization
            {{{', '.join(str(toCode(varInfo[u].lb) if varInfo[u].lb != -float('inf') else "MINUS_INFINITY") for u in self.uVars)}}},  // lb u
            {{{', '.join(str(toCode(varInfo[u].ub) if varInfo[u].ub != float('inf') else "PLUS_INFINITY") for u in self.uVars)}}},  // ub u
            {{{', '.join(str(toCode(varInfo[p].initialGuess)) for p in self.pVars)}}},  // p0 initial guesses for optimization
            {{{', '.join(str(toCode(varInfo[p].lb) if varInfo[p].lb != -float('inf') else "MINUS_INFINITY") for p in self.pVars)}}},  // lb p
            {{{', '.join(str(toCode(varInfo[p].ub) if varInfo[p].ub != float('inf') else "PLUS_INFINITY") for p in self.pVars)}}},  // ub p
            {"Mayer" + self.name + "::create()" if self.M else "{}"},
            {"Lagrange" + self.name + "::create()" if self.L else "{}"},
            std::move(F),
            std::move(G),
            std::move(R),
            std::move(A),
            "{self.name}");
    
    #ifdef USER_SCALING
    problem.nominalsX = {{{', '.join(str(toCode("1" if varInfo[x].nominal is None else varInfo[x].nominal)) for x in self.xVars)}}};  // nominal states
    problem.nominalsU = {{{', '.join(str(toCode("1" if varInfo[u].nominal is None else varInfo[u].nominal)) for u in self.uVars)}}};  // nominal inputs
    problem.nominalsP = {{{', '.join(str(toCode("1" if varInfo[p].nominal is None else varInfo[p].nominal)) for p in self.pVars)}}};  // nominal inputs
    
    problem.nominalObjective = {self.objectiveNominal()};  // nominal objective
    problem.nominalsF = {{{', '.join(str(toCode("1" if f.nominal is None else f.nominal)) for f in self.F)}}};  // nominal dynamic
    problem.nominalsG = {{{', '.join(str(toCode("1" if g.nominal is None else g.nominal)) for g in self.G)}}};  // nominal path
    problem.nominalsR = {{{', '.join(str(toCode("1" if r.nominal is None else r.nominal)) for r in self.R)}}};  // nominal final
    problem.nominalsA = {{{', '.join(str(toCode("1" if a.nominal is None else a.nominal)) for a in self.A)}}};  // nominal parametric
    #endif
    
    #ifdef INITIAL_STATES_PATH
    problem.initialStatesPath = INITIAL_STATES_PATH "/initialValues.csv";
    #endif
    {self.generateProblemCondition()}
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
    
    #ifdef USER_SCALING
    solver.userScaling = USER_SCALING;
    #endif
    
    // optimize
    int status = solver.solve();
    return status;
}}        
        """
        print(".cpp: codegen done.\n")

        os.makedirs(f".generated/{self.name}", exist_ok=True)

        with open(f".generated/{self.name}/{filename}.cpp", "w") as file:
            file.write(OUTPUT)

        print(f"Generated model to .generated/{self.name}/{filename}.cpp.\n")
        print(f"Model creation, derivative calculations, and code generation took {round(timer.process_time() - self.creationTime, 4)} seconds.\n")
        return 0

    def optimize(self, tf=1, steps=1, rksteps=1, flags={}, meshFlags={}, resimulate=False):

        # generate corresponding main function with flags, mesh, refinement
        # set runtime parameter file from map
        # run the code

        if not resimulate:  # use everything from the previous optimization
            if not self.addedDummy:  # use provided values / standard values
                self.setFinalTime(tf)
                self.setSteps(steps)
                self.setRkSteps(rksteps)
            else:  # purely parametric
                print("\nSetting tf = 0, steps = 1, rksteps = 1, since the model is purely parametric.")
                self.setFinalTime(0)
                self.setSteps(1)
                self.setRkSteps(1)

            self.setFlags(flags)

            self.setMeshFlags(meshFlags)
        else:
            self.resultHistory = {}  # clear result history for new optimization

        if self.initVars == InitVars.SOLVE:
            print("Solving IVP for initial state guesses...")

            solveStart = timer.time()
            timeVals, stateVals = self.solveDynamic()
            print(f"\nSolving the IVP took {round(timer.time() - solveStart, 4)} seconds.")

            print(f"\nWriting guesses to {self.initialStatesPath + '/initialValues.csv'}...")
            with open(self.initialStatesPath + "/initialValues.csv", "w") as file:
                for i in range(len(timeVals)):
                    row = [str(timeVals[i])]
                    for dim in range(len(self.xVars)):
                        row.append(str(stateVals[dim][i]))

                    file.write(",".join(row) + "\n")
            print("\nInitial guesses done.\n")

        ### main codegen
        filename = self.name + "Generated"
        OUTPUT = "//defines\n\n"
        OUTPUT += f"#define INIT_VARS InitVars::{self.initVars.name}\n"
        OUTPUT += f"#define RADAU_INTEGRATOR IntegratorSteps::Steps{self.rksteps}\n"
        OUTPUT += f"#define INTERVALS {self.steps}\n"
        OUTPUT += f"#define FINAL_TIME {self.tf}\n"
        OUTPUT += f"#define LINEAR_SOLVER LinearSolver::{self.linearSolver.name}\n"
        OUTPUT += f"#define MESH_ALGORITHM MeshAlgorithm::{self.meshAlgorithm.name}\n"
        OUTPUT += f"#define MESH_ITERATIONS {self.meshIterations}\n"
        OUTPUT += f"#define TOLERANCE {self.tolerance}\n"
        OUTPUT += f"#define USER_SCALING {'true' if self.userScaling else 'false'}\n"

        if self.outputFilePath:
            OUTPUT += f'#define EXPORT_OPTIMUM_PATH "{self.outputFilePath}"\n'
        if self.exportHessianPath:
            OUTPUT += f'#define EXPORT_HESSIAN_PATH "{self.exportHessianPath}"\n'
        if self.exportJacobianPath:
            OUTPUT += f'#define EXPORT_JACOBIAN_PATH "{self.exportJacobianPath}"\n'
        if self.initVars == InitVars.SOLVE and self.initialStatesPath:
            OUTPUT += f'#define INITIAL_STATES_PATH "{self.initialStatesPath}"\n'
        if self.meshSigma:
            OUTPUT += f"#define SIGMA {self.meshSigma}\n"
        if self.meshLevel:
            OUTPUT += f"#define LEVEL {self.meshLevel}\n"
        if self.meshCTol:
            OUTPUT += f"#define C_TOL {self.meshCTol}\n"
        if self.rpVars:
            OUTPUT += "\n// values for runtime parameters\n"

        for rp in self.rpVars:
            info = varInfo[rp]
            OUTPUT += f'#define {str(info.symbol).upper().strip(" ")}_VALUE {info.value}\n'

        with open(f".generated/{self.name}/{filename}Params.h", "w") as file:
            file.write(OUTPUT)

        print("Compiling generated code...\n")
        compileStart = timer.time()
        # TODO: investigate if -ffast-math is save here
        os.system(
            f"g++ .generated/{self.name}/{filename}.cpp -O3 -ffast-math -I../src/ -L../cmake-build-release/src -lipopt_do -o.generated/{self.name}/{self.name}"
        )
        print(f"Compiling to C++ took {round(timer.time() - compileStart, 4)} seconds.")

        os.system(f"LD_LIBRARY_PATH=../cmake-build-release/src/ ./.generated/{self.name}/{self.name}")

        self.initAnalysis()

        return 0

    def initVarNames(self):
        self.xVarNames = [info.symbol for variable, info in varInfo.items() if isinstance(info, StateStruct)]
        self.uVarNames = [info.symbol for variable, info in varInfo.items() if isinstance(info, InputStruct)]
        self.pVarNames = [info.symbol for variable, info in varInfo.items() if isinstance(info, ParameterStruct)]
        self.rpVarNames = [info.symbol for variable, info in varInfo.items() if isinstance(info, RuntimeParameterStruct)]

    def getResults(self, meshIteration=None):
        meshIteration = self.checkMeshIteration(meshIteration)

        if meshIteration not in self.resultHistory:
            try:
                results = pd.read_csv(self.outputFilePath + "/" + self.name + str(meshIteration) + ".csv", delimiter=",")
            except:
                raise Exception("meshIteration out of range. Set Model.maxMeshIteration to the maximum mesh iteration!")

            # remove dummy column for purely parametric models
            if self.addedDummy:
                results = results.drop(columns=["x[0]"])

            alias = {variable.name: varInfo[variable].symbol for variable in varInfo}
            for col in results.columns:
                if col not in alias:
                    alias[col] = col
            results.rename(columns=alias, inplace=True)
            self.resultHistory[meshIteration] = results
        self.initVarNames()
        return self.resultHistory[meshIteration]

    def printResults(self, meshIteration=None):
        if meshIteration is None:
            meshIteration = self.modelInfo["maxMeshIteration"]
        meshIteration = self.checkMeshIteration(meshIteration)
        self.getResults(meshIteration)
        meshIteration = self.modelInfo["maxMeshIteration"]
        print("")
        print(self.resultHistory[meshIteration])

    def printResultParameters(self, meshIteration=None):
        meshIteration = self.checkMeshIteration(meshIteration)
        self.getResults(meshIteration)
        print("")
        print(self.resultHistory[meshIteration][self.pVarNames].iloc[0].to_string())
        print("")

    def checkMeshIteration(self, meshIteration):
        maxMeshIteration = self.modelInfo["maxMeshIteration"]
        if meshIteration is None:
            return maxMeshIteration
        if type(meshIteration) == int:
            if meshIteration > maxMeshIteration:
                print(f"meshIteration too large. Setting meshIteration to maximum value of {maxMeshIteration}.")
                meshIteration = maxMeshIteration
        return meshIteration

    def setPlotDefaults(self):
        plt.rcParams.update(
            {
                "font.serif": ["Times New Roman"],
                "axes.labelsize": 14,
                "axes.titlesize": 16,
                "xtick.labelsize": 12,
                "ytick.labelsize": 12,
                "legend.fontsize": 12,
                "legend.frameon": True,
                "legend.loc": "best",
                "grid.alpha": 0.3,
                "grid.linestyle": "--",
                "grid.linewidth": 0.5,
                "figure.figsize": (12, 8),
                "axes.titlepad": 20,
                "axes.labelpad": 10,
            }
        )

    def plotStates(self, meshIteration=None, interval=None, dots=Dots.OFF):
        self.initVarNames()
        self.plotGeneral(meshIteration=meshIteration, interval=interval, dots=dots, specifCols=self.xVarNames)

    def plot(self, specifCols=None, meshIteration=None, interval=None, dots=Dots.OFF):
        self.initVarNames()
        self.plotGeneral(meshIteration=meshIteration, interval=interval, dots=dots, specifCols=specifCols)

    def plotInputs(self, meshIteration=None, interval=None, dots=Dots.OFF):
        self.initVarNames()
        self.plotGeneral(meshIteration=meshIteration, interval=interval, dots=dots, specifCols=self.uVarNames)

    def parametricPlot(self, varX, varY, meshIteration=None, interval=None, dots=Dots.OFF):
        self.initVarNames()
        if interval is None:
            interval = [0, self.tf]
        self.plotParametric(varX=varX, varY=varY, meshIteration=meshIteration, interval=interval, dots=dots)

    def plotInputsAndRefinement(self, meshIteration=None, interval=None, markerSize=30, dotsMesh=Dots.BASE, dotsGraph=Dots.OFF, epsilon=1e-14):
        self.initVarNames()
        self.plotColsAndRefinement(meshIteration=meshIteration, interval=interval, specifCols=self.uVarNames, markerSize=markerSize, dotsMesh=dotsMesh, dotsGraph=dotsGraph, epsilon=epsilon)

    def plotColsAndRefinement(self, meshIteration=None, interval=None, specifCols=None, markerSize=30, dotsMesh=Dots.BASE, dotsGraph=Dots.OFF, epsilon=1e-14):
        if interval is None:
            interval = [0, self.tf]

        figMesh, axMesh = self._plotMeshRefinement(interval=interval, markerSize=markerSize, dots=dotsMesh, epsilon=epsilon)
        figGraph, axsGraph = self._plotGeneral(meshIteration=meshIteration, interval=interval, specifCols=specifCols, dots=dotsGraph)

        plt.close(figGraph)
        plt.close(figMesh)

        fig, axs = plt.subplots(len(axsGraph) + 1, 1, figsize=(12, 10))

        for i in range(len(axsGraph)):
            axs[i].set_xticks([])
            axs[i].plot(axsGraph[i].lines[0].get_xdata(), axsGraph[i].lines[0].get_ydata())
            axs[i].set_ylabel(axsGraph[i].get_ylabel())
            axs[i].set_xlim(interval)
            if axsGraph[i].collections:
                for coll in axsGraph[i].collections:
                    offsets = coll.get_offsets()
                    axs[i].scatter(offsets[:, 0], offsets[:, 1], color="red", s=30, edgecolor="black", alpha=0.8, zorder=5)

        if axMesh.collections:
            for coll in axMesh.collections:
                offsets = coll.get_offsets()
                axs[len(axs) - 1].scatter(offsets[:, 0], offsets[:, 1], color="red", s=markerSize, edgecolor="black", alpha=0.8)
            axs[len(axs) - 1].set_ylabel(axMesh.get_ylabel())
            axs[len(axs) - 1].set_xlim(interval)

        plt.tight_layout()
        plt.subplots_adjust(left=0.075, right=0.95, top=0.925, bottom=0.075, hspace=0.1)
        plt.show()

    def plotMeshRefinement(self, interval=None, markerSize=30, dots=Dots.BASE, epsilon=1e-14):
        fig, ax = self._plotMeshRefinement(interval=interval, markerSize=markerSize, dots=dots, epsilon=epsilon)
        ax.set_xlabel("Time")
        ax.set_title("Inserted Mesh Points Over Time")
        plt.tight_layout()
        plt.subplots_adjust(left=0.075, right=0.95, top=0.925, bottom=0.075, hspace=0.1)
        plt.show()

    def plotGeneral(self, meshIteration=None, interval=None, specifCols=None, dots=Dots.OFF):
        fig, axs = self._plotGeneral(meshIteration=meshIteration, interval=interval, specifCols=specifCols, dots=dots)
        axs[-1].set_xlabel("Time")
        axs[0].set_title(f"Optimal Solution: {self.name}")
        plt.tight_layout()
        plt.subplots_adjust(left=0.075, right=0.95, top=0.925, bottom=0.075, hspace=0.1)
        plt.show()

    def _plotGeneral(self, meshIteration=None, interval=None, specifCols=None, dots=Dots.OFF):
        self.setPlotDefaults()

        meshIteration = self.checkMeshIteration(meshIteration)
        if interval is None:
            interval = [0, self.tf]
        self.getResults(meshIteration=meshIteration)

        columns_to_plot = specifCols or self.resultHistory[meshIteration].columns[1:]

        fig, axs = plt.subplots(len(columns_to_plot), 1, figsize=(12, 8 * len(columns_to_plot)), sharex=True)

        if len(columns_to_plot) == 1:
            axs = [axs]

        for idx, column in enumerate(columns_to_plot):
            ax = axs[idx]
            ax.plot(
                self.resultHistory[meshIteration]["time"],
                self.resultHistory[meshIteration][column],
                label=column,
                linewidth=2,
                linestyle="-",
                color="steelblue",
            )
            self._applyDots(ax, dots, meshIteration, column)
            ax.set_ylabel(column)
            ax.set_xlim(interval)
            ax.legend(frameon=True, loc="best")

        return fig, axs

    def plotParametric(self, varX, varY, meshIteration=None, interval=None, dots=Dots.OFF):
        self.setPlotDefaults()

        meshIteration = self.checkMeshIteration(meshIteration)
        if interval is None:
            interval = [0, self.tf]

        self.getResults(meshIteration=meshIteration)

        x_data = self.resultHistory[meshIteration][varInfo[varX].symbol]
        y_data = self.resultHistory[meshIteration][varInfo[varY].symbol]
        time_data = self.resultHistory[meshIteration]["time"]
        time_filtered = (time_data >= interval[0]) & (time_data <= interval[1])

        plt.plot(x_data[time_filtered], y_data[time_filtered])

        if dots != Dots.OFF:
            self._applyDotsParametric(x_data[time_filtered], y_data[time_filtered], dots)

        plt.title(f"{varInfo[varX].symbol} vs {varInfo[varY].symbol}")
        plt.xlabel(varInfo[varX].symbol)
        plt.ylabel(varInfo[varY].symbol)
        plt.tight_layout()
        plt.show()

    # Helper function to apply dots in regular plots
    def _applyDots(self, ax, dots, meshIteration, column):
        if dots == Dots.ALL:
            ax.scatter(
                self.resultHistory[meshIteration]["time"],
                self.resultHistory[meshIteration][column],
                color="red",
                s=30,
                edgecolor="black",
                alpha=0.8,
                zorder=5,
            )
        elif dots == Dots.BASE:
            ax.scatter(
                [x for i, x in enumerate(self.resultHistory[meshIteration]["time"]) if i % self.rksteps == 0],
                [x for i, x in enumerate(self.resultHistory[meshIteration][column]) if i % self.rksteps == 0],
                color="red",
                s=30,
                edgecolor="black",
                alpha=0.8,
                zorder=5,
            )

    # Helper function to apply dots in parametric plots
    def _applyDotsParametric(self, x_data, y_data, dots):
        if dots == Dots.ALL:
            plt.scatter(x_data, y_data, color="red", s=30, edgecolor="black", alpha=0.8, zorder=5)
        elif dots == Dots.BASE:
            plt.scatter(
                [x for i, x in enumerate(x_data) if i % self.rksteps == 0],
                [y for i, y in enumerate(y_data) if i % self.rksteps == 0],
                color="red",
                s=30,
                edgecolor="black",
                alpha=0.8,
                zorder=5,
            )

    # Private mesh refinement plot function
    def _plotMeshRefinement(self, ax=None, interval=None, markerSize=30, dots=Dots.BASE, epsilon=1e-14):
        self.setPlotDefaults()
        from matplotlib.ticker import MaxNLocator

        if interval is None:
            interval = [0, self.tf]

        prev_points = np.array([])

        fig, ax = plt.subplots(figsize=(12, 8)) if ax is None else (None, ax)

        for m in range(self.modelInfo["maxMeshIteration"] + 1):
            self.getResults(m)
            arr = self.resultHistory[m]["time"].to_numpy()[:: self._getModulo(dots)]

            if m == 0:
                ax.scatter(arr, np.ones_like(arr) * m, color="red", s=markerSize, edgecolor="black", alpha=0.8, zorder=5)
            else:
                new_points = np.array([pt for pt in arr if prev_points.size == 0 or not np.any(np.abs(prev_points - pt) < epsilon)])
                ax.scatter(new_points, np.ones_like(new_points) * m, color="red", s=markerSize, edgecolor="black", alpha=0.8, zorder=5)

            prev_points = arr

        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_xlim(interval)
        ax.set_xlabel("Time")
        ax.set_ylabel("Iteration")
        ax.set_title("Inserted Mesh Points Over Time")
        return fig, ax

    # Helper function to get modulo value based on dots argument
    def _getModulo(self, dots):
        return 1 if dots == Dots.ALL else self.rksteps

    def plotSparseMatrix(self, matrixType):
        from matplotlib.patches import Rectangle
        import scipy.sparse as sp

        if matrixType == MatrixType.JACOBIAN:
            file_path = self.exportJacobianPath + f"/{self.name}_jacobian.csv"
        elif matrixType == MatrixType.HESSIAN:
            file_path = self.exportHessianPath + f"/{self.name}_hessian.csv"
        else:
            raise InvalidMatrix("Plotting is only possible for matrixTypes JACOBIAN or HESSIAN.")

        with open(file_path, "r") as f:
            lines = f.readlines()

        dim_row, dim_col = map(int, lines[1].strip().split(","))

        rows, cols = [], []
        for line in lines[3:]:
            row, col = map(int, line.strip().split(","))
            rows.append(row)
            cols.append(col)

        data = np.ones(len(rows))
        m = sp.coo_matrix((data, (rows, cols)), shape=(dim_row, dim_col))

        fig, ax = plt.subplots()

        for x, y, data in zip(m.col, m.row, m.data):
            ax.add_patch(Rectangle(xy=(x, y), width=1, height=1, edgecolor="black", facecolor="blue", alpha=0.6))

        ax.set_xlim(0, m.shape[1])
        ax.set_ylim(0, m.shape[0])
        ax.invert_yaxis()

        if matrixType == MatrixType.JACOBIAN:
            ax.set_xlabel("Variables")
            ax.set_ylabel("Equations")
            ax.set_title("Jacobian Sparsity")
        elif matrixType == MatrixType.HESSIAN:
            ax.set_xlabel("Variables")
            ax.set_ylabel("Variables")
            ax.set_title("Hessian Sparsity")

        plt.show()


# time symbol
t = Symbol("t")
time = t

### GLOBAL ALIAS AND GLOBAL VAR DEFINITIONS

# variables
Model.addControl = Model.addInput
Model.addContinuous = Model.addInput
Model.addX = Model.addState
Model.addU = Model.addInput
Model.addP = Model.addParameter
Model.addRP = Model.addRuntimeParameter

# objective
Model.addM = Model.addMayer
Model.addL = Model.addLagrange

# constraints
Model.addOde = Model.addDynamic
Model.addF = Model.addDynamic
Model.addG = Model.addPath
Model.addR = Model.addFinal
Model.addA = Model.addParametric
