import os
from optimization.variables import *
from optimization.expressions import *
from optimization.structures import *
import time as timer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# set global precision
pd.set_option('display.precision', 16)

# TODO: only diff if first diff != 0
# add vectorized eval of RHS = [f, g]^T, vectorized evalDiff, evalDiff2?
# or with colored jacobian
# TODO: make framework more robust, cleaner

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
        self.alias = {}
        self.addedDummy = False
        
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

        # stuff for analyzing the results
        self.resultHistory = {}
        self.modelInfo = {}

    def addState(self, start, symbol=None, lb=-float("inf"), ub=float("inf")):
        
        # adds a state x for optimization, must a have starting value
        # specify lb or ub if needed
        
        variable = State(start, symbol=symbol, lb=lb, ub=ub)
        self.alias[variable.symbol] = variable.name
        self.xVars.append(variable)
        return variable
    
    def addInput(self, symbol=None, lb=-float("inf"), ub=float("inf")):
        
        # adds an input / control u for optimization
        # specify lb or ub if needed
        
        variable = Input(symbol=symbol, lb=lb, ub=ub)
        self.alias[variable.symbol] = variable.name
        self.uVars.append(variable)
        return variable
        
    def addParameter(self, symbol=None, lb=-float("inf"), ub=float("inf")):
        
        # adds a parameter p for optimization
        # specify lb or ub if needed

        variable = Parameter(symbol=symbol, lb=lb, ub=ub)
        self.alias[variable.symbol] = variable.name
        self.pVars.append(variable)
        return variable
    
    def addRuntimeParameter(self, default, symbol=None):
        
        # adds a runtime parameter: can be changed for optimization
        # will be handled symbolically -> thus quite slow
        # for actual constants simply write pythonic var = value
        
        variable = RuntimeParameter(default, symbol=symbol)
        self.alias[variable.symbol] = variable.name
        self.rpVars.append(variable)
        return variable

    def addObjective(self, mayer, lagrange, obj=Objective.MINIMIZE):
        self.addMayer(mayer, obj)
        self.addLagrange(lagrange, obj)

    def addMayer(self, expr, obj=Objective.MINIMIZE):
        
        # adds the mayer term: min/max expr(tf)
        
        if self.M:
            raise InvalidModel("Mayer term already set")
        else:
            if obj == Objective.MAXIMIZE:
                self.M = Expression(-1 * expr)
                print("Setting Mayer term as -1 * mayer, since maximization is chosen.\n")
            else:
                self.M = Expression(expr)
                
    def addLagrange(self, expr, obj=Objective.MINIMIZE):
        
        # adds the lagrange term: min/max integral_0_tf expr dt
        
        if self.L:
            raise InvalidModel("Lagrange term already set")
        else:
            if obj == Objective.MAXIMIZE:
                self.L = Expression(-1 * expr)
                print("Setting Lagrange term as -1 * lagrange, since maximization is chosen.\n")
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
        self.addDynamic(x, 0)
        self.addedDummy = True
    
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
        
        OUTPUT += "// runtime parameters and global constants\n"
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
            {{{', '.join(str(toCode(x.start)) for x in self.xVars)}}},  // x0
            {{{', '.join(str(toCode(x.lb) if x.lb != -float('inf') else "MINUS_INFINITY") for x in self.xVars)}}},  // lb x
            {{{', '.join(str(toCode(x.ub) if x.ub != float('inf') else "PLUS_INFINITY") for x in self.xVars)}}},  // ub x
            {{{', '.join(str(toCode(u.lb) if u.lb != -float('inf') else "MINUS_INFINITY") for u in self.uVars)}}},  // lb u
            {{{', '.join(str(toCode(u.ub) if u.ub != float('inf') else "PLUS_INFINITY") for u in self.uVars)}}},  // ub u
            {{{', '.join(str(toCode(p.lb) if p.lb != -float('inf') else "MINUS_INFINITY") for p in self.pVars)}}},  // lb p
            {{{', '.join(str(toCode(p.ub) if p.ub != float('inf') else "PLUS_INFINITY") for p in self.pVars)}}},  // ub p
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


        os.makedirs(f".generated/{self.name}",exist_ok=True)

        with open(f'.generated/{self.name}/{filename}.cpp', 'w') as file:
            file.write(OUTPUT)

        print(f"Generated model to .generated/{self.name}/{filename}.cpp\n")
        print(f"Model creation, derivative calculations, and code generation took {round(timer.process_time() - self.creationTime, 4)} seconds.")
        return 0
        
    def optimize(self, tf=0, steps=1, rksteps=1, flags={}, meshFlags={}):
        
        # generate corresponding main function with flags, mesh, refinement
        # set runtime parameter file from map
        # run the code

        # always with setter to ensure some security
        if not self.addedDummy:
            self.setFinalTime(tf)
            self.setSteps(steps)
            self.setRkSteps(rksteps)
        else:
            print("\nSetting tf = 0, steps = 1, rksteps = 1, since the model is purely parametric.")
            self.setFinalTime(0)
            self.setSteps(1)
            self.setRkSteps(1)

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

        with open(f'.generated/{self.name}/{filename}Params.h', 'w') as file:
            file.write(OUTPUT)

        os.system(f"g++ .generated/{self.name}/{filename}.cpp -O3 -I../src/ -L../cmake-build-release/src -lipopt_do -o.generated/{self.name}/{self.name}") # vorher src lipopt_do

        os.system(f"LD_LIBRARY_PATH=../cmake-build-release/src/ ./.generated/{self.name}/{self.name}")

        self.initAnalysis()

        return 0

    def initAnalysis(self):
        with open('/tmp/modelinfo.txt', 'r') as file:
            for line in file:
                line = line.strip()
                key, value = line.split(',')
                key = key.strip()
                value = value.strip()
                self.modelInfo[key] = int(value)

    def checkMeshIteration(self, meshIteration):
        maxMeshIteration = self.modelInfo["maxMeshIteration"]
        if meshIteration is None:
            return maxMeshIteration
        if type(meshIteration) == int:
            if meshIteration > maxMeshIteration:
                print(f"meshIteration too large. Setting meshIteration to maximum value of {maxMeshIteration}.")
                meshIteration = maxMeshIteration
        return meshIteration

    def plotStates(self, meshIteration=None, interval=None, dots=False):
        self.plot(meshIteration=meshIteration, interval=interval, dots=dots, specifCols=[v.name for v in self.xVars])

    def plotInputs(self, meshIteration=None, interval=None, dots=False):
        self.plot(meshIteration=meshIteration, interval=interval, dots=dots, specifCols=[v.name for v in self.uVars])

    def plot(self, meshIteration=None, interval=None, specifCols=None, dots=False):
        meshIteration = self.checkMeshIteration(meshIteration)
        if interval is None:
            interval = [0, self.tf]
        self.getResults(meshIteration=meshIteration)
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
            columns_to_plot = self.resultHistory[meshIteration].columns[1:]
        else:
            columns_to_plot = specifCols

        num_plots = len(columns_to_plot)
        fig, axs = plt.subplots(num_plots, 1, figsize=(12, 8 * num_plots), sharex=True)

        if num_plots == 1:
            axs = [axs]

        for idx, column in enumerate(columns_to_plot):
            ax = axs[idx]
            ax.plot(self.resultHistory[meshIteration]['time'], self.resultHistory[meshIteration][column], label=column, linewidth=2, linestyle='-', color='steelblue')
            if dots:
                ax.scatter(self.resultHistory[meshIteration]['time'], self.resultHistory[meshIteration][column], color='red', s=30, edgecolor='black', alpha=0.8, zorder=5)
            ax.set_xlabel('time')
            ax.set_ylabel(column)
            ax.set_xlim(interval[0], interval[1])
            ax.legend(frameon=True, loc='best')
            ax.grid(True)
            ax.title.set_fontsize(16)

        plt.tight_layout()
        plt.show()

    def getResults(self, meshIteration=None):
        meshIteration = self.checkMeshIteration(meshIteration)

        if meshIteration not in self.resultHistory:
            try:
                results = pd.read_csv(self.outputFilePath + "/" + self.name + str(meshIteration) + ".csv", sep=",")
            except:
                raise Exception("meshIteration out of range. Set Model.maxMeshIteration to the maximum mesh iteration!")
            # remove dummy column for purely parametric models
            if self.addedDummy:
                results = results.drop(columns=['x[0]'])

            results.rename(columns=self.alias, inplace=True)
            self.resultHistory[meshIteration] = results
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
        print(self.resultHistory[meshIteration][[v.name for v in self.pVars]].iloc[0].to_string())

    def plotPointDensity(self,  meshIteration=None):

        # test version of a density function to visualize the thickening at jumps and steep sections

        meshIteration = self.checkMeshIteration(meshIteration)

        if meshIteration == "all":    # all iterations
            meshIteration = range(0, self.modelInfo["maxMeshIteration"] + 1)
        elif type(meshIteration) == int:    # specific iteration
            meshIteration = [meshIteration]

        # else given iterable of iterations

        for m in meshIteration[::-1]:   # reverse the list -> s.t. you can view at which iteration some detections stopped
            self.getResults(m)
            arr = self.resultHistory[m]["time"].to_numpy()

            window_size = self.tf /self.steps
            step_size = self.tf / (self.steps * self.rksteps)

            arr_min, arr_max = min(arr), max(arr)
            arr_range = np.arange(arr_min, arr_max, step_size)

            local_density = []
            for center in arr_range:
                count = np.sum((arr >= center - window_size / 2) & (arr <= center + window_size / 2))
                density = count / window_size
                local_density.append(density)

            plt.plot(arr_range, local_density, label=f'Mesh Iteration {m}')

            if len(meshIteration) == 1:
                plt.scatter(arr, np.zeros_like(arr) - 0.01, color='red', s=30, edgecolor='black', alpha=0.8, zorder=5)

        plt.xlabel('time')
        plt.ylabel('Local Density')
        plt.title('Local Mesh Point Density')
        plt.legend()

        plt.show()

    def plotSparseMatrix(self, matrixType):
        from matplotlib.patches import Rectangle
        import scipy.sparse as sp
        if matrixType == MatrixType.JACOBIAN:
            file_path = self.exportJacobianPath + f"/{self.name}_jacobian.csv"
        elif matrixType == MatrixType.HESSIAN:
            file_path = self.exportHessianPath  + f"/{self.name}_hessian.csv"
        else:
            raise InvalidMatrix("Plotting is only possible for matrixTypes JACOBIAN or HESSIAN.")

        with open(file_path, 'r') as f:
            lines = f.readlines()

        dim_row, dim_col = map(int, lines[1].strip().split(','))

        rows, cols = [], []
        for line in lines[3:]:
            row, col = map(int, line.strip().split(','))
            rows.append(row)
            cols.append(col)

        data = np.ones(len(rows))
        m = sp.coo_matrix((data, (rows, cols)), shape=(dim_row, dim_col))

        fig, ax = plt.subplots()

        for (x, y, data) in zip(m.col, m.row, m.data):
            ax.add_patch(Rectangle(
                xy=(x, y), width=1, height=1, edgecolor='black', facecolor='blue', alpha=0.6))

        ax.set_xlim(0, m.shape[1])
        ax.set_ylim(0, m.shape[0])
        ax.invert_yaxis()

        if matrixType == MatrixType.JACOBIAN:
            ax.set_xlabel('Variables')
            ax.set_ylabel('Equations')
            ax.set_title('Jacobian Sparsity')
        elif matrixType == MatrixType.HESSIAN:
            ax.set_xlabel('Variables')
            ax.set_ylabel('Variables')
            ax.set_title('Hessian Sparsity')

        plt.show()

"""
### GLOBAL ALIAS AND GLOBAL VAR DEFINITIONS

# variables
Continuous = Input
Control = Input

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
"""
# time symbol
t = Symbol("t")
time = t
