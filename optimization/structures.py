from enum import Enum


class InvalidModel(Exception):
    pass


class InvalidMatrix(Exception):
    pass


class InitVars(Enum):
    CONST = 1
    SOLVE = 2
    SOLVE_EXPLICIT = 3
    SOLVE_EXPLICIT_EULER = 4


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


class RefinementMethod(Enum):
    LINEAR_SPLINE = 1
    POLYNOMIAL = 2


class MeshAlgorithm(Enum):
    NONE = 1
    BASIC = 2
    L2_BOUNDARY_NORM = 3


class MatrixType(Enum):
    JACOBIAN = 1
    HESSIAN = 2


class Dots(Enum):
    OFF = 1
    ALL = 2
    BASE = 3


class IVPSolver(Enum):
    Radau = 1
    BDF = 2
    LSODA = 3
    RK45 = 4
    DOP853 = 5
    RK23 = 6
    RADAU = Radau
