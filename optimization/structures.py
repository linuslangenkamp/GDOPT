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


class MeshAlgorithm(Enum):
    NONE = 1
    BASIC = 2
    L2_BOUNDARY_NORM = 3


class MatrixType(Enum):
    DEFAULT = 1
    JACOBIAN = 2
    HESSIAN = 3
