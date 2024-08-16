from sympy import *


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
