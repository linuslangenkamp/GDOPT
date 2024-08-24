class VariableStruct:
    id_counter = 0

    def __init__(self, symbol, lb=-float("inf"), ub=float("inf")):
        self.symbol = symbol
        self.lb = lb
        self.ub = ub
        self.id = VariableStruct.id_counter
        VariableStruct.id_counter += 1

class StateStruct(VariableStruct):
    id_counter = 0

    def __init__(self, start, symbol=None, lb=-float("inf"), ub=float("inf")):
        super().__init__(symbol, lb, ub)
        self.start = start
        self.symbol = symbol if symbol is not None else f'x[{self.id}]'
        self.id = StateStruct.id_counter
        StateStruct.id_counter += 1

class InputStruct(VariableStruct):
    id_counter = 0

    def __init__(self, symbol=None, lb=-float("inf"), ub=float("inf")):
        super().__init__(symbol, lb, ub)
        self.symbol = symbol if symbol is not None else f'u[{self.id}]'
        self.id = InputStruct.id_counter
        InputStruct.id_counter += 1


class ParameterStruct(VariableStruct):
    id_counter = 0

    def __init__(self, symbol=None, lb=-float("inf"), ub=float("inf")):
        super().__init__(symbol, lb, ub)
        self.symbol = symbol if symbol is not None else f'p[{self.id}]'
        self.id = ParameterStruct.id_counter
        ParameterStruct.id_counter += 1


class RuntimeParameterStruct(VariableStruct):
    id_counter = 0

    def __init__(self, default, symbol=None, lb=-float("inf"), ub=float("inf")):
        super().__init__(symbol, lb, ub)
        self.symbol = symbol if symbol is not None else f'rp[{self.id}]'
        self.id = ParameterStruct.id_counter
        self.value = default
        ParameterStruct.id_counter += 1

class_order = {
    'State': 0,
    'Input': 1,
    'Parameter': 2
}


def get_sort_key(symbol, obj_map):
    obj = obj_map[symbol]  # Get the object that corresponds to the symbol
    class_name = obj.__class__.__name__.replace('Struct', '')  # Get the class name
    order = class_order.get(class_name, float('inf'))  # Determine the type order
    return (order, obj.id)  # Return a tuple with type order and id for sorting


def sort_symbols(symbols, obj_map):
    return symbols.sort(key=lambda s: get_sort_key(s, obj_map))