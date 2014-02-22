from sympy import symbols, Function, Matrix, Symbol, zeros, oo
from sympy.core.function import AppliedUndef, Derivative

t = symbols('t')

def time_variables(var_names, time=t):
    out_list = []
    for v in var_names:
        x = symbols(v, cls=Function)
        out_list += [x(time)]
    return out_list


