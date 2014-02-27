from sympy import (symbols, Function, Matrix, Symbol, zeros, oo, sympify,
                   DeferredVector as DV, lambdify)
from sympy.core.function import AppliedUndef, Derivative
import numpy as np


t = symbols('t')

def time_variables(var_names, time=t):
    out_list = []
    for v in var_names:
        x = symbols(v, cls=Function)
        out_list += [x(time)]
    return out_list

# These are just to make the translation process easier, so that the parts of
# the DeferredVector can be picked out.
class DeferredVector(DV):
    def __getitem__(self, i):
        if i == -0:
            i = 0
        if i < 0:
            raise IndexError('DeferredVector index out of range')
        component_name = '%s[%d]' % (self.name, i)
        DS = DeferredSymbol(component_name)
        DS._num = i
        return DS

class DeferredSymbol(Symbol):
    _num = 0
