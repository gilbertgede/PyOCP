from sympy import (symbols, Function, Matrix, Symbol, zeros, oo, sympify,
                   DeferredVector as DV, lambdify)
from sympy.core.function import AppliedUndef, Derivative
import numpy as np
from numba import autojit

eps = np.finfo(np.float64).eps
t = symbols('t')

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


def time_variables(var_names, time=t):
    out_list = []
    for v in var_names:
        x = symbols(v, cls=Function)
        out_list += [x(time)]
    return out_list


def num_diff(f_in, x0, var=0):
    """
    A numerical differencing function.

    Assumes x0 is a numpy array.
    """
    f = lambda x: f_in(*x)
    def f3(h):
        return 0.5 / (sum(h)**3) * (-f(x0 - 2*h) + 2 * f(x0 - h) -
                                    2 * f(x0 + h) + f(x0 + 2*h))
    dx = np.zeros(len(x0))
    h0 = 1.e-2
    dx[var] = h0
    print(x0)
    print(dx)
    print(f(x0))
    for i in range(5):
        dx[var] = (3 * eps / abs(f3(dx)))**(1./3.)
        if (h0 - dx[var]) < eps:
            h0 = dx[var]
            break
        h0 = dx[var]

    return 0.5 / h0 * (f(x0 + dx) - f(x0 - dx))







