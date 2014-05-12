from sympy import (symbols, Function, Matrix, Symbol, zeros, oo, sympify,
                   DeferredVector as DV, lambdify, Subs, Dummy, ccode, S)
from sympy.core.function import AppliedUndef, Derivative
import numpy
import os
from numpy.linalg import norm
from numba import autojit

eps = numpy.finfo(numpy.float64).eps
big = 1 / eps**2
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


def time_variables(names, time=t):
    esses = symbols(names, cls=Function)
    if hasattr(esses, '__iter__'):
        return [e(time) for e in esses]
    else:
        return esses(t)


def num_diff(f_in, x0, var=0, num_args=1):
    """
    A numerical differencing function.

    Assumes x0 is a numpy array.
    """
    f = f_in
    if num_args > 1:
        f = lambda x: f_in(*x)

    f3 = lambda h: 0.5 / (sum(h)**3) * (-f(x0 - 2*h) + 2 * f(x0 - h) -
                                        2 * f(x0 + h) + f(x0 + 2*h))
    dx = numpy.zeros(len(x0))
    dx[var] = 1.e-4
    for i in range(5):
        h0 = dx[var]
        temp = norm(f3(dx))
        if temp < eps:     # This is for functions where f'''(x) = 0
            dx[var] = 1    # If they're polynomial and this happens, the
            break          # step size doesn't really matter?
        dx[var] = (3 * eps / temp)**(1./3.)
        if abs(h0 - dx[var]) < eps:
            break
    h0 = dx[var]

    return 0.5 / h0 * (f(x0 + dx) - f(x0 - dx))


def pull_out_derivatives(D_vector, name, dv_ind_var, func_dict):
    """
    Takes in a list of quantities, which have present Derivative terms -
    meaning that there are derivatives of numerical functions in them.

    Returns a list, dict, and deferred vector to evaluate derivatives.
    """
    ders = dict()
    varss = dict()

    # This is assuming that Subs are only have a Derivative object inside them
    for d in D_vector:
        s_atoms = list(d.atoms(Subs))
        temp = []
        for s in s_atoms:
            d = s.args[0].copy()
            if d.__class__ != Derivative:
                raise Exception("Unexpected Subs in a numerical diff conversion")
            d._args = (d.args[0].subs(s.args[1][0], s.args[2][0]), s.args[2][0])
            ders.update({s : d.expr})
            varss.update({s : d.variables[0]})

    for d in D_vector:
        d_atoms = list(d.atoms(Derivative))
        d_atoms = [i for i in d_atoms if i.atoms(Dummy) == set()]
        ders.update(zip(d_atoms, [da.expr for da in d_atoms]))
        varss.update(zip(d_atoms, [d.variables[0] for d in d_atoms]))

    key_list = ders.keys()
    f = [func_dict[str(ders[d].func)] for d in key_list]
    x0_func = [lambdify(dv_ind_var, ders[d].args, func_dict) for d in key_list]
    index = [ders[d].args.index(varss[d]) for d in key_list]
    num_args = [len(ders[d].args) for d in key_list]

    # Needed due to the scope issues w/ lambda functions
    def cbf(f, xf, i, na):
        return lambda x: num_diff(f, xf(x), i, na)

    D_func = DeferredVector(name)
    # Order is: D_func, D_func_list, D_func_dict
    return (D_func,
            [cbf(f[i], x0_func[i], index[i], num_args[i]) for i in range(len(f))],
            dict(zip(key_list, D_func)))


