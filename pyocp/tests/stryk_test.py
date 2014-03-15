from pyocp.common import *
from pyocp.arc import arc
from sympy import sin, cos, sqrt
import numpy as np
from scipy.sparse import coo_matrix
import ipopt

# States
x, v, w = time_variables('x, v, w', t)
x_d, v_d, w_d = [i.diff(t) for i in [x, v, w]]
# Controls
u = time_variables('u', t)
# Params
l = symbols('l')
# ODEs
odes = [x_d - v,
        v_d - u,
        w_d - u**2 / 2]

arc1 = arc(1, t)

arc1.objective(term=w)
arc1.odes(odes)
arc1.bounds()
arc1.constraints(initial=[x, v - 1, w], terminal=[x, v + 1], g=Matrix([l-x]), g_l=Matrix([0]), g_u=Matrix([oo]))
arc1.problem_definition()
arc1.numerical_defs(m=11, args={l:1/9})

f, f_x, c, c_x = arc1.generate_funcs()
# numbers from analytic solution
xn = np.array([0.        ,  0.073     ,  0.104     , 0.111     ,  0.11111111,
            0.11111111,  0.11111111,  0.111     , 0.104     ,  0.073     ,  0.])
vn = np.array([ 1.  ,  0.49,  0.16,  0.01,  0.  ,  0.  ,  0.  , -0.01, -0.16,
            -0.49, -1.  ])
wn = np.array([0.   ,  1.314,  1.872,  1.998,  2.   ,  2.   ,  2.   ,  2.002,
            2.128,  2.686,  4.   ])
un = np.array([-6. , -4.2, -2.4, -0.6,  0. ,  0. ,  0. , -0.6, -2.4, -4.2, -6. ])
nlp_x = np.hstack([xn,vn,wn,un])
for i in range(len(xn)):
    nlp_x[0 + 4 * i] = xn[i]
    nlp_x[1 + 4 * i] = vn[i]
    nlp_x[2 + 4 * i] = wn[i]
    nlp_x[3 + 4 * i] = un[i]



class stryk(object):
    def objective(self, xi):
        return np.array([f(xi)])
    def gradient(self, xi):
        return f_x(xi)
    def constraints(self, xi):
        return c(xi)
    def jacobian(self, xi):
        r, c, d = c_x(xi)
        return coo_matrix((d, (r, c))).toarray()


x0 = nlp_x
lb = [-1e19] * len(x0)
ub = [ 1e19] * len(x0)

cl = [0, 0, 0] + [0, 0, 0, 0] * 10 + [0] + [0, 0]
cu = [0, 0, 0] + [1e19, 0, 0, 0] * 10 + [1e19] + [0, 0]

nlp = ipopt.problem(
            n=len(x0),
            m=len(cl),
            problem_obj=stryk(),
            lb=lb,
            ub=ub,
            cl=cl,
            cu=cu)

#
# Set solver options
#
#nlp.addOption('derivative_test', 'second-order')
nlp.addOption('mu_strategy', 'adaptive')
nlp.addOption('tol', 1e-7)

#
# Scale the problem (Just for demonstration purposes)
#

#
# Solve the problem
#
x, info = nlp.solve(x0)

print("Solution of the primal variables: x=%s\n" % repr(x))

print("Solution of the dual variables: lambda=%s\n" % repr(info['mult_g']))

print("Objective=%s\n" % repr(info['obj_val']))







