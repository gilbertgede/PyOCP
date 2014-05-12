from pyocp.common import *
from pyocp.arc_c import arc
from pyocp.problem import problem
from sympy import sin, cos, sqrt

# States
x, v, w = time_variables('x, v, w', t)
x_d, v_d, w_d = [i.diff(t) for i in [x, v, w]]
# Controls
u = time_variables('u', t)
T = Function('T')
# Params
l = symbols('l')
# ODEs
odes = [x_d - v,
        v_d - u,#T(x, v, w) * u,
        w_d - u**2/2]#T(x, v, w) * u**2 / 2]

DI = problem('Double Integrator', 3, t)
for i in range(3):
    DI[i].odes(odes)
    DI[i].path_constraints(g=Matrix([l-x]), g_l=Matrix([0]), g_u=Matrix([oo]))
    DI[i].continuous_bounds()
    DI[i].parameter_bounds(arguments=[l])

DI[2].objective(term=w)

DI.link_arcs(left=None, right=0, right_cons=[x, v - 1, w], fixed_time=0)
DI.link_arcs(left=0, right=1, left_cons=[x, v, w], right_cons=[x, v, w], free_time=0.3)
DI.link_arcs(left=1, right=2, left_cons=[x, v, w], right_cons=[x, v, w], free_time=0.7)
DI.link_arcs(left=2, right=None, left_cons=[x, v + 1], fixed_time=1)
DI.complete_links()

s_g = {x : [0, 1./9, 1./9, 0],
        v : [1, 0, 0, -1],
        w : [0, 2, 2, 4]}

c_g = {u : [-6, 0, 0, -6]}
tg = [0, .3, .7, 1]


DI.numerical_defs(m=[10, 10, 10], states_guess=s_g, controls_guess=c_g,
                  time_guess_grid=tg, parameters_guess={}, arguments={l : 1. / 9})

DI.write_c_files()

#for i in range(3):
#    DI[i].numerical_defs(m=10, states_guess=s_g[i], controls_guess=c_g[i],
#                         parameters_guess={}, arguments={l : 1. / 9})










"""

arc1 = arc(1, t)

arc1.objective(term=(w))#*T(x,v,w)))
arc1.odes(odes)
arc1.boundary_constraints(initial=[x, v - 1, w], terminal=[x, v + 1])
arc1.path_constraints(g=Matrix([l-x]), g_l=Matrix([0]), g_u=Matrix([oo]))
arc1.continuous_bounds()
arc1.parameter_bounds(arguments=[l])
print(arc1)

# numbers from analytic solution
xn = array([0.        ,  0.073     ,  0.104     , 0.111     ,  0.11111111,
            0.11111111,  0.11111111,  0.111     , 0.104     ,  0.073     ,  0.])
vn = array([ 1.  ,  0.49,  0.16,  0.01,  0.  ,  0.  ,  0.  , -0.01, -0.16,
            -0.49, -1.  ])
wn = array([0.   ,  1.314,  1.872,  1.998,  2.   ,  2.   ,  2.   ,  2.002,
            2.128,  2.686,  4.   ])
states_guess = {x : xn, v : vn, w : wn}
un = array([-6. , -4.2, -2.4, -0.6,  0. ,  0. ,  0. , -0.6, -2.4, -4.2, -6. ])
controls_guess = {u : un}
parameters_guess = {}

arc1.numerical_defs(m=11, states_guess=states_guess,
                    controls_guess=controls_guess,
                    parameters_guess=parameters_guess, arguments={l:1/9})#,
#                    external_include='#include "t_func.h"')
#f, f_x, c, c_x = arc1.generate_c_files()
arc1.generate_c_files()


nlp_x = hstack([xn,vn,wn,un])
for i in range(len(xn)):
    nlp_x[0 + 4 * i] = xn[i]
    nlp_x[1 + 4 * i] = vn[i]
    nlp_x[2 + 4 * i] = wn[i]
    nlp_x[3 + 4 * i] = un[i]

from scipy.sparse import coo_matrix
rows, cols, data = c_x(nlp_x)
c_x0 = coo_matrix((data, (rows, cols))).toarray()

ll = []
for i in range(len(nlp_x)):
    temp = num_diff(c, nlp_x, i)
    ll += [max(abs(temp - c_x0[:, i]))]

print(ll)
"""


