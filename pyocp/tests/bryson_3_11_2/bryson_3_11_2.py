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
DI.run_and_read()
