from pyocp.common import *
from pyocp.arc import arc
from sympy import sin, cos, sqrt

# States
r, u, v, phi = time_variables(['r', 'u', 'v', 'phi'], t)
t, r_0, m_0, m_dot, mu, T = symbols('t, r_0, m_0, m_dot, mu, T')
t_i, t_f = symbols('t_i, t_f')

r_dot, u_dot, v_dot = [v.diff(t) for v in [r, u, v]]

m = m_0 - m_dot * t

odes = [-r_dot + u,
        -u_dot + v**2 / r - mu / r**2 + T * sin(phi) / m,
        -v_dot + -u * v / r + T * cos(phi) / m]

arc1 = arc(1, t)

arc1.objective(term=r)
arc1.odes(odes)
arc1.bounds(t_f=(t_f, t_f))
arc1.constraints(initial=[r - r_0, u, sqrt(mu / r_0)],
                 terminal=[u, v - sqrt(mu / r)])
arc1.problem_definition()
arc1.numerical_defs(m=10, args={T:0.85, m_0:10000, t_f:193*24*60*60, m_dot:12.9/24/60/60, mu:0.00060498220})

f, f_x, c = arc1.generate_funcs()

