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
arc1.odes(odes)
arc1.bounds()
arc1.constraints(initial=[r - r_0, u, sqrt(mu / r_0)],
                 terminal=[u, v - sqrt(mu / r)])
arc1.problem_definition()

