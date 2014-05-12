from .common import *

class arc(object):
    """
    The main building block of an PyOCP problem.

    Each arc has it's own state variables, ODEs, constraints, and time grid.

    Glossary of terms:
        states | state variables
            Symbolic quantities which are only functions of time and represent
            the states (configuration) of the system and whose behavior is
            governed by ODEs.
        controls | control variables
            Symbolic quantities which are only functions of time and represent
            the controls (inputs) to the system - the time history of the
            controls is what PyOCP helps find...
        parameters
            Symbolic quantities which are not functions, but only symbols. Will
            be considered "free" to be optimized as part of the optimization
            process.
        arguments
            Symbolic quantities which are only symbols. Will _NOT_ be
            considered "free" to be optimized. Must instead have a numerical
            value supplied at some point.
            controls, specified, parameters, arguments, etc. Needs to be
            either later supplied as: numerical functions, or sympy
            expressions.

    This is still a work in progress...
    """

    def __init__(self, number, time_symbol, name=None, subarc=False):
        """
        Initializer for arc class.

        Needs to be supplied with a number, and optionally a name (e.g.,
        "atmospheric_phase"). If no name is given, the arc number is used.
        """

        self.number = number
        if name:
            self.name = name
        else:
            self.name = "Arc " + str(number)

        self.time = sympify(time_symbol)

        self._obj_init = sympify(0)
        self._obj_term = sympify(0)

        self._initial = Matrix([[]])
        self._terminal = Matrix([[]])
        self._g = Matrix([[]])

        self._t_i = (0, 0, 0)
        self._t_f = (1, 1, 1)
        self._initial_guesses = {}
        self._param_limits = {}

        self.subarc = subarc
        self._left_internal = False
        self._right_internal = False
        self._subarc_jgr = []
        self._subarc_jgc = []

        self.template_names = {}

    def _sanitize_bounds(self, bounds, n, name, cons=False):
        """
        Santitizing the bounds the user provides.

        Valid entries in the bounds Matrix are: SymPy Symbols, constant
        numerical scalars, SymPy functions of the arc's time variable which
        will be later substituted for, or valid SymPy expressions with no
        dependencies on states (e.g., a combination of a function of time and a
        parametric scaling factor).
        """
        if (not isinstance(bounds, Matrix)) or (bounds.shape != (n, 1)):
            raise Exception('Bounds for ' + name + ' need to be ' +
                            'supplied as a SymPy matrix with size: ' +
                            str(n) + 'x1.')

        # first checking to make sure no states/controls are in the bounds
        pararg_list = [i for i in self._params_args]
        func_list = [i for i in self._funcs]

        for v in bounds:
            syms = v.atoms(Symbol) - {self.time,}
            pararg_list += [i for i in syms if i not in pararg_list]

            bad = v.atoms(AppliedUndef)
            for b in bad:
                if not cons:
                    if (b in self._states) or (b in self._controls) or (b in self._funcs):
                        raise Exception('Bounds cannot contain states or ' +
                                        'controls or functions.')
                else:
                    if not ((b in self._states) or (b in self._controls) or
                            (b.func in func_list)):
                        func_list += [b.func]

        self._params_args = pararg_list
        self._funcs = Matrix(func_list)


    def _print_bounds(self, l, x, u, lead):
        """
        Just used for some prettier printing.

        e.g.

          0  <=  x(t)  <=      1
         10  <=  y(t)  <=  v * 2
          1  <=  u(t)  <=      3

        """

        if len(x) == 0:
            return lead + "n/a"

        ls = repr(l).split('\n')
        xs = repr(x).split('\n')
        us = repr(u).split('\n')

        if len(xs) == 1:
            outstr = lead
            for v in [ls, xs, us]:
                # start='Matrix([['=len 9 - end=')]]'=len 3
                outstr += v[0][9:-3] + ' ≤ '
            return outstr[:-3]
        else:
            for v in [ls, xs, us]:
                v.pop(0)

        ls = [i[1:-2] for i in ls[:-1]] + [ls[-1][1:-3]]
        xs = [i[1:-2] for i in xs[:-1]] + [xs[-1][1:-3]]
        us = [i[1:-2] for i in us[:-1]] + [us[-1][1:-3]]
        le = [' ≤ '] * len(xs)

        out_list = list(zip(ls, le, xs, le, us))
        outstr = ''
        out_list = [outstr.join(s) for s in out_list]
        for s in out_list:
            outstr += lead + s + '\n'
        outstr = outstr[:-1]
        return outstr


    def _rename_funcs(self, x, append):
        """
        Appends names of functions; assumes a matrix is passed in for x.
        """
        c = x.T.tolist()[0]
        for i, v in enumerate(c):
            t = v.args[0]
            name = v.func.__name__
            new_var = time_variables([name + append], t)
            c[i] = new_var
        return Matrix(c)


    def objective(self, **kwargs):
        """
        Provide the objective function for the problem.

        Needs to be in form:
            J = term + init

        An integral term is not yet supported, and instead needs to be added as
        a state. This is the Mayer form.

        kwargs to supply are:
            init : SymPy Expression
                A SymPy Expression which is made of the states/controls/time at
                the initial time; e.g., init = f(x(t_i), u(t_i)).
            term : SymPy Expression
                A SymPy Expression which is made of the states/controls/time at
                the final time; e.g., term = f(x(t_f), u(t_f)).

        """

        if 'init' in kwargs:
            self._obj_init += kwargs['init']
        if 'term' in kwargs:
            self._obj_term += kwargs['term']


    def odes(self, ode_list=None, specified_list=None):
        """
        Class where symbolic ODEs are provided.

        Provide a list of SymPy objects which are first order ODEs, in the form:
            x' - f(x, u, t)

        Also, if there are quantities which are specified functions of time
        (e.g., planetary orbits), that show up symbolically in the ODEs, they
        need to be identified and supplied in the specified_list argument.
        """

        t = self.time

        allvars_set = set()
        state_set = set()
        state_list = []
        der_list = []
        params_set = set()

        if specified_list:
            funcs_set = list(specified_list)
        else:
            funcs_set = set()

        for ode in ode_list:
            varss = [i for i in ode.atoms(AppliedUndef) if i.args==(t,)]
            allvars_set.update(varss)
            ders = [i for i in ode.atoms(Derivative) if i.expr in allvars_set]
            if len(ders) > 1:
                raise Exception('ODEs are not in first order form')
            der_list += ders
        for v in der_list:
            if v.expr in allvars_set:
                state_list += [v.expr]
                state_set.update([v.expr])
        # sorting list of states to be in order ode's were provided in
        for ode in ode_list:
            ders = [i for i in ode.atoms(Derivative) if i.expr in allvars_set]

        self._states = Matrix(state_list)
        self._state_ders = self._states.diff(t)

        odes = []
        for der in self._state_ders:
            ode = [o for o in ode_list if der in o][0]
            sign = ode.diff(der)
            if sign == 1:
                odes += [-1 * (ode - der)]
            elif sign == -1:
                odes += [ode + der]
            else:
                raise Exception("Invalid format for ODE")
        self._odes = Matrix(odes)

        controls_set = allvars_set - state_set - funcs_set
        self._controls = Matrix(list(controls_set))

        for ode in ode_list:
            funcs = [i for i in ode.atoms(AppliedUndef) if i.args!=(t,)]
            funcs_set.update(funcs)
        self._funcs = Matrix(list(funcs_set))

        for ode in ode_list:
            pars = ode.atoms(Symbol)
            params_set.update(pars)
        params_set -= {t,}
        self._params_args = list(params_set)


    def boundary_constraints(self, **kwargs):
        """
        Function to supply more constraints at the beginning and end of an arc.

        Valid keyword arguments are:
            initial : list
                A list of constraint equations for t=t_i. Can be function of
                states, controls, parameters, and functions.
            terminal : list
                A list of constraint equations for t=t_f. Can be function of
                states, controls, parameters, and functions.
            t_i : 1-tuple or 3-tuple
                If a fixed initial, supply a numeric value. If not fixed,
                supply the initial guess for the value and lower and upper
                boundaries on the initial time.
            t_f : 1-tuple or 3-tuple
                Same requirements for t_i, except for the final time for the arc.
        """

        if 't_i' in kwargs:
            v = kwargs['t_i']
            if len(v) == 1:
                self._t_i = (v[0], v[0], v[0])
            else:
                t_i = symbols(self.time.name + '_i')
                self._t_i = (t_i, v[1], v[2])
                self._params_args += [t_i]
                self._params = Matrix(list(self._params) + [t_i])
                self._param_limits[t_i] = (v[1], v[2])
                self._initial_guesses[t_i] = v[0]

        if 't_f' in kwargs:
            v = kwargs['t_f']
            if len(v) == 1:
                self._t_f = (v[0], v[0], v[0])
            else:
                t_f = symbols(self.time.name + '_f')
                self._t_f = (t_f, v[1], v[2])
                self._params_args += [t_f]
                self._params = Matrix(list(self._params) + [t_f])
                self._param_limits[t_f] = (v[1], v[2])
                self._initial_guesses[t_f] = v[0]

        i_to_add = []
        if 'left_internal' in kwargs and kwargs['left_internal'] is True:
            if self._t_i[1] != self._t_i[2]:
                i_to_add = [-t_i]
            self._left_internal = True
        t_to_add = []
        if 'right_internal' in kwargs and kwargs['right_internal'] is True:
            if self._t_f[1] != self._t_f[2]:
                t_to_add = [t_f]
            self._right_internal = True

        if 'initial' in kwargs:
            self._initial = Matrix(kwargs['initial'] + i_to_add)
            self._sanitize_bounds(self._initial, len(self._initial), 'initial', True)
        if 'terminal' in kwargs:
            self._terminal = Matrix(kwargs['terminal'] + t_to_add)
            self._sanitize_bounds(self._terminal, len(self._terminal), 'terminal', True)


    def path_constraints(self, **kwargs):
        """
        Function to supply path equality and inequaltiy constraints.

        Valid keyword arguments are:
            g : SymPy Matrix
                A SymPy matrix containing equality constraint expressions
                (requires g_l, g_u to be supplied). For equality constraints,
                set the appropriate entries in g_l, g_u to be equal.  Can be
                functions of parameters, time, states, controls, and functions.
            g_l : SymPy Matrix
                Represents lower bounds on path constraints. Should be a SymPy
                Matrix of m x 1, where m is the number of inequality constraint
                functions. Valid entries in the Matrix are: SymPy Symbols,
                constant numerical scalars, or valid SymPy expressions with no
                dependencies on states (e.g., a combination of a function of
                time and an argument).
            g_u : SymPy Matrix
                Represents upper bounds on path constraints. Follows same
                format as g_l.
        """

        if 'g' in kwargs:
            self._g = kwargs['g']
            self._sanitize_bounds(self._g, len(self._g), 'g', True)

        ng = len(self._g)
        self._g_l = Matrix([-oo] * ng)
        self._g_u = Matrix([ oo] * ng)
        if 'g_l' in kwargs:
            self._g_l = Matrix(kwargs['g_l'])
            self._sanitize_bounds(self._g_l, ng, 'g_l')
        if 'g_u' in kwargs:
            self._g_u = Matrix(kwargs['g_u'])
            self._sanitize_bounds(self._g_u, ng, 'g_u')


    def continuous_bounds(self, **kwargs):
        """
        Function to define state and control boundaries.

        Valid keyword arguments are:
            x_limits : dict
                Represents the lower/upper bounds on state variables. If not
                provided bounds are assumed to be (-oo, oo). Format is a
                dict: {state_1 : (lower, upper), state_2 : (lower, upper), ...}.
                The user is responsible for making sure all desired states are
                included.  Valid entries in the Matrix are: SymPy Symbols,
                constant numerical scalars, or valid SymPy expressions with no
                dependencies on "free quantities" (states, controls,
                parameters, etc.).
            u_limits : dict
                Represents lower/upper bounds on control variables; same input
                requirements as x_limits. User is responsible for identifying
                all desired controls are included.
        """

        nx = len(self._states)
        nu = len(self._controls)

        x_lower = []
        x_upper = []
        xlims = {}
        if 'x_limits' in kwargs:
            xlims = kwargs['x_limits']
        for s in self._states:
            if s in xlims:
                x_lower += [xlims[s][0]]
                x_upper += [xlims[s][1]]
            else:
                x_lower += [-oo]
                x_upper += [ oo]
        self._x_l = Matrix(x_lower)
        self._x_u = Matrix(x_upper)
        self._sanitize_bounds(self._x_l, nx, 'x_l')
        self._sanitize_bounds(self._x_u, nx, 'x_u')

        u_lower = []
        u_upper = []
        ulims = {}
        if 'u_limits' in kwargs:
            ulims = kwargs['u_limits']
        for c in self._controls:
            if c in ulims:
                u_lower += [ulims[c][0]]
                u_upper += [ulims[c][1]]
            else:
                u_lower += [-oo]
                u_upper += [ oo]
        self._u_l = Matrix(u_lower)
        self._u_u = Matrix(u_upper)
        self._sanitize_bounds(self._u_l, nu, 'u_l')
        self._sanitize_bounds(self._u_u, nu, 'u_u')


    def parameter_bounds(self, arguments=[], param_limits={}):
        """
        Function to define bounds on parameters and which symbols are fixed
        constants.

        Arguments are:
            arguments : list
                A list of symbolic quantities which are fixed needs to be
                supplied, to allow for proper identificiation of "free"
                parameters.
            param_limits : dict
                Represents the lower/upper bounds on parameters; same input If
                not provided bounds are assumed to be (-oo, oo). Format is a
                dict: {param_1 : (lower, upper), param_2 : (lower, upper),
                ...}. Valid entries in the Matrix are: SymPy Symbols, constant
                numerical scalars, or valid SymPy expressions with no
                dependencies on "free quantities" (states, controls,
                parameters, etc.). The user is responsible for identifying
                parameters.
        """

        self._param_limits.update(param_limits)
        temp_set = set(self._params_args)
        temp_set -= set(arguments)
        self._params = Matrix(list(temp_set))
        self._arguments = Matrix(arguments)


    """
    End of functions used to _symbolically_ define the problem. Functions after
    this break are used to define the relevant computational parts of the arc.

    At this point, the following have been defined:
        _states
        _controls
        _funcs
        _odes
        _t_i
        _t_f
        _x_l
        _x_u
        _u_l
        _u_u

    There has NOT been a distinction between parameters and arguments within
        _params_args

    """

    def __str__(self):
        """
        Prints out the current problem, as well as how all the symbolic
        quantities have been partitioned.
        """

        out_str = ""
        lead = '      '
        t = self.time
        ti = symbols(t.name + '_i')
        tf = symbols(t.name + '_f')
        out_str += '\nDefinition for %s\n\n' % self.name
        out_str += 'Objective Function\n'
        out_str += '  min:\n'
        obj_i = self._obj_init.subs(t, ti)
        obj_f = self._obj_term.subs(t, tf)
        out_str += lead + str(obj_i + obj_f) + '\n\n'

        out_str += 'Quantities\n'
        text = ['states', 'controls', 'functions', 'parameters', 'arguments']
        vals = [self._states, self._controls, self._funcs, self._params, self._arguments]
        for i in range(len(text)):
            out_str += '  %s:\n' % text[i]
            out_str += lead + str([v for v in vals[i]]) + '\n'
        out_str += '\n'

        out_str += 'Constraints\n'
        out_str += '  times:\n'
        tees = [(self._t_i, 't_i'), (self._t_f, 't_f')]
        for t in tees:
            #t[0][1] == t[0][2]
            if t[0][0] == 'fixed':
                out_str += lead + '%s = ' % t[1] + str(t[0][0]) + '\n'
            else:
                out_str += lead + str(t[0][1]) + ' ≤ %s ≤ ' % t[1] + str(t[0][2]) + '\n'
        out_str += '  state bounds:\n'
        out_str += self._print_bounds(self._x_l, self._states, self._x_u, lead) + '\n'
        out_str += '  control bounds:\n'
        out_str += self._print_bounds(self._u_l, self._controls, self._u_u, lead) + '\n'
        # TODO fix this part
        #if len(self._params) > 0:
        #    out_str += '  parameter bounds:\n'
        #    out_str += self._print_bounds(self._p_l, self._params, self._p_u, lead) + '\n'

        if self._initial != Matrix([[]]) and self._left_internal is False:
            out_str += '  initial conditions:\n'
            for i in self._initial:
                out_str += lead + '0 = %s\n' % str(i)

        if self._terminal != Matrix([[]]) and self._right_internal is False:
            out_str += '  terminal conditions:\n'
            for i in self._terminal:
                out_str += lead + '0 = %s\n' % str(i)

        if self._g != Matrix([[]]):
            out_str += '  path constraints:\n'
            out_str += self._print_bounds(self._g_l, self._g, self._g_u, lead) + '\n'
        return out_str

    __repr__ = __str__


    def numerical_defs(self, m, states_guess, controls_guess, parameters_guess,
                       arguments={}, time_grid=None, external_include=None):
        """
        Place to provide some numerical solution definitions of the arc.

        Arguments are:
            m : int
                The number of time points in the arc.
            states_guess : dict
                The initial guess for state values. Provided as a dictionary
                with each key being a symbolic variable and the value an
                interable of length m.
            controls_guess : dict
                The initial guess for control values. Follows the same format
                as states_guess.
            parameters_guess : dict
                The initial guess for parameter values. Provided as a
                dictionary with each key being a symbolic parameter and the
                value being the initial guess.
            arguments : dict
                A mapping of symbols (not variables) in the equations
                from symbol to number.
            time_grid : iterable
                An iterable over the interval (0, 1) which is used for a
                non-uniform time grid. Needs to have a length of m. If not
                provided, a uniform time grid is used.
            external_include : str
                A C include statement for external functions. These functions
                are expected to return a float, and their symbolic name and
                call signature must match the C name and call signature - they
                can only take in values, no pointers.
        """

        if not time_grid:
            self._time_grid = []
            for i in range(m):
                self._time_grid += [i / (m - 1)]
        else:
            if len(time_grid) != m:
                raise Exception("Length of time_grid does not equal m")
            self._time_grid = []
            i = time_grid[0]
            f = time_grid[-1]
            s = f - i
            for j in time_grid:
                self._time_grid += [j / s - i]

        self._m = m

        self._argument_values = {}
        for a in self._arguments:
            if a not in arguments:
                raise Exception("A defined argument did not have a value assigned")
            self._argument_values[a] = arguments[a]

        ti = self._t_i[0]
        tf = self._t_f[0]

        t_points = self._time_grid.copy()
        scale = tf - ti
        for i in range(m):
            t_points[i] = ti + scale * t_points[i]

        self._t_points = list(Matrix(t_points))

        if external_include:
            self.template_names['user_h_include'] = external_include
        else:
            self.template_names['user_h_include'] = ""

        funcs = self._funcs
        temp_wrap = dict()
        for f in funcs:
            temp_args = []
            for i in range(len(f.args)):
                temp_args += ['x[%d]' % i]
            base = "double UW_%s(double *x)" % str(type(f))
            temp_wrap[base] = (base + "\n{\n" +
             "  return %s(%s);\n" % (str(type(f)), ', '.join(temp_args)) +
             "}\n")
        if self.subarc is False:
            temp_str = "\n".join(temp_wrap.values())
            self.template_names['wrap_user_funcs'] = temp_str
        else:
            self.template_names['wrap_user_funcs'] = '\n'
            self._wrapped_funcs = temp_wrap

        n = len(self._states)
        guesses = []
        in_dict = self._initial_guesses
        in_dict.update(states_guess)
        in_dict.update(controls_guess)
        in_dict.update(parameters_guess)
        temp_dict = {}
        temp_list = list(self._states) + list(self._controls)

        for i in temp_list:
            if len(in_dict[i]) != len(self._time_grid):
                tg = [min(self._time_grid), max(self._time_grid)]
                temp_dict[i] = numpy.interp(self._time_grid, tg, in_dict[i])
            else:
                temp_dict[i] = in_dict[i]

        idx = 0
        for i in range(m):
            for j in temp_list:
                guesses += ["  x[%d] " + "= %s;\n" % ccode(S(temp_dict[j][i]))]
                idx += 1
        for p in self._params:
            guesses += ["  x[%d] " + "= %s;\n" % ccode(S(in_dict[p]))]
            idx += 1
        self._guesses = guesses


        g_lower = [0] * len(self._initial)
        g_upper = [0] * len(self._initial)
        for i in range(m - 1):
            g_lower += list(self._g_l)
            g_upper += list(self._g_u)
            g_lower += [0] * n
            g_upper += [0] * n
        g_lower += list(self._g_l)
        g_upper += list(self._g_u)
        g_lower += [0] * len(self._terminal)
        g_upper += [0] * len(self._terminal)
        cbd = []
        for i in range(len(g_lower)):
            temp = "  g_L[%d" + " + %d] = %s;\n" % (i, ccode(S(g_lower[i])))
            temp += "  g_U[%d" + " + %d] = %s;\n" % (i, ccode(S(g_upper[i])))
            cbd += [temp]
        self._constraints_bounds_definition = cbd

        p_lower = []
        p_upper = []
        for p in self._params:
            if p in self._param_limits:
                p_lower += [self._param_limits[p][0]]
                p_upper += [self._param_limits[p][1]]
            else:
                p_lower += [-oo]
                p_upper += [oo]
        self._p_l = Matrix(p_lower)
        self._p_u = Matrix(p_upper)


        lower = []
        upper = []
        for i in range(m):
            lower += [j for j in self._x_l]
            lower += [j for j in self._u_l]
            upper += [j for j in self._x_u]
            upper += [j for j in self._u_u]
        lower += [l for l in self._p_l]
        upper += [u for u in self._p_u]
        nlp_b_def = []
        for i in range(len(lower)):
            temp = "  x_L[%d" + " + %d] = %s;\n" % (i, ccode(S(lower[i])))
            temp += "  x_U[%d" + " + %d] = %s;\n" % (i, ccode(S(upper[i])))
            nlp_b_def += [temp]
        self._nlp_bounds_def = nlp_b_def


    def generate_c_files(self, **kwargs):
        """
        Generates functions for the NLP solver.

        Right now, only works for IPOPT.
        """

        t = self.time
        states = self._states
        controls = self._controls
        params = self._params
        args = self._arguments
        odes = self._odes
        funcs = self._funcs
        g = self._g

        t_points = self._t_points

        n = len(states)
        nu = len(controls)
        m = self._m
        ng = len(g)

        # Filling the initial values

        self.template_names['initial_nlp_values'] = ""
        for i in range(len(self._guesses)):
            self.template_names['initial_nlp_values'] += self._guesses[i] % i

        # the values for arguments (user controllable constants)
        arg_defs = ["  Number %s;\n" % str(a) for a in self._arguments]
        arg_struct_def = "struct MyUserData\n{\n%s};\n" % "".join(arg_defs)

        arg_struct_fill = ""
        for a in self._arguments:
            arg_struct_fill += "  Number %s = %s;\n" % (str(a), ccode(S(self._argument_values[a])))
            arg_struct_fill += "  user_data.%s = %s;\n" % (str(a), str(a))

        arg_struct_unpack = "".join(["  Number %s = my_data->%s;\n" % (str(a), str(a)) for a in self._arguments])


        self.template_names['arg_struct_def'] = arg_struct_def
        self.template_names['arg_struct_fill'] = arg_struct_fill

        if 'arg_struct_unpack' in kwargs:
            self.template_names['arg_struct_unpack'] = kwargs['arg_struct_unpack']
        else:
            self.template_names['arg_struct_unpack'] = arg_struct_unpack

        # The mapping from sympy to numeric
        X = DeferredVector('X')
        x = []
        for i in range(m):
            x += [j.subs(t, t_points[i]) for j in states]
            x += [j.subs(t, t_points[i]) for j in controls]
        x += [i for i in params]
        to_num = dict(zip(x, X))
        self.template_names['num_nlp_vars'] = len(x)

        nbd = ""
        for bd in self._nlp_bounds_def:
            nbd += bd % (0, 0)
        self.template_names['nlp_bounds_definition'] = nbd


        self.template_names['num_nlp_constraints'] = (len(self._initial) +
                                                      len(self._terminal) +
                                                      m * len(self._g) +
                                                      (m - 1) * n)

        self.template_names['constraints_bounds_definition'] = ""
        for cbd in self._constraints_bounds_definition:
            self.template_names['constraints_bounds_definition'] += cbd % (0, 0)


        t_points_num = list(Matrix(t_points).subs(to_num))
        # This gives the numerical values for t at each time point
        t_ccode =  ("    Number t[%d];\n" % len(t_points_num))
        t_ccode += ("    Number grid[%d];\n" % len(self._time_grid))
        idx = 0
        for i in t_points_num:
            t_ccode += "    t[%d] = %s;\n" % (idx, ccode(i))
            idx += 1
        idx = 0
        for i in self._time_grid:
            t_ccode += "    grid[%d] = %s;\n" % (idx, ccode(S(i)))
            idx += 1

        self.template_names['t_ccode'] = t_ccode




        """
        # Getting the objective function
        """
        obj = sympify(0)
        for i in [(self._obj_init, t_points[0]), (self._obj_term, t_points[-1])]:
            z = [(v, v.subs(t, i[1])) for v in list(states) + list(controls)]
            obj += i[0].subs(z)

        obj_num = obj.subs(to_num)
        objective_definition = "*obj_value = %s;\n" % ccode(obj_num)

        self.template_names['objective_definition'] = objective_definition








        """
        # Getting the gradient of the objective function
        """
        present = sorted([i._num for i in obj_num.atoms(DeferredSymbol)])
        obj_diffed = []
        for i in present:
            obj_diffed += [obj_num.diff(X[i])]
        temp_ders, out_string, max_args = pull_out_derivatives(obj_diffed, "  ")

        grad_f_definition = ""
        if len(temp_ders) > 0:
            grad_f_definition += ("  Number diffs[%d];\n" % len(temp_ders) +
                                  "  Number xtemp[%d];\n" % max_args)

            idx = 0
            for ostr in out_string:
                for o in ostr[:-1]:
                    grad_f_definition += o
                grad_f_definition += "  diffs[%d] = %s;\n" % (idx, ostr[-1])
                idx += 1

        for i in range(len(x)):
            if i not in present:
                grad_f_definition += "  grad_f[%d] = 0.0;\n" % i

        for i in range(len(present)):
            g_f_s = obj_diffed[i].subs(dict(zip(temp_ders, DeferredVector('diffs'))))
            grad_f_definition += "  grad_f[%d] = %s;\n" % (present[i], ccode(g_f_s))

        self.template_names['grad_f_definition'] = grad_f_definition










        """
        # Getting the constraint function
        """

        # c has structure init, g, d, g, d, g, d, g, term - where init is
        # intial constraints, term is terminal constraints, g are the path
        # constraints, and d are the defects

        num_nlp_g = self.template_names['num_nlp_constraints']
        # Initial/Terminal constraints
        g_def = ""
        initial  = self._initial.subs(t, t_points[0])
        terminal = self._terminal.subs(t, t_points[-1])
        init_ccode = [ccode(i.subs(to_num)) for i in initial]
        term_ccode = [ccode(i.subs(to_num)) for i in terminal]
        temp = ""
        if self._left_internal is True:
            temp = "+"
        for i in range(len(initial)):
            g_def += "  g[%d] %s= %s;\n" % (i, temp, init_ccode[i])
        if self._right_internal is True:
            temp = "+"
        for i in range(len(terminal)):
            g_def += "  g[%d] %s= %s;\n" % (num_nlp_g + i - len(terminal), temp, term_ccode[i])

        # Used for the path constraints
        L = DeferredVector('local')
        local_num_list = list(states) + list(controls) + list(params) + [t]
        local_num = dict(zip(local_num_list, L))

        g_ccode = [ccode(i.subs(local_num)) for i in g]
        g_def += "  Number local[%d];\n" % len(local_num_list)

        g_def += "  for (int iii = 0; iii < %d; iii++)\n  {\n" % m
        idx = 0
        for i in list(states) + list(controls):
            g_def += "    local[%d] = X[iii * %d + %d];\n" % (idx, n + nu, idx)
            idx += 1
        for i in range(len(params)):
            g_def += "    local[%d] = X[%d];\n" % (idx, len(x) - len(params) + i)
            idx += 1
        g_def += "    local[%d] = t[iii];\n" % idx
        for i in range(len(g_ccode)):
            g_def += "    g[iii * %d + %d + %d] = %s;\n" % (ng + n, i, len(init_ccode), g_ccode[i])
        g_def += "  }\n"



        # Used for the defects
        t_i = self._t_i[0]
        t_f = self._t_f[0]
        grid_l = symbols('grid_l')
        grid_r = symbols('grid_r')
        t_l = grid_l * (t_f - t_i) + t_i
        t_r = grid_r * (t_f - t_i) + t_i
        dt = t_r - t_l
        y_0 = self._rename_funcs(states, '_0').subs(t, t_l)
        y_1 = self._rename_funcs(states, '_1').subs(t, t_r)
        u_0 = self._rename_funcs(controls, '_0').subs(t, t_l)
        u_1 = self._rename_funcs(controls, '_1').subs(t, t_r)
        f_0 = odes.subs(list(zip(states, y_0)) + list(zip(controls, u_0))).subs(t, t_l)
        f_1 = odes.subs(list(zip(states, y_1)) + list(zip(controls, u_1))).subs(t, t_r)
        y_c = (y_0 + y_1) / 2 + dt / 8 * (f_0 - f_1)
        u_c = (u_0 + u_1) / 2
        t_c = (t_l + t_r) / 2
        defect = y_1 - y_0 - dt / 6 * (f_1 + f_0 + 4 *
                                       odes.subs(list(zip(states, y_c)) +
                                                 list(zip(controls, u_c))).subs(t, t_c))

        D = DeferredVector('d_local')
        defect_num_list = list(y_0) + list(u_0) + list(y_1) + list(u_1) + list(params) + [grid_l, grid_r]
        defect_num = dict(zip(defect_num_list, D))
        d_ccode = [ccode(i.subs(defect_num)) for i in defect]

        #d_locals = [(list(x[j * (n + nu) : (j + 2) * (n + nu)]) + pars +
        #                t_local[j : j + 2]) for j in range(m - 1)]

        g_def += "  Number d_local[%d];\n" % len(defect_num_list)

        g_def += "  for (int iii = 0; iii < %d; iii++)\n  {\n" % (m - 1)
        idx = 0
        for i in range(2 * len(list(states) + list(controls))):
            g_def += "    d_local[%d] = X[iii * %d + %d];\n" % (idx, n + nu, idx)
            idx += 1
        for i in range(len(params)):
            g_def += "    d_local[%d] = X[%d];\n" % (idx, len(x) - len(params) + i)
            idx += 1
        g_def += "    d_local[%d] = grid[iii];\n" % idx
        idx += 1
        g_def += "    d_local[%d] = grid[iii + 1];\n" % idx
        for i in range(len(d_ccode)):
            g_def += "    g[iii * %d + %d + %d] = %s;\n" % (n + ng, i, ng + len(init_ccode), d_ccode[i])
        g_def += "  }\n"

        self.template_names['g_definition'] = g_def













        """
        # Getting the Jacobian of the constraint function
        """

        jac_g_rows_cols = ""
        jac_g_values_only = ""
        jac_g_vals = ""
        jac_g_rc_idx = 0
        der_list = []
        diff_idx = 0
        jac_g_diffs = ""

        # Initial/Terminal constraints
        initial = self._initial.subs(t, t_points[0]).subs(to_num)
        init_ccode_j = []
        present_i_diffed = []

        for i in range(initial.rows):
            temp = sorted([j._num for j in initial[i].atoms(DeferredSymbol)])
            for j in range(len(temp)):
                jac_g_rows_cols += "    iRow[%d] = %d;\n" % (jac_g_rc_idx, i)
                self._subarc_jgr += ["    iRow[%d" % jac_g_rc_idx + " + %d] = " + "%d" % i + " + %d;\n"]
                jac_g_rows_cols += "    jCol[%d] = %d;\n" % (jac_g_rc_idx, temp[j])
                self._subarc_jgc += ["    jCol[%d" % jac_g_rc_idx + " + %d] = " + "%d" % temp[j] + " + %d;\n"]
                jac_g_rc_idx += 1
                present_i_diffed += [initial[i].diff(X[temp[j]])]
        init_ders, out_strs, max_args_i = pull_out_derivatives(present_i_diffed, "    ")
        der_list += init_ders
        for ostr in out_strs:
            for o in ostr[:-1]:
                jac_g_diffs += o
            jac_g_diffs += "    diffs[%d] = %s;\n" % (diff_idx, ostr[-1])
            diff_idx += 1
        diff_dict = dict(zip(der_list, DeferredVector('diffs')))
        for j in range(len(present_i_diffed)):
            init_ccode_j += [ccode(present_i_diffed[j].subs(diff_dict))]

        terminal = self._terminal.subs(t, t_points[-1]).subs(to_num)
        term_ccode_j = []
        present_t_diffed = []
        for i in range(terminal.rows):
            temp = sorted([j._num for j in terminal[i].atoms(DeferredSymbol)])
            for j in range(len(temp)):
                r = len(initial) + (m - 1) * n + m * len(g) + i
                jac_g_rows_cols += "    iRow[%d] = %d;\n" % (jac_g_rc_idx, r)
                self._subarc_jgr += ["    iRow[%d" % jac_g_rc_idx + " + %d] = " + "%d" % r + " + %d;\n"]
                c = temp[j]
                jac_g_rows_cols += "    jCol[%d] = %d;\n" % (jac_g_rc_idx, c)
                self._subarc_jgc += ["    jCol[%d" % jac_g_rc_idx + " + %d] = " + "%d" % c + " + %d;\n"]
                jac_g_rc_idx += 1
                present_t_diffed += [terminal[i].diff(X[temp[j]])]
        term_ders, out_strs, max_args_t = pull_out_derivatives(present_t_diffed)
        der_list += term_ders
        for ostr in out_strs:
            for o in ostr[:-1]:
                jac_g_diffs += o
            jac_g_diffs += "    diffs[%d] = %s;\n" % (diff_idx, ostr[-1], "    ")
            diff_idx += 1
        diff_dict = dict(zip(der_list, DeferredVector('diffs')))
        for j in range(len(present_t_diffed)):
            term_ccode_j += [ccode(present_t_diffed[j].subs(diff_dict))]



        # Path constraints
        path_list = list(states) + list(controls) + list(params)
        path_def_syms = [local_num[i] for i in path_list]
        D_g = g.subs(local_num).jacobian(path_def_syms)
        g_ders, out_strs, max_args_g = pull_out_derivatives(D_g, "      ")
        der_list += g_ders

        dg_ccode = []
        diff_dict = dict(zip(g_ders, DeferredVector('l_diffs')))
        for i in range(D_g.rows):
            temp = []
            for j in range(D_g.cols):
                temp_subs = D_g[i, j].subs(diff_dict)
                temp += [ccode(temp_subs)]
            dg_ccode += [temp]


        for j in range(m):
            for k in range(D_g.rows):
                for ii in range(D_g.cols):
                    r = len(initial) + j * (D_g.rows + n) + k
                    jac_g_rows_cols += "    iRow[%d] = %d;\n" % (jac_g_rc_idx, r)
                    self._subarc_jgr += ["    iRow[%d" % jac_g_rc_idx + " + %d] = " + "%d" % r + " + %d;\n"]
                    c = j * (n + 1) + ii
                    if ii >= n + nu:
                        c = (n + nu) * (m-1) + ii
                    jac_g_rows_cols += "    jCol[%d] = %d;\n" % (jac_g_rc_idx, c)
                    self._subarc_jgc += ["    jCol[%d" % jac_g_rc_idx + " + %d] = " + "%d" % c + " + %d;\n"]
                    jac_g_rc_idx += 1

        # only do this step if it's useful
        if len(g_ders) > 0:
            jac_g_diffs += "    for (int iii = 0; iii < %d; iii++);\n" % (m)
            jac_g_diffs += "    {\n"
            idx = 0
            for i in list(states) + list(controls):
                jac_g_diffs += "      local[%d] = X[iii * %d + %d];\n" % (idx, n + nu, idx)
                idx += 1
            for i in range(len(params)):
                jac_g_diffs += "      local[%d] = X[%d];\n" % (idx, len(x) - len(params) + i)
                idx += 1
            jac_g_diffs += "      local[%d] = t[iii];\n" % idx

            for ostr in out_strs:
                for o in ostr[:-1]:
                    jac_g_diffs += o
                jac_g_diffs += "      diffs[%d + iii * %d] = %s;\n" % (diff_idx, len(out_strs), ostr[-1])
                diff_idx += 1
            jac_g_diffs += "    }\n"
            diff_idx += (m - 1) * len(out_strs)



        # Defect parts
        defect_list = list(y_0) + list(u_0) + list(y_1) + list(u_1) + list(params)
        defect_def_syms = [defect_num[i] for i in defect_list]
        D_defect = defect.subs(defect_num).jacobian(defect_def_syms)


        d_ders, out_strs, max_args_d = pull_out_derivatives(D_defect, "      ")
        der_list += d_ders

        dd_ccode = []
        diff_dict = dict(zip(d_ders, DeferredVector('d_diffs')))
        for i in range(D_defect.rows):
            temp = []
            for j in range(D_defect.cols):
                temp_subs = D_defect[i, j].subs(diff_dict)
                temp += [ccode(temp_subs)]
            dd_ccode += [temp]

        for j in range(m - 1):
            for k in range(D_defect.rows):
                for ii in range(D_defect.cols):
                    r = len(initial) + D_g.rows + j * (D_g.rows + n) + k
                    jac_g_rows_cols += "    iRow[%d] = %d;\n" % (jac_g_rc_idx, r)
                    self._subarc_jgr += ["    iRow[%d" % jac_g_rc_idx + " + %d] = " + "%d" % r + " + %d;\n"]
                    c = j * (n + 1) + ii
                    if ii >= 2 * (n + nu):
                        c = (n + nu) * (m-2) + ii
                    jac_g_rows_cols += "    jCol[%d] = %d;\n" % (jac_g_rc_idx, c)
                    self._subarc_jgc += ["    jCol[%d" % jac_g_rc_idx + " + %d] = " + "%d" % c + " + %d;\n"]
                    jac_g_rc_idx += 1



        # only do this step if it's useful
        if len(d_ders) > 0:
            jac_g_diffs += "    for (int iii = 0; iii < %d; iii++)\n    {\n" % (m - 1)
            idx = 0
            for i in range(2 * len(list(states) + list(controls))):
                jac_g_diffs += "      d_local[%d] = X[iii * %d + %d];\n" % (idx, n + nu, idx)
                idx += 1
            for i in range(len(params)):
                jac_g_diffs += "      d_local[%d] = X[%d];\n" % (idx, len(x) - len(params) + i)
                idx += 1
            jac_g_diffs += "      d_local[%d] = grid[iii];\n" % idx
            idx += 1
            jac_g_diffs += "      d_local[%d] = grid[iii + 1];\n" % idx

            for ostr in out_strs:
                for o in ostr[:-1]:
                    jac_g_diffs += o
                jac_g_diffs += "      diffs[%d + iii * %d] = %s;\n" % (diff_idx, len(out_strs), ostr[-1])
                diff_idx += 1
            jac_g_diffs += "    }\n"
            diff_idx += (m - 2) * len(out_strs)


        num_nonzero_jac = jac_g_rc_idx
        self.template_names['jac_g_rows_cols'] = jac_g_rows_cols
        self.template_names['num_nonzero_jac'] = num_nonzero_jac

        temp_str =  ("    Number l_diffs[%d];\n" % len(g_ders) +
                     "    Number d_diffs[%d];\n" % len(d_ders) +
                     "    Number xtemp[%d];\n" %
                     max(max_args_i, max_args_t, max_args_g, max_args_d) +
                     "    Number local[%d];\n" % len(local_num_list) +
                     "    Number d_local[%d];\n" % len(defect_num_list) +
                     "    Number diffs[%d];\n" % diff_idx)

        jac_g_diffs = temp_str + "\n" + jac_g_diffs

        v_idx = 0
        for i in init_ccode_j:
            jac_g_values_only += "    values[%d] = %s;\n" % (v_idx, i)
            v_idx += 1
        for t in term_ccode_j:
            jac_g_values_only += "    values[%d] = %s;\n" % (v_idx, t)
            v_idx += 1
        jac_g_values_only += ("    for (int iii = 0; iii < %d; iii++)\n" % m +
                              "    {\n")
        idx = 0
        for i in list(states) + list(controls):
            jac_g_values_only += "      local[%d] = X[iii * %d + %d];\n" % (idx, n + nu, idx)
            idx += 1
        for i in range(len(params)):
            jac_g_values_only += "      local[%d] = X[%d];\n" % (idx, len(x) - len(params) + i)
            idx += 1
        jac_g_values_only += "      local[%d] = t[iii];\n" % idx
        for i in range(len(g_ders)):
            jac_g_values_only += "      l_diffs[%d] = diffs[%d + iii * %d];\n" % (i, len(init_ders) + len(term_ders) + i, len(g_ders))

        for i in range(len(dg_ccode)):
            for j in dg_ccode[i]:
                jac_g_values_only += "      values[%d + iii * %d] = %s;\n" % (v_idx, len(dg_ccode) * len(dg_ccode[0]), j)
                v_idx += 1
        jac_g_values_only += "    }\n"


        jac_g_values_only += "    for (int iii = 0; iii < %d; iii++)\n    {\n" % (m - 1)

        idx = 0
        for i in range(2 * len(list(states) + list(controls))):
            jac_g_values_only += "      d_local[%d] = X[iii * %d + %d];\n" % (idx, n + nu, idx)
            idx += 1
        for i in range(len(params)):
            jac_g_values_only += "      d_local[%d] = X[%d];\n" % (idx, len(x) - len(params) + i)
            idx += 1
        jac_g_values_only += "      d_local[%d] = grid[iii];\n" % idx
        idx += 1
        jac_g_values_only += "      d_local[%d] = grid[iii + 1];\n" % idx

        for i in range(len(d_ders)):
            jac_g_values_only += "      d_diffs[%d] = diffs[%d + iii * %d];\n" % (i, len(init_ders) + len(term_ders) + m * len(g_ders) + i, len(d_ders))

        for i in range(len(dd_ccode)):
            for j in dd_ccode[i]:
                jac_g_values_only += "      values[%d + iii * %d] = %s;\n" % (v_idx + len(dg_ccode) * len(dg_ccode[0]) * (m - 1), len(dd_ccode) * len(dd_ccode[0]), j)
                v_idx += 1
        jac_g_values_only += "    }\n"

        jac_g_values_only = jac_g_diffs + "\n\n" + jac_g_values_only
        self.template_names['jac_g_values_only'] = jac_g_values_only

        self.template_names['base_filename'] = self.name.replace(" ", "_")

        self.template_names['num_nonzero_hess'] = 0

        if self.subarc is True:
            self.template_names['hide_open'] = '/*'
            self.template_names['hide_close'] = '*/'
        else:
            self.template_names['hide_open'] = ''
            self.template_names['hide_close'] = ''


        from pyocp.ipopt_templates.makefile import IPOPT_makefile_template
        from pyocp.ipopt_templates.arc import IPOPT_arc_c_template, IPOPT_arc_h_template
        from pyocp.ipopt_templates.common import PyOCP_common_c, PyOCP_common_h


        c_code = IPOPT_arc_c_template.format(**self.template_names)
        h_code = IPOPT_arc_h_template.format(**self.template_names)
        m_code = IPOPT_makefile_template.format(**self.template_names)
        with open(self.template_names['base_filename'] + '.c', 'w') as f:
            f.write(c_code)
        with open(self.template_names['base_filename'] + '.h', 'w') as f:
            f.write(h_code)
        if self.subarc is False:
            with open('Makefile', 'w') as f:
                f.write(m_code)
            with open('PyOCP_common.h', 'w') as f:
                f.write(PyOCP_common_h)
            with open('PyOCP_common.c', 'w') as f:
                f.write(PyOCP_common_c)


    def run_and_read(self):
        os.system('make clean')
        os.system('make')
        os.system('./' + self.template_names['base_filename'])
        self.read_results()

    def read_results(self, filename=None, results=None):
        """
        File to read results and put them into a dictionary of lists.

        Takes in a filename (if non-auto-generated) and optionally, a list of
        results in the order output by the NLP solver.
        """

        if not results:
            temp = []
            if not filename:
                filename = self.template_names['base_filename'] + '_results.csv'
            with open(filename, 'r') as f:
                temp += f.readlines()

            nlp_vars_results = temp[0].split(',')[:-1]
        else:
            nlp_vars_results = results

        states = self._states
        controls = self._controls
        params = self._params
        m = self._m
        n = len(states)
        nu = len(controls)
        result_list = (list(states) + list(controls)) * m + list(params)
        store_list = {}
        for r in (list(states) + list(controls) + list(params)):
            store_list[r] = []
        for i in range(len(result_list)):
            store_list[result_list[i]] += [float(nlp_vars_results[i])]

        t = Matrix(self._t_points)
        if len(params) > 0:
            t = t.subs(dict(zip(params, nlp_vars_results[-len(params):])))

        store_list[self.time] = list(t)
        self.results = store_list





def pull_out_derivatives(has_derivatives, spacing=""):
    """
    Takes in an interable expressions with Derivatives, returns lists of
    relevant information.
    """

    temp_func = []
    temp_args = []
    temp_vars = []

    dont_use = []

    temp_ders = []
    for d in has_derivatives:
        temp_ders += list(d.atoms(Subs))
    temp_ders = list(set(temp_ders))
    for s in temp_ders:
        dont_use += [s.args[0]]
        d = s.args[0].copy()
        if d.__class__ != Derivative:
            raise Exception("Unexpected Subs in a numerical diff conversion")
        d._args = (d.args[0].subs(s.args[1][0], s.args[2][0]), s.args[2][0])
        temp_func += [d.expr.func]
        temp_args += [d.expr.args]
        temp_vars += [d.expr.args.index(d.args[1])]

    all_ders = temp_ders

    temp_ders = []
    for d in has_derivatives:
        temp_ders += list(d.atoms(Derivative))
    temp_ders = list(set(temp_ders))

    inds_to_remove = []
    for i in range(len(temp_ders)):
        td = temp_ders[i]
        if td in dont_use:
            inds_to_remove += [i]
            continue
        temp_func += [td.expr.func]
        temp_args += [td.expr.args]
        temp_vars += [td.expr.args.index(td.args[1])]
    for i in reversed(inds_to_remove):
        temp_ders.remove(temp_ders[i])

    all_ders += temp_ders

    out_string = []
    for i in range(len(all_ders)):
        temp_str = []
        for j in range(len(temp_args[i])):
            temp_str += [spacing + "xtemp[%d] = %s;\n" % (j, ccode(temp_args[i][j]))]
        temp_str += ["num_diff(&UW_%s, xtemp, %d, %d)" % (str(temp_func[i]), temp_vars[i], len(temp_args[i]))]
        out_string += [temp_str]
    max_args = max([len(i) for i in temp_args] + [0])
    return all_ders, out_string, max_args


