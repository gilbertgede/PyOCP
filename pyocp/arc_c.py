from .common import *
import re

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

    template_names = {}

    def __init__(self, number, time_symbol, name=None):
        """
        Initializer for arc class.

        Needs to be supplied with a number, and optionally a name (e.g.,
        "atmospheric_phase"). If no name is given, the arc number is used.
        """

        self._number = number
        self._time = sympify(time_symbol)
        if name:
            self._name = name
        else:
            self._name = "Arc " + str(number)
        self._constraints_set = False


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
            syms = v.atoms(Symbol) - {self._time,}
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

        self._params_args = Matrix(pararg_list)
        self._funcs = Matrix(func_list)


    def _print_bounds(self, l, x, u, lead):
        """
        Just used for some prettier printing.

        e.g.

          0  <=  x(t)  <=      1
         10  <=  y(t)  <=  v * 2
          1  <=  u(t)  <=      3

        """

        ls = repr(l).split('\n')
        xs = repr(x).split('\n')
        us = repr(u).split('\n')

        if len(xs) == 1:
            outstr = lead
            for v in [ls, xs, us]:
                start = 'Matrix([['
                end = ')]]'
                res = v[0][len(start):-len(end)]
                outstr += res + ' ≤ '
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

        self._obj_init = sympify(0)
        if 'init' in kwargs:
            self._obj_init += kwargs['init']
        self._obj_term = sympify(0)
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

        t = self._time

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
        self._params_args = Matrix(list(params_set))


    def constraints(self, **kwargs):
        """
        Function to supply more complex constraints.

        Allows initial, terminal, and path equality and inequaltiy constraints.

        Valid keyword arguments are:
            initial : list
                A list of constraint equations for t=t_i. Can be function of
                states, controls, parameters, and functions.
            terminal : list
                A list of constraint equations for t=t_f. Can be function of
                states, controls, parameters, and functions.
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

        times = ('initial', 'terminal')
        for v in times:
            if v in kwargs:
                val = Matrix(kwargs[v])
                self._sanitize_bounds(val, len(val), v, True)
                self.__setattr__('_' + v, val)
            else:
                self.__setattr__('_' + v, Matrix([[]]))

        if 'g' in kwargs:
            val = kwargs['g']
            self._sanitize_bounds(val, len(val), 'g', True)
            self.__setattr__('_g', val)
        else:
            self.__setattr__('_g', Matrix([[]]))

        ng = len(self._g)
        bounds = (('g_l', -oo, ng), ('g_u', oo, ng))
        for v in bounds:
            if v[0] in kwargs:
                val = kwargs[v[0]]
                self._sanitize_bounds(val, v[2], v[0])
                self.__setattr__('_' + v[0], val)
            else:
                self.__setattr__('_' + v[0], Matrix([v[1]] * v[2]))

        self._constraints_set = True


    def bounds(self, arguments, **kwargs):
        """
        Function to define time, state, control, and parameter boundaries.

        Arguments are:
            arguments : list
                A list of symbolic quantities which are fixed needs to be
                supplied, to allow for proper identificiation of "free"
                parameters.

        Valid keyword arguments are:
            t_i : 2-tuple
                Represents upper and lower bounds on initial time. If the 2
                values are not equal, value is free. Tuple entries need to be
                sympifiable, and not contain states, controls, specified, or
                parameters. If not provided, assumed to be fixed at 0.
            t_f : 2-tuple
                Represents upper and lower bounds on final time. If the 2
                values are not equal, value is free. Tuple entries need to be
                sympifiable, and not contain states, controls, specified, or
                parameters. If not provided, assumed to be fixed at 1.
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
            param_limits : dict
                Represents the lower/upper bounds on parameters; same input
                requirements as x_limits.  The user is responsible for
                identifying parameters.
        """

        if self._constraints_set is not True:
            raise("Constraints need to be defined before bounds are defined.")

        times = (('t_i', 0), ('t_f', 1))
        for v in times:
            if v[0] in kwargs:
                val = kwargs[v[0]]
                flag = 'free'
                if val[0] == val[1]:
                    flag = 'fixed'
                self._sanitize_bounds(Matrix(val), 2, v[0])
                self.__setattr__('_' + v[0], (flag, val[0], val[1]))
            else:
                self.__setattr__('_' + v[0], ('fixed', v[1], v[1]))

        nx = len(self._states)
        nu = len(self._controls)

        x_lower = []
        x_upper = []
        try:
            xlims = kwargs.pop('x_limits')
        except:
            xlims = {}
        for s in self._states:
            if s in xlims:
                x_lower += [xlims[s][0]]
                x_upper += [xlims[s][1]]
            else:
                x_lower += [-oo]
                x_upper += [oo]
        self._x_l = Matrix(x_lower)
        self._x_u = Matrix(x_upper)
        self._sanitize_bounds(self._x_l, nx, 'x_l')
        self._sanitize_bounds(self._x_u, nx, 'x_u')

        u_lower = []
        u_upper = []
        try:
            ulims = kwargs.pop('u_limits')
        except:
            ulims = {}
        for c in self._controls:
            if c in ulims:
                u_lower += [ulims[c][0]]
                u_upper += [ulims[c][1]]
            else:
                u_lower += [-oo]
                u_upper += [oo]
        self._u_l = Matrix(u_lower)
        self._u_u = Matrix(u_upper)
        self._sanitize_bounds(self._u_l, nu, 'u_l')
        self._sanitize_bounds(self._u_u, nu, 'u_u')


        temp_set = set(list(self._params_args))
        temp_set -= set(arguments)
        self._params = Matrix(list(temp_set))
        self._arguments = Matrix(arguments)

        p_lower = []
        p_upper = []
        try:
            plims = kwargs.pop('param_limits')
        except:
            plims = {}
        for p in self._params:
            if p in plims:
                p_lower += [plims[p][0]]
                p_upper += [plims[p][1]]
            else:
                p_lower += [-oo]
                p_upper += [oo]
        self._p_l = Matrix(p_lower)
        self._p_u = Matrix(p_upper)


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

    def problem_definition(self, **kwargs):
        """
        Prints out the current problem, as well as how all the symbolic
        quantities have been partitioned.

        If no keyword arguments are provided, simply prints out the problem.
        If keyword args are provided, checks the changes between symbolic
        classification, and regroups the variables as instructed. Use the
        following kwargs:

            to be developed...

        """

        if len(kwargs) == 0:
            lead = '      '
            t = self._time
            ti = symbols(t.name + '_i')
            tf = symbols(t.name + '_f')
            print('Definition for ' + self._name + '\n')
            print('Objective Function')
            print('  min:')
            obj_i = self._obj_init.subs(t, ti)
            obj_f = self._obj_term.subs(t, tf)
            print(lead + str(obj_i + obj_f))
            print('')

            print('Quantities')
            text = ['states', 'controls', 'functions', 'parameters', 'arguments']
            vals = [self._states, self._controls, self._funcs, self._params, self._arguments]
            for i in range(len(text)):
                print('  %s:' % text[i])
                print(lead + str([v for v in vals[i]]))
            print('')

            print('Constraints')
            print('  times:')
            tees = [(self._t_i, 't_i'), (self._t_f, 't_f')]
            for t in tees:
                if t[0][0] == 'fixed':
                    print(lead + '%s = ' % t[1] + str(t[0][1]))
                else:
                    print(lead + str(t[0][1]) + ' ≤ %s ≤ ' % t[1] + str(t_i[2]))
            print('  state bounds:')
            print(self._print_bounds(self._x_l, self._states, self._x_u, lead))
            print('  control bounds:')
            print(self._print_bounds(self._u_l, self._controls, self._u_u, lead))
            print('  parameter bounds:')
            print(self._print_bounds(self._p_l, self._params, self._p_u, lead))
            init_term = [(self._initial, 'initial conditions'), (self._terminal, 'terminal conditions')]
            for it in init_term:
                if it != Matrix([[]]):
                    print('  %s:' % it[1])
                    strs = [str(i) for i in it[0]]
                    for s in strs:
                        print(lead + '0 = %s' % s)
            if self._g != Matrix([[]]):
                print('  path constraints:')
                print(self._print_bounds(self._g_l, self._g, self._g_u, lead))
            print('')
        else:
            pass


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
                raise Exception("A defined arguement did not have a value assigned")
            self._argument_values[a] = arguments[a]

        if self._t_i[0] == 'fixed':
            ti = self._t_i[1]
        else:
            ti = symbols(self._time.name + '_i')
        if self._t_f[0] == 'fixed':
            tf = self._t_f[1]
        else:
            tf = symbols(self._time.name + '_f')

        t_points = self._time_grid
        scale = tf - ti
        for i in range(m):
            t_points[i] = ti + scale * t_points[i]

        self._t_points = t_points

        if external_include:
            self.template_names['user_h_include'] = external_include
        else:
            self.template_names['user_h_include'] = ""

        funcs = self._funcs
        temp_wrap = ""
        for f in funcs:
            temp_args = []
            for i in range(len(f.args)):
                temp_args += ['*x[%d]' % i]
            temp_wrap += (
             "double UW_%s(double *x)\n" % str(type(f))+
             "{\n" +
             "  return %s(%s);\n" % (str(type(f)), ', '.join(temp_args)) +
             "}\n")
        self.template_names['wrap_user_funcs'] = temp_wrap

        guesses = ""
        temp_dict = states_guess.copy()
        temp_dict.update(controls_guess)
        temp_list = list(self._states) + list(self._controls)
        idx = 0
        for i in range(m):
            for j in temp_list:
                guesses += "  x[%d] = %s;\n" % (idx, ccode(S(temp_dict[j][i])))
                idx += 1
        for p in self._params:
            guesses += "  x[%d] = %s;\n" % (idx, ccode(S(parameters_guess[p])))
            idx += 1
        self.template_names['initial_nlp_values'] = guesses


    def generate_c_files(self):
        """
        Generates functions for the NLP solver.

        Right now, only works for IPOPT, producing sparse matrices, and does so
        in a slow fashion (Python, not C).
        """

        t = self._time
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

        # the values for arguments (user controllable constants)
        arg_struct_def = "struct MyUserData\n{\n"
        for i in self._arguments:
            arg_struct_def += "  Number %s;\n" % str(i)
        arg_struct_def += "  };\n"

        arg_struct_fill = ""
        for i in self._arguments:
            arg_struct_fill += "  Number %s = %s;\n" % (str(i), ccode(S(self._argument_values[i])))
            arg_struct_fill += "  user_data.%s = %s;\n" % (str(i), str(i))

        arg_struct_unpack = ""
        for i in self._arguments:
            arg_struct_unpack += "  Number %s = my_data->%s;\n" % (str(i), str(i))

        self.template_names['arg_struct_def'] = arg_struct_def
        self.template_names['arg_struct_fill'] = arg_struct_fill
        self.template_names['arg_struct_unpack'] = arg_struct_unpack


        # The mapping from sympy to numeric
        lower = []
        upper = []
        X = DeferredVector('X')
        x = []
        for i in range(m):
            x += [j.subs(t, t_points[i]) for j in states]
            x += [j.subs(t, t_points[i]) for j in controls]
            lower += [j for j in self._x_l]
            lower += [j for j in self._u_l]
            upper += [j for j in self._x_u]
            upper += [j for j in self._u_u]
        x += [i for i in params]
        lower += [l for l in self._p_l]
        upper += [u for u in self._u_u]
        to_num = dict(zip(x, X))

        nlp_bounds_definition = ""
        for i in range(len(lower)):
            nlp_bounds_definition += "  x_L[%d] = %s;\n" % (i, ccode(S(lower[i])))
            nlp_bounds_definition += "  x_U[%d] = %s;\n" % (i, ccode(S(upper[i])))

        self.template_names['num_nlp_vars'] = len(x)
        self.template_names['nlp_bounds_definition'] = nlp_bounds_definition


        self.template_names['num_nlp_constraints'] = (len(self._initial) +
                                                      len(self._terminal) +
                                                      m * len(self._g) +
                                                      (m - 1) * n)
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
        cbd = ""
        for i in range(len(g_lower)):
            cbd += "  g_L[%d] = %s;\n" % (i, ccode(S(g_lower[i])))
            cbd += "  g_U[%d] = %s;\n" % (i, ccode(S(g_upper[i])))
        self.template_names['constraints_bound_definition'] = cbd


        t_points_num = list(Matrix(t_points).subs(to_num))
        # This gives the numerical values for t at each time point
        t_ccode = ("    Number *t;\n" +
                   "    t = (Number*)malloc(sizeof(Number)*n);\n")
        idx = 0
        for i in t_points_num:
            t_ccode += "    t[%d] = " % idx + ccode(i) + ";\n"
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
                                  "  Number *xtemp;\n" +
                                  "  xtemp = (Number*)malloc(sizeof(Number) * %d);\n"
                                  % max_args)

            idx = 0
            for os in out_string:
                for o in os[:-1]:
                    grad_f_definition += o
                grad_f_definition += "  diffs[%d] = %s;\n" % (idx, os[-1])
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
        for i in range(len(initial)):
            g_def += "  g[%d] = %s;\n" % (i, init_ccode[i])
        for i in range(len(terminal)):
            g_def += "  g[%d] = %s;\n" % (num_nlp_g + i - len(terminal), term_ccode[i])

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
        t_l = symbols(t.name + '_l')
        t_r = symbols(t.name + '_r')
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
        defect_num_list = list(y_0) + list(u_0) + list(y_1) + list(u_1) + list(params) + [t_l, t_r]
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
        g_def += "    d_local[%d] = t[iii];\n" % idx
        idx += 1
        g_def += "    d_local[%d] = t[iii + 1];\n" % idx
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
                jac_g_rows_cols += "    jCol[%d] = %d;\n" % (jac_g_rc_idx, temp[j])
                jac_g_rc_idx += 1
                present_i_diffed += [initial[i].diff(X[temp[j]])]
        init_ders, out_strs, max_args_i = pull_out_derivatives(present_i_diffed, "    ")
        der_list += init_ders
        for os in out_strs:
            for o in os[:-1]:
                grad_f_definition += o
            grad_f_definition += "    diffs[%d] = %s;\n" % (diff_idx, os[-1])
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
                jac_g_rows_cols += "    iRow[%d] = %d;\n" % (jac_g_rc_idx, len(initial) + (m - 1) * n + m * len(g) + i)
                jac_g_rows_cols += "    jCol[%d] = %d;\n" % (jac_g_rc_idx, temp[j])
                jac_g_rc_idx += 1
                present_t_diffed += [terminal[i].diff(X[temp[j]])]
        term_ders, out_strs, max_args_t = pull_out_derivatives(present_t_diffed)
        der_list += term_ders
        for os in out_strs:
            for o in os[:-1]:
                grad_f_definition += o
            grad_f_definition += "    diffs[%d] = %s;\n" % (diff_idx, os[-1], "    ")
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
        # TODO This isn't right - as der_list isn't the full length of diffs
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
                    c = j * (n + 1) + ii
                    if ii >= n + nu:
                        c = (n + nu) * m + ii
                    jac_g_rows_cols += "    jCol[%d] = %d;\n" % (jac_g_rc_idx, c)
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

            for os in out_strs:
                for o in os[:-1]:
                    grad_f_definition += o
                grad_f_definition += "      diffs[%d] = %s;\n" % (diff_idx, os[-1])
                diff_idx += 1
            jac_g_diffs += "    }\n"



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

        for j in range(m):
            for k in range(D_defect.rows):
                for ii in range(D_defect.cols):
                    r = len(initial) + D_g.rows + j * (D_g.rows + n) + k
                    jac_g_rows_cols += "    iRow[%d] = %d;\n" % (jac_g_rc_idx, r)
                    c = j * (n + 1) + ii
                    if ii >= 2 * (n + nu):
                        c = (n + nu) * m + ii
                    jac_g_rows_cols += "    jCol[%d] = %d;\n" % (jac_g_rc_idx, c)
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
            jac_g_diffs += "      d_local[%d] = t[iii];\n" % idx
            idx += 1
            jac_g_diffs += "      d_local[%d] = t[iii + 1];\n" % idx

            for os in out_strs:
                for o in os[:-1]:
                    grad_f_definition += o
                grad_f_definition += "      diffs[%d] = %s;\n" % (diff_idx, os[-1])
                diff_idx += 1
            jac_g_diffs += "    }\n"


        num_nonzero_jac = jac_g_rc_idx
        self.template_names['jac_g_rows_cols'] = jac_g_rows_cols
        self.template_names['num_nonzero_jac'] = num_nonzero_jac

        temp_str =  ("    Number l_diffs[%d];\n" % len(g_ders) +
                     "    Number d_diffs[%d];\n" % len(d_ders) +
                     "    Number *xtemp;\n" +
                     "    xtemp = (Number*)malloc(sizeof(Number) * %d);\n" %
                     max(max_args_i, max_args_t, max_args_g, max_args_d) +
                     "    Number local[%d];\n" % len(local_num_list) +
                     "    Number d_local[%d];\n" % len(defect_num_list))

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
            jac_g_values_only += "      l_diffs[%d] = diffs[%d + iii * %d];\n" % (i, len(init_ders) + len(term_ders), len(g_ders))

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
        jac_g_values_only += "      d_local[%d] = t[iii];\n" % idx
        idx += 1
        jac_g_values_only += "      d_local[%d] = t[iii + 1];\n" % idx

        for i in range(len(d_ders)):
            jac_g_values_only += "      d_diffs[%d] = diffs[%d + iii * %d];\n" % (i, len(init_ders) + len(term_ders) + m * len(d_ders), len(d_ders))

        for i in range(len(dd_ccode)):
            for j in dd_ccode[i]:
                jac_g_values_only += "      values[%d + iii * %d] = %s;\n" % (v_idx + len(dg_ccode) * len(dg_ccode[0]) * (m - 1), len(dd_ccode) * len(dd_ccode[0]), j)
                v_idx += 1
        jac_g_values_only += "    }\n"

        jac_g_values_only = jac_g_diffs + "\n\n" + jac_g_values_only
        self.template_names['jac_g_values_only'] = jac_g_values_only

        self.template_names['base_filename'] = self._name.replace(" ", "_")

        self.template_names['num_nonzero_hess'] = 0


        from .ipopt_template import IPOPT_c_template, IPOPT_h_template, IPOPT_makefile_template

        c_code = IPOPT_c_template.format(**self.template_names)
        h_code = IPOPT_h_template.format(**self.template_names)
        m_code = IPOPT_makefile_template.format(**self.template_names)
        with open(self.template_names['base_filename'] + '.c', 'w') as f:
            f.write(c_code)
        with open(self.template_names['base_filename'] + '.h', 'w') as f:
            f.write(h_code)
        with open('Makefile', 'w') as f:
            f.write(m_code)








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
        temp_strs = []
        for j in range(len(temp_args[i])):
            temp_strs += [spacing + "xtemp[%d] = %s;\n" % (j, ccode(temp_args[i][j]))]
        temp_str = ["num_diff(&UW_%s, xtemp, %d, %d)" % (str(temp_func[i]), temp_vars[i], len(temp_args[i]))]
        out_string += [temp_str]
    max_args = max([len(i) for i in temp_args] + [0])
    return all_ders, out_string, max_args


