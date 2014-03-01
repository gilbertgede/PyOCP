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
        specified
            Symbolic quantities which are only functions of time and represent
            quantities which are not states, not controls, and are prescribed
            ahead of time; e.g. position of planetary bodies as a function of
            time.  It is assumed that these will be later supplied as:
            numerical functions, or sympy expressions.
        functions
            Symbolic quantities which can be functions of time, states,
            controls, specified, parameters, arguments, etc. Needs to be
            either later supplied as: numerical functions, or sympy
            expressions.

    This is still a work in progress...
    """

    def __init__(self, number, time_symbol, name=None):
        """
        Initializer for arc class.

        Needs to be supplied with a number, and optionally a name (e.g.,
        atmospheric phase). If no name is given, the arc number is used.
        """

        self._number = number
        if isinstance(time_symbol, Symbol):
            self._time = time_symbol
        else:
            raise Exception("Invalid time symbol supplied (use SymPy's symbol)")
        if name:
            if isinstance(name, str):
                self._name = name
            else:
                raise Exception("Invalid name supplied")
        else:
            self._name = "Arc " + str(number)


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
        spec_list = [i for i in self._specified]
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
                            (b in self._specs)):
                        if b.free_symbols == set([self._time]):
                            if b not in spec_list:
                                spec_list += [b]
                        else:
                            if b not in func_list:
                                func_list += [b]

        self._specified = Matrix(spec_list)
        self._params_args = Matrix(pararg_list)
        self._funcs = Matrix(func_list)


    def _print_bounds(self, l, x, u, lead):
        """
        Just used for some prettier printing.

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

        an integral term is not yet supported, and instead needs to be added as
        a state.

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
        funcs_set = set()
        params_set = set()

        if specified_list:
            spec_set = list(specified_list)
        else:
            spec_set = set()
        self._specified = Matrix(list(spec_set))

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

        controls_set = allvars_set - state_set - spec_set
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


    def bounds(self, **kwargs):
        """
        Function to define time, state, and control boundaries.

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
            x_l : SymPy Matrix
                Represents lower bounds on state variables. If not provided,
                assumed to be -oo. Should be a SymPy Matrix of n x 1, where n
                is the number of states. Valid entries in the Matrix are: SymPy
                Symbols, constant numerical scalars, SymPy functions of the
                arc's time variable which will be later substituted for, or
                valid SymPy expressions with no dependencies on states (e.g., a
                combination of a function of time and a parametric scaling
                factor).
            x_u : SymPy Matrix
                Represents upper bounds on state variables. If not provided,
                assumed to be oo. Same input requirements as x_l.
            u_l : SymPy Matrix
                Represents lower bounds on state variables. If not provided,
                assumed to be -oo. Same input requirements as x_l, except the
                shape neeeds to be len(u) x 1.
            u_u : SymPy Matrix
                Represents upper bounds on state variables. If not provided,
                assumed to be oo. Same input requirements as x_l, except the
                shape neeeds to be len(u) x 1.
        """

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
        bounds = (('x_l', -oo, nx), ('x_u', oo, nx),
                  ('u_l', -oo, nu), ('u_u', oo, nu))
        for v in bounds:
            if v[0] in kwargs:
                val = kwargs[v[0]]
                self._sanitize_bounds(val, v[2], v[0])
                self.__setattr__('_' + v[0], val)
            else:
                self.__setattr__('_' + v[0],
                                 Matrix([v[1] for i in range(v[2])]))


    def constraints(self, **kwargs):
        """
        Function to supply more complex constraints.

        Allows initial, terminal, and path equality and inequaltiy constraints.

        Valid keyword arguments are:
            initial : list
                A list of constraint equations for t=t_i. Can be function of
                states, controls, parameters, and specifieds.
            terminal : list
                A list of constraint equations for t=t_f. Can be function of
                states, controls, parameters, and specifieds.
            g : SymPy Matrix
                A SymPy matrix containing equality constraint expressions
                (requires g_l, g_u to be supplied). For equality constraints,
                set the appropriate entries in g_l, g_u to be equal.  Can be
                functions of parameters, time, states, controls, specified
                functions, and external functions.
            g_l : SymPy Matrix
                Represents lower bounds on path constraints. Should be a SymPy
                Matrix of m x 1, where m is the number of inequality constraint
                functions. Valid entries in the Matrix are: SymPy Symbols,
                constant numerical scalars, SymPy functions of the arc's time
                variable which will be later substituted for, or valid SymPy
                expressions with no dependencies on states (e.g., a combination
                of a function of time and a parametric scaling factor).
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


    """
    End of functions used to _symbolically_ define the problem. Functions after
    this break are used to define the relevant computational parts of the arc.

    At this point, the following have been defined:
        _states
        _controls
        _specified
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
            print('  states:')
            print(lead + str([i for i in self._states]))
            print('  controls:')
            print(lead + str([i for i in self._controls]))
            print('  specifieds:')
            print(lead + str([i for i in self._specified]))
            print('  functions:')
            print(lead + str([i for i in self._funcs]))
            print('  parameters and arguments:')
            print(lead + str([i for i in self._params_args]))
            print('')
            print('Constraints')
            print('  times:')
            t_i = self._t_i
            if t_i[0] == 'fixed':
                print(lead + 't_i = ' + str(t_i[1]))
            else:
                print(lead + str(t_i[1]) + ' ≤ t_i ≤ ' + str(t_i[2]))
            t_f = self._t_f
            if t_f[0] == 'fixed':
                print(lead + 't_f = ' + str(t_f[1]))
            else:
                print(lead + str(t_f[1]) + ' ≤ t_f ≤ ' + str(t_f[2]))
            print('  state bounds:')
            print(self._print_bounds(self._x_l, self._states, self._x_u, lead))
            print('  control bounds:')
            print(self._print_bounds(self._u_l, self._controls, self._u_u, lead))
            if self._initial != Matrix([[]]):
                print('  initial conditions:')
                strs = [str(i) for i in self._initial]
                l = max([len(i) for i in strs])
                for s in strs:
                    print(lead + '0 = ' + '%ls' % s)
            if self._terminal != Matrix([[]]):
                print('  terminal conditions:')
                strs = [str(i) for i in self._terminal]
                l = max([len(i) for i in strs])
                for s in strs:
                    print(lead + '0 = ' + '%ls' % s)
            if self._g != Matrix([[]]):
                print('  path constraints:')
                print(self._print_bounds(self._g_l, self._g, self._g_u, lead))
        else:
            pass


    def numerical_defs(self, m, args=None, time_grid=None):
        """
        Place to provide some numerical solution defitions of the arc.

        Arguements are:
            m : int
                The number of time points in the arc.
            args : dict
                A mapping of which symbols (not variables) in the equations
                from symbol to number.
            time_grid : iteratable
                An iterable over the interval (0, 1) which is used for a
                non-uniform time grid. Needs to have a length of m. If not
                provided, a uniform time grid is used.
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

        if args:
            for i in args:
                if float(args[i]) != args[i]:
                    raise Exception("Non-float values supplied as an argument")
            params_args = set(self._params_args)
            arg_keys = set(args.keys())
            arg_vals = set(args.values())
            self._params = Matrix(list(params_args - arg_keys))
            self._args = Matrix(list(arg_keys))
            self._args_vals = args
        else:
            self._params = self._params_args
            self._args = Matrix([])
            self._args_vals = dict()


    def generate_funcs(self):
        """
        Generates functions for the NLP solver.

        Right now, only works for IPOPT, producing sparse matrices, and does so
        in a slow fashion (Python, not C).
        """

        t = self._time
        states = self._states
        controls = self._controls
        params = self._params
        args = self._args
        odes = self._odes
        g = self._g

        n = len(states)
        nu = len(controls)
        m = self._m
        ng = len(self._g)

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
        dt_points = []
        for i in range(m - 1):
            dt_points += [t_points[i + 1] - t_points[i]]

        x = []
        for i in range(m):
            for j in states:
                x += [j.subs(t, t_points[i])]
            for j in controls:
                x += [j.subs(t, t_points[i])]
        for i in params:
            x += [i]

        X = DeferredVector('X')
        to_num = dict(zip(x, X))
        to_num.update(self._args_vals)
        t_lam = lambdify(X, list(Matrix(t_points).subs(to_num)))
        dt_lam = lambdify(X, list(Matrix(dt_points).subs(to_num)))

        """
        # Getting the objective function
        """
        obj = sympify(0)
        for i in [(self._obj_init, t_points[0]), (self._obj_term, t_points[-1])]:
            z = []
            for v in [states, controls]:
                z += list(zip(v, v.subs(t, i[1])))
            obj += i[0].subs(z)

        obj_num = obj.subs(to_num)
        f = lambdify(X, obj_num)

        """
        # Getting the gradient of the objective function
        """
        present = [i._num for i in obj_num.atoms(DeferredSymbol)]
        ret = []
        for i in present:
            ret += [lambdify(X, obj_num.diff(X[i]))]
        def f_x(x):
            f_x = np.zeros(len(x))
            for i, v in enumerate(present):
                f_x[v] = ret[i](x)
            return f_x

        """
        # Getting the constraint function
        """
        # Initial/Terminal constraints
        initial = self._initial.subs(t, t_points[0])
        terminal =self._terminal.subs(t, t_points[-1])
        init_lam = [lambdify(X, i.subs(to_num)) for i in initial]
        term_lam = [lambdify(X, i.subs(to_num)) for i in terminal]

        # Used for the path constraints
        L = DeferredVector('L')
        local_num = dict(zip(list(states) + list(controls) + list(params) + [t], L))
        local_num.update(self._args_vals)
        g_lam = [lambdify(L, i.subs(local_num)) for i in g]

        # Used for the defects
        dt = symbols('dt')
        t_l = symbols(t.name + '_l')
        t_r = symbols(t.name + '_r')
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
        D = DeferredVector('D')
        defect_num = dict(zip(list(y_0) + list(u_0) + list(y_1) + list(u_1) +
                              list(params) + [t_l, t_r] + [dt], D))
        defect_num.update(self._args_vals)
        d_lam = [lambdify(D, i.subs(defect_num)) for i in defect]

        # c has structure init, g, d, g, d, g, d, g, term - where init is
        # intial constraints, term is terminal constraints, g are the path
        # constraints, and d are the defects

        def c(x):
            c = np.zeros(len(initial) + m * ng + n * (m - 1) + len(terminal))
            # make time grid, dt grid - done here in case there is any time
            # dependence on parameters (e.g. free final time)
            t_local = t_lam(x)
            dt_local = dt_lam(x)
            if len(params) == 0:
                pars = []
            else:
                pars = list(x[-len(params):])
            # Fill the return vector
            i = 0
            for j in init_lam:
                c[i] = j(x)
                i += 1
            for j in range(m):
                local = (list(x[j * (n + nu) : (j + 1) * (n + nu)]) +
                         pars + [t_local[j]])
                for k in g_lam:
                    c[i] = k(local)
                    i += 1
                if j == m - 1:      # if we're at the last time point, no
                    break           # defects are computed...
                d_local = (list(x[j * (n + nu) : (j + 2) * (n + nu)]) +
                           pars + t_local[j : j + 2] + [dt_local[j]])
                for k in d_lam:
                    c[i] = k(d_local)
                    i += 1
            for j in term_lam:
                c[i] = j(x)
                i += 1
            return c


        return f, f_x, c





