r"""
Mixed integer linear programming
"""

include "../ext/stdsage.pxi"

class MixedIntegerLinearProgram:
    r"""
    The ``MixedIntegerLinearProgram`` class is the link between Sage, linear
    programming (LP) and  mixed integer programming (MIP) solvers. See the
    Wikipedia article on
    `linear programming <http://en.wikipedia.org/wiki/Linear_programming>`_
    for further information. A mixed integer program consists of variables,
    linear constraints on these variables, and an objective function which is
    to be maximised or minimised under these constraints. An instance of
    ``MixedIntegerLinearProgram`` also requires the information on the
    direction of the optimization.

    A ``MixedIntegerLinearProgram`` (or ``LP``) is defined as a maximization
    if ``maximization=True`` and is a minimization if ``maximization=False``.

    INPUT:

    - ``maximization``

      - When set to ``True`` (default), the ``MixedIntegerLinearProgram`` is
        defined as a maximization.
      - When set to ``False``, the ``MixedIntegerLinearProgram`` is defined as
        a minimization.

    EXAMPLES::

         sage: ### Computation of a maximum stable set in Petersen's graph
         sage: g = graphs.PetersenGraph()
         sage: p = MixedIntegerLinearProgram(maximization=True)
         sage: b = p.new_variable()
         sage: p.set_objective(sum([b[v] for v in g]))
         sage: for (u,v) in g.edges(labels=None):
         ...       p.add_constraint(b[u] + b[v], max=1)
         sage: p.set_binary(b)
         sage: p.solve(objective_only=True)     # optional - requires Glpk or COIN-OR/CBC
         4.0
    """

    def __init__(self, maximization=True):
        r"""
        Constructor for the ``MixedIntegerLinearProgram`` class.

        INPUT:

        - ``maximization``

          - When set to ``True`` (default), the ``MixedIntegerLinearProgram``
            is defined as a maximization.
          - When set to ``False``, the ``MixedIntegerLinearProgram`` is
            defined as a minimization.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True)
        """
        try:
            from sage.numerical.mipCoin import solveCoin
            self.default_solver = "Coin"
        except:
            try:
                from sage.numerical.mipGlpk import solveGlpk
                self.default_solver = "GLPK"
            except:
                self.default_solver = None

        from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing
        from sage.rings.real_double import RealDoubleField as RR
        P = InfinitePolynomialRing(RR(), names=('x',))
        (self.x,) = P._first_ngens(1)

        # number of variables
        self.count = 0
        self.maximization = maximization
        self.objective = None

        self.variables = {}
        self.constraints = []

        #constains the min and max bounds on the variables
        self.min = {}
        self.max = {}

        # constains the variables types
        self.types = {}

        #constains the variables' values when
        # solve(objective_only=False) is called
        self.values = {}

        # Several constants
        self.__BINARY = 1
        self.__REAL = -1
        self.__INTEGER = 0

    def __repr__(self):
         r"""
         Returns a short description of the ``MixedIntegerLinearProgram``.

         EXAMPLE::

             sage: p = MixedIntegerLinearProgram()
             sage: v = p.new_variable()
             sage: p.add_constraint(v[1] + v[2], max=2)
             sage: print p
             Mixed Integer Program ( maximization, 2 variables, 1 constraints )
         """
         return "Mixed Integer Program ( " + \
             ( "maximization" if self.maximization else "minimization" ) + \
             ", " + str(len(self.variables)) + " variables, " +  \
             str(len(self.constraints)) + " constraints )"

    def new_variable(self, vtype=-1, dim=1):
        r"""
        Returns an instance of ``MIPVariable`` associated
        to the current instance of ``MixedIntegerLinearProgram``.

        A new variable ``x`` is defined by::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()

        It behaves exactly as a usual dictionary would. It can use any key
        argument you may like, as ``x[5]`` or ``x["b"]``, and has methods
        ``items()`` and ``keys()``.

        Any of its fields exists, and is uniquely defined.

        INPUT:

        - ``dim`` (integer) -- Defines the dimension of the dictionary.
          If ``x`` has dimension `2`, its fields will be of the form
          ``x[key1][key2]``.
        - ``vtype`` (integer) -- Defines the type of the variables
          (default is ``REAL``).

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: # available types are p.__REAL, p.__INTEGER and p.__BINARY
            sage: x = p.new_variable(vtype=p.__REAL)
            sage: y = p.new_variable(dim=2)
            sage: p.add_constraint(x[2] + y[3][5], max=2)
        """
        return MIPVariable(self, vtype, dim=dim)

    def export(self, format="text"):
        r"""
        Exports the ``MixedIntegerLinearProgram`` to a string in
        different formats.

        INPUT:

        - ``format``

          - ``"text"`` -- (default) human-readable format

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2)
            sage: print p.export(format="text")
            Maximization:
              x2 + x1
            Constraints:
              2.0*x2 - 3.0*x1
            Variables:
              x2 is a real variable (min=0.0, max=+oo)
              x1 is a real variable (min=0.0, max=+oo)
        """
        if format == "text":
            value = ( "Maximization:\n"
                      if self.maximization
                      else "Minimization:\n" )
            value += "  " + ( str(self.objective)
                              if self.objective != None
                              else "Undefined" )
            value += "\nConstraints:"
            for c in self.constraints:
                value += "\n  " + str(c["function"])
            value += "\nVariables:"
            for v in self.variables.keys():
                value += "\n  " + str(v) + " is"
                if self.is_integer(v):
                    value += " an integer variable"
                elif self.is_binary(v):
                    value += " an boolean variable"
                else:
                    value += " a real variable"
                value += " (min=" + \
                    ( str(self.get_min(v))
                      if self.get_min(v) != None
                      else "-oo" ) + \
                    ", max=" + \
                    ( str(self.get_max(v))
                      if self.get_max(v) != None
                      else "+oo" ) + \
                    ")"
            return value
        else:
            raise ValueError("Only human-readable format is currently defined.")

    def get_values(self, *lists):
        r"""
        Return values found by the previous call to ``solve()``.

        INPUT:

        - Any instance of ``MIPVariable`` (or one of its elements),
          or lists of them.

        OUTPUT:

        - Each instance of ``MIPVariable`` is replaced by a dictionary
          containing the numerical values found for each
          corresponding variable in the instance.
        - Each element of an instance of a ``MIPVariable`` is replaced
          by its corresponding numerical value.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: y = p.new_variable(dim=2)
            sage: p.set_objective(x[3] + y[2][9] + x[5])
            sage: p.add_constraint(x[3] + y[2][9] + 2*x[5], max=2)
            sage: p.solve() # optional - requires Glpk or COIN-OR/CBC
            2.0
            sage: #
            sage: # Returns the optimal value of x[3]
            sage: p.get_values(x[3]) # optional - requires Glpk or COIN-OR/CBC
            2.0
            sage: #
            sage: # Returns a dictionary identical to x
            sage: # containing values for the corresponding
            sage: # variables
            sage: x_sol = p.get_values(x)
            sage: x_sol.keys()
            [3, 5]
            sage: #
            sage: # Obviously, it also works with
            sage: # variables of higher dimension
            sage: y_sol = p.get_values(y)
            sage: #
            sage: # We could also have tried :
            sage: [x_sol, y_sol] = p.get_values(x, y)
            sage: # Or
            sage: [x_sol, y_sol] = p.get_values([x, y])
        """
        val = []
        for l in lists:
            if isinstance(l, MIPVariable):
                if l.depth() == 1:
                    c = {}
                    for (k,v) in l.items():
                        c[k] = self.values[v] if self.values.has_key(v) else None
                    val.append(c)
                else:
                    c = {}
                    for (k,v) in l.items():
                        c[k] = self.get_values(v)
                    val.append(c)
            elif isinstance(l, list):
                if len(l) == 1:
                    val.append([self.get_values(l[0])])
                else:
                    c = []
                    [c.append(self.get_values(ll)) for ll in l]
                    val.append(c)
            elif self.variables.has_key(l):
                val.append(self.values[l])
        if len(lists) == 1:
            return val[0]
        else:
            return val

    def show(self):
        r"""
        Prints the ``MixedIntegerLinearProgram`` in a human-readable format.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2)
            sage: p.show()
            Maximization:
              x2 + x1
            Constraints:
              2.0*x2 - 3.0*x1
            Variables:
              x2 is a real variable (min=0.0, max=+oo)
              x1 is a real variable (min=0.0, max=+oo)
        """
        print self.export(format="text")

    def set_objective(self,obj):
        r"""
        Sets the objective of the ``MixedIntegerLinearProgram``.

        INPUT:

        - ``obj`` -- A linear function to be optimized.

        EXAMPLE:

        Let's solve the following linear program::

            Maximize:
              x + 5 * y
            Constraints:
              x + 0.2 y       <= 4
              1.5 * x + 3 * y <= 4
            Variables:
              x is Real (min = 0, max = None)
              y is Real (min = 0, max = None)

        This linear program can be solved as follows::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 0.2*x[2], max=4)
            sage: p.add_constraint(1.5*x[1]+3*x[2], max=4)
            sage: p.solve()     # optional - requires Glpk or COIN-OR/CBC
            6.6666666666666661
        """
        self.objective=obj

    def add_constraint(self, linear_function, max=None, min=None):
        r"""
        Adds a constraint to the ``MixedIntegerLinearProgram``.

        INPUT:

        - ``consraint`` -- A linear function.
        - ``max`` -- An upper bound on the constraint (set to ``None``
          by default).
        - ``min`` -- A lower bound on the constraint.

        EXAMPLE:

        Consider the following linear program::

            Maximize:
              x + 5 * y
            Constraints:
              x + 0.2 y       <= 4
              1.5 * x + 3 * y <= 4
            Variables:
              x is Real (min = 0, max = None)
              y is Real (min = 0, max = None)

        This linear program can be solved as follows::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 0.2*x[2], max=4)
            sage: p.add_constraint(1.5*x[1] + 3*x[2], max=4)
            sage: p.solve()     # optional - requires Glpk or COIN-OR/CBC
            6.6666666666666661
        """
        max = float(max) if max != None else None
        min = float(min) if min != None else None
        self.constraints.append({
                "function": linear_function,
                "min": min,
                "max": max,
                "card": len(linear_function.variables()) })

    def set_binary(self, e):
        r"""
        Sets a variable or a ``MIPVariable`` as binary.

        INPUT:

        - ``e`` -- An instance of ``MIPVariable`` or one of
          its elements.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: # With the following instruction, all the variables
            sage: # from x will be binary.
            sage: p.set_binary(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)
            sage: #
            sage: # It is still possible, though, to set one of these
            sage: # variables as real while keeping the others as they are.
            sage: p.set_real(x[3])
        """
        if isinstance(e, MIPVariable):
            e.vtype = self.__BINARY
            if e.depth() == 1:
                for v in e.values():
                    self.types[v] = self.__BINARY
            else:
                for v in e.keys():
                    self.set_binary(e[v])
        elif self.variables.has_key(e):
            self.types[e] = self.__BINARY
        else:
            raise ValueError("e must be an instance of MIPVariable or one of its elements.")

    def is_binary(self, e):
        r"""
        Tests whether the variable ``e`` is binary. Variables are real by
        default.

        INPUT:

        - ``e`` -- A variable (not a ``MIPVariable``, but one of its elements.)

        OUTPUT:

        ``True`` if the variable ``e`` is binary; ``False`` otherwise.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[1])
            sage: p.is_binary(v[1])
            False
            sage: p.set_binary(v[1])
            sage: p.is_binary(v[1])
            True
        """
        # Returns an exception if the variable does not exist.
        # For example if the user tries to find out the type of
        # a MIPVariable or anything else.
        self.variables[e]
        if self.types.has_key(e) and self.types[e] == self.__BINARY:
            return True
        return False

    def set_integer(self, e):
        r"""
        Sets a variable or a ``MIPVariable`` as integer.

        INPUT:

        - ``e`` -- An instance of ``MIPVariable`` or one of
          its elements.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: # With the following instruction, all the variables
            sage: # from x will be integers
            sage: p.set_integer(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)
            sage: #
            sage: # It is still possible, though, to set one of these
            sage: # variables as real while keeping the others as they are.
            sage: p.set_real(x[3])
        """
        if isinstance(e, MIPVariable):
            e.vtype = self.__INTEGER
            if e.depth() == 1:
                for v in e.values():
                    self.types[v] = self.__INTEGER
            else:
                for v in e.keys():
                    self.set_integer(e[v])
        elif self.variables.has_key(e):
            self.types[e] = self.__INTEGER
        else:
            raise ValueError("e must be an instance of MIPVariable or one of its elements.")

    def is_integer(self, e):
        r"""
        Tests whether the variable is an integer. Variables are real by
        default.

        INPUT:

        - ``e`` -- A variable (not a ``MIPVariable``, but one of its elements.)

        OUTPUT:

        ``True`` if the variable ``e`` is an integer; ``False`` otherwise.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[1])
            sage: p.is_integer(v[1])
            False
            sage: p.set_integer(v[1])
            sage: p.is_integer(v[1])
            True
        """
        # Returns an exception if the variable does not exist.
        # For example if the user tries to find out the type of
        # a MIPVariable or anything else.
        self.variables[e]
        if self.types.has_key(e) and self.types[e] == self.__INTEGER:
            return True
        return False

    def set_real(self,e):
        r"""
        Sets a variable or a ``MIPVariable`` as real.

        INPUT:

        - ``e`` -- An instance of ``MIPVariable`` or one of
          its elements.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: # With the following instruction, all the variables
            sage: # from x will be real (they are by default, though).
            sage: p.set_real(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)
            sage: #
            sage: # It is still possible, though, to set one of these
            sage: # variables as binary while keeping the others as they are.
            sage: p.set_binary(x[3])
        """
        if isinstance(e, MIPVariable):
            e.vtype = self.__REAL
            if e.depth() == 1:
                for v in e.values():
                    self.types[v] = self.__REAL
            else:
                for v in e.keys():
                    self.set_real(e[v])
        elif self.variables.has_key(e):
            self.types[e] = self.__REAL
        else:
            raise ValueError("e must be an instance of MIPVariable or one of its elements.")

    def is_real(self, e):
        r"""
        Tests whether the variable is real. Variables are real by default.

        INPUT:

        - ``e`` -- A variable (not a ``MIPVariable``, but one of its elements.)

        OUTPUT:

        ``True`` if the variable is real; ``False`` otherwise.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[1])
            sage: p.is_real(v[1])
            True
            sage: p.set_binary(v[1])
            sage: p.is_real(v[1])
            False
            sage: p.set_real(v[1])
            sage: p.is_real(v[1])
            True
        """
        # Returns an exception if the variable does not exist.
        # For example if the user tries to find out the type of
        # a MIPVariable or anything else.
        self.variables[e]
        if (not self.types.has_key(e)) or self.types[e] == self.__REAL:
            return True
        return False

    def solve(self, solver=None, log=False, objective_only=False):
        r"""
        Solves the ``MixedIntegerLinearProgram``.

        INPUT:

        - ``solver`` -- 3 solvers should be available through this class:

          - GLPK (``solver="GLPK"``). See the
            `GLPK <http://www.gnu.org/software/glpk/>`_ web site.

          - COIN Branch and Cut  (``solver="Coin"``). See the
            `COIN-OR <http://www.coin-or.org>`_ web site.

          - CPLEX (``solver="CPLEX"``). See the
            `CPLEX <http://www.ilog.com/products/cplex/>`_ web site.
            An interface to CPLEX is not yet implemented.

          ``solver`` should then be equal to one of ``"GLPK"``, ``"Coin"``,
          ``"CPLEX"``, or ``None``. If ``solver=None`` (default), the default
          solver is used (COIN if available, GLPK otherwise).

        - ``log`` -- This boolean variable indicates whether progress should
          be printed during the computations.

        - ``objective_only`` -- Boolean variable.

          - When set to ``True``, only the objective function is returned.
          - When set to ``False`` (default), the optimal numerical values
            are stored (takes computational time).

        OUTPUT:

        The optimal value taken by the objective function.

        EXAMPLES:

        Consider the following linear program::

            Maximize:
              x + 5 * y
            Constraints:
              x + 0.2 y       <= 4
              1.5 * x + 3 * y <= 4
            Variables:
              x is Real (min = 0, max = None)
              y is Real (min = 0, max = None)

        This linear program can be solved as follows::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 0.2*x[2], max=4)
            sage: p.add_constraint(1.5*x[1] + 3*x[2], max=4)
            sage: p.solve()           # optional - requires Glpk or COIN-OR/CBC
            6.6666666666666661
            sage: p.get_values(x)     # optional - requires Glpk or COIN-OR/CBC
            {1: 0.0, 2: 1.3333333333333333}

            sage: ### Computation of a maximum stable set in Petersen's graph
            sage: g = graphs.PetersenGraph()
            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: b = p.new_variable()
            sage: p.set_objective(sum([b[v] for v in g]))
            sage: for (u,v) in g.edges(labels=None):
            ...       p.add_constraint(b[u] + b[v], max=1)
            sage: p.set_binary(b)
            sage: p.solve(objective_only=True)     # optional - requires Glpk or COIN-OR/CBC
            4.0
        """
        if self.objective == None:
            raise ValueError("No objective function has been defined.")

        if solver == None:
            solver = self.default_solver

        if solver == None:
            raise ValueError("There does not seem to be any solver installed. Please visit http://www.sagemath.org/doc/tutorial/tour_LP.html for more informations.")
        elif solver == "Coin":
            try:
                from sage.numerical.mipCoin import solveCoin
            except:
                raise NotImplementedError("Coin/CBC is not installed and cannot be used to solve this MixedIntegerLinearProgram. To install it, you can type in Sage: install_package('cbc')")
            return solveCoin(self, log=log, objective_only=objective_only)
        elif solver == "GLPK":
            try:
                from sage.numerical.mipGlpk import solveGlpk
            except:
                raise NotImplementedError("GLPK is not installed and cannot be used to solve this MixedIntegerLinearProgram. To install it, you can type in Sage: install_package('glpk')")
            return solveGlpk(self, log=log, objective_only=objective_only)
        elif solver == "CPLEX":
            raise NotImplementedError("The support for CPLEX is not implemented yet.")
        else:
            raise NotImplementedError("'solver' should be set to 'GLPK', 'Coin', 'CPLEX' or None (in which case the default one is used).")

    def _NormalForm(self, exp):
        r"""
        Returns a dictionary built from the linear function.

        INPUT:

        - ``exp`` -- The expression representing a linear function.

        OUTPUT:

        A dictionary whose keys are the IDs of the variables and whose
        values are their coefficients. The value corresponding to key
        `-1` is the constant coefficient.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p._NormalForm(v[0] + v[1])
            {1: 1.0, 2: 1.0, -1: 0.0}
        """
        d = dict( zip([self.variables[v] for v in exp.variables()],
                      exp.coefficients()) )
        d[-1] = exp.constant_coefficient()
        return d

    def _add_element_to_ring(self, vtype):
        r"""
        Creates a new variable from the main ``InfinitePolynomialRing``.

        INPUT:

        - ``vtype`` (integer) -- Defines the type of the variables
          (default is ``REAL``).

        OUTPUT:

        - The newly created variable.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.count
            0
            sage: p._add_element_to_ring(p.__REAL)
            x1
            sage: p.count
            1
        """
        self.count += 1
        v = self.x[self.count]
        self.variables[v] = self.count
        self.types[v] = vtype
        self.min[v] = 0.0
        return v

    def set_min(self, v, min):
        r"""
        Sets the minimum value of a variable.

        INPUT:

        - ``v`` -- a variable (not a ``MIPVariable``, but one of its
          elements).
        - ``min`` -- the minimum value the variable can take.
          When ``min=None``, the variable has no lower bound.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[1])
            sage: p.get_min(v[1])
            0.0
            sage: p.set_min(v[1],6)
            sage: p.get_min(v[1])
            6.0
        """
        self.min[v] = min

    def set_max(self, v, max):
        r"""
        Sets the maximum value of a variable.

        INPUT

        - ``v`` -- a variable (not a ``MIPVariable``, but one of its
          elements).
        - ``max`` -- the maximum value the variable can take.
          When ``max=None``, the variable has no upper bound.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[1])
            sage: p.get_max(v[1])
            sage: p.set_max(v[1],6)
            sage: p.get_max(v[1])
            6.0
        """
        self.max[v] = max

    def get_min(self, v):
        r"""
        Returns the minimum value of a variable.

        INPUT:

        - ``v`` -- a variable (not a ``MIPVariable``, but one of its elements).

        OUTPUT:

        Minimum value of the variable, or ``None`` if
        the variable has no lower bound.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[1])
            sage: p.get_min(v[1])
            0.0
            sage: p.set_min(v[1],6)
            sage: p.get_min(v[1])
            6.0
        """
        return float(self.min[v]) if self.min.has_key(v) else None

    def get_max(self, v):
        r"""
        Returns the maximum value of a variable.

        INPUT:

        - ``v`` -- a variable (not a ``MIPVariable``, but one of its elements).

        OUTPUT:

        Maximum value of the variable, or ``None`` if
        the variable has no upper bound.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[1])
            sage: p.get_max(v[1])
            sage: p.set_max(v[1],6)
            sage: p.get_max(v[1])
            6.0
        """
        return float(self.max[v]) if self.max.has_key(v) else None

class MIPSolverException(Exception):
    r"""
    Exception raised when the solver fails.
    """

    def __init__(self, value):
        r"""
        Constructor for ``MIPSolverException``.

        ``MIPSolverException`` is the exception raised when the solver fails.

        EXAMPLE::

            sage: MIPSolverException("Error")
            MIPSolverException()
        """
        self.value = value

    def __str__(self):
        r"""
        Returns the value of the instance of ``MIPSolverException``.

        EXAMPLE::

            sage: e = MIPSolverException("Error")
            sage: print e
            'Error'
        """
        return repr(self.value)

class MIPVariable:
    r"""
    ``MIPVariable`` is a variable used by the class
    ``MixedIntegerLinearProgram``.
    """

    def __init__(self, p, vtype, dim=1):
        r"""
        Constructor for ``MIPVariable``.

        INPUT:

        - ``p`` -- the instance of ``MixedIntegerLinearProgram`` to which the
          variable is to be linked.
        - ``vtype`` (integer) -- Defines the type of the variables
          (default is ``REAL``).
        - ``dim`` -- the integer defining the definition of the variable.

        For more informations, see the method
        ``MixedIntegerLinearProgram.new_variable``.

        EXAMPLE::

            sage: p=MixedIntegerLinearProgram()
            sage: v=p.new_variable()
        """
        self.dim = dim
        self.dict = {}
        self.p = p
        self.vtype = vtype

    def __getitem__(self, i):
        r"""
        Returns the symbolic variable corresponding to the key.

        Returns the element asked, otherwise creates it.
        (When depth>1, recursively creates the variables).

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: v[0]
            x1
        """
        if self.dict.has_key(i):
            return self.dict[i]
        elif self.dim == 1:
            self.dict[i] = self.p._add_element_to_ring(self.vtype)
            return self.dict[i]
        else:
            self.dict[i] = MIPVariable(self.p, self.vtype, dim=self.dim-1)
            return self.dict[i]

    def keys(self):
        r"""
        Returns the keys already defined in the dictionary.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: v.keys()
            [0, 1]
        """
        return self.dict.keys()

    def items(self):
        r"""
        Returns the pairs (keys,value) contained in the dictionary.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: v.items()
            [(0, x1), (1, x2)]
        """
        return self.dict.items()

    def depth(self):
        r"""
        Returns the current variable's depth.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: v.depth()
            1
        """
        return self.dim

    def values(self):
        r"""
        Returns the symbolic variables associated to the current dictionary.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: v.values()
            [x1, x2]
        """
        return self.dict.values()
