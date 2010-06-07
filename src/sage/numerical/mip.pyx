r"""
Mixed integer linear programming
"""

include "../ext/stdsage.pxi"
include "../ext/interrupt.pxi"
from copy import deepcopy

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

    EXAMPLES:

    Computation of a maximum stable set in Petersen's graph::

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

        self._default_solver = None

        try:
            if self._default_solver == None:
                from sage.numerical.mip_cplex import solve_cplex
                self._default_solver = "CPLEX"
        except ImportError:
            pass

        try:
            if self._default_solver == None:
                from sage.numerical.mip_coin import solve_coin
                self._default_solver = "Coin"
        except ImportError:
            pass

        try:
            if self._default_solver == None:
                from sage.numerical.mip_glpk import solve_glpk
                self._default_solver = "GLPK"
        except ImportError:
            pass


        # List of all the MIPVariables linked to this instance of
        # MixedIntegerLinearProgram
        self._mipvariables = []

        # Associates an index to the variables
        self._variables = {}

        # contains the variables' values when
        # solve(objective_only=False) is called
        self._values = {}

        # Several constants
        self.__BINARY = 1
        self.__REAL = -1
        self.__INTEGER = 0

        # ######################################################
        # The information data of a Linear Program
        #
        # - name
        # - maximization
        # - objective
        #    --> name
        #    --> i
        #    --> values
        # - variables
        #   --> names
        #   --> type
        #   --> bounds
        #       --> min
        #       --> max
        # - constraints
        #   --> names
        #   --> matrix
        #       --> i
        #       --> j
        #       --> values
        #   --> bounds
        #       --> min
        #       --> max
        #
        # The Constraint matrix being almost always sparse, it is stored
        # as a list of positions (i,j) in the matrix with an associated value.
        #
        # This is how matrices are exchanged in GLPK's or Cbc's libraries
        # By storing the data this way, we do no have to convert them
        # ( too often ) and this process is a bit faster.
        #
        # ######################################################


        self._name = None
        self._maximization = maximization
        self._objective_i = None
        self._objective_values = None
        self._objective_name = None
        self._variables_name = []
        self._variables_type = []
        self._variables_bounds_min = []
        self._variables_bounds_max = []
        self._constraints_name = []
        self._constraints_matrix_i = []
        self._constraints_matrix_j = []
        self._constraints_matrix_values = []
        self._constraints_bounds_max = []
        self._constraints_bounds_min = []



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
         return "Mixed Integer Program "+("\""+self._name+"\"" if self._name!=None else "")+" ( " + \
             ( "maximization" if self._maximization else "minimization" ) + \
             ", " + str(len(self._variables)) + " variables, " +  \
             str(len(self._constraints_bounds_min)) + " constraints )"

    def __eq__(self,p):
        r"""
        Test of equality.

        INPUT:

        - ``p`` -- an instance of ``MixedIntegerLinearProgram`` to be tested
          against ``self``.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.add_constraint(v[1] + v[2], max=2)
            sage: p == loads(dumps(p))
            True
            sage: p2 = loads(dumps(p))
            sage: p2.add_constraint(2*v[1] + 3*v[2], max=1)
            sage: p == p2
            False
        """

        return (
            self._name == p._name and
            self._maximization == p._maximization and
            self._objective_i == p._objective_i and
            self._objective_values == p._objective_values and
            self._objective_name == p._objective_name and
            self._variables_name == p._variables_name and
            self._variables_type == p._variables_type and
            self._variables_bounds_min == p._variables_bounds_min and
            self._variables_bounds_max == p._variables_bounds_max and
            self._constraints_name == p._constraints_name and
            self._constraints_matrix_i == p._constraints_matrix_i and
            self._constraints_matrix_j == p._constraints_matrix_j and
            self._constraints_matrix_values == p._constraints_matrix_values and
            self._constraints_bounds_max == p._constraints_bounds_max
            )

    def __getitem__(self, v):
        r"""
        Returns the symbolic variable corresponding to the key
        from a default dictionary.

        It returns the element asked, and otherwise creates it.
        If necessary, it also creates the default dictionary.

        This method lets the user define LinearProgram without having to
        define independent dictionaries when it is not necessary for him.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: p.set_objective(p['x'] + p['z'])
            sage: p['x']
            x_0
        """

        try:
            return self._default_mipvariable[v]
        except AttributeError:
            self._default_mipvariable = self.new_variable()
            return self._default_mipvariable[v]

    def set_problem_name(self,name):
        r"""
        Sets the name of the ``MixedIntegerLinearProgram``.

        INPUT:

        - ``name`` -- A string representing the name of the
          ``MixedIntegerLinearProgram``.

        EXAMPLE::

            sage: p=MixedIntegerLinearProgram()
            sage: p.set_problem_name("Test program")
            sage: p
            Mixed Integer Program "Test program" ( maximization, 0 variables, 0 constraints )
        """

        self._name=name

    def set_objective_name(self,name):
        r"""
        Sets the name of the objective function.

        INPUT:

        - ``name`` -- A string representing the name of the
          objective function.

        EXAMPLE::

            sage: p=MixedIntegerLinearProgram()
            sage: p.set_objective_name("Objective function")
        """

        self._objective_name=name

    def _update_variables_name(self):
        r"""
        Updates the names of the variables.

        Only called before writing the Problem to a MPS or LP file.

        EXAMPLE::

            sage: p=MixedIntegerLinearProgram()
            sage: v=p.new_variable(name="Test")
            sage: v[5]+v[99]
            x_0 +x_1
            sage: p._update_variables_name()
        """

        self._variables_name=['']*len(self._variables)
        for v in self._mipvariables:
            v._update_variables_name()


    def new_variable(self, real=False, binary=False, integer=False, dim=1,name=None):
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

        - ``binary, integer, real`` (boolean) -- Set one of these arguments
          to ``True`` to ensure that the variable gets the corresponding
          type. The default type is ``real``.

        - ``name`` (string) -- Associates a name to the variable. This is
          only useful when exporting the linear program to a file using
          ``write_mps`` or ``write_lp``, and has no other effect.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()

         To define two dictionaries of variables, the first being
         of real type, and the second of integer type ::

            sage: x = p.new_variable(real=True)
            sage: y = p.new_variable(dim=2, integer=True)
            sage: p.add_constraint(x[2] + y[3][5], max=2)
            sage: p.is_integer(x[2])
            False
            sage: p.is_integer(y[3][5])
            True

        An exception is raised when two types are supplied ::

            sage: z = p.new_variable(real = True, integer = True)
            Traceback (most recent call last):
            ...
            ValueError: Exactly one of the available types has to be True
        """
        if name==None:
            name="V"+str(len(self._mipvariables))

        if sum([real, binary, integer]) >= 2:
            raise ValueError("Exactly one of the available types has to be True")

        if binary:
            vtype = self.__BINARY
        elif integer:
            vtype = self.__INTEGER
        else:
            vtype = self.__REAL

        v=MIPVariable(self, vtype, dim=dim,name=name)
        self._mipvariables.append(v)
        return v

    def constraints(self):
        r"""
        Returns the list of constraints.

        This function returns the constraints as a list
        of tuples ``(linear_function,min_bound,max_bound)``, representing
        the constraint:

        .. MATH::

            \text{min\_bound}
            \leq
            \text{linear\_function}
            \leq
            \text{max\_bound}

        Variables ``min_bound`` (respectively ``max_bound``) is set
        to ``None`` when the function has no lower ( respectively upper )
        bound.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 2/10*x[2], max=4)
            sage: p.add_constraint(15/10*x[1]+3*x[2], max=4)
            sage: p.constraints()
            [(x_0 +1/5 x_1, None, 4), (3/2 x_0 +3 x_1, None, 4)]
        """

        d = [0]*len(self._variables)
        for (v,id) in self._variables.iteritems():
            d[id]=v

        constraints=[0]*len(self._constraints_bounds_min)
        for (i,j,value) in zip(self._constraints_matrix_i,self._constraints_matrix_j,self._constraints_matrix_values):
            constraints[i]+=value*d[j]
        return zip(constraints,self._constraints_bounds_min,self._constraints_bounds_max)


    def show(self):
        r"""
        Displays the ``MixedIntegerLinearProgram`` in a human-readable
        way.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2, name="Constraint_1")
            sage: p.show()
            Maximization:
              x_0 +x_1
            Constraints:
              Constraint_1: -3 x_0 +2 x_1 <= 2
            Variables:
              x_0 is a real variable (min=0.0, max=+oo)
              x_1 is a real variable (min=0.0, max=+oo)
        """

        self._update_variables_name()

        inv_variables = [0]*len(self._variables)
        for (v,id) in self._variables.iteritems():
            inv_variables[id]=v


        value = ( "Maximization:\n"
                  if self._maximization
                  else "Minimization:\n" )
        value+="  "
        if self._objective_i==None:
            value+="Undefined"
        else:
            value+=str(sum([inv_variables[i]*c for (i,c) in zip(self._objective_i, self._objective_values)]))

        value += "\nConstraints:"
        for (c,min,max), name in zip(self.constraints(), self._constraints_name):
            value += "\n  "+(name+":" if name is not None else "")+" " + (str(min)+" <= " if min!=None else "")+str(c)+(" <= "+str(max) if max!=None else "")
        value += "\nVariables:"
        for _,v in sorted([(str(x),x) for x in self._variables.keys()]):

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
        print value

    def write_mps(self,filename,modern=True):
        r"""
        Write the linear program as a MPS file.

        This function export the problem as a MPS file.

        INPUT:

        - ``filename`` -- The file in which you want the problem
          to be written.

        - ``modern`` -- Lets you choose between Fixed MPS and Free MPS

            - ``True`` -- Outputs the problem in Free MPS
            - ``False`` -- Outputs the problem in Fixed MPS

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2,name="OneConstraint")
            sage: p.write_mps(SAGE_TMP+"/lp_problem.mps") # optional - requires GLPK

        For information about the MPS file format :
        http://en.wikipedia.org/wiki/MPS_%28format%29
        """

        try:
            from sage.numerical.mip_glpk import write_mps
        except:
            raise NotImplementedError("You need GLPK installed to use this function. To install it, you can type in Sage: install_package('glpk')")

        self._update_variables_name()
        write_mps(self, filename, modern)


    def write_lp(self,filename):
        r"""
        Write the linear program as a LP file.

        This function export the problem as a LP file.

        INPUT:

        - ``filename`` -- The file in which you want the problem
          to be written.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2)
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp") # optional - requires GLPK

        For more information about the LP file format :
        http://lpsolve.sourceforge.net/5.5/lp-format.htm
        """
        try:
            from sage.numerical.mip_glpk import write_lp
        except:
            raise NotImplementedError("You need GLPK installed to use this function. To install it, you can type in Sage: install_package('glpk')")

        self._update_variables_name()
        write_lp(self, filename)


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
            sage: p.set_objective(x[3] + 3*y[2][9] + x[5])
            sage: p.add_constraint(x[3] + y[2][9] + 2*x[5], max=2)
            sage: p.solve() # optional - requires Glpk or COIN-OR/CBC
            6.0

        To return  the optimal value of ``y[2][9]``::

            sage: p.get_values(y[2][9]) # optional - requires Glpk or COIN-OR/CBC
            2.0

        To get a dictionary identical to ``x`` containing optimal
        values for the corresponding variables ::

            sage: x_sol = p.get_values(x)
            sage: x_sol.keys()
            [3, 5]

        Obviously, it also works with variables of higher dimension::

            sage: y_sol = p.get_values(y)

        We could also have tried ::

            sage: [x_sol, y_sol] = p.get_values(x, y)

        Or::

            sage: [x_sol, y_sol] = p.get_values([x, y])
        """
        val = []
        for l in lists:
            if isinstance(l, MIPVariable):
                if l.depth() == 1:
                    c = {}
                    for (k,v) in l.items():
                        c[k] = self._values[v] if self._values.has_key(v) else None
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
            elif self._variables.has_key(l):
                val.append(self._values[l])
        if len(lists) == 1:
            return val[0]
        else:
            return val

    def set_objective(self,obj):
        r"""
        Sets the objective of the ``MixedIntegerLinearProgram``.

        INPUT:

        - ``obj`` -- A linear function to be optimized.
          ( can also be set to ``None`` or ``0`` when just
          looking for a feasible solution )

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
            sage: p.add_constraint(x[1] + 2/10*x[2], max=4)
            sage: p.add_constraint(1.5*x[1]+3*x[2], max=4)
            sage: round(p.solve(),5)     # optional - requires Glpk or COIN-OR/CBC
            6.66667
            sage: p.set_objective(None)
            sage: p.solve() #optional - requires Glpk or COIN-OR/CBC
            0.0
        """


        self._objective_i = []
        self._objective_values = []

        # If the objective is None, or a constant, we want to remember
        # that the objective function has been defined ( the user did not
        # forget it ). In some LP problems, you just want a feasible solution
        # and do not care about any function being optimal.

        if obj != None:
            f = obj.dict()
        else:
            return None

        f.pop(-1,0)

        for (v,coeff) in f.iteritems():
            self._objective_i.append(v)
            self._objective_values.append(coeff)

    def add_constraint(self, linear_function, max=None, min=None, name=None):
        r"""
        Adds a constraint to self.

        INPUT:

        - ``linear_function`` -- Two different types of arguments are possible:
            - A linear function. In this case, arguments ``min`` or ``max``
              have to be specified.
            - A linear constraint of the form ``A <= B``, ``A >= B``,
              ``A <= B <= C``, ``A >= B >= C`` or ``A == B``. In this
              case, arguments ``min`` and ``max`` will be ignored.
        - ``max`` -- An upper bound on the constraint (set to ``None``
          by default). This must be a numerical value.
        - ``min`` -- A lower bound on the constraint.  This must be a
          numerical value
        - ``name`` -- A name for the constraint.

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

        It can be solved as follows::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 0.2*x[2], max=4)
            sage: p.add_constraint(1.5*x[1] + 3*x[2], max=4)
            sage: round(p.solve(),6)     # optional - requires Glpk or COIN-OR/CBC
            6.666667

        There are two different ways to add the constraint
        ``x[5] + 3*x[7] <= x[6] + 3`` to self.

        The first one consists in giving ``add_constraint`` this
        very expression::

            sage: p.add_constraint( x[5] + 3*x[7] <= x[6] + 3 )

        The second (slightly more efficient) one is to use the
        arguments ``min`` or ``max``, which can only be numerical
        values::

            sage: p.add_constraint( x[5] + 3*x[7] - x[6], max = 3 )

        One can also define double-bounds or equality using symbols
        ``<=``, ``>=`` and ``==``::

            sage: p.add_constraint( x[5] + 3*x[7] == x[6] + 3 )
            sage: p.add_constraint( x[5] + 3*x[7] <= x[6] + 3 <= x[8] + 27 )

        The previous program can be rewritten::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 0.2*x[2] <= 4)
            sage: p.add_constraint(1.5*x[1] + 3*x[2] <= 4)
            sage: p.solve()     # optional - requires Glpk or COIN-OR/CBC
            6.6666666666666661


        TESTS::

            sage: p=MixedIntegerLinearProgram()
            sage: p.add_constraint(sum([]),min=2)
        """
        if linear_function is None or linear_function is 0:
            return None

        if isinstance(linear_function, LinearFunction):

            f = linear_function.dict()

            self._constraints_name.append(name)

            constant_coefficient = f.get(-1,0)

            # We do not want to ignore the constant coefficient
            max = (max - constant_coefficient) if max != None else None
            min = (min - constant_coefficient) if min != None else None

            c = len(self._constraints_bounds_min)

            for (v,coeff) in f.iteritems():
                if v != -1:
                    self._constraints_matrix_i.append(c)
                    self._constraints_matrix_j.append(v)
                    self._constraints_matrix_values.append(coeff)

            self._constraints_bounds_max.append(max)
            self._constraints_bounds_min.append(min)

        elif isinstance(linear_function,LinearConstraint):
            functions = linear_function.constraints

            if linear_function.equality:
                self.add_constraint(functions[0] - functions[1], min=0, max=0, name=name)

            elif len(functions) == 2:
                self.add_constraint(functions[0] - functions[1], max=0, name=name)

            else:
                self.add_constraint(functions[0] - functions[1], max=0, name=name)
                self.add_constraint(functions[1] - functions[2], max=0, name=name)

    def set_binary(self, e):
        r"""
        Sets a variable or a ``MIPVariable`` as binary.

        INPUT:

        - ``e`` -- An instance of ``MIPVariable`` or one of
          its elements.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()

        With the following instruction, all the variables
        from x will be binary::

            sage: p.set_binary(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)

        It is still possible, though, to set one of these
        variables as real while keeping the others as they are::

            sage: p.set_real(x[3])
        """
        if isinstance(e, MIPVariable):
            e._vtype = self.__BINARY
            if e.depth() == 1:
                for v in e.values():
                    self._variables_type[self._variables[v]] = self.__BINARY
            else:
                for v in e.keys():
                    self.set_binary(e[v])
        elif self._variables.has_key(e):
            self._variables_type[self._variables[e]] = self.__BINARY
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

        if self._variables_type[self._variables[e]] == self.__BINARY:
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

        With the following instruction, all the variables
        from x will be integers::

            sage: p.set_integer(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)

        It is still possible, though, to set one of these
        variables as real while keeping the others as they are::

            sage: p.set_real(x[3])
        """
        if isinstance(e, MIPVariable):
            e._vtype = self.__INTEGER
            if e.depth() == 1:
                for v in e.values():
                    self._variables_type[self._variables[v]] = self.__INTEGER
            else:
                for v in e.keys():
                    self.set_integer(e[v])
        elif self._variables.has_key(e):
            self._variables_type[self._variables[e]] = self.__INTEGER
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

        if self._variables_type[self._variables[e]] == self.__INTEGER:
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

        With the following instruction, all the variables
        from x will be real (they are by default, though)::

            sage: p.set_real(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)

         It is still possible, though, to set one of these
         variables as binary while keeping the others as they are::

            sage: p.set_binary(x[3])
        """
        if isinstance(e, MIPVariable):
            e._vtype = self.__REAL
            if e.depth() == 1:
                for v in e.values():
                    self._variables_type[self._variables[v]] = self.__REAL
            else:
                for v in e.keys():
                    self.set_real(e[v])
        elif self._variables.has_key(e):
            self._variables_type[self._variables[e]] = self.__REAL
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

        if self._variables_type[self._variables[e]] == self.__REAL:
            return True
        return False

    def solve(self, solver=None, log=0, objective_only=False, threads = 0):
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

        - ``log`` -- integer (default: ``0``) The verbosity level. Indicates
          whether progress should be printed during computation.

        - ``threads`` -- Number of threads to use. This option is only useful
          when Coin is used to solve the problem. Set ``threads`` to 0 (its
          default value) to use one thread per core available.

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
            sage: round(p.solve(),6)  # optional - requires Glpk or COIN-OR/CBC
            6.666667
            sage: p.get_values(x)     # optional random - requires Glpk or COIN-OR/CBC
            {0: 0.0, 1: 1.3333333333333333}

         Computation of a maximum stable set in Petersen's graph::

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


        if self._objective_i == None:
            raise ValueError("No objective function has been defined.")

        if solver == None:
            solver = self._default_solver


        if solver == None:
            raise ValueError("There does not seem to be any (Mixed) Integer Linear Program solver installed. Please visit http://www.sagemath.org/doc/constructions/linear_programming.html to learn more about the solvers available.")


        try:
            if solver == "Coin":
                from sage.numerical.mip_coin import solve_coin as solve
            elif solver == "GLPK":
                from sage.numerical.mip_glpk import solve_glpk as solve
            elif solver == "CPLEX":
                from sage.numerical.mip_cplex import solve_cplex as solve
            else:
                NotImplementedError("'solver' should be set to 'GLPK', 'Coin', 'CPLEX' or None (in which case the default one is used).")

        except ImportError:
            raise NotImplementedError("The required solver is not installed and cannot be used to solve this (Mixed) Integer Linear Program. To install it, follow the instructions given at http://www.sagemath.org/doc/constructions/linear_programming.html")


        if solver=="Coin":
            return solve(self, log=log, objective_only=objective_only, threads=threads)
        else:
            return solve(self, log=log, objective_only=objective_only)


    def _add_element_to_ring(self, vtype):
        r"""
        Creates a new variable in the (Mixed) Integer Linear Program.

        INPUT:

        - ``vtype`` (integer) -- Defines the type of the variables
          (default is ``REAL``).

        OUTPUT:

        - The newly created variable.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: len(p._variables_type)
            0
            sage: p._add_element_to_ring(p.__REAL)
            x_0
            sage: len(p._variables_type)
            1
        """

        v = LinearFunction({len(self._variables) : 1})

        self._variables[v] = len(self._variables)
        self._variables_type.append(vtype)
        self._variables_bounds_min.append(0)
        self._variables_bounds_max.append(None)
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
        self._variables_bounds_min[self._variables[v]] = min

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

        self._variables_bounds_max[self._variables[v]] = max

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
        return float(self._variables_bounds_min[self._variables[v]]) if self._variables_bounds_min[self._variables[v]] != None else None

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
        return float(self._variables_bounds_max[self._variables[v]])  if self._variables_bounds_max[self._variables[v]] != None else None

class MIPSolverException(Exception):
    r"""
    Exception raised when the solver fails.
    """

    def __init__(self, value):
        r"""
        Constructor for ``MIPSolverException``.

        ``MIPSolverException`` is the exception raised when the solver fails.

        EXAMPLE::

            sage: from sage.numerical.mip import MIPSolverException
            sage: MIPSolverException("Error")
            MIPSolverException()

        TESTS:

         No continuous solution::

            sage: p=MixedIntegerLinearProgram()
            sage: v=p.new_variable()
            sage: p.add_constraint(v[0],max=5.5)
            sage: p.add_constraint(v[0],min=7.6)
            sage: p.set_objective(v[0])

         Tests of GLPK's Exceptions::

            sage: p.solve(solver="GLPK") # optional - requires GLPK
            Traceback (most recent call last):
            ...
            MIPSolverException: 'GLPK : Solution is undefined'

         No integer solution::

            sage: p=MixedIntegerLinearProgram()
            sage: v=p.new_variable()
            sage: p.add_constraint(v[0],max=5.6)
            sage: p.add_constraint(v[0],min=5.2)
            sage: p.set_objective(v[0])
            sage: p.set_integer(v)

         Tests of GLPK's Exceptions::

            sage: p.solve(solver="GLPK") # optional - requires GLPK
            Traceback (most recent call last):
            ...
            MIPSolverException: 'GLPK : Solution is undefined'


        """
        self.value = value

    def __str__(self):
        r"""
        Returns the value of the instance of ``MIPSolverException``.

        EXAMPLE::

            sage: from sage.numerical.mip import MIPSolverException
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

    def __init__(self, p, vtype, dim=1, name=None):
        r"""
        Constructor for ``MIPVariable``.

        INPUT:

        - ``p`` -- the instance of ``MixedIntegerLinearProgram`` to which the
          variable is to be linked.
        - ``vtype`` (integer) -- Defines the type of the variables
          (default is ``REAL``).
        - ``dim`` -- the integer defining the definition of the variable.
        - ``name`` -- A name for the ``MIPVariable``.

        For more informations, see the method
        ``MixedIntegerLinearProgram.new_variable``.

        EXAMPLE::

            sage: p=MixedIntegerLinearProgram()
            sage: v=p.new_variable()
        """
        self._dim = dim
        self._dict = {}
        self._p = p
        self._vtype = vtype
        self._name=name


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
            x_0
        """
        if self._dict.has_key(i):
            return self._dict[i]
        elif self._dim == 1:
            self._dict[i] = self._p._add_element_to_ring(self._vtype)
            return self._dict[i]
        else:
            self._dict[i] = MIPVariable(self._p, self._vtype, dim=self._dim-1)
            return self._dict[i]

    def _update_variables_name(self, prefix=None):
        r"""
        Updates the names of the variables in the parent instant of ``MixedIntegerLinearProgram``.

        Only called before writing the Problem to a MPS or LP file.

        EXAMPLE::

            sage: p=MixedIntegerLinearProgram()
            sage: v=p.new_variable(name="Test")
            sage: v[5]+v[99]
            x_0 +x_1
            sage: p._variables_name=['']*2
            sage: v._update_variables_name()
        """

        if prefix == None:
            prefix = self._name

        if self._dim == 1:
            for (k,v) in self._dict.iteritems():
                self._p._variables_name[self._p._variables[v]]=prefix + "[" + str(k) + "]"
        else:
            for v in self._dict.itervalues():
                v._update_variables_name(prefix=prefix + "[" + str(k) + "]")



    def __repr__(self):
        r"""
        Returns a representation of self.

        EXAMPLE::

            sage: p=MixedIntegerLinearProgram()
            sage: v=p.new_variable(dim=3)
            sage: v
            MIPVariable of dimension 3.
            sage: v[2][5][9]
            x_0
            sage: v
            MIPVariable of dimension 3.
        """
        return "MIPVariable of dimension " + str(self._dim) + "."

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
        return self._dict.keys()

    def items(self):
        r"""
        Returns the pairs (keys,value) contained in the dictionary.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: v.items()
            [(0, x_0), (1, x_1)]
        """
        return self._dict.items()

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
        return self._dim

    def values(self):
        r"""
        Returns the symbolic variables associated to the current dictionary.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[0] + v[1])
            sage: v.values()
            [x_0, x_1]
        """
        return self._dict.values()


class LinearFunction:
    r"""
    An elementary algebra to represent symbolic linear functions.
    """

    def __init__(self,f):
        r"""
        Constructor taking a dictionary or a numerical value as its argument.

        A linear function is represented as a dictionary. The
        values are the coefficient of the variable represented
        by the keys ( which are integers ). The key ``-1``
        corresponds to the constant term.

        EXAMPLES:

        With a dictionary::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({0 : 1, 3 : -8})
            x_0 -8 x_3

        Using the constructor with a numerical value::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction(25)
            25
        """
        if isinstance(f, dict):
            self._f = f
        else:
            self._f = {-1:f}

    def dict(self):
        r"""
        Returns the dictionary corresponding to the Linear Function.

        A linear function is represented as a dictionary. The
        value are the coefficient of the variable represented
        by the keys ( which are integers ). The key ``-1``
        corresponds to the constant term.

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: lf = LinearFunction({0 : 1, 3 : -8})
            sage: lf.dict()
            {0: 1, 3: -8}
        """
        return self._f

    def __add__(self,b):
        r"""
        Defining the + operator

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({0 : 1, 3 : -8}) + LinearFunction({2 : 5, 3 : 2}) - 16
            -16 +x_0 +5 x_2 -6 x_3
        """
        if isinstance(b, LinearFunction):
            e = deepcopy(self._f)
            for (id,coeff) in b.dict().iteritems():
                e[id] = self._f.get(id,0) + coeff
            return LinearFunction(e)
        else:
            el = deepcopy(self)
            el.dict()[-1] = el.dict().get(-1,0) + b
            return el

    def __neg__(self):
        r"""
        Defining the - operator (opposite).

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: -LinearFunction({0 : 1, 3 : -8})
            -1 x_0 +8 x_3
        """
        return LinearFunction(dict([(id,-coeff) for (id, coeff) in self._f.iteritems()]))

    def __sub__(self,b):
        r"""
        Defining the - operator (substraction).

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({2 : 5, 3 : 2}) - 3
            -3 +5 x_2 +2 x_3
            sage: LinearFunction({0 : 1, 3 : -8}) - LinearFunction({2 : 5, 3 : 2}) - 16
            -16 +x_0 -5 x_2 -10 x_3
        """
        if isinstance(b, LinearFunction):
            e = deepcopy(self._f)
            for (id,coeff) in b.dict().iteritems():
                e[id] = self._f.get(id,0) - coeff
            return LinearFunction(e)
        else:
            el = deepcopy(self)
            el.dict()[-1] = self._f.get(-1,0) - b
            return el

    def __radd__(self,b):
        r"""
        Defining the + operator (right side).

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: 3 + LinearFunction({2 : 5, 3 : 2})
            3 +5 x_2 +2 x_3
        """
        if isinstance(self,LinearFunction):
            return self.__add__(b)
        else:
            return b.__add__(self)

    def __rsub__(self,b):
        r"""
        Defining the - operator (right side).

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: 3 - LinearFunction({2 : 5, 3 : 2})
            3 -5 x_2 -2 x_3
        """
        if isinstance(self,LinearFunction):
            return (-self).__add__(b)
        else:
            return b.__sub__(self)

    def __mul__(self,b):
        r"""
        Defining the * operator.

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({2 : 5, 3 : 2}) * 3
            15 x_2 +6 x_3
        """
        return LinearFunction(dict([(id,b*coeff) for (id, coeff) in self._f.iteritems()]))

    def __rmul__(self,b):
        r"""
        Defining the * operator (right side).

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: 3 * LinearFunction({2 : 5, 3 : 2})
            15 x_2 +6 x_3
        """
        return self.__mul__(b)

    def __repr__(self):
        r"""
        Returns a string version of the linear function.

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({2 : 5, 3 : 2})
            5 x_2 +2 x_3
        """
        cdef dict d = deepcopy(self._f)
        cdef bool first = True
        t = ""

        if d.has_key(-1):
            coeff = d.pop(-1)
            if coeff!=0:
                t = str(coeff)
                first = False

        cdef list l = sorted(d.items())
        for id,coeff in l:
            if coeff!=0:
                if not first:
                    t += " "
                t += ("+" if (not first and coeff >= 0) else "") + (str(coeff) + " " if coeff != 1 else "") + "x_" + str(id)
                first = False
        return t

    def __le__(self,other):
        r"""
        Defines the <= operator

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({2 : 5, 3 : 2}) <= LinearFunction({2 : 3, 9 : 2})
            5 x_2 +2 x_3 <= 3 x_2 +2 x_9
        """
        return LinearConstraint(self).__le__(other)

    def __lt__(self,other):
        r"""
        Defines the < operator

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({2 : 5, 3 : 2}) < LinearFunction({2 : 3, 9 : 2})
            Traceback (most recent call last):
            ...
            ValueError: The strict operators are not defined. Use <= and >= instead.
        """
        return LinearConstraint(self).__lt__(other)

    def __gt__(self,other):
        r"""
        Defines the > operator

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({2 : 5, 3 : 2}) > LinearFunction({2 : 3, 9 : 2})
            Traceback (most recent call last):
            ...
            ValueError: The strict operators are not defined. Use <= and >= instead.
        """
        return LinearConstraint(self).__gt__(other)


    def __ge__(self,other):
        r"""
        Defines the >= operator

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({2 : 5, 3 : 2}) >= LinearFunction({2 : 3, 9 : 2})
            3 x_2 +2 x_9 <= 5 x_2 +2 x_3
        """
        return LinearConstraint(self).__ge__(other)

    def __hash__(self):
        r"""
        Defines a ``__hash__`` function

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: d = {}
            sage: d[LinearFunction({2 : 5, 3 : 2})] = 3
        """
        return id(self)

    def __eq__(self,other):
        r"""
        Defines the == operator

        EXAMPLE::

            sage: from sage.numerical.mip import LinearFunction
            sage: LinearFunction({2 : 5, 3 : 2}) == LinearFunction({2 : 3, 9 : 2})
            5 x_2 +2 x_3 = 3 x_2 +2 x_9
        """
        return LinearConstraint(self).__eq__(other)


class LinearConstraint:
    """
    A class to represent formal Linear Constraints.

    A Linear Constraint being an inequality between
    two linear functions, this class lets the user
    write ``LinearFunction1 <= LinearFunction2``
    to define the corresponding constraint, which
    can potentially involve several layers of such
    inequalities (``(A <= B <= C``), or even equalities
    like ``A == B``.

    This class has no reason to be instanciated by the
    user, and is meant to be used by instances of
    MixedIntegerLinearProgram.

    INPUT:

    - ``c`` -- A ``LinearFunction``

    EXAMPLE::

        sage: p = MixedIntegerLinearProgram()
        sage: b = p.new_variable()
        sage: b[2]+2*b[3] <= b[8]-5
        x_0 +2 x_1 <= -5 +x_2
    """

    def __init__(self, c):
        r"""
        Constructor for ``LinearConstraint``

        INPUT:

        - ``c`` -- A linear function (see ``LinearFunction``).

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: b[2]+2*b[3] <= b[8]-5
            x_0 +2 x_1 <= -5 +x_2
        """
        self.equality = False
        self.constraints = []
        if isinstance(c, LinearFunction):
            self.constraints.append(c)
        else:
            self.constraints.append(LinearFunction(c))

    def __repr__(self):
        r"""
        Returns a string representation of the constraint.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: print b[3] <= b[8] + 9
            x_0 <= 9 +x_1
        """
        if self.equality:
            return str(self.constraints[0]) + " = " + str(self.constraints[1])
        else:
            first = True
            s = ""
            for c in self.constraints:
                s += (" <= " if not first else "") + c.__repr__()
                first = False
            return s

    def __eq__(self,other):
        r"""
        Defines the == operator

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: print b[3] == b[8] + 9
            x_0 = 9 +x_1
        """
        if not isinstance(other, LinearConstraint):
            other = LinearConstraint(other)

        if len(self.constraints) == 1 and len(other.constraints) == 1:
            self.constraints.extend(other.constraints)
            self.equality = True
            return self
        else:
            raise ValueError("Impossible to mix equality and inequality in the same equation")

    def __le__(self,other):
        r"""
        Defines the <= operator

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: print b[3] <= b[8] + 9
            x_0 <= 9 +x_1
        """

        if not isinstance(other, LinearConstraint):
            other = LinearConstraint(other)

        if self.equality or other.equality:
            raise ValueError("Impossible to mix equality and inequality in the same equation")

        self.constraints.extend(other.constraints)
        return self

    def __lt__(self, other):
        r"""
        Prevents the use of the stricts operators ``<`` and ``>``

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: print b[3] < b[8] + 9
            Traceback (most recent call last):
            ...
            ValueError: The strict operators are not defined. Use <= and >= instead.
            sage: print b[3] > b[8] + 9
            Traceback (most recent call last):
            ...
            ValueError: The strict operators are not defined. Use <= and >= instead.
        """
        raise ValueError("The strict operators are not defined. Use <= and >= instead.")

    __gt__ = __lt__

    def __ge__(self,other):
        r"""
        Defines the >= operator

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: print b[3] >= b[8] + 9
            9 +x_1 <= x_0
        """
        if not isinstance(other, LinearConstraint):
            other = LinearConstraint(other)

        if self.equality or other.equality:
            raise ValueError("Impossible to mix equality and inequality in the same equation")

        self.constraints = other.constraints + self.constraints
        return self
