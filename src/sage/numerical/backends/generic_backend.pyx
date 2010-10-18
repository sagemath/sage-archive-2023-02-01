r"""
Generic Backend for LP solvers

This class only lists the methods that should be defined by any
interface with a LP Solver. All these methods immediately raise
``NotImplementedError`` exceptions when called, and are obviously
meant to be replaced by the solver-specific method. This file can also
be used as a template to create a new interface : one would only need
to replace the occurences of ``"Nonexistent_LP_solver"`` by the
solver's name, and replace ``GenericBackend`` by
``SolverName(GenericBackend)`` so that the new solver extends this
class.
"""

cdef class GenericBackend:

    cpdef int add_variable(self):
        r"""
        Adds a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")    # optional - Nonexistent_LP_solver
            sage: p.ncols()                                         # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                   # optional - Nonexistent_LP_solver
            1
            sage: p.ncols()                                         # optional - Nonexistent_LP_solver
            1
        """
        raise NotImplementedError()

    cpdef int add_variables(self, int n):
        r"""
        Adds ``number`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")    # optional - Nonexistent_LP_solver
            sage: p.ncols()                                         # optional - Nonexistent_LP_solver
            0
            sage: p.add_variables(5)                                 # optional - Nonexistent_LP_solver
            5
            sage: p.ncols()                                         # optional - Nonexistent_LP_solver
            5
        """

        raise NotImplementedError()

    cpdef  set_variable_type(self, int variable, int vtype):
        r"""
        Sets the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1  Integer
            *  0  Binary
            * -1 Real

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")   # optional - Nonexistent_LP_solver
            sage: p.ncols()                                        # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                  # optional - Nonexistent_LP_solver
            1
            sage: p.set_variable_type(0,1)                          # optional - Nonexistent_LP_solver
            sage: p.is_variable_integer(0)                          # optional - Nonexistent_LP_solver
            True
        """

        raise NotImplementedError()

    cpdef set_sense(self, int sense):
        r"""
        Sets the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.is_maximization()                              # optional - Nonexistent_LP_solver
            True
            sage: p.set_sense(-1)                              # optional - Nonexistent_LP_solver
            sage: p.is_maximization()                              # optional - Nonexistent_LP_solver
            False
        """

        raise NotImplementedError()

    cpdef  set_objective_coefficient(self, int variable, double coeff):
        r"""
        Sets the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.get_objective_coeff(0)                         # optional - Nonexistent_LP_solver
            0.0
            sage: p.set_objective_coefficient(0,2)                       # optional - Nonexistent_LP_solver
            sage: p.get_objective_coeff(0)                         # optional - Nonexistent_LP_solver
            2.0
        """

        raise NotImplementedError()

    cpdef  set_objective(self, list coeff):
        r"""
        Sets the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")    # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                                 # optional - Nonexistent_LP_solver
            5
            sage: p.set_objective([1, 1, 2, 1, 3])                   # optional - Nonexistent_LP_solver
            sage: map(lambda x :p.get_objective_coeff(x), range(5))  # optional - Nonexistent_LP_solver
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """

        raise NotImplementedError()

    cpdef set_verbosity(self, int level):
        r"""
        Sets the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")   # optional - Nonexistent_LP_solver
            sage: p.set_verbosity(2)                                # optional - Nonexistent_LP_solver

        """

        raise NotImplementedError()

    cpdef add_constraint(self, list indices, list coeffs, int direction, double bound):
        r"""
        Adds a constraint.

        INPUT:

        - ``indices`` (list of integers) -- this list constains the
          indices of the variables whose coefficient is nonzero in the
          constraint.

        - ``coeffs`` (list of real values) -- associates a coefficient
          to the variables listed by ``indices``. Namely, the ith
          entry of ``coeffs`` corresponds to the coefficient of the
          variable represented by the ith entry in ``indices``.

        - ``direction`` (integer) -- the direction of the constraint,
          where :

              * +1 indicates : function `\leq` ``bound``
              *  0 indicates : function `=` ``bound``
              * -1 indicates : function `\geq` ``bound``

        - ``bound`` (double) -- value of the right-hand side (as
          illustrated immediately above).

        .. NOTE::

            ``indices`` and ``coeffs`` are expected to be of the same
            length.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                              # optional - Nonexistent_LP_solver
            5
            sage: p.add_constraint(range(5), range(5), 0, 2)      # optional - Nonexistent_LP_solver
            sage: p.row(0)                                    # optional - Nonexistent_LP_solver
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)                             # optional - Nonexistent_LP_solver
            (2.0, 2.0)
        """
        raise NotImplementedError()

    cpdef add_col(self, list indices, list coeffs):
        r"""
        Adds a column.

        INPUT:

        - ``indices`` (list of integers) -- this list constains the
          indices of the constraints in which the variable's
          coefficient is nonzero

        - ``coeffs`` (list of real values) -- associates a coefficient
          to the variable in each of the constraints in which it
          appears. Namely, the ith entry of ``coeffs`` corresponds to
          the coefficient of the variable in the constraint
          represented by the ith entry in ``indices``.

        .. NOTE::

            ``indices`` and ``coeffs`` are expected to be of the same
            length.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.nrows()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_constraints(5, -1, 0)                      # optional - Nonexistent_LP_solver
            sage: p.add_col(range(5), range(5))                    # optional - Nonexistent_LP_solver
            sage: p.nrows()                                       # optional - Nonexistent_LP_solver
            5
        """

        raise NotImplementedError()

    cpdef add_constraints(self, int number, int direction, double bound):
        r"""
        Adds constraints.

        INPUT:

        - ``number`` (integer) -- the number of constraints to add.

        - ``direction`` (integer) -- the direction of the constraint,
          where :

              * +1 indicates : function `\leq` ``bound``
              *  0 indicates : function `=` ``bound``
              * -1 indicates : function `\geq` ``bound``

        - ``bound`` (double) -- value of the right-hand side (as
          illustrated immediately above).

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")   # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                                # optional - Nonexistent_LP_solver
            5
            sage: p.add_constraints(5, +1, 2)                       # optional - Nonexistent_LP_solver
            sage: p.row(4)                                      # optional - Nonexistent_LP_solver
            ([], [])
            sage: p.row_bounds(4)                               # optional - Nonexistent_LP_solver
            (None, 2.0)
        """

        raise NotImplementedError()

    cpdef int solve(self) except -1:
        r"""
        Solves the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_constraints(5, -1, 0)                     # optional - Nonexistent_LP_solver
            sage: p.add_col(range(5), range(5))                   # optional - Nonexistent_LP_solver
            sage: p.solve()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.set_objective_coefficient(0,1)                      # optional - Nonexistent_LP_solver
            sage: p.solve()                                       # optional - Nonexistent_LP_solver
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """
        raise NotImplementedError()

    cpdef double get_objective_value(self):
        r"""
        Returns the value of the objective function.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)          # optional - Nonexistent_LP_solver
            sage: p.set_objective([2, 5])                          # optional - Nonexistent_LP_solver
            sage: p.solve()                                        # optional - Nonexistent_LP_solver
            0
            sage: p.get_objective_value()                          # optional - Nonexistent_LP_solver
            7.5
            sage: p.get_variable_value(0)                          # optional - Nonexistent_LP_solver
            0.0
            sage: p.get_variable_value(1)                          # optional - Nonexistent_LP_solver
            1.5
        """

        raise NotImplementedError()

    cpdef double get_variable_value(self, int variable):
        r"""
        Returns the value of a variable given by the solver.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                              # optional - Nonexistent_LP_solver
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)         # optional - Nonexistent_LP_solver
            sage: p.set_objective([2, 5])                         # optional - Nonexistent_LP_solver
            sage: p.solve()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.get_objective_value()                         # optional - Nonexistent_LP_solver
            7.5
            sage: p.get_variable_value(0)                         # optional - Nonexistent_LP_solver
            0.0
            sage: p.get_variable_value(1)                         # optional - Nonexistent_LP_solver
            1.5
        """

        raise NotImplementedError()

    cpdef int ncols(self):
        r"""
        Returns the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            2
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            2
        """

        raise NotImplementedError()

    cpdef int nrows(self):
        r"""
        Returns the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.nrows()                                      # optional - Nonexistent_LP_solver
            0
            sage: p.add_constraints(2, -1, 2)                     # optional - Nonexistent_LP_solver
            sage: p.nrows()                                      # optional - Nonexistent_LP_solver
            2
        """

        raise NotImplementedError()

    cpdef name(self):
        raise NotImplementedError()

    cpdef bint is_maximization(self):
        r"""
        Tests whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.is_maximization()                             # optional - Nonexistent_LP_solver
            True
            sage: p.set_sense(-1)                             # optional - Nonexistent_LP_solver
            sage: p.is_maximization()                             # optional - Nonexistent_LP_solver
            False
        """
        raise NotImplementedError()

    cpdef problem_name(self, char * name = NULL):
        r"""
        Returns or defines the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")   # optional - Nonexistent_LP_solver
            sage: p.problem_name("There once was a french fry") # optional - Nonexistent_LP_solver
            sage: print p.get_problem_name()                        # optional - Nonexistent_LP_solver
            There once was a french fry
        """

        raise NotImplementedError()

    cpdef write_lp(self, char * name):
        r"""
        Writes the problem to a .lp file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)          # optional - Nonexistent_LP_solver
            sage: p.set_objective([2, 5])                          # optional - Nonexistent_LP_solver
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp")            # optional - Nonexistent_LP_solver
        """
        raise NotImplementedError()

    cpdef write_mps(self, char * name, int modern):
        r"""
        Writes the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)          # optional - Nonexistent_LP_solver
            sage: p.set_objective([2, 5])                          # optional - Nonexistent_LP_solver
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp")            # optional - Nonexistent_LP_solver
        """
        raise NotImplementedError()

    cpdef row(self, int i):
        r"""
        Returns a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_constraint`` method.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                               # optional - Nonexistent_LP_solver
            5
            sage: p.add_constraint(range(5), range(5), 0, 2)       # optional - Nonexistent_LP_solver
            sage: p.row(0)                                     # optional - Nonexistent_LP_solver
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)                              # optional - Nonexistent_LP_solver
            (2.0, 2.0)
        """

        raise NotImplementedError()

    cpdef double get_objective_coeff(self, int i):
        r"""
        Sets the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.get_objective_coeff(0)                         # optional - Nonexistent_LP_solver
            0.0
            sage: p.set_objective_coefficient(0,2)                       # optional - Nonexistent_LP_solver
            sage: p.get_objective_coeff(0)                         # optional - Nonexistent_LP_solver
            2.0
        """

        raise NotImplementedError()

    cpdef row_bounds(self, int index):
        r"""
        Returns the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                               # optional - Nonexistent_LP_solver
            5
            sage: p.add_constraint(range(5), range(5), 0, 2)       # optional - Nonexistent_LP_solver
            sage: p.row(0)                                     # optional - Nonexistent_LP_solver
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)                              # optional - Nonexistent_LP_solver
            (2.0, 2.0)
        """
        raise NotImplementedError()

    cpdef col_bounds(self, int index):
        r"""
        Returns the bounds of a specific variable.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, None)
            sage: p.variable_max(0, 5)                         # optional - Nonexistent_LP_solver
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, 5.0)
        """

        raise NotImplementedError()

    cpdef bint is_variable_binary(self, int index):
        r"""
        Tests whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.set_variable_type(0,0)                         # optional - Nonexistent_LP_solver
            sage: p.is_variable_binary(0)                          # optional - Nonexistent_LP_solver
            True

        """

        raise NotImplementedError()

    cpdef bint is_variable_integer(self, int index):
        r"""
        Tests whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.set_variable_type(0,1)                         # optional - Nonexistent_LP_solver
            sage: p.is_variable_integer(0)                         # optional - Nonexistent_LP_solver
            True
        """

        raise NotImplementedError()

    cpdef bint is_variable_continuous(self, int index):
        r"""
        Tests whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.is_variable_continuous(0)                      # optional - Nonexistent_LP_solver
            True
            sage: p.set_variable_type(0,1)                         # optional - Nonexistent_LP_solver
            sage: p.is_variable_continuous(0)                      # optional - Nonexistent_LP_solver
            False

        """

        raise NotImplementedError()

    cpdef row_name(self, int index, char * name = NULL):
        r"""
        Returns or defines the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_constraints(1, -1, 2)                      # optional - Nonexistent_LP_solver
            sage: p.row_name(0, "Empty constraint 1")          # optional - Nonexistent_LP_solver
            sage: p.row_name(0)                                # optional - Nonexistent_LP_solver
            'Empty constraint 1'

        """

        raise NotImplementedError()

    cpdef col_name(self, int index, char * name = NULL):
        r"""
        Returns or defines the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.col_name(0, "I am a variable")             # optional - Nonexistent_LP_solver
            sage: p.col_name(0)                                # optional - Nonexistent_LP_solver
            'I am a variable'
        """

        raise NotImplementedError()

    cpdef  variable_max(self, int index, value = None):
        r"""
        Returns or defines the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``None``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, None)
            sage: p.variable_max(0, 5)                         # optional - Nonexistent_LP_solver
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, 5.0)
        """

        raise NotImplementedError()

    cpdef  variable_min(self, int index, value = None):
        r"""
        Returns or defines the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``None``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            1
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, None)
            sage: p.variable_min(0, 5)                         # optional - Nonexistent_LP_solver
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (5.0, None)
        """

        raise NotImplementedError()


cpdef GenericBackend get_solver(constraint_generation = False, solver = None):
    r"""
    Returns a solver according to the given preferences

    INPUT:

    - ``solver`` -- 3 solvers should be available through this class:

        - GLPK (``solver="GLPK"``). See the `GLPK
          <http://www.gnu.org/software/glpk/>`_ web site.

        - COIN Branch and Cut (``solver="Coin"``). See the `COIN-OR
          <http://www.coin-or.org>`_ web site.

        - CPLEX (``solver="CPLEX"``). See the
          `CPLEX <http://www.ilog.com/products/cplex/>`_ web site.
          An interface to CPLEX is not yet implemented.

        ``solver`` should then be equal to one of ``"GLPK"``,
        ``"Coin"``, ``"CPLEX"``, or ``None``. If ``solver=None``
        (default), the solvers are tried in this order :

            * CPLEX
            * Coin
            * GLPK

        A backend corresponding to the first solver available is then
        returned

    - ``constraint_generation`` (boolean) -- whether the solver
      returned is to be used for constraint/variable generation. As
      the interface with Coin does not support constraint/variable
      generation, setting ``constraint_generation`` to ``False``
      ensures that the backend to Coin is not returned when ``solver =
      None``. This is set to ``False`` by default.

    EXAMPLE::

        sage: from sage.numerical.backends.generic_backend import get_solver
        sage: p = get_solver()

    """

    if solver is None:

        try:
            from sage.numerical.backends.cplex_backend import CPLEXBackend
            return CPLEXBackend()
        except ImportError:
            pass

        try:
            if not constraint_generation:
                from sage.numerical.backends.coin_backend import CoinBackend
                return CoinBackend()
        except ImportError:
            pass

        from sage.numerical.backends.glpk_backend import GLPKBackend
        return GLPKBackend()

    elif solver == "Coin":
        from sage.numerical.backends.coin_backend import CoinBackend
        return CoinBackend()

    elif solver == "GLPK":
        from sage.numerical.backends.glpk_backend import GLPKBackend
        return GLPKBackend()

    elif solver == "CPLEX":
        from sage.numerical.backends.cplex_backend import CPLEXBackend
        return CPLEXBackend()

    else:
        raise ValueError("'solver' should be set to 'GLPK', 'Coin', 'CPLEX' or None (in which case the default one is used).")


