r"""
Mixed integer linear programming

A linear program (`LP <http://en.wikipedia.org/wiki/Linear_programming>`_)
is an `optimization problem <http://en.wikipedia.org/wiki/Optimization_%28mathematics%29>`_
in the following form

.. MATH::
     \max \{ c^T x \;|\; A x \leq b, x \geq 0 \}

with given `A \in \mathbb{R}^{m,n}`, `b \in \mathbb{R}^m`,
`c \in \mathbb{R}^n` and unknown `x \in \mathbb{R}^{n}`.
If some or all variables in the vector `x` are restricted over
the integers `\mathbb{Z}`, the problem is called mixed integer
linear program (`MILP <http://en.wikipedia.org/wiki/Mixed_integer_linear_programming>`_).
A wide variety of problems in optimization
can be formulated in this standard form. Then, solvers are
able to calculate a solution.

Imagine you want to solve the following linear system of three equations:

 - `w_0 + w_1 + w_2 - 14 w_3 = 0`
 - `w_1 + 2 w_2 - 8 w_3 = 0`
 - `2 w_2 - 3 w_3 = 0`

and this additional inequality:

 - `w_0 - w_1 - w_2 \geq 0`

where all `w_i \in \mathbb{Z}^+`. You know that the trivial solution is `w_i=0`,
but what is the first non-trivial one with `w_3 \geq 1`?

A mixed integer linear program can give you an answer:

  #. You have to create an instance of :class:`MixedIntegerLinearProgram` and
     -- in our case -- specify that it is a minimization.
  #. Create an dictionary ``w`` of integer variables ``w`` via ``w =
     p.new_variable(integer=True)`` (note that **by default all variables are
     non-negative**, cf :meth:`~MixedIntegerLinearProgram.new_variable`).
  #. Add those three equations as equality constraints via
     :meth:`add_constraint <sage.numerical.mip.MixedIntegerLinearProgram.add_constraint>`.
  #. Also add the inequality constraint.
  #. Add an inequality constraint `w_3 \geq 1` to exclude the trivial solution.
  #. By default, all variables are non-negative. We remove that constraint
     via ``p.set_min(variable, None)``, see :meth:`set_min <sage.numerical.mip.MixedIntegerLinearProgram.set_min>`.
  #. Specify the objective function via :meth:`set_objective <sage.numerical.mip.MixedIntegerLinearProgram.set_objective>`.
     In our case that is just `w_3`. If it
     is a pure constraint satisfaction problem, specify it as ``None``.
  #. To check if everything is set up correctly, you can print the problem via
     :meth:`show <sage.numerical.mip.MixedIntegerLinearProgram.show>`.
  #. :meth:`Solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>` it and print the solution.

The following example shows all these steps::

    sage: p = MixedIntegerLinearProgram(maximization=False, solver = "GLPK")
    sage: w = p.new_variable(integer=True) # all variables are non-negative by default
    sage: p.add_constraint(w[0] + w[1] + w[2] - 14*w[3] == 0)
    sage: p.add_constraint(w[1] + 2*w[2] - 8*w[3] == 0)
    sage: p.add_constraint(2*w[2] - 3*w[3] == 0)
    sage: p.add_constraint(w[0] - w[1] - w[2] >= 0)
    sage: p.add_constraint(w[3] >= 1)
    sage: _ = [ p.set_min(w[i], None) for i in range(1,4) ]
    sage: p.set_objective(w[3])
    sage: p.show()
    Minimization:
       x_3
    Constraints:
      0.0 <= x_0 + x_1 + x_2 - 14.0 x_3 <= 0.0
      0.0 <= x_1 + 2.0 x_2 - 8.0 x_3 <= 0.0
      0.0 <= 2.0 x_2 - 3.0 x_3 <= 0.0
      - x_0 + x_1 + x_2 <= 0.0
      - x_3 <= -1.0
    Variables:
      x_0 is an integer variable (min=0.0, max=+oo)
      x_1 is an integer variable (min=-oo, max=+oo)
      x_2 is an integer variable (min=-oo, max=+oo)
      x_3 is an integer variable (min=-oo, max=+oo)
    sage: print 'Objective Value:', p.solve()
    Objective Value: 2.0
    sage: for i, v in p.get_values(w).iteritems():\
              print 'w_%s = %s' % (i, int(round(v)))
    w_0 = 15
    w_1 = 10
    w_2 = 3
    w_3 = 2

Different backends compute with different base fields, for example::

    sage: p = MixedIntegerLinearProgram(solver = 'GLPK')
    sage: p.base_ring()
    Real Double Field
    sage: x = p.new_variable()
    doctest:839: DeprecationWarning: The default behaviour of new_variable() will soon change ! It will return 'real' variables instead of nonnegative ones. Please be explicit and call new_variable(nonnegative=True) instead.
    See http://trac.sagemath.org/15521 for details.
    sage: 0.5 + 3/2*x[1]
    0.5 + 1.5*x_0

    sage: p = MixedIntegerLinearProgram(solver = 'ppl')
    sage: p.base_ring()
    Rational Field
    sage: x = p.new_variable()
    sage: 0.5 + 3/2*x[1]
    1/2 + 3/2*x_0

AUTHORS:

- Risan (2012/02): added extension for exact computation

Index of functions and methods
------------------------------

Below are listed the methods of :class:`MixedIntegerLinearProgram`. This module
also implements the :class:`MIPSolverException` exception, as well as the
:class:`MIPVariable` class.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~MixedIntegerLinearProgram.add_constraint`            | Adds a constraint to the ``MixedIntegerLinearProgram``
    :meth:`~MixedIntegerLinearProgram.base_ring`                 | Return the base ring
    :meth:`~MixedIntegerLinearProgram.constraints`               | Returns a list of constraints, as 3-tuples
    :meth:`~MixedIntegerLinearProgram.get_backend`               | Returns the backend instance used
    :meth:`~MixedIntegerLinearProgram.get_max`                   | Returns the maximum value of a variable
    :meth:`~MixedIntegerLinearProgram.get_min`                   | Returns the minimum value of a variable
    :meth:`~MixedIntegerLinearProgram.get_values`                | Return values found by the previous call to ``solve()``
    :meth:`~MixedIntegerLinearProgram.is_binary`                 | Tests whether the variable ``e`` is binary
    :meth:`~MixedIntegerLinearProgram.is_integer`                | Tests whether the variable is an integer
    :meth:`~MixedIntegerLinearProgram.is_real`                   | Tests whether the variable is real
    :meth:`~MixedIntegerLinearProgram.linear_constraints_parent` | Return the parent for all linear constraints
    :meth:`~MixedIntegerLinearProgram.linear_function`           | Construct a new linear function
    :meth:`~MixedIntegerLinearProgram.linear_functions_parent`   | Return the parent for all linear functions
    :meth:`~MixedIntegerLinearProgram.new_variable`              | Returns an instance of ``MIPVariable`` associated
    :meth:`~MixedIntegerLinearProgram.number_of_constraints`     | Returns the number of constraints assigned so far
    :meth:`~MixedIntegerLinearProgram.number_of_variables`       | Returns the number of variables used so far
    :meth:`~MixedIntegerLinearProgram.polyhedron`                | Returns the polyhedron defined by the Linear Program
    :meth:`~MixedIntegerLinearProgram.remove_constraint`         | Removes a constraint from self
    :meth:`~MixedIntegerLinearProgram.remove_constraints`        | Remove several constraints
    :meth:`~MixedIntegerLinearProgram.set_binary`                | Sets a variable or a ``MIPVariable`` as binary
    :meth:`~MixedIntegerLinearProgram.set_integer`               | Sets a variable or a ``MIPVariable`` as integer
    :meth:`~MixedIntegerLinearProgram.set_max`                   | Sets the maximum value of a variable
    :meth:`~MixedIntegerLinearProgram.set_min`                   | Sets the minimum value of a variable
    :meth:`~MixedIntegerLinearProgram.set_objective`             | Sets the objective of the ``MixedIntegerLinearProgram``
    :meth:`~MixedIntegerLinearProgram.set_problem_name`          | Sets the name of the ``MixedIntegerLinearProgram``
    :meth:`~MixedIntegerLinearProgram.set_real`                  | Sets a variable or a ``MIPVariable`` as real
    :meth:`~MixedIntegerLinearProgram.show`                      | Displays the ``MixedIntegerLinearProgram`` in a human-readable
    :meth:`~MixedIntegerLinearProgram.solve`                     | Solves the ``MixedIntegerLinearProgram``
    :meth:`~MixedIntegerLinearProgram.solver_parameter`          | Return or define a solver parameter
    :meth:`~MixedIntegerLinearProgram.sum`                       | Efficiently computes the sum of a sequence of LinearFunction elements
    :meth:`~MixedIntegerLinearProgram.write_lp`                  | Write the linear program as a LP file
    :meth:`~MixedIntegerLinearProgram.write_mps`                 | Write the linear program as a MPS file

Classes and methods
-------------------
"""

#*****************************************************************************
#       Copyright (C) 2012 Nathann Cohen <nathann.cohen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/stdsage.pxi"
include "sage/ext/interrupt.pxi"
include "sage/ext/cdefs.pxi"

from sage.structure.sage_object cimport SageObject
from sage.misc.cachefunc import cached_method
from sage.numerical.linear_functions import is_LinearFunction, is_LinearConstraint
from sage.misc.superseded import deprecated_function_alias, deprecation
from sage.misc.superseded import deprecated_function_alias

cdef class MixedIntegerLinearProgram(SageObject):
    r"""
    The ``MixedIntegerLinearProgram`` class is the link between Sage, linear
    programming (LP) and mixed integer programming (MIP) solvers.

    A Mixed Integer Linear Program (MILP) consists of variables, linear
    constraints on these variables, and an objective function which is to be
    maximised or minimised under these constraints.

    See the :wikipedia:`Linear_programming` for further information on linear
    programming, and the :mod:`MILP module <sage.numerical.mip>` for its use in
    Sage.

    INPUT:

    - ``solver`` -- selects a solver:

      - GLPK (``solver="GLPK"``). See the `GLPK
        <http://www.gnu.org/software/glpk/>`_ web site.

      - COIN Branch and Cut (``solver="Coin"``). See the `COIN-OR
        <http://www.coin-or.org>`_ web site.

      - CPLEX (``solver="CPLEX"``). See the `CPLEX
        <http://www.ilog.com/products/cplex/>`_ web site.

      - Gurobi (``solver="Gurobi"``). See the `Gurobi <http://www.gurobi.com/>`_
          web site.

      - PPL (``solver="PPL"``). See the `PPL <http://bugseng.com/products/ppl>`_
          web site.

      - If ``solver=None`` (default), the default solver is used (see
        :func:`default_mip_solver`)

    - ``maximization``

      - When set to ``True`` (default), the ``MixedIntegerLinearProgram``
        is defined as a maximization.

      - When set to ``False``, the ``MixedIntegerLinearProgram`` is
        defined as a minimization.

    - ``constraint_generation`` -- Only used when ``solver=None``.

      - When set to ``True``, after solving the ``MixedIntegerLinearProgram``,
        it is possible to add a constraint, and then solve it again.
        The effect is that solvers that do not support this feature will not be
        used.

      - Defaults to ``False``.

    .. WARNING::

        All LP variables are non-negative by default (see :meth:`new_variable`
        and :meth:`set_min`).

    .. SEEALSO::

     - :func:`default_mip_solver` -- Returns/Sets the default MIP solver.

    EXAMPLES:

    Computation of a maximum stable set in Petersen's graph::

         sage: g = graphs.PetersenGraph()
         sage: p = MixedIntegerLinearProgram(maximization=True)
         sage: b = p.new_variable()
         sage: p.set_objective(sum([b[v] for v in g]))
         sage: for (u,v) in g.edges(labels=None):
         ...       p.add_constraint(b[u] + b[v], max=1)
         sage: p.set_binary(b)
         sage: p.solve(objective_only=True)
         4.0
    """

    def __init__(self, solver=None, maximization=True,
                 constraint_generation=False, check_redundant=False):
        r"""
        Constructor for the ``MixedIntegerLinearProgram`` class.

        INPUT:

        - ``solver`` -- the following solvers should be available through this class:

          - GLPK (``solver="GLPK"``). See the `GLPK
            <http://www.gnu.org/software/glpk/>`_ web site.

          - COIN Branch and Cut (``solver="Coin"``). See the `COIN-OR
            <http://www.coin-or.org>`_ web site.

          - CPLEX (``solver="CPLEX"``). See the `CPLEX
            <http://www.ilog.com/products/cplex/>`_ web site.  An interface to
            CPLEX is not yet implemented.

          - Gurobi (``solver="Gurobi"``). See the `Gurobi
            <http://www.gurobi.com/>`_ web site.

          - PPL (``solver="PPL"``). See the `PPL
            <http://bugseng.com/products/ppl>`_ web site.

          -If ``solver=None`` (default), the default solver is used (see
           ``default_mip_solver`` method.

        - ``maximization``

          - When set to ``True`` (default), the ``MixedIntegerLinearProgram``
            is defined as a maximization.
          - When set to ``False``, the ``MixedIntegerLinearProgram`` is
            defined as a minimization.

        - ``constraint_generation`` -- Only used when ``solver=None``.

          - When set to ``True``, after solving the 
            ``MixedIntegerLinearProgram``, it is possible to add a constraint,
            and then solve it again. The effect is that solvers that do not
            support this feature will not be used.

          - Defaults to ``False``.

        - ``check_redundant`` -- whether to check that constraints added to the
          program are redundant with constraints already in the program.
          Only obvious redundancies are checked: to be considered redundant,
          either a constraint is equal to another constraint in the program,
          or it is a constant multiple of the other. To make this search
          effective and efficient, constraints are normalized; thus, the
          constraint `-x_1 < 0` will be stored as `x_1 > 0`.

        .. SEEALSO::

        - :meth:`default_mip_solver` -- Returns/Sets the default MIP solver.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True)

        TESTS:

        Checks that the objects are deallocated without invoking the cyclic garbage
        collector (cf. :trac:`12616`)::

            sage: del p
            sage: def just_create_variables():
            ...       p = MixedIntegerLinearProgram()
            ...       b = p.new_variable(nonnegative=True)
            ...       p.add_constraint(b[3]+b[6] <= 2)
            ...       p.solve()
            sage: C = sage.numerical.mip.MixedIntegerLinearProgram
            sage: import gc
            sage: _ = gc.collect()  # avoid side effects of other doc tests
            sage: sum([1 for x in gc.get_objects() if isinstance(x,C)])
            0

        We now disable the cyclic garbage collector. Since :trac:`12616` avoids
        a reference cycle, the mixed integer linear program created in
        ``just_create_variables()`` is removed even without the cyclic garbage
        collection::

            sage: gc.disable()
            sage: just_create_variables()
            sage: sum([1 for x in gc.get_objects() if isinstance(x,C)])
            0
            sage: gc.enable()
        """
        self.__BINARY = 0
        self.__REAL = -1
        self.__INTEGER = 1

        from sage.numerical.backends.generic_backend import get_solver
        self._backend = get_solver(solver=solver,
                                   constraint_generation=constraint_generation)
        if not maximization:
            self._backend.set_sense(-1)

        # Associates an index to the variables
        self._variables = {}

        # Check for redundant constraints
        self._check_redundant = check_redundant
        if check_redundant:
            self._constraints = list()

    def linear_functions_parent(self):
        """
        Return the parent for all linear functions

        EXAMPLES::

             sage: p = MixedIntegerLinearProgram()
             sage: p.linear_functions_parent()
             Linear functions over Real Double Field
        """
        if self._linear_functions_parent is None:
            base_ring = self._backend.base_ring()
            from sage.numerical.linear_functions import LinearFunctionsParent
            self._linear_functions_parent = LinearFunctionsParent(base_ring)
        return self._linear_functions_parent

    def linear_constraints_parent(self):
        """
        Return the parent for all linear constraints

        See :mod:`~sage.numerical.linear_functions` for more
        details.

        EXAMPLES::

             sage: p = MixedIntegerLinearProgram()
             sage: p.linear_constraints_parent()
             Linear constraints over Real Double Field
        """
        if self._linear_constraints_parent is None:
            from sage.numerical.linear_functions import LinearConstraintsParent
            LF = self.linear_functions_parent()
            self._linear_constraints_parent = LinearConstraintsParent(LF)
        return self._linear_constraints_parent

    def __call__(self, x):
        """
        Construct a new linear function

        EXAMPLES::

             sage: p = MixedIntegerLinearProgram()
             sage: p.linear_function({1:3, 4:5})
             3*x_1 + 5*x_4

        This is equivalent to::

            sage: p({1:3, 4:5})
            3*x_1 + 5*x_4
        """
        parent = self.linear_functions_parent()
        return parent(x)

    linear_function = __call__

    def _repr_(self):
         r"""
         Returns a short description of the ``MixedIntegerLinearProgram``.

         EXAMPLE::

             sage: p = MixedIntegerLinearProgram()
             sage: v = p.new_variable()
             sage: p.add_constraint(v[1] + v[2], max=2)
             sage: print p
             Mixed Integer Program ( maximization, 2 variables, 1 constraints )
         """
         cdef GenericBackend b = self._backend

         return ("Mixed Integer Program "+

                 ( "\"" +self._backend.problem_name()+ "\""
                   if (str(self._backend.problem_name()) != "") else "")+

                 " ( " + ("maximization" if b.is_maximization() else "minimization" ) +

                 ", " + str(b.ncols()) + " variables, " +
                 str(b.nrows()) + " constraints )")

    def __copy__(self):
        r"""
        Returns a copy of the current ``MixedIntegerLinearProgram`` instance.

        EXAMPLE::

        sage: p = MixedIntegerLinearProgram()
        sage: p.add_constraint(p[0] + p[1], max = 10)
        sage: q = copy(p)
        sage: q.number_of_constraints()
        1
        """
        cdef MixedIntegerLinearProgram p = \
            MixedIntegerLinearProgram(solver="GLPK")
        from copy import copy
        try:
            p._variables = copy(self._variables)
        except AttributeError:
            pass

        try:
            p._default_mipvariable = self._default_mipvariable
        except AttributeError:
            pass

        try:
            p._check_redundant = self._check_redundant
            p._constraints = copy(self._constraints)
        except AttributeError:
            pass

        p._backend = (<GenericBackend> self._backend).copy()
        return p

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
        except TypeError:
            self._default_mipvariable = self.new_variable()
            return self._default_mipvariable[v]

    def base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        A ring. The coefficients that the chosen solver supports.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver='GLPK')
            sage: p.base_ring()
            Real Double Field
            sage: p = MixedIntegerLinearProgram(solver='ppl')
            sage: p.base_ring()
            Rational Field
        """
        return self._backend.base_ring()

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
        self._backend.problem_name(name)

    def new_variable(self, real=False, nonnegative=False, binary=False, integer=False, dim=1,name=""):
        r"""
        Returns an instance of ``MIPVariable`` associated
        to the current instance of ``MixedIntegerLinearProgram``.

        A new variable ``x`` is defined by::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()

        It behaves exactly as an usual dictionary would. It can use any key
        argument you may like, as ``x[5]`` or ``x["b"]``, and has methods
        ``items()`` and ``keys()``.

        .. WARNING::

            By default, all ``x[i]`` are assumed to be non-negative. See
            :meth:`set_min` to set a different lower bound.

        INPUT:

        - ``dim`` (integer) -- Defines the dimension of the dictionary.
          If ``x`` has dimension `2`, its fields will be of the form
          ``x[key1][key2]``. Deprecated.

        - ``binary, integer, nonnegative`` (boolean) -- Set one of these
          arguments to ``True`` to ensure that the variable gets the
          corresponding type. The default type is ``nonnegative``, which
          represents nonnegative real variables.

        - ``name`` (string) -- Associates a name to the variable. This is
          only useful when exporting the linear program to a file using
          ``write_mps`` or ``write_lp``, and has no other effect.

        .. SEEALSO::

            - :meth:`set_min`,:meth:`get_min` -- set/get the lower bound of a
              variable. Note that by default, all variables are non-negative.

            - :meth:`set_max`,:meth:`get_max` -- set/get the upper bound of a
              variable.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()

         To define two dictionaries of variables, the first being
         of real type, and the second of integer type ::

            sage: x = p.new_variable(real=True)
            doctest:839: DeprecationWarning: The meaning of 'real' will change, to represent real variables instead of nonnegative ones. Please use the new 'nonnegative' variable type.
            See http://trac.sagemath.org/15521 for details.
            sage: y = p.new_variable(integer=True)
            sage: p.add_constraint(x[2] + y[3,5], max=2)
            sage: p.is_integer(x[2])
            False
            sage: p.is_integer(y[3,5])
            True

        An exception is raised when two types are supplied ::

            sage: z = p.new_variable(real = True, integer = True)
            Traceback (most recent call last):
            ...
            ValueError: Exactly one of the available types has to be True

        TESTS:

        Default behaviour (:trac:`15521`)::

            sage: x = p.new_variable()
            sage: p.get_min(x[0])
            0.0
        """
        if sum([real, binary, integer]) >= 2:
            raise ValueError("Exactly one of the available types has to be True")

        if binary:
            vtype = self.__BINARY
        elif integer:
            vtype = self.__INTEGER
        elif real:
            deprecation(15521, "The meaning of 'real' will change, to "+
                        "represent real variables instead of nonnegative "+
                        "ones. Please use the new 'nonnegative' variable type.")
            vtype = self.__REAL
        elif nonnegative:
            vtype = self.__REAL
        else:
            deprecation(15521, "The default behaviour of new_variable() will "+
                        "soon change ! It will return 'real' variables instead "+
                        "of nonnegative ones. Please be explicit and call "+
                        "new_variable(nonnegative=True) instead.")
            vtype = self.__REAL

        v=MIPVariable(self, vtype, dim=dim,name=name)
        return v

    cpdef int number_of_constraints(self):
      r"""
      Returns the number of constraints assigned so far.

      EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)
            sage: p.add_constraint(p[0] - 2*p[1], min = 1)
            sage: p.number_of_constraints()
            2
      """
      return self._backend.nrows()

    cpdef int number_of_variables(self):
      r"""
      Returns the number of variables used so far.

      Note that this is backend-dependent, i.e. we count solver's
      variables rather than user's variables. An example of the latter
      can be seen below: Gurobi converts double inequalities,
      i.e. inequalities like `m <= c^T x <= M`, with `m<M`, into
      equations, by adding extra variables: `c^T x + y = M`, `0 <= y
      <= M-m`.

      EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(p[0] - p[2], max = 4)
            sage: p.number_of_variables()
            2
            sage: p.add_constraint(p[0] - 2*p[1], min = 1)
            sage: p.number_of_variables()
            3
            sage: p = MixedIntegerLinearProgram(solver="glpk")
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)
            sage: p.number_of_variables()
            2
            sage: p = MixedIntegerLinearProgram(solver="gurobi")   # optional - Gurobi
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)  # optional - Gurobi
            sage: p.number_of_variables()                          # optional - Gurobi
            3
      """
      return self._backend.ncols()

    def constraints(self, indices = None):
        r"""
        Returns a list of constraints, as 3-tuples.

        INPUT:

        - ``indices`` -- select which constraint(s) to return

            - If ``indices = None``, the method returns the list of all the
              constraints.

            - If ``indices`` is an integer `i`, the method returns constraint
              `i`.

            - If ``indices`` is a list of integers, the method returns the list
              of the corresponding constraints.

        OUTPUT:

        Each constraint is returned as a triple ``lower_bound, (indices,
        coefficients), upper_bound``.  For each of those entries, the
        corresponding linear function is the one associating to variable
        ``indices[i]`` the coefficient ``coefficients[i]``, and `0` to all the
        others.

        ``lower_bound`` and ``upper_bound`` are numerical values.

        EXAMPLE:

        First, let us define a small LP::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)
            sage: p.add_constraint(p[0] - 2*p[1], min = 1)

        To obtain the list of all constraints::

            sage: p.constraints()          # not tested
            [(1.0, ([1, 0], [-1.0, 1.0]), 4.0), (1.0, ([2, 0], [-2.0, 1.0]), None)]

        Or constraint `0` only::

            sage: p.constraints(0)         # not tested
            (1.0, ([1, 0], [-1.0, 1.0]), 4.0)

        A list of constraints containing only `1`::

            sage: p.constraints([1])       # not tested
            [(1.0, ([2, 0], [-2.0, 1.0]), None)]

        TESTS:

        As the ordering of the variables in each constraint depends on the
        solver used, we define a short function reordering it before it is
        printed. The output would look the same without this function applied::

            sage: def reorder_constraint((lb,(ind,coef),ub)):
            ...     d = dict(zip(ind, coef))
            ...     ind.sort()
            ...     return (lb, (ind, [d[i] for i in ind]), ub)

        Running the examples from above, reordering applied::

            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: p.add_constraint(p[0] - p[2], min = 1, max = 4)
            sage: p.add_constraint(p[0] - 2*p[1], min = 1)
            sage: sorted(map(reorder_constraint,p.constraints()))
            [(1.0, ([0, 1], [1.0, -1.0]), 4.0), (1.0, ([0, 2], [1.0, -2.0]), None)]
            sage: reorder_constraint(p.constraints(0))
            (1.0, ([0, 1], [1.0, -1.0]), 4.0)
            sage: sorted(map(reorder_constraint,p.constraints([1])))
            [(1.0, ([0, 2], [1.0, -2.0]), None)]

        """
        from sage.rings.integer import Integer as Integer
        cdef int i
        cdef str s
        cdef GenericBackend b = self._backend

        result = list()

        # If indices is None, we actually want to return all constraints
        if indices is None:
          indices = range(b.nrows())

        # Only one constraint
        if isinstance(indices, int) or isinstance(indices, Integer):
            lb, ub = b.row_bounds(indices)
            return (lb, b.row(indices), ub)

        # List of constraints
        elif isinstance(indices, list):
            for i in indices:
                lb, ub = b.row_bounds(i)
                result.append((lb, b.row(i), ub))

            return result

        # Weird Input
        else:
          raise ValueError, "constraints() requires a list of integers, though it will accommodate None or an integer."

    def polyhedron(self, **kwds):
        r"""
        Returns the polyhedron defined by the Linear Program.

        INPUT:

        All arguments given to this method are forwarded to the constructor of
        the :func:`Polyhedron` class.

        OUTPUT:

        A :func:`Polyhedron` object whose `i`-th variable represents the `i`-th
        variable of ``self``.

        .. warning::

            The polyhedron is built from the variables stored by the LP solver
            (i.e. the output of :meth:`show`). While they usually match the ones
            created explicitely when defining the LP, a solver like Gurobi has
            been known to introduce additional variables to store constraints of
            the type ``lower_bound <= linear_function <= upper bound``. You
            should be fine if you did not install Gurobi or if you do not use it
            as a solver, but keep an eye on the number of variables in the
            polyhedron, or on the output of :meth:`show`. Just in case.

        EXAMPLES:

        A LP on two variables::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(2*p['x'] + p['y'] <= 1)
            sage: p.add_constraint(3*p['y'] + p['x'] <= 2)
            sage: P = p.polyhedron(); P
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        3-D Polyhedron::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(2*p['x'] + p['y'] + 3*p['z'] <= 1)
            sage: p.add_constraint(2*p['y'] + p['z'] + 3*p['x'] <= 1)
            sage: p.add_constraint(2*p['z'] + p['x'] + 3*p['y'] <= 1)
            sage: P = p.polyhedron(); P
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices

        An empty polyhedron::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(2*p['x'] + p['y'] + 3*p['z'] <= 1)
            sage: p.add_constraint(2*p['y'] + p['z'] + 3*p['x'] <= 1)
            sage: p.add_constraint(2*p['z'] + p['x'] + 3*p['y'] >= 2)
            sage: P = p.polyhedron(); P
            The empty polyhedron in QQ^3

        An unbounded polyhedron::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(2*p['x'] + p['y'] - p['z'] <= 1)
            sage: P = p.polyhedron(); P
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 3 vertices and 3 rays

        A square (see :trac:`14395`) ::

            sage: p = MixedIntegerLinearProgram()
            sage: x,y = p['x'], p['y']
            sage: p.set_min(x,None)
            sage: p.set_min(y,None)
            sage: p.add_constraint( x <= 1 )
            sage: p.add_constraint( x >= -1 )
            sage: p.add_constraint( y <= 1 )
            sage: p.add_constraint( y >= -1 )
            sage: p.polyhedron()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        from copy import copy

        cdef GenericBackend b = self._backend
        cdef int i

        inequalities = []
        equalities = []
        nvar = self.number_of_variables()

        # Constraints
        for lb, (indices, values), ub in self.constraints():

            coeffs = dict(zip(indices, values))

            # Equalities
            if (not lb is None) and lb == ub:
                linear_function = []
                linear_function = [coeffs.get(i,0) for i in range(nvar)]
                linear_function.insert(0,-lb)
                equalities.append(linear_function)
                continue

            # Lower Bound
            if not lb is None:
                linear_function = []
                linear_function = [coeffs.get(i,0) for i in range(nvar)]
                linear_function.insert(0,-lb)
                inequalities.append(linear_function)

            # Upper Bound
            if not ub is None:
                linear_function = []
                linear_function = [-coeffs.get(i,0) for i in range(nvar)]
                linear_function.insert(0,ub)
                inequalities.append(linear_function)

        zero = [0] * nvar

        # Variable bounds
        for 0<= i < nvar:
            lb, ub = b.col_bounds(i)

            # Fixed variable
            if (not lb is None) and lb == ub:
                linear_function = copy(zero)
                linear_function[i] = 1
                linear_function.insert(0,-lb)
                equalities.append(linear_function)
                continue

            # Lower bound
            if not lb is None:
                linear_function = copy(zero)
                linear_function[i] = 1
                linear_function.insert(0,-lb)
                inequalities.append(linear_function)

            # Upper bound
            if not ub is None:
                linear_function = copy(zero)
                linear_function[i] = -1
                linear_function.insert(0,ub)
                inequalities.append(linear_function)

        return Polyhedron(ieqs = inequalities, eqns = equalities)

    def show(self):
        r"""
        Displays the ``MixedIntegerLinearProgram`` in a human-readable
        way.

        EXAMPLES:

        When constraints and variables have names ::

            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: x = p.new_variable(name="Hey")
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2, name="Constraint_1")
            sage: p.show()
            Maximization:
              Hey[1] + Hey[2]
            Constraints:
              Constraint_1: -3.0 Hey[1] + 2.0 Hey[2] <= 2.0
            Variables:
              Hey[1] is a continuous variable (min=0.0, max=+oo)
              Hey[2] is a continuous variable (min=0.0, max=+oo)

        Without any names ::

            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2)
            sage: p.show()
            Maximization:
              x_0 + x_1
            Constraints:
              -3.0 x_0 + 2.0 x_1 <= 2.0
            Variables:
              x_0 is a continuous variable (min=0.0, max=+oo)
              x_1 is a continuous variable (min=0.0, max=+oo)

        With `\QQ` coefficients::

            sage: p = MixedIntegerLinearProgram(solver= 'ppl')
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + 1/2*x[2])
            sage: p.add_constraint(-3/5*x[1] + 2/7*x[2], max=2/5)
            sage: p.show()
            Maximization:
              x_0 + 1/2 x_1
            Constraints:
              constraint_0: -3/5 x_0 + 2/7 x_1 <= 2/5
            Variables:
              x_0 is a continuous variable (min=0, max=+oo)
              x_1 is a continuous variable (min=0, max=+oo)
        """

        cdef int i, j
        cdef GenericBackend b = self._backend

        # inv_variables associates a MIPVariable object to an id
        inv_variables = {}

        for (v, id) in self._variables.iteritems():
            inv_variables[id]=v

        # varid_name associates variables id to names
        varid_name = {}

        for 0<= i < b.ncols():
            s = b.col_name(i)
            varid_name[i] = ("x_"+str(i)) if s == "" else s

        ##### Sense and objective function

        print ("Maximization:" if b.is_maximization() else "Minimization:")
        print " ",

        first = True
        for 0<= i< b.ncols():
            c = b.objective_coefficient(i)
            if c == 0:
                continue

            print (("+ " if (not first and c>0) else "") +
                   ("" if c == 1 else ("- " if c == -1 else str(c)+" "))+varid_name[i]
                   ),
            first = False

        if b.obj_constant_term > self._backend.zero(): print "+", b.obj_constant_term
        elif b.obj_constant_term < self._backend.zero(): print "-", -b.obj_constant_term

        print

        ##### Constraints
        print "Constraints:"

        for 0<= i < b.nrows():

            indices, values = b.row(i)

            lb, ub = b.row_bounds(i)

            print " ",


            # Constraint's name
            if b.row_name(i):
                print b.row_name(i)+":",

            # Lower bound
            if lb is not None:
                print str(lb)+" <=",

            first = True

            for j,c in sorted(zip(indices, values)):

                if c == 0:
                    continue

                print (("+ " if (not first and c>0) else "") +
                       ("" if c == 1 else ("- " if c == -1 else (str(c) + " " if first and c < 0 else ("- " + str(abs(c)) + " " if c < 0 else str(c) + " "))))+varid_name[j]
                       ),
                first = False

            # Upper bound
            print ("<= "+str(ub) if ub is not None else "")


        ##### Variables
        print "Variables:"

        for 0<= i < b.ncols():
            print "  " + varid_name[i] + " is",

            if b.is_variable_integer(i):
                print "an integer variable",
            elif b.is_variable_binary(i):
                print "a boolean variable",
            else:
                print "a continuous variable",

            lb, ub = b.col_bounds(i)

            print "(min=" + ( str(lb) if lb is not None else "-oo" )+",",
            print "max=" + ( str(ub) if ub is not None else "+oo" )+")"


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

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2,name="OneConstraint")
            sage: p.write_mps(os.path.join(SAGE_TMP, "lp_problem.mps"))
            Writing problem data to ...
            17 records were written

        For information about the MPS file format :
        http://en.wikipedia.org/wiki/MPS_%28format%29
        """
        self._backend.write_mps(filename, modern)

    def write_lp(self,filename):
        r"""
        Write the linear program as a LP file.

        This function export the problem as a LP file.

        INPUT:

        - ``filename`` -- The file in which you want the problem
          to be written.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
            sage: x = p.new_variable()
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2)
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))
            Writing problem data to ...
            9 lines were written

        For more information about the LP file format :
        http://lpsolve.sourceforge.net/5.5/lp-format.htm
        """

        self._backend.write_lp(filename)

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

        .. NOTE::

            While a variable may be declared as binary or integer, its value as
            returned by the solver is of type ``float``.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: y = p.new_variable()
            sage: p.set_objective(x[3] + 3*y[2,9] + x[5])
            sage: p.add_constraint(x[3] + y[2,9] + 2*x[5], max=2)
            sage: p.solve()
            6.0

        To return  the optimal value of ``y[2,9]``::

            sage: p.get_values(y[2,9])
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

        TESTS:

        When 'dim' will be removed, also remove from this function the code that
        uses it::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable(dim=2)
            doctest:839: DeprecationWarning: The 'dim' argument will soon disappear. Fortunately variable[1,2] is easier to use than variable[1][2]
            See http://trac.sagemath.org/15489 for details.
            sage: p.add_constraint(b[1][2] +  b[2][3] == 0)
            sage: _ = p.solve()
            sage: _ = p.get_values(b)
        """
        val = []
        for l in lists:
            if isinstance(l, MIPVariable):
                if l.depth() == 1:
                    c = {}
                    for (k,v) in l.items():
                        #c[k] = self._values[v] if self._values.has_key(v) else None
                        c[k] = self._backend.get_variable_value(self._variables[v])
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
            elif l in self._variables:
                #val.append(self._values[l])
                val.append(self._backend.get_variable_value(self._variables[l]))

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
            sage: round(p.solve(),5)
            6.66667
            sage: p.set_objective(None)
            sage: _ = p.solve()
        """
        cdef list values = []

        # If the objective is None, or a constant, we want to remember
        # that the objective function has been defined ( the user did not
        # forget it ). In some LP problems, you just want a feasible solution
        # and do not care about any function being optimal.

        cdef int i

        if obj is not None:
            f = obj.dict()
        else:
            f = {-1 : 0}

        d = f.pop(-1,self._backend.zero())

        for i in range(self._backend.ncols()):
            values.append(f.get(i,self._backend.zero()))

        self._backend.set_objective(values, d)

    def add_constraint(self, linear_function, max=None, min=None, name=None):
        r"""
        Adds a constraint to the ``MixedIntegerLinearProgram``.

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
          numerical value.
        - ``name`` -- A name for the constraint.

        To set a lower and/or upper bound on the variables use the methods
        ``set_min`` and/or ``set_max`` of ``MixedIntegerLinearProgram``.

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
            sage: round(p.solve(),6)
            6.666667

        There are two different ways to add the constraint
        ``x[5] + 3*x[7] <= x[6] + 3`` to a ``MixedIntegerLinearProgram``.

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
            sage: round(p.solve(), 5)
            6.66667

        TESTS:

        Complex constraints::

            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: b = p.new_variable()
            sage: p.add_constraint( b[8] - b[15] <= 3*b[8] + 9)
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              -2.0 x_0 - x_1 <= 9.0
            Variables:
              x_0 is a continuous variable (min=0.0, max=+oo)
              x_1 is a continuous variable (min=0.0, max=+oo)

        Empty constraint::

            sage: p=MixedIntegerLinearProgram()
            sage: p.add_constraint(sum([]),min=2)

        Min/Max are numerical ::

            sage: v = p.new_variable()
            sage: p.add_constraint(v[3] + v[5], min = v[6])
            Traceback (most recent call last):
            ...
            ValueError: min and max arguments are required to be numerical
            sage: p.add_constraint(v[3] + v[5], max = v[6])
            Traceback (most recent call last):
            ...
            ValueError: min and max arguments are required to be numerical

        Do not add redundant elements (notice only one copy of each constraint is added)::

            sage: lp = MixedIntegerLinearProgram(solver = "GLPK", check_redundant=True)
            sage: for each in xrange(10): lp.add_constraint(lp[0]-lp[1],min=1)
            sage: lp.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              1.0 <= x_0 - x_1
            Variables:
              x_0 is a continuous variable (min=0.0, max=+oo)
              x_1 is a continuous variable (min=0.0, max=+oo)

        We check for constant multiples of constraints as well::

            sage: for each in xrange(10): lp.add_constraint(2*lp[0]-2*lp[1],min=2)
            sage: lp.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              1.0 <= x_0 - x_1
            Variables:
              x_0 is a continuous variable (min=0.0, max=+oo)
              x_1 is a continuous variable (min=0.0, max=+oo)

        But if the constant multiple is negative, we should add it anyway (once)::

              sage: for each in xrange(10): lp.add_constraint(-2*lp[0]+2*lp[1],min=-2)
              sage: lp.show()
              Maximization:
              <BLANKLINE>
              Constraints:
                1.0 <= x_0 - x_1
                x_0 - x_1 <= 1.0
              Variables:
                x_0 is a continuous variable (min=0.0, max=+oo)
                x_1 is a continuous variable (min=0.0, max=+oo)

        TESTS:

        Catch ``True`` / ``False`` as INPUT (:trac:`13646`)::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: p.add_constraint(True)
            Traceback (most recent call last):
            ...
            ValueError: argument must be a linear function or constraint, got True
        """
        if linear_function is 0:
            return

        # Raising an exception when min/max are not as expected
        from sage.rings.all import RR
        if ((min is not None and min not in RR)
            or (max is not None and max not in RR)):
            raise ValueError("min and max arguments are required to be numerical")

        if is_LinearFunction(linear_function):
            f = linear_function.dict()
            constant_coefficient = f.get(-1,0)

            # We do not want to ignore the constant coefficient
            max = (max - constant_coefficient) if max is not None else None
            min = (min - constant_coefficient) if min is not None else None

            indices = []
            values = []

            if self._check_redundant:
              b = self._backend
              from __builtin__ import min as min_function
              i = min_function([v for (v,coeff) in f.iteritems() if coeff != 0])
              c = f[i]
              C = [(v,coeff/c) for (v,coeff) in f.iteritems() if v != -1]
              if c > 0:
                min = min/c if min is not None else None
                max = max/c if max is not None else None
              else:
                tempmin = max/c if max is not None else None
                tempmax = min/c if min is not None else None
                min, max = tempmin, tempmax
              if (tuple(C),min,max) in self._constraints:
                return None
              else:
                self._constraints.append((tuple(C),min,max))
            else:
              C = [(v,coeff) for (v,coeff) in f.iteritems() if v != -1]

            if min is None and max is None:
                raise ValueError("Both max and min are set to None ? Weird!")

            self._backend.add_linear_constraint(C, min, max, name)

        elif is_LinearConstraint(linear_function):
            constraint = linear_function
            for lhs, rhs in constraint.equations():
                self.add_constraint(lhs-rhs, min=0, max=0, name=name)
            for lhs, rhs in constraint.inequalities():
                self.add_constraint(lhs-rhs, max=0, name=name)
        else:
            raise ValueError('argument must be a linear function or constraint, got '+str(linear_function))

    def remove_constraint(self, int i):
        r"""
        Removes a constraint from self.

        INPUT:

        - ``i`` -- Index of the constraint to remove.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x, y = p[0], p[1]
            sage: p.add_constraint(x + y, max = 10)
            sage: p.add_constraint(x - y, max = 0)
            sage: p.add_constraint(x, max = 4)
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              x_0 + x_1 <= 10.0
              x_0 - x_1 <= 0.0
              x_0 <= 4.0
            ...
            sage: p.remove_constraint(1)
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              x_0 + x_1 <= 10.0
              x_0 <= 4.0
            ...
            sage: p.number_of_constraints()
            2
        """
        if self._check_redundant: self._constraints.pop(i)
        self._backend.remove_constraint(i)

    def remove_constraints(self, constraints):
        r"""
        Remove several constraints.

        INPUT:

        - ``constraints`` -- an iterable containing the indices of the rows to remove.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x, y = p[0], p[1]
            sage: p.add_constraint(x + y, max = 10)
            sage: p.add_constraint(x - y, max = 0)
            sage: p.add_constraint(x, max = 4)
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              x_0 + x_1 <= 10.0
              x_0 - x_1 <= 0.0
              x_0 <= 4.0
            ...
            sage: p.remove_constraints([0, 1])
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              x_0 <= 4.0
            ...
            sage: p.number_of_constraints()
            1

        When checking for redundant constraints, make sure you remove only
        the constraints that were actually added. Problems could arise if
        you have a function that builds lps non-interactively, but it fails
        to check whether adding a constraint actually increases the number of
        constraints. The function might later try to remove constraints that
        are not actually there::

            sage: p = MixedIntegerLinearProgram(check_redundant=True)
            sage: x, y = p[0], p[1]
            sage: p.add_constraint(x + y, max = 10)
            sage: for each in xrange(10): p.add_constraint(x - y, max = 10)
            sage: p.add_constraint(x, max = 4)
            sage: p.number_of_constraints()
            3
            sage: p.remove_constraints(range(1,9))
            Traceback (most recent call last):
            ...
            IndexError: pop index out of range
            sage: p.remove_constraint(1)
            sage: p.number_of_constraints()
            2

        We should now be able to add the old constraint back in::

            sage: for each in xrange(10): p.add_constraint(x - y, max = 10)
            sage: p.number_of_constraints()
            3
        """
        if self._check_redundant:
          for i in sorted(constraints,reverse=True):
            self._constraints.pop(i)
        self._backend.remove_constraints(constraints)

    def set_binary(self, ee):
        r"""
        Sets a variable or a ``MIPVariable`` as binary.

        INPUT:

        - ``ee`` -- An instance of ``MIPVariable`` or one of
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
        variables as nonnegative while keeping the others as they are::

            sage: p.set_nonnegative(x[3])

        TESTS:

        When 'dim' will be removed, also remove all the ``is_*`` and ``set_*``
        functions the code that uses it::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable(dim=2)
            sage: p.add_constraint(b[1][2] +  b[2][3] == 0)
            sage: p.set_binary(b)
        """
        cdef MIPVariable e
        e = <MIPVariable> ee

        if isinstance(e, MIPVariable):
            e._vtype = self.__BINARY
            if e.depth() == 1:
                for v in e.values():
                    self._backend.set_variable_type(self._variables[v],self.__BINARY)
            else:
                for v in e.keys():
                    self.set_binary(e[v])
        elif e in self._variables:
            self._backend.set_variable_type(self._variables[e],self.__BINARY)
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

        return self._backend.is_variable_binary(self._variables[e])

    def set_integer(self, ee):
        r"""
        Sets a variable or a ``MIPVariable`` as integer.

        INPUT:

        - ``ee`` -- An instance of ``MIPVariable`` or one of
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
        variables as nonnegative while keeping the others as they are::

            sage: p.set_nonnegative(x[3])
        """
        cdef MIPVariable e
        e = <MIPVariable> ee

        if isinstance(e, MIPVariable):
            e._vtype = self.__INTEGER
            if e.depth() == 1:
                for v in e.values():
                    self._backend.set_variable_type(self._variables[v],self.__INTEGER)
            else:
                for v in e.keys():
                    self.set_integer(e[v])
        elif e in self._variables:
            self._backend.set_variable_type(self._variables[e],self.__INTEGER)
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
        return self._backend.is_variable_integer(self._variables[e])

    set_real = deprecated_function_alias(15521, set_nonnegative)
    def set_nonnegative(self,ee):
        r"""
        Sets a variable or a ``MIPVariable`` as nonnegative.

        INPUT:

        - ``ee`` -- An instance of ``MIPVariable`` or one of
          its elements.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()

        With the following instruction, all the variables
        from x will be nonnegative::

            sage: p.set_nonnegative(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)

         It is still possible, though, to set one of these
         variables as binary while keeping the others as they are::

            sage: p.set_binary(x[3])

        TESTS:

        :trac:`15521`::

            sage: p.set_real(x[3])
            doctest:1: DeprecationWarning: set_real is deprecated. Please use set_nonnegative instead.
            See http://trac.sagemath.org/15521 for details.
        """

        cdef MIPVariable e
        e = <MIPVariable> ee

        if isinstance(e, MIPVariable):
            e._vtype = self.__REAL
            if e.depth() == 1:
                for v in e.values():
                    self._backend.set_variable_type(self._variables[v],self.__REAL)
                    self._backend.variable_lower_bound(self._variables[v], 0)
            else:
                for v in e.keys():
                    self.set_real(e[v])
        elif e in self._variables:
            self._backend.set_variable_type(self._variables[e],self.__REAL)
            self._backend.variable_lower_bound(self._variables[e], 0)
        else:
            raise ValueError("e must be an instance of MIPVariable or one of its elements.")

    is_real = deprecated_function_alias(15521,is_nonnegative)
    def is_nonnegative(self, e):
        r"""
        Tests whether the variable is nonnegative.

        INPUT:

        - ``e`` -- A variable (not a ``MIPVariable``, but one of its elements.)

        OUTPUT:

        ``True`` if the variable is nonnegative; ``False`` otherwise.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[1])
            sage: p.is_nonnegative(v[1])
            True
            sage: p.set_binary(v[1])
            sage: p.is_nonnegative(v[1])
            False
            sage: p.set_nonnegative(v[1])
            sage: p.is_nonnegative(v[1])
            True
        """
        return (self._backend.is_variable_continuous(self._variables[e]) and
                self._backend.variable_lower_bound(self._variables[e]) == 0)

    def solve(self, log=None, objective_only=False):
        r"""
        Solves the ``MixedIntegerLinearProgram``.

        INPUT:

        - ``log`` -- integer (default: ``None``) The verbosity level. Indicates
          whether progress should be printed during computation. The solver is
          initialized to report no progress.

        - ``objective_only`` -- Boolean variable.

          - When set to ``True``, only the objective function is returned.
          - When set to ``False`` (default), the optimal numerical values
            are stored (takes computational time).

        OUTPUT:

        The optimal value taken by the objective function.

        .. WARNING::

            By default, all variables of a LP are assumed to be
            non-negative. See :meth:`set_min` to change it.

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
            sage: round(p.solve(),6)
            6.666667
            sage: x = p.get_values(x)
            sage: round(x[1],6)
            0.0
            sage: round(x[2],6)
            1.333333

         Computation of a maximum stable set in Petersen's graph::

            sage: g = graphs.PetersenGraph()
            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: b = p.new_variable()
            sage: p.set_objective(sum([b[v] for v in g]))
            sage: for (u,v) in g.edges(labels=None):
            ...       p.add_constraint(b[u] + b[v], max=1)
            sage: p.set_binary(b)
            sage: p.solve(objective_only=True)
            4.0

        Constraints in the objective function are respected::

            sage: p = MixedIntegerLinearProgram()
            sage: x, y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.set_integer(x); p.set_integer(y)
            sage: p.solve()
            9.0
        """
        if log is not None: self._backend.set_verbosity(log)

        self._backend.solve()

        return self._backend.get_objective_value()

    def set_min(self, v, min):
        r"""
        Sets the minimum value of a variable.

        .. WARNING::

            By default, all variables are defined to be non-negative.

        INPUT:

        - ``v`` -- a variable (not a ``MIPVariable``, but one of its
          elements).
        - ``min`` -- the minimum value the variable can take.
          When ``min=None``, the variable has no lower bound.

        .. SEEALSO::

            - :meth:`get_min` -- get the minimum value of a variable.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: p.set_objective(v[1])
            sage: p.get_min(v[1])
            0.0
            sage: p.set_min(v[1],6)
            sage: p.get_min(v[1])
            6.0
            sage: p.set_min(v[1], None)
            sage: p.get_min(v[1])

        """
        self._backend.variable_lower_bound(self._variables[v], min)

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
        self._backend.variable_upper_bound(self._variables[v], max)

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
            sage: p.set_min(v[1], None)
            sage: p.get_min(v[1])
        """
        return self._backend.variable_lower_bound(self._variables[v])

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
        return self._backend.variable_upper_bound(self._variables[v])

    def solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        The solver parameters are by essence solver-specific, which
        means their meaning heavily depends on the solver used.

        (If you do not know which solver you are using, then you use
        use GLPK).

        Aliases:

        Very common parameters have aliases making them
        solver-independent. For example, the following::

            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: p.solver_parameter("timelimit", 60)

        Sets the solver to stop its computations after 60 seconds, and
        works with GLPK, CPLEX and Gurobi.

            - ``"timelimit"`` -- defines the maximum time spent on a
              computation. Measured in seconds.

        Solver-specific parameters:

            - GLPK : We have implemented very close to comprehensive coverage of
              the GLPK solver parameters for the simplex and integer
              optimization methods. For details, see the documentation of
              :meth:`GLPKBackend.solver_parameter
              <sage.numerical.backends.glpk_backend.GLPKBackend.solver_parameter>`.

            - CPLEX's parameters are identified by a string. Their
              list is available `on ILOG's website
              <http://publib.boulder.ibm.com/infocenter/odmeinfo/v3r4/index.jsp?topic=/ilog.odms.ide.odme.help/Content/Optimization/Documentation/ODME/_pubskel/ODME_pubskels/startall_ODME34_Eclipse1590.html>`_.

              The command ::

                  sage: p = MixedIntegerLinearProgram(solver = "CPLEX") # optional - CPLEX
                  sage: p.solver_parameter("CPX_PARAM_TILIM", 60)       # optional - CPLEX

              works as intended.

            - Gurobi's parameters should all be available through this
              method. Their list is available on Gurobi's website
              `<http://www.gurobi.com/documentation/5.5/reference-manual/node798>`_.

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: p.solver_parameter("timelimit", 60)
            sage: p.solver_parameter("timelimit")
            60.0
        """
        if value is None:
            return self._backend.solver_parameter(name)
        else:
            self._backend.solver_parameter(name, value)

    cpdef sum(self, L):
        r"""
        Efficiently computes the sum of a sequence of
        :class:`~sage.numerical.linear_functions.LinearFunction` elements

        INPUT:

        - ``mip`` -- the :class:`MixedIntegerLinearProgram` parent.

        - ``L`` -- list of
          :class:`~sage.numerical.linear_functions.LinearFunction` instances.

        .. NOTE::

            The use of the regular ``sum`` function is not recommended
            as it is much less efficient than this one

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()

        The following command::

            sage: s = p.sum([v[i] for i in xrange(90)])

        is much more efficient than::

            sage: s = sum([v[i] for i in xrange(90)])
        """
        d = {}
        for v in L:
            for id,coeff  in v.iteritems():
                d[id] = coeff + d.get(id,0)
        return self.linear_functions_parent()(d)

    def get_backend(self):
        r"""
        Returns the backend instance used.

        This might be useful when acces to additional functions provided by
        the backend is needed.

        EXAMPLE:

        This example uses the simplex algorthm and prints information::

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
            sage: x, y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: b = p.get_backend()
            sage: b.solver_parameter("simplex_or_intopt", "simplex_only")
            sage: b.solver_parameter("verbosity_simplex", "GLP_MSG_ALL")
            sage: p.solve()  # tol 0.00001
            GLPK Simplex Optimizer, v4.44
            2 rows, 2 columns, 4 non-zeros
            *     0: obj =   7.000000000e+00  infeas =  0.000e+00 (0)
            *     2: obj =   9.400000000e+00  infeas =  0.000e+00 (0)
            OPTIMAL SOLUTION FOUND
            9.4
        """
        return self._backend

class MIPSolverException(RuntimeError):
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

            sage: p=MixedIntegerLinearProgram(solver="GLPK")
            sage: v=p.new_variable()
            sage: p.add_constraint(v[0],max=5.5)
            sage: p.add_constraint(v[0],min=7.6)
            sage: p.set_objective(v[0])

        Tests of GLPK's Exceptions::

            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: 'GLPK : Solution is undefined'

        No integer solution::

            sage: p=MixedIntegerLinearProgram(solver="GLPK")
            sage: v=p.new_variable()
            sage: p.add_constraint(v[0],max=5.6)
            sage: p.add_constraint(v[0],min=5.2)
            sage: p.set_objective(v[0])
            sage: p.set_integer(v)

        Tests of GLPK's Exceptions::

            sage: p.solve()
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

cdef class MIPVariable(SageObject):
    r"""
    ``MIPVariable`` is a variable used by the class
    ``MixedIntegerLinearProgram``.
    """

    def __cinit__(self, p, vtype, dim=1, name=""):
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

        self._hasname = (len(name) >0)

        if dim > 1:
            from sage.misc.superseded import deprecation
            deprecation(15489, "The 'dim' argument will soon disappear. "+
                        "Fortunately variable[1,2] is easier to use than "+
                        "variable[1][2]")

        # create a temporary char *
        cdef char *name_c = name
        # and copy it over
        self._name = <char*>sage_malloc(len(name)+1)
        strcpy(self._name, name_c)

    def __dealloc__(self):
        if self._name:
            sage_free(self._name)

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
        cdef MIPVariable s = self

        cdef int j

        if i in self._dict:
            return self._dict[i]
        elif self._dim == 1:
            zero = self._p._backend.zero()
            j = self._p._backend.add_variable(zero , None, False, True, False, zero,
                                              (str(self._name) + "[" + str(i) + "]")
                                               if self._hasname else None)

            v = self._p.linear_function({j : 1})
            self._p._variables[v] = j
            self._p._backend.set_variable_type(j,self._vtype)
            self._dict[i] = v

            return v

        else:
            self._dict[i] = MIPVariable(
                self._p,
                self._vtype,
                dim=self._dim-1,
                name = ("" if not self._hasname
                        else (str(self._name) + "[" + str(i) + "]")))

            return self._dict[i]

    def _repr_(self):
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

def Sum(x):
    """
    Only for legacy support, use :meth:`MixedIntegerLinearProgram.sum` instead.

    EXAMPLES::

        sage: from sage.numerical.mip import Sum
        sage: Sum([])
        doctest:...: DeprecationWarning: use MixedIntegerLinearProgram.sum() instead
        See http://trac.sagemath.org/13646 for details.

        sage: p = MixedIntegerLinearProgram()
        sage: x = p.new_variable()
        sage: Sum([ x[0]+x[1], x[1]+x[2], x[2]+x[3] ])   # deprecation is only shown once
        x_0 + 2*x_1 + 2*x_2 + x_3
    """
    from sage.misc.superseded import deprecation
    deprecation(13646, 'use MixedIntegerLinearProgram.sum() instead')
    if not x:
        return None
    parent = x[0].parent()
    return parent.sum(x)
