r"""
Mixed Integer Linear Programming

This module implements classes and methods for the efficient solving of Linear
Programs (`LP <http://en.wikipedia.org/wiki/Linear_programming>`_) and Mixed
Integer Linear Programs (`MILP
<http://en.wikipedia.org/wiki/Mixed_integer_linear_programming>`_).

*Do you want to understand how the simplex method works?* See the
:mod:`~sage.numerical.interactive_simplex_method` module (educational purposes
only)

Definition
----------

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

Example
-------

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
  #. Create a dictionary ``w`` of integer variables ``w`` via ``w =
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
    sage: w = p.new_variable(integer=True, nonnegative=True)
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
    sage: for i, v in p.get_values(w).iteritems():
    ....:     print 'w_%s = %s' % (i, int(round(v)))
    w_0 = 15
    w_1 = 10
    w_2 = 3
    w_3 = 2

Different backends compute with different base fields, for example::

    sage: p = MixedIntegerLinearProgram(solver='GLPK')
    sage: p.base_ring()
    Real Double Field
    sage: x = p.new_variable(real=True, nonnegative=True)
    sage: 0.5 + 3/2*x[1]
    0.5 + 1.5*x_0

    sage: p = MixedIntegerLinearProgram(solver='ppl')
    sage: p.base_ring()
    Rational Field
    sage: x = p.new_variable(nonnegative=True)
    sage: 0.5 + 3/2*x[1]
    1/2 + 3/2*x_0


Linear Variables and Expressions
--------------------------------

The underlying linear programming backends always work with matrices
where each column corresponds to a linear variable. These variables
can be accessed using the :meth:`MixedIntegerLinearProgram.gen` method
or by calling with a dictionary variable index to coefficient::

    sage: mip = MixedIntegerLinearProgram()
    sage: 5 + mip.gen(0) + 2*mip.gen(1)
    5 + x_0 + 2*x_1
    sage: mip({-1:5, 0:1, 1:2})
    5 + x_0 + 2*x_1

However, this alone is often not convenient to construct a linear
program. To make your code more readable, you can construct
:class:`MIPVariable` objects that can be arbitrarily named and
indexed. Internally, this is then translated back to the `x_i`
variables. For example::

    sage: mip.<a,b> = MixedIntegerLinearProgram()
    sage: a
    MIPVariable of dimension 1
    sage: 5 + a[1] + 2*b[3]
    5 + x_0 + 2*x_1

Indices can be any object, not necessarily integers. Multi-indices are
also allowed::

    sage: a[4, 'string', QQ]
    x_2
    sage: a[4, 'string', QQ] - 7*b[2]
    x_2 - 7*x_3
    sage: mip.show()
    Maximization:
    <BLANKLINE>
    Constraints:
    Variables:
      a[1] = x_0 is a continuous variable (min=-oo, max=+oo)
      b[3] = x_1 is a continuous variable (min=-oo, max=+oo)
      a[(4, 'string', Rational Field)] = x_2 is a continuous variable (min=-oo, max=+oo)
      b[2] = x_3 is a continuous variable (min=-oo, max=+oo)


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
    :meth:`~MixedIntegerLinearProgram.best_known_objective_bound`| Return the value of the currently best known bound
    :meth:`~MixedIntegerLinearProgram.constraints`               | Returns a list of constraints, as 3-tuples
    :meth:`~MixedIntegerLinearProgram.get_backend`               | Returns the backend instance used
    :meth:`~MixedIntegerLinearProgram.get_max`                   | Returns the maximum value of a variable
    :meth:`~MixedIntegerLinearProgram.get_min`                   | Returns the minimum value of a variable
    :meth:`~MixedIntegerLinearProgram.get_objective_value`       | Return the value of the objective function
    :meth:`~MixedIntegerLinearProgram.get_relative_objective_gap`| Return the relative objective gap of the best known solution
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

AUTHORS:

- Risan (2012/02): added extension for exact computation
"""

#*****************************************************************************
#       Copyright (C) 2012 Nathann Cohen <nathann.cohen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from copy import copy
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.misc.cachefunc import cached_method
from sage.misc.superseded import deprecation

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

      - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT <http://www.cvxopt.org/>`_
          web site.

      - PPL (``solver="PPL"``). See the `PPL <http://bugseng.com/products/ppl>`_
          web site.

      - If ``solver=None`` (default), the default solver is used (see
        :func:`default_mip_solver`)

      - ``solver`` can also be a callable,
        see :func:`sage.numerical.backends.generic_backend.get_solver` for
        examples.

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
         sage: b = p.new_variable(binary=True)
         sage: p.set_objective(sum([b[v] for v in g]))
         sage: for (u,v) in g.edges(labels=None):
         ....:     p.add_constraint(b[u] + b[v], max=1)
         sage: p.solve(objective_only=True)
         4.0

    TESTS:

    Check that :trac:`16497` is fixed::

        sage: from sage.numerical.mip import MixedIntegerLinearProgram
        sage: for type in ["binary", "integer"]:
        ....:     k = 3
        ....:     items = [1/5, 1/3, 2/3, 3/4, 5/7]
        ....:     maximum=1
        ....:     p=MixedIntegerLinearProgram()
        ....:     box=p.new_variable(nonnegative=True, **{type:True})
        ....:     for b in range(k):
        ....:          p.add_constraint(p.sum([items[i]*box[i,b] for i in range(len(items))]) <= maximum)
        ....:     for i in range(len(items)):
        ....:         p.add_constraint(p.sum([box[i,b] for b in range(k)]) == 1)
        ....:     p.set_objective(None)
        ....:     _ = p.solve()
        ....:     box=p.get_values(box)
        ....:     print(all(v in ZZ for v in box.values()))
        True
        True
    """

    def __init__(self, solver=None, maximization=True,
                 constraint_generation=False, check_redundant=False,
                 names=tuple(), base_ring=None):
        r"""
        Constructor for the ``MixedIntegerLinearProgram`` class.

        INPUT:

        - ``solver`` -- the following solvers should be available through this class:

          - GLPK (``solver="GLPK"``). See the `GLPK
            <http://www.gnu.org/software/glpk/>`_ web site.

          - GLPK's implementation of an exact rational simplex
            method (``solver="GLPK/exact"``).

          - COIN Branch and Cut (``solver="Coin"``). See the `COIN-OR
            <http://www.coin-or.org>`_ web site.

          - CPLEX (``solver="CPLEX"``). See the `CPLEX
            <http://www.ilog.com/products/cplex/>`_ web site.  An interface to
            CPLEX is not yet implemented.

          - Gurobi (``solver="Gurobi"``). See the `Gurobi
            <http://www.gurobi.com/>`_ web site.

          - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT <http://www.cvxopt.org/>`_
              web site.

          - PPL (``solver="PPL"``). See the `PPL
            <http://bugseng.com/products/ppl>`_ web site.

          - If ``solver=None`` (default), the default solver is used, see
            :func:`default_mip_solver`.

          - ``solver`` can also be a callable,
            see :func:`sage.numerical.backends.generic_backend.get_solver` for
            examples.

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

        - ``names`` -- list/tuple/iterable of string. Default names of
          the MIP variables. Used to enable the ``MIP.<x> =
          MixedIntegerLinearProgram()`` syntax.

        .. SEEALSO::

        - :meth:`default_mip_solver` -- Returns/Sets the default MIP solver.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True)

        TESTS:

        Checks that the objects are deallocated without invoking the cyclic garbage
        collector (cf. :trac:`12616`)::

            sage: del p
            sage: def just_create_variables():
            ....:     p = MixedIntegerLinearProgram()
            ....:     b = p.new_variable(nonnegative=True)
            ....:     p.add_constraint(b[3]+b[6] <= 2)
            ....:     p.solve()
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
        self._first_variable_names = list(names)
        from sage.numerical.backends.generic_backend import get_solver
        self._backend = get_solver(solver=solver,
                                   constraint_generation=constraint_generation,
                                   base_ring=base_ring)
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
             sage: v = p.new_variable(nonnegative=True)
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
            sage: v = p.new_variable(nonnegative=True)
            sage: p.add_constraint(v[0] + v[1], max = 10)
            sage: q = copy(p)
            sage: q.number_of_constraints()
            1
        """
        def copying_solver(**kwdargs):
            return (<GenericBackend> self._backend).copy()

        cdef MixedIntegerLinearProgram p = \
            MixedIntegerLinearProgram(solver=copying_solver)
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

        return p

    def __deepcopy__(self, memo={}):
        """
        Return a deep copy of ``self``.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: p.add_constraint(b[1] + b[2] <= 6)
            sage: p.set_objective(b[1] + b[2])
            sage: cp = deepcopy(p)
            sage: cp.solve()
            6.0

        TEST:

        Test that `deepcopy` makes actual copies but preserves identities::

            sage: mip = MixedIntegerLinearProgram()
            sage: ll = [mip, mip]
            sage: dcll=deepcopy(ll)
            sage: ll[0] is dcll[0]
            False
            sage: dcll[0] is dcll[1]
            True

        """
        return self.__copy__()

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
            sage: from sage.rings.all import AA
            sage: p = MixedIntegerLinearProgram(solver='InteractiveLP', base_ring=AA)
            sage: p.base_ring()
            Algebraic Real Field
            sage: d = polytopes.dodecahedron()
            sage: p = MixedIntegerLinearProgram(base_ring=d.base_ring())
            sage: p.base_ring()
            Number Field in sqrt5 with defining polynomial x^2 - 5
        """
        return self._backend.base_ring()

    def set_problem_name(self,name):
        r"""
        Sets the name of the ``MixedIntegerLinearProgram``.

        INPUT:

        - ``name`` -- A string representing the name of the
          ``MixedIntegerLinearProgram``.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: p.set_problem_name("Test program")
            sage: p
            Mixed Integer Program "Test program" ( maximization, 0 variables, 0 constraints )
        """
        self._backend.problem_name(name)

    def new_variable(self, real=False, binary=False, integer=False, nonnegative=False, name=""):
        r"""
        Return a new MIPVariable

        A new variable ``x`` is defined by::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable(nonnegative=True)

        It behaves exactly as a usual dictionary would. It can use any key
        argument you may like, as ``x[5]`` or ``x["b"]``, and has methods
        ``items()`` and ``keys()``.

        .. SEEALSO::

            - :meth:`set_min`, :meth:`get_min` -- set/get the lower bound of a
              variable. Note that by default, all variables are non-negative.

            - :meth:`set_max`, :meth:`get_max` -- set/get the upper bound of a
              variable.

        INPUT:

        - ``binary, integer, real`` -- boolean. Set one of these
          arguments to ``True`` to ensure that the variable gets the
          corresponding type.

        - ``nonnegative`` -- boolean, default ``False``. Whether the
          variable should be assumed to be nonnegative. Rather useless
          for the binary type.

        - ``name`` -- string. Associates a name to the variable. This
          is only useful when exporting the linear program to a file
          using ``write_mps`` or ``write_lp``, and has no other
          effect.

        OUTPUT:

        A new instance of :class:`MIPVariable` associated to the
        current :class:`MixedIntegerLinearProgram`.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable(); x
            MIPVariable of dimension 1
            sage: x0 = x[0]; x0
            x_0

        By default, variables are unbounded::

            sage: print p.get_min(x0)
            None
            sage: print p.get_max(x0)
            None

         To define two dictionaries of variables, the first being
         of real type, and the second of integer type ::

            sage: x = p.new_variable(real=True, nonnegative=True)
            sage: y = p.new_variable(integer=True, nonnegative=True)
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

        Unbounded variables::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable(real=True)
            sage: y = p.new_variable(integer=True)
            sage: p.add_constraint(x[0]+x[3] <= 8)
            sage: p.add_constraint(y[0] >= y[1])
            sage: p.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              x_0 + x_1 <= 8.0
              - x_2 + x_3 <= 0.0
            Variables:
              x_0 is a continuous variable (min=-oo, max=+oo)
              x_1 is a continuous variable (min=-oo, max=+oo)
              x_2 is an integer variable (min=-oo, max=+oo)
              x_3 is an integer variable (min=-oo, max=+oo)

        On the Sage command line, generator syntax is accepted as a
        shorthand for generating new variables with default settings::

            sage: mip.<x, y, z> = MixedIntegerLinearProgram()
            sage: mip.add_constraint(x[0] + y[1] + z[2] <= 10) 
            sage: mip.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              x[0] + y[1] + z[2] <= 10.0
            Variables:
              x[0] = x_0 is a continuous variable (min=-oo, max=+oo)
              y[1] = x_1 is a continuous variable (min=-oo, max=+oo)
              z[2] = x_2 is a continuous variable (min=-oo, max=+oo)
        """
        if sum([real, binary, integer]) > 1:
            raise ValueError("Exactly one of the available types has to be True")

        if binary:
            vtype = self.__BINARY
        elif integer:
            vtype = self.__INTEGER
        else:
            vtype = self.__REAL

        if not name and self._first_variable_names:
            name = self._first_variable_names.pop(0)

        return mip_variable_parent(self,
                      vtype,
                      name=name,
                      lower_bound=0 if (nonnegative or binary) else None,
                      upper_bound=1 if binary else None)

    def _first_ngens(self, n):
        """
        Construct the first `n` MIPVariables.

        This method is used for the generater syntax (see below). You
        probably shouldn't use it for anything else.

        INPUT:

        - ``n`` -- integer. The number of variables to construct.

        OUTPUT:

        A tuple of not necessarily positive :class:`MIPVariable`
        instances.

        EXAMPLES::

            sage: mip.<a,b> = MixedIntegerLinearProgram()
            sage: a[0] + b[2]
            x_0 + x_1
            sage: mip.show()
            Maximization:
            <BLANKLINE>
            Constraints:
            Variables:
              a[0] = x_0 is a continuous variable (min=-oo, max=+oo)
              b[2] = x_1 is a continuous variable (min=-oo, max=+oo)
        """
        return tuple(self.new_variable() for i in range(n))

    def gen(self, i):
        """
        Return the linear variable `x_i`.
        
        EXAMPLE::

            sage: mip = MixedIntegerLinearProgram()
            sage: mip.gen(0)
            x_0
            sage: [mip.gen(i) for i in range(10)]
            [x_0, x_1, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9]
        """
        return self.linear_functions_parent().gen(i)

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
            ....:    d = dict(zip(ind, coef))
            ....:    ind.sort()
            ....:    return (lb, (ind, [d[i] for i in ind]), ub)

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
          raise ValueError("constraints() requires a list of integers, though it will accommodate None or an integer.")

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
            created explicitly when defining the LP, a solver like Gurobi has
            been known to introduce additional variables to store constraints of
            the type ``lower_bound <= linear_function <= upper bound``. You
            should be fine if you did not install Gurobi or if you do not use it
            as a solver, but keep an eye on the number of variables in the
            polyhedron, or on the output of :meth:`show`. Just in case.

        .. SEEALSO::

            :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.to_linear_program`
            -- return the :class:`MixedIntegerLinearProgram` object associated
            with a :func:`Polyhedron` object.

        EXAMPLES:

        A LP on two variables::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(0 <= 2*p['x'] + p['y'] <= 1)
            sage: p.add_constraint(0 <= 3*p['y'] + p['x'] <= 2)
            sage: P = p.polyhedron(); P
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        3-D Polyhedron::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(0 <= 2*p['x'] + p['y'] + 3*p['z'] <= 1)
            sage: p.add_constraint(0 <= 2*p['y'] + p['z'] + 3*p['x'] <= 1)
            sage: p.add_constraint(0 <= 2*p['z'] + p['x'] + 3*p['y'] <= 1)
            sage: P = p.polyhedron(); P
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices

        An empty polyhedron::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.add_constraint(2*v['x'] + v['y'] + 3*v['z'] <= 1)
            sage: p.add_constraint(2*v['y'] + v['z'] + 3*v['x'] <= 1)
            sage: p.add_constraint(2*v['z'] + v['x'] + 3*v['y'] >= 2)
            sage: P = p.polyhedron(); P
            The empty polyhedron in QQ^3

        An unbounded polyhedron::

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(2*p['x'] + p['y'] - p['z'] <= 1)
            sage: P = p.polyhedron(); P
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex, 1 ray, 2 lines

        A square (see :trac:`14395`) ::

            sage: p = MixedIntegerLinearProgram()
            sage: x,y = p['x'], p['y']
            sage: p.add_constraint( x <= 1 )
            sage: p.add_constraint( x >= -1 )
            sage: p.add_constraint( y <= 1 )
            sage: p.add_constraint( y >= -1 )
            sage: p.polyhedron()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
        """
        from sage.geometry.polyhedron.constructor import Polyhedron
        cdef GenericBackend b = self._backend
        cdef int i

        # Constraints
        inequalities = []
        equalities = []
        nvar = self.number_of_variables()
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

        # Variable bounds
        zero = [0] * nvar
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

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
            sage: x = p.new_variable(name="Hey")
            sage: p.set_objective(x[1] + x[2])
            sage: p.add_constraint(-3*x[1] + 2*x[2], max=2, name="Constraint_1")
            sage: p.show()
            Maximization:
              Hey[1] + Hey[2]
            Constraints:
              Constraint_1: -3.0 Hey[1] + 2.0 Hey[2] <= 2.0
            Variables:
              Hey[1] = x_0 is a continuous variable (min=-oo, max=+oo)
              Hey[2] = x_1 is a continuous variable (min=-oo, max=+oo)

        Without any names ::

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
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

            sage: p = MixedIntegerLinearProgram(solver='ppl')
            sage: x = p.new_variable(nonnegative=True)
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

        With a constant term in the objective::

            sage: p = MixedIntegerLinearProgram(solver='ppl')
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[0] + 42)
            sage: p.show()
            Maximization:
              x_0 + 42
            Constraints:
            Variables:
              x_0 is a continuous variable (min=0, max=+oo)
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
            varid_name[i] = s if s else 'x_'+str(i)

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
        d = b.objective_constant_term()
        if d > self._backend.zero(): print "+", d,
        elif d < self._backend.zero(): print "-", -d,
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
            for j, c in sorted(zip(indices, values)):
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
            if b.is_variable_integer(i):
                var_type = 'an integer'
            elif b.is_variable_binary(i):
                var_type = 'a boolean'
            else:
                var_type = 'a continuous'
            if varid_name[i] == str(self.gen(i)):
                name = varid_name[i]
            else:
                name = '{0} = {1}'.format(varid_name[i], self.gen(i))
            lb, ub = b.col_bounds(i)
            print('  {0} is {1} variable (min={2}, max={3})'.format(
                name, var_type, 
                lb if lb is not None else "-oo",
                ub if ub is not None else "+oo"))

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
            sage: x = p.new_variable(nonnegative=True)
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
            sage: x = p.new_variable(nonnegative=True)
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
            sage: x = p.new_variable(nonnegative=True)
            sage: y = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[3] + 3*y[2,9] + x[5])
            sage: p.add_constraint(x[3] + y[2,9] + 2*x[5], max=2)
            sage: p.solve()
            6.0

        To return the optimal value of ``y[2,9]``::

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

        Test that an error is reported when we try to get the value
        of something that is not a variable in this problem::

            sage: p.get_values("Something strange")
            Traceback (most recent call last):
            ...
            TypeError: Not a MIPVariable: ...
            sage: p.get_values("Something stranger", 4711)
            Traceback (most recent call last):
            ...
            TypeError: Not a MIPVariable: ...
            sage: M1 = MixedIntegerLinearProgram()
            sage: M2 = MixedIntegerLinearProgram()
            sage: x = M1.new_variable()
            sage: y = M1.new_variable()
            sage: z = M2.new_variable()
            sage: M2.add_constraint(z[0] <= 5)
            sage: M2.solve()
            0.0
            sage: M2.get_values(x)
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: M2.get_values(y)
            Traceback (most recent call last):
            ...
            ValueError: ...
        """
        val = []
        for l in lists:
            if isinstance(l, MIPVariable):
                if self != l.mip():
                    raise ValueError("Variable {!r} is a variable from a different problem".format(l))
                c = {}
                for (k,v) in l.items():
                    c[k] = self._backend.get_variable_value(self._variables[v])
                val.append(c)
            elif isinstance(l, list):
                if len(l) == 1:
                    val.append([self.get_values(l[0])])
                else:
                    c = []
                    [c.append(self.get_values(ll)) for ll in l]
                    val.append(c)
            elif l in self._variables:
                val.append(self._backend.get_variable_value(self._variables[l]))
            else:
                raise TypeError("Not a MIPVariable: {!r}".format(l))

        if len(lists) == 1:
            return val[0]
        else:
            return val

    def set_objective(self,obj):
        r"""
        Sets the objective of the ``MixedIntegerLinearProgram``.

        INPUT:

        - ``obj`` -- A linear function to be optimized.
          ( can also be set to ``None`` or ``0`` or any number when just
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
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 2/10*x[2], max=4)
            sage: p.add_constraint(1.5*x[1]+3*x[2], max=4)
            sage: round(p.solve(),5)
            6.66667
            sage: p.set_objective(None)
            sage: _ = p.solve()

        TESTS:

        Test whether numbers as constant objective functions are accepted::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(42)
            sage: p.solve() # tol 1e-8
            42

        """
        cdef list values = []

        # If the objective is None, or a constant, we want to remember
        # that the objective function has been defined ( the user did not
        # forget it ). In some LP problems, you just want a feasible solution
        # and do not care about any function being optimal.
        cdef int i

        if obj is None:
            f = {-1 : 0}
        else:
            # See if it is a constant
            R = self.base_ring()
            try:
                f = {-1: R(obj)}
            except TypeError:
                # Should be a linear function
                f = obj.dict()
        d = f.pop(-1,self._backend.zero())

        for i in range(self._backend.ncols()):
            values.append(f.get(i,self._backend.zero()))
        self._backend.set_objective(values, d)

    def add_constraint(self, linear_function, max=None, min=None, name=None):
        r"""
        Adds a constraint to the ``MixedIntegerLinearProgram``.

        INPUT:

        - ``linear_function`` -- Four different types of arguments are
          admissible:

            - A linear function. In this case, one of the arguments
              ``min`` or ``max`` has to be specified.

            - A linear constraint of the form ``A <= B``, ``A >= B``,
              ``A <= B <= C``, ``A >= B >= C`` or ``A == B``.

            - A vector-valued linear function, see
              :mod:`~sage.numerical.linear_tensor`. In this case, one
              of the arguments ``min`` or ``max`` has to be specified.

            - An (in)equality of vector-valued linear functions, that
              is, elements of the space of linear functions tensored
              with a vector space. See
              :mod:`~sage.numerical.linear_tensor_constraints` for
              details.

        - ``max`` -- constant or ``None`` (default). An upper bound on
          the linear function. This must be a numerical value for
          scalar linear functions, or a vector for vector-valued
          linear functions. Not allowed if the ``linear_function``
          argument is a symbolic (in)-equality.

        - ``min`` -- constant or ``None`` (default). A lower bound on
          the linear function. This must be a numerical value for
          scalar linear functions, or a vector for vector-valued
          linear functions. Not allowed if the ``linear_function``
          argument is a symbolic (in)-equality.

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
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[0] + 5*x[1])
            sage: p.add_constraint(x[0] + 0.2*x[1], max=4)
            sage: p.add_constraint(1.5*x[0] + 3*x[1], max=4)
            sage: p.solve()     # rel tol 1e-15
            6.666666666666666

        There are two different ways to add the constraint
        ``x[5] + 3*x[7] <= x[6] + 3`` to a ``MixedIntegerLinearProgram``.

        The first one consists in giving ``add_constraint`` this
        very expression::

            sage: p.add_constraint(x[5] + 3*x[7] <= x[6] + 3)

        The second (slightly more efficient) one is to use the
        arguments ``min`` or ``max``, which can only be numerical
        values::

            sage: p.add_constraint(x[5] + 3*x[7] - x[6], max=3)

        One can also define double-bounds or equality using symbols
        ``<=``, ``>=`` and ``==``::

            sage: p.add_constraint(x[5] + 3*x[7] == x[6] + 3)
            sage: p.add_constraint(x[5] + 3*x[7] <= x[6] + 3 <= x[8] + 27)
        
        Using this notation, the previous program can be written as::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[0] + 5*x[1])
            sage: p.add_constraint(x[0] + 0.2*x[1] <= 4)
            sage: p.add_constraint(1.5*x[0] + 3*x[1] <= 4)
            sage: p.solve()     # rel tol 1e-15
            6.666666666666666

        The two constraints can alse be combined into a single
        vector-valued constraint::

            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[0] + 5*x[1])
            sage: f_vec = vector([1, 1.5]) * x[0] + vector([0.2, 3]) * x[1];  f_vec
            (1.0, 1.5)*x_0 + (0.2, 3.0)*x_1
            sage: p.add_constraint(f_vec, max=vector([4, 4]))
            sage: p.solve()     # rel tol 1e-15
            6.666666666666666

        Instead of specifying the maximum in the optional ``max``
        argument, we can also use (in)equality notation for
        vector-valued linear functions::

            sage: f_vec <= 4    # constant rhs becomes vector
            (1.0, 1.5)*x_0 + (0.2, 3.0)*x_1 <= (4.0, 4.0)
            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[0] + 5*x[1])
            sage: p.add_constraint(f_vec <= 4)
            sage: p.solve()     # rel tol 1e-15
            6.666666666666666

        Finally, one can use the matrix * :class:`MIPVariable`
        notation to write vector-valued linear functions::

            sage: m = matrix([[1.0, 0.2], [1.5, 3.0]]);  m
            [ 1.00000000000000 0.200000000000000]
            [ 1.50000000000000  3.00000000000000]
            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[0] + 5*x[1])
            sage: p.add_constraint(m * x <= 4)
            sage: p.solve()     # rel tol 1e-15
            6.666666666666666

        TESTS:

        Complex constraints::

            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: b = p.new_variable(nonnegative=True)
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

            sage: p = MixedIntegerLinearProgram()
            sage: p.add_constraint(sum([]),min=2)

        Min/Max are numerical ::

            sage: v = p.new_variable(nonnegative=True)
            sage: p.add_constraint(v[3] + v[5], min = v[6])
            Traceback (most recent call last):
            ...
            ValueError: min and max arguments are required to be constants
            sage: p.add_constraint(v[3] + v[5], max = v[6])
            Traceback (most recent call last):
            ...
            ValueError: min and max arguments are required to be constants

        Do not add redundant elements (notice only one copy of each constraint is added)::

            sage: lp = MixedIntegerLinearProgram(solver="GLPK", check_redundant=True)
            sage: for each in xrange(10): lp.add_constraint(lp[0]-lp[1],min=1)
            sage: lp.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              1.0 <= x_0 - x_1
            Variables:
              x_0 is a continuous variable (min=-oo, max=+oo)
              x_1 is a continuous variable (min=-oo, max=+oo)

        We check for constant multiples of constraints as well::

            sage: for each in xrange(10): lp.add_constraint(2*lp[0]-2*lp[1],min=2)
            sage: lp.show()
            Maximization:
            <BLANKLINE>
            Constraints:
              1.0 <= x_0 - x_1
            Variables:
              x_0 is a continuous variable (min=-oo, max=+oo)
              x_1 is a continuous variable (min=-oo, max=+oo)

        But if the constant multiple is negative, we should add it anyway (once)::

              sage: for each in xrange(10): lp.add_constraint(-2*lp[0]+2*lp[1],min=-2)
              sage: lp.show()
              Maximization:
              <BLANKLINE>
              Constraints:
                1.0 <= x_0 - x_1
                -2.0 <= -2.0 x_0 + 2.0 x_1 
              Variables:
                x_0 is a continuous variable (min=-oo, max=+oo)
                x_1 is a continuous variable (min=-oo, max=+oo)

        Catch ``True`` / ``False`` as INPUT (:trac:`13646`)::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(True)
            Traceback (most recent call last):
            ...
            ValueError: argument must be a linear function or constraint, got True
        """
        if linear_function is 0:
            return

        from sage.numerical.linear_functions import is_LinearFunction, is_LinearConstraint
        from sage.numerical.linear_tensor import is_LinearTensor
        from sage.numerical.linear_tensor_constraints import  is_LinearTensorConstraint
        if is_LinearFunction(linear_function) or is_LinearTensor(linear_function):
            # Find the parent for the coefficients
            if is_LinearFunction(linear_function):
                M = linear_function.parent().base_ring()
            elif is_LinearTensor(linear_function):
                if not linear_function.parent().is_vector_space():
                    raise ValueError('the linear function must be vector-valued')
                M = linear_function.parent().free_module()
            else:
                assert False, 'unreachable'
            # Normalize min/max
            try:
                min = None if min is None else M(min)
                max = None if max is None else M(max)
            except (ValueError, TypeError):
                raise ValueError("min and max arguments are required to be constants")
            if min is None and max is None:
                raise ValueError('at least one of "max" or "min" must be set')
            # Shift constant away
            constraint = copy(linear_function.dict())
            try:
                constant_coefficient = constraint.pop(-1)
                max = (max - constant_coefficient) if max is not None else None
                min = (min - constant_coefficient) if min is not None else None
            except KeyError:
                pass
            # Send to backend
            if is_LinearFunction(linear_function):
                if self._check_redundant and self._is_redundant_constraint(constraint, min, max):
                    return
                self._backend.add_linear_constraint(constraint.items(), min, max, name)
            elif is_LinearTensor(linear_function):
                self._backend.add_linear_constraint_vector(M.degree(), constraint.items(), min, max, name)
            else:
                assert False, 'unreachable'
        elif is_LinearConstraint(linear_function):
            if not(min is None and max is None):
                raise ValueError('min and max must not be specified for (in)equalities')
            relation = linear_function
            for lhs, rhs in relation.equations():
                self.add_constraint(lhs-rhs, min=0, max=0, name=name)
            for lhs, rhs in relation.inequalities():
                self.add_constraint(lhs-rhs, max=0, name=name)
        elif is_LinearTensorConstraint(linear_function):
            if not(min is None and max is None):
                raise ValueError('min and max must not be specified for (in)equalities')
            relation = linear_function
            M = relation.parent().linear_tensors().free_module()
            self.add_constraint(relation.lhs() - relation.rhs(), 
                                min=M(0) if relation.is_equation() else None, 
                                max=M(0), name=name)
        else:
            raise ValueError('argument must be a linear function or constraint, got '+str(linear_function))

    def _is_redundant_constraint(self, constraint, min_bound, max_bound):
        """
        Check whether the (scalar) constraint is redundant.

        INPUT:

        - ``constraint`` -- dictionary of a non-zero linear function
          without constant term.

        - ``min_bound``, ``max_bound`` -- base ring elements or
          ``None``. The lower and upper bound.

        OUTPUT:

        Boolean. Whether the (normalized) constraint has already been added.

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram(check_redundant=True)
            sage: mip.add_constraint(x[0], min=1)
            sage: mip._is_redundant_constraint((x[0]).dict(), 1, None)
            True
            sage: mip._is_redundant_constraint((-2*x[0]).dict(), None, -2)
            True
            sage: mip._is_redundant_constraint((x[1]).dict(), 1, None)
            False
        """
        assert self._constraints is not None, 'must be initialized with check_redundant=True'
        assert -1 not in constraint, 'no constant term allowed'
        i0 = min([i for i, c in constraint.items() if c != 0])
        rescale = constraint[i0]
        constraint = tuple((i, c/rescale) for i, c in constraint.items())
        if rescale > 0:
            min_scaled = min_bound/rescale if min_bound is not None else None
            max_scaled = max_bound/rescale if max_bound is not None else None
        else:
            min_scaled = max_bound/rescale if max_bound is not None else None
            max_scaled = min_bound/rescale if min_bound is not None else None
        key = (constraint, min_scaled, max_scaled)
        if key in self._constraints:
            return True
        else:
            self._constraints.append(key)
            return False

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
            sage: x = p.new_variable(nonnegative=True)

        With the following instruction, all the variables
        from x will be binary::

            sage: p.set_binary(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)

        It is still possible, though, to set one of these
        variables as integer while keeping the others as they are::

            sage: p.set_integer(x[3])

        """
        cdef MIPVariable e
        e = <MIPVariable> ee

        if isinstance(e, MIPVariable):
            e._vtype = self.__BINARY
            for v in e.values():
                self._backend.set_variable_type(self._variables[v],self.__BINARY)
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
            sage: v = p.new_variable(nonnegative=True)
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
            sage: x = p.new_variable(nonnegative=True)

        With the following instruction, all the variables
        from x will be integers::

            sage: p.set_integer(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)

        It is still possible, though, to set one of these
        variables as binary while keeping the others as they are::

            sage: p.set_binary(x[3])
        """
        cdef MIPVariable e
        e = <MIPVariable> ee

        if isinstance(e, MIPVariable):
            e._vtype = self.__INTEGER
            for v in e.values():
                self._backend.set_variable_type(self._variables[v],self.__INTEGER)
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
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[1])
            sage: p.is_integer(v[1])
            False
            sage: p.set_integer(v[1])
            sage: p.is_integer(v[1])
            True
        """
        return self._backend.is_variable_integer(self._variables[e])

    def set_real(self,ee):
        r"""
        Sets a variable or a ``MIPVariable`` as real.

        INPUT:

        - ``ee`` -- An instance of ``MIPVariable`` or one of
          its elements.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable(nonnegative=True)

        With the following instruction, all the variables
        from x will be real::

            sage: p.set_real(x)
            sage: p.set_objective(x[0] + x[1])
            sage: p.add_constraint(-3*x[0] + 2*x[1], max=2)

         It is still possible, though, to set one of these
         variables as binary while keeping the others as they are::

            sage: p.set_binary(x[3])

        """
        cdef MIPVariable e
        e = <MIPVariable> ee

        if isinstance(e, MIPVariable):
            e._vtype = self.__REAL
            for v in e.values():
                self._backend.set_variable_type(self._variables[v],self.__REAL)
                self._backend.variable_lower_bound(self._variables[v], 0)
        elif e in self._variables:
            self._backend.set_variable_type(self._variables[e],self.__REAL)
            self._backend.variable_lower_bound(self._variables[e], 0)
        else:
            raise ValueError("e must be an instance of MIPVariable or one of its elements.")

    def is_real(self, e):
        r"""
        Tests whether the variable is real.

        INPUT:

        - ``e`` -- A variable (not a ``MIPVariable``, but one of its elements.)

        OUTPUT:

        ``True`` if the variable is real; ``False`` otherwise.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
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
        return self._backend.is_variable_continuous(self._variables[e])

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
            sage: x = p.new_variable(nonnegative=True)
            sage: p.set_objective(x[1] + 5*x[2])
            sage: p.add_constraint(x[1] + 0.2*x[2], max=4)
            sage: p.add_constraint(1.5*x[1] + 3*x[2], max=4)
            sage: round(p.solve(),6)
            6.666667
            sage: x = p.get_values(x)
            sage: round(x[1],6) # abs tol 1e-15
            0.0
            sage: round(x[2],6)
            1.333333

         Computation of a maximum stable set in Petersen's graph::

            sage: g = graphs.PetersenGraph()
            sage: p = MixedIntegerLinearProgram(maximization=True)
            sage: b = p.new_variable(nonnegative=True)
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

        INPUT:

        - ``v`` -- a variable.

        - ``min`` -- the minimum value the variable can take.  When
          ``min=None``, the variable has no lower bound.

        .. SEEALSO::

            - :meth:`get_min` -- get the minimum value of a variable.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[1])
            sage: p.get_min(v[1])
            0.0
            sage: p.set_min(v[1],6)
            sage: p.get_min(v[1])
            6.0
            sage: p.set_min(v[1], None)
            sage: p.get_min(v[1])

        With a :class:`MIPVariable` as an argument::

            sage: vv = p.new_variable(real=True)
            sage: p.get_min(vv)
            sage: p.get_min(vv[0])
            sage: p.set_min(vv,5)
            sage: p.get_min(vv[0])
            5.0
            sage: p.get_min(vv[9])
            5.0
        """
        try:
            v.set_min(min)
        except AttributeError:
            self._backend.variable_lower_bound(self._variables[v], min)

    def set_max(self, v, max):
        r"""
        Sets the maximum value of a variable.

        INPUT

        - ``v`` -- a variable.

        - ``max`` -- the maximum value the variable can take.  When
          ``max=None``, the variable has no upper bound.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[1])
            sage: p.get_max(v[1])
            sage: p.set_max(v[1],6)
            sage: p.get_max(v[1])
            6.0

        With a :class:`MIPVariable` as an argument::

            sage: vv = p.new_variable(real=True)
            sage: p.get_max(vv)
            sage: p.get_max(vv[0])
            sage: p.set_max(vv,5)
            sage: p.get_max(vv[0])
            5.0
            sage: p.get_max(vv[9])
            5.0
        """
        try:
            v.set_max(max)
        except AttributeError:
            self._backend.variable_upper_bound(self._variables[v], max)

    def get_min(self, v):
        r"""
        Returns the minimum value of a variable.

        INPUT:

        - ``v`` -- a variable

        OUTPUT:

        Minimum value of the variable, or ``None`` if the variable has no lower
        bound.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[1])
            sage: p.get_min(v[1])
            0.0
            sage: p.set_min(v[1],6)
            sage: p.get_min(v[1])
            6.0
            sage: p.set_min(v[1], None)
            sage: p.get_min(v[1])
        """
        try:
            return (<MIPVariable?>v)._lower_bound
        except TypeError:
            return self._backend.variable_lower_bound(self._variables[v])

    def get_max(self, v):
        r"""
        Returns the maximum value of a variable.

        INPUT:

        - ``v`` -- a variable.

        OUTPUT:

        Maximum value of the variable, or ``None`` if the variable has no upper
        bound.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[1])
            sage: p.get_max(v[1])
            sage: p.set_max(v[1],6)
            sage: p.get_max(v[1])
            6.0
        """
        try:
            return (<MIPVariable?>v)._upper_bound
        except TypeError:
            return self._backend.variable_upper_bound(self._variables[v])

    def solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        The solver parameters are by essence solver-specific, which means their
        meaning heavily depends on the solver used.

        (If you do not know which solver you are using, then you use GLPK).

        Aliases:

        Very common parameters have aliases making them solver-independent. For
        example, the following::

            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: p.solver_parameter("timelimit", 60)

        Sets the solver to stop its computations after 60 seconds, and works
        with GLPK, CPLEX and Gurobi.

            - ``"timelimit"`` -- defines the maximum time spent on a
              computation. Measured in seconds.

        Another example is the ``"logfile"`` parameter, which is used to specify
        the file in which computation logs are recorded. By default, the logs
        are not recorded, and we can disable this feature providing an empty
        filename. This is currently working with CPLEX and Gurobi::

            sage: p = MixedIntegerLinearProgram(solver = "CPLEX") # optional - CPLEX
            sage: p.solver_parameter("logfile")                   # optional - CPLEX
            ''
            sage: p.solver_parameter("logfile", "/dev/null")      # optional - CPLEX
            sage: p.solver_parameter("logfile")                   # optional - CPLEX
            '/dev/null'
            sage: p.solver_parameter("logfile", '')               # optional - CPLEX
            sage: p.solver_parameter("logfile")                   # optional - CPLEX
            ''

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
            sage: v = p.new_variable(nonnegative=True)

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

        This might be useful when access to additional functions provided by
        the backend is needed.

        EXAMPLE:

        This example uses the simplex algorithm and prints information::

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
            sage: x, y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: b = p.get_backend()
            sage: b.solver_parameter("simplex_or_intopt", "simplex_only")
            sage: b.solver_parameter("verbosity_simplex", "GLP_MSG_ALL")
            sage: p.solve()  # rel tol 1e-5
            GLPK Simplex Optimizer, v4.55
            2 rows, 2 columns, 4 non-zeros
            *     0: obj =   7.000000000e+00  infeas =  0.000e+00 (0)
            *     2: obj =   9.400000000e+00  infeas =  0.000e+00 (0)
            OPTIMAL LP SOLUTION FOUND
            9.4
        """
        return self._backend

    def get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
            sage: x, y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.solve()  # rel tol 1e-5
            9.4
            sage: p.get_objective_value()  # rel tol 1e-5
            9.4
        """
        return self._backend.get_objective_value()

    def best_known_objective_bound(self):
        r"""
        Return the value of the currently best known bound.

        This method returns the current best upper (resp. lower) bound on the
        optimal value of the objective function in a maximization
        (resp. minimization) problem. It is equal to the output of
        :meth:get_objective_value if the MILP found an optimal solution, but it
        can differ if it was interrupted manually or after a time limit (cf
        :meth:solver_parameter).

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: g = graphs.CubeGraph(9)
            sage: p = MixedIntegerLinearProgram(solver="GLPK")
            sage: p.solver_parameter("mip_gap_tolerance",100)
            sage: b = p.new_variable(binary=True)
            sage: p.set_objective(p.sum(b[v] for v in g))
            sage: for v in g:
            ....:     p.add_constraint(b[v]+p.sum(b[u] for u in g.neighbors(v)) <= 1)
            sage: p.add_constraint(b[v] == 1) # Force an easy non-0 solution
            sage: p.solve() # rel tol 100
            1.0
            sage: p.best_known_objective_bound() # random
            48.0
        """
        return self._backend.best_known_objective_bound()

    def get_relative_objective_gap(self):
        r"""
        Return the relative objective gap of the best known solution.

        For a minimization problem, this value is computed by
        `(\texttt{bestinteger} - \texttt{bestobjective}) / (1e-10 +
        |\texttt{bestobjective}|)`, where ``bestinteger`` is the value returned
        by :meth:`~MixedIntegerLinearProgram.get_objective_value` and
        ``bestobjective`` is the value returned by
        :meth:`~MixedIntegerLinearProgram.best_known_objective_bound`. For a
        maximization problem, the value is computed by `(\texttt{bestobjective}
        - \texttt{bestinteger}) / (1e-10 + |\texttt{bestobjective}|)`.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: g = graphs.CubeGraph(9)
            sage: p = MixedIntegerLinearProgram(solver="GLPK")
            sage: p.solver_parameter("mip_gap_tolerance",100)
            sage: b = p.new_variable(binary=True)
            sage: p.set_objective(p.sum(b[v] for v in g))
            sage: for v in g:
            ....:     p.add_constraint(b[v]+p.sum(b[u] for u in g.neighbors(v)) <= 1)
            sage: p.add_constraint(b[v] == 1) # Force an easy non-0 solution
            sage: p.solve() # rel tol 100
            1.0
            sage: p.get_relative_objective_gap() # random
            46.99999999999999

        TESTS:

        Just make sure that the variable *has* been defined, and is not just
        undefined::

            sage: p.get_relative_objective_gap() > 1
            True
        """
        return self._backend.get_relative_objective_gap()

    def interactive_lp_problem(self,form='standard'):
        r"""
        Returns an InteractiveLPProblem and, if available, a basis.

        INPUT:

        - ``form`` -- (default: ``"standard"``) a string specifying return type: either
          ``None``, or ``"std"`` or ``"standard"``, respectively returns an instance of
          :class:`InteractiveLPProblem` or of :class:`InteractiveLPProblemStandardForm`

        OUTPUT:

        A 2-tuple consists of an instance of class :class:`InteractiveLPProblem` or
        :class:`InteractiveLPProblemStandardForm` that is constructed based on a given
        :class:`MixedIntegerLinearProgram`, and a list of basic
        variables (the basis) if standard form is chosen (by default), otherwise ``None``.

        All variables must have 0 as lower bound and no upper bound.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(names=['m'], solver="GLPK")
            sage: x = p.new_variable(nonnegative=True)
            sage: y = p.new_variable(nonnegative=True, name='n')
            sage: v = p.new_variable(nonnegative=True)
            sage: p.add_constraint( x[0] + x[1] - 7*y[0] + v[0]<= 2, name='K' )
            sage: p.add_constraint( x[1] + 2*y[0] - v[0] <= 3 )
            sage: p.add_constraint( 5*x[0] + y[0] <= 21, name='L' )
            sage: p.set_objective( 2*x[0] + 3*x[1] + 4*y[0] + 5*v[0])
            sage: lp, basis = p.interactive_lp_problem()
            sage: basis
            ['K', 'w_1', 'L']
            sage: lp.constraint_coefficients()
            [ 1.0  1.0 -7.0  1.0]
            [ 0.0  1.0  2.0 -1.0]
            [ 5.0  0.0  1.0  0.0]
            sage: lp.b()
            (2.0, 3.0, 21.0)
            sage: lp.objective_coefficients()
            (2.0, 3.0, 4.0, 5.0)
            sage: lp.decision_variables()
            (m_0, m_1, n_0, x_3)
            sage: view(lp) #not tested
            sage: d = lp.dictionary(*basis)
            sage: view(d) #not tested

        TESTS::

            sage: b = p.get_backend()
            sage: import sage.numerical.backends.glpk_backend as backend
            sage: b.solver_parameter(backend.glp_simplex_or_intopt, backend.glp_simplex_only)
            sage: b.solve()
            0
            sage: lp2, basis2 = p.interactive_lp_problem()
            sage: set(basis2)
            {'n_0', 'w_1', 'x_3'}
            sage: d2 = lp2.dictionary(*basis2)
            sage: d2.is_optimal()
            True
            sage: view(d2) #not tested

            sage: lp3, _ = p.interactive_lp_problem(form=None)
            sage: lp3.constraint_coefficients()
            [ 1.0  1.0 -7.0  1.0]
            [ 0.0  1.0  2.0 -1.0]
            [ 5.0  0.0  1.0  0.0]
            sage: lp3.b()
            (2.0, 3.0, 21.0)
            sage: lp3.objective_coefficients()
            (2.0, 3.0, 4.0, 5.0)
            sage: lp3.decision_variables()
            (m_0, m_1, n_0, x_3)
            sage: view(lp3) #not tested
        """
        back_end = self.get_backend()
        for i in range(self.number_of_variables()):
            if back_end.variable_lower_bound(i) != 0:
                raise ValueError('Problem variables must have 0 as lower bound')
            if back_end.variable_upper_bound(i) is not None:
                raise ValueError('Problem variables must not have upper bound') 

        # Construct 'A'
        coef_matrix = []
        for constraint in self.constraints():
            coef_row = [0] * self.number_of_variables()
            for index, value in zip(constraint[1][0],constraint[1][1]):
                coef_row[index] = value
            coef_matrix.append(coef_row)

        # Construct 'b'
        upper_bound_vector = [c[2] for c in self.constraints()]

        # Raise exception if exist lower bound
        for constraint in self.constraints():
            if constraint[0] != None:
                raise ValueError('Problem constraints cannot have lower bounds')

        # Construct 'c'
        def get_obj_coef(i):
            return back_end.objective_coefficient(i)
        objective_coefs_vector = [get_obj_coef(i) for i in range(self.number_of_variables())]

        def format(name, prefix, index):
            if name:
                return name.replace('[','_').strip(']')
            else:
                return prefix + '_' + str(index)

        # Construct 'x'
        var_names = [format(back_end.col_name(i), 'x', i) for i in range(back_end.ncols())]

        A = coef_matrix
        b = upper_bound_vector
        c = objective_coefs_vector
        x = var_names

        if form is None:
            from sage.numerical.interactive_simplex_method import InteractiveLPProblem
            return InteractiveLPProblem(A, b, c, x), None
        elif form == 'standard' or form == 'std':
            # Construct slack names
            slack_names = [format(back_end.row_name(i), 'w', i) for i in range(back_end.nrows())]
            w = slack_names
            from sage.numerical.interactive_simplex_method import InteractiveLPProblemStandardForm
            lp = InteractiveLPProblemStandardForm(A, b, c, x, slack_variables=w)
            basic_variables = []
            for i, e in enumerate(lp.x()):
                if back_end.is_variable_basic(i):
                    basic_variables.append(str(e))
                elif not back_end.is_variable_nonbasic_at_lower_bound(i):
                    raise ValueError('Invalid column status')
            for i, e in enumerate(lp.slack_variables()):
                if back_end.is_slack_variable_basic(i):
                    basic_variables.append(str(e))
                elif not back_end.is_slack_variable_nonbasic_at_lower_bound(i):
                    raise ValueError('Invalid row status')
            return lp, basic_variables
        else:
            raise ValueError('Form of interactive_lp_problem must be either None or \'standard\'')

class MIPSolverException(RuntimeError):
    r"""
    Exception raised when the solver fails.

    EXAMPLE::

        sage: from sage.numerical.mip import MIPSolverException
        sage: e = MIPSolverException("Error")
        sage: e
        MIPSolverException('Error',)
        sage: print e
        Error

    TESTS:

    No continuous solution::

        sage: p = MixedIntegerLinearProgram(solver="GLPK")
        sage: v = p.new_variable(nonnegative=True)
        sage: p.add_constraint(v[0],max=5.5)
        sage: p.add_constraint(v[0],min=7.6)
        sage: p.set_objective(v[0])

    Tests of GLPK's Exceptions::

        sage: p.solve()
        Traceback (most recent call last):
        ...
        MIPSolverException: GLPK: Problem has no feasible solution

    No integer solution::

        sage: p = MixedIntegerLinearProgram(solver="GLPK")
        sage: v = p.new_variable(nonnegative=True)
        sage: p.add_constraint(v[0],max=5.6)
        sage: p.add_constraint(v[0],min=5.2)
        sage: p.set_objective(v[0])
        sage: p.set_integer(v)

    Tests of GLPK's Exceptions::

        sage: p.solve()
        Traceback (most recent call last):
        ...
        MIPSolverException: GLPK: Problem has no feasible solution
    """
    pass


cdef class MIPVariable(Element):
    r"""
    ``MIPVariable`` is a variable used by the class
    ``MixedIntegerLinearProgram``.

    .. warning::

        You should not instantiate this class directly. Instead, use
        :meth:`MixedIntegerLinearProgram.new_variable`.
    """

    def __init__(self, parent, mip, vtype, name, lower_bound, upper_bound):
        r"""
        Constructor for ``MIPVariable``.

        INPUT:

        - ``parent`` -- :class:`MIPVariableParent`. The parent of the
          MIP variable.

        - ``mip`` -- :class:`MixedIntegerLinearProgram`. The
          underlying linear program.

        - ``vtype`` (integer) -- Defines the type of the variables
          (default is ``REAL``).

        - ``name`` -- A name for the ``MIPVariable``.

        - ``lower_bound``, ``upper_bound`` -- lower bound and upper
          bound on the variable. Set to ``None`` to indicate that the
          variable is unbounded.

        For more informations, see the method
        ``MixedIntegerLinearProgram.new_variable``.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: p.new_variable(nonnegative=True)
            MIPVariable of dimension 1
        """
        super(MIPVariable, self).__init__(parent)
        self._dict = {}
        self._p = mip
        self._vtype = vtype
        self._lower_bound = lower_bound
        self._upper_bound = upper_bound
        self._name = name

    def __getitem__(self, i):
        r"""
        Returns the symbolic variable corresponding to the key.

        Returns the element asked, otherwise creates it.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[0] + v[1])
            sage: v[0]
            x_0
        """
        cdef int j
        if i in self._dict:
            return self._dict[i]
        zero = self._p._backend.zero()
        name = self._name + "[" + str(i) + "]" if self._name else None
        j = self._p._backend.add_variable(
            lower_bound=self._lower_bound,
            upper_bound=self._upper_bound,
            binary=False,
            continuous=True,
            integer=False,
            obj=zero,
            name=name)
        v = self._p.linear_function({j : 1})
        self._p._variables[v] = j
        self._p._backend.set_variable_type(j, self._vtype)
        self._dict[i] = v
        return v

    def set_min(self, min):
        r"""
        Sets a lower bound on the variable.

        INPUT:

        - ``min`` -- a lower bound, or ``None`` to mean that the variable is
          unbounded.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(real=True, nonnegative=True)
            sage: p.get_min(v)
            0
            sage: p.get_min(v[0])
            0.0
            sage: p.set_min(v,4)
            sage: p.get_min(v)
            4
            sage: p.get_min(v[0])
            4.0
        """
        self._lower_bound = min
        for v in self._p._variables:
            self._p.set_min(v,min)

    def set_max(self, max):
        r"""
        Sets an upper bound on the variable.

        INPUT:

        - ``max`` -- an upper bound, or ``None`` to mean that the variable is
          unbounded.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(real=True, nonnegative=True)
            sage: p.get_max(v)
            sage: p.get_max(v[0])
            sage: p.set_max(v,4)
            sage: p.get_max(v)
            4
            sage: p.get_max(v[0])
            4.0
        """
        self._upper_bound = max
        for v in self._p._variables:
            self._p.set_max(v,max)

    def _repr_(self):
        r"""
        Returns a representation of self.

        EXAMPLE::

            sage: p=MixedIntegerLinearProgram()
            sage: v=p.new_variable()
            sage: v
            MIPVariable of dimension 1
        """
        return "MIPVariable of dimension 1"

    def keys(self):
        r"""
        Returns the keys already defined in the dictionary.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
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
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[0] + v[1])
            sage: v.items()
            [(0, x_0), (1, x_1)]
        """
        return self._dict.items()

    def values(self):
        r"""
        Returns the symbolic variables associated to the current dictionary.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p.set_objective(v[0] + v[1])
            sage: v.values()
            [x_0, x_1]
        """
        return self._dict.values()

    def mip(self):
        r"""
        Returns the :class:`MixedIntegerLinearProgram` in which ``self`` is a variable.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable(nonnegative=True)
            sage: p == v.mip()
            True
        """
        return self._p

    cdef _matrix_rmul_impl(self, m):
        """
        Implement the action of a matrix multiplying from the right.
        """
        result = dict()
        for i, row in enumerate(m.rows()):
            x = self[i]
            assert len(x.dict()) == 1
            x_index = x.dict().keys()[0]
            result[x_index] = row
        from sage.modules.free_module import FreeModule
        V = FreeModule(self._p.base_ring(), m.ncols())
        T = self._p.linear_functions_parent().tensor(V)
        return T(result)

    cdef _matrix_lmul_impl(self, m):
        """
        Implement the action of a matrix multiplying from the left.
        """
        result = dict()
        for i, col in enumerate(m.columns()):
            x = self[i]
            assert len(x.dict()) == 1
            x_index = x.dict().keys()[0]
            result[x_index] = col
        from sage.modules.free_module import FreeModule
        V = FreeModule(self._p.base_ring(), m.nrows())
        T = self._p.linear_functions_parent().tensor(V)
        return T(result)

    cpdef _acted_upon_(self, mat, bint self_on_left):
        """
        Act with matrices on MIPVariables.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: v = p.new_variable()
            sage: m = matrix([[1,2], [3,4]])
            sage: v * m
            (1.0, 2.0)*x_0 + (3.0, 4.0)*x_1
            sage: m * v
            (1.0, 3.0)*x_0 + (2.0, 4.0)*x_1
        """
        from sage.matrix.matrix import is_Matrix
        if is_Matrix(mat):
            return self._matrix_rmul_impl(mat) if self_on_left else self._matrix_lmul_impl(mat)


cdef class MIPVariableParent(Parent):
    """
    Parent for :class:`MIPVariable`.

    .. warning::

        This class is for internal use. You should not instantiate it
        yourself. Use :meth:`MixedIntegerLinearProgram.new_variable`
        to generate mip variables.
    """

    Element = MIPVariable

    def _repr_(self):
        r"""
        Return representation of self.

        OUTPUT:

        String.

        EXAMPLES::

            sage: mip.<v> = MixedIntegerLinearProgram()
            sage: v.parent()
            Parent of MIPVariables
        """
        return 'Parent of MIPVariables'

    def _an_element_(self):
        """
        Construct a MIP variable.

        OUTPUT:

        This is required for the coercion framework. We raise a
        ``TypeError`` to abort search for any coercion to another
        parent for binary operations. The only interesting operations
        involving :class:`MIPVariable` elements are actions by
        matrices.

        EXAMPLES::

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: parent = x.parent()
            sage: parent.an_element()    # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: disallow coercion
        """
        raise TypeError('disallow coercion')

    def _element_constructor_(self, mip, vtype, name="", lower_bound=0, upper_bound=None):
        """
        The Element constructor

        INPUT/OUTPUT:

        See :meth:`MIPVariable.__init__`.

        EXAMPLES::

            sage: mip = MixedIntegerLinearProgram()
            sage: mip.new_variable()    # indirect doctest
            MIPVariable of dimension 1
        """
        return self.element_class(self, mip, vtype, name, lower_bound, upper_bound)


mip_variable_parent = MIPVariableParent()

