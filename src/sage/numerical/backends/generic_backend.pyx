r"""
Generic Backend for LP solvers

This class only lists the methods that should be defined by any
interface with a LP Solver. All these methods immediately raise
``NotImplementedError`` exceptions when called, and are obviously
meant to be replaced by the solver-specific method. This file can also
be used as a template to create a new interface : one would only need
to replace the occurrences of ``"Nonexistent_LP_solver"`` by the
solver's name, and replace ``GenericBackend`` by
``SolverName(GenericBackend)`` so that the new solver extends this
class.

AUTHORS:

- Nathann Cohen (2010-10)      : initial implementation
- Risan (2012-02)              : extension for PPL backend
- Ingolfur Edvardsson (2014-06): extension for CVXOPT backend

"""

#*****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy

cdef class GenericBackend:

    cpdef base_ring(self):
        from sage.rings.all import RDF
        return RDF

    cpdef zero(self):
        return self.base_ring()(0)

    cpdef int add_variable(self, lower_bound=0, upper_bound=None,
                           binary=False, continuous=True, integer=False,
                           obj=None, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        INPUT:

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default: ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of this variable in the objective function (default: 0.0)

        - ``name`` - an optional name for the newly added variable (default: ``None``).

        OUTPUT: The index of the newly created variable

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")    # optional - Nonexistent_LP_solver
            sage: p.ncols()                                           # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                    # optional - Nonexistent_LP_solver
            0
            sage: p.ncols()                                           # optional - Nonexistent_LP_solver
            1
            sage: p.add_variable(binary=True)                         # optional - Nonexistent_LP_solver
            1
            sage: p.add_variable(lower_bound=-2.0, integer=True)      # optional - Nonexistent_LP_solver
            2
            sage: p.add_variable(continuous=True, integer=True)       # optional - Nonexistent_LP_solver
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x',obj=1.0)                    # optional - Nonexistent_LP_solver
            3
            sage: p.col_name(3)                                       # optional - Nonexistent_LP_solver
            'x'
            sage: p.objective_coefficient(3)                          # optional - Nonexistent_LP_solver
            1.0
        """
        raise NotImplementedError()

    cpdef int add_variables(self, int n, lower_bound=0, upper_bound=None, binary=False, continuous=True, integer=False, obj=None, names=None) except -1:
        """
        Add ``n`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        INPUT:

        - ``n`` - the number of new variables (must be > 0)

        - ``lower_bound`` - the lower bound of the variable (default: 0)

        - ``upper_bound`` - the upper bound of the variable (default: ``None``)

        - ``binary`` - ``True`` if the variable is binary (default: ``False``).

        - ``continuous`` - ``True`` if the variable is binary (default: ``True``).

        - ``integer`` - ``True`` if the variable is binary (default: ``False``).

        - ``obj`` - (optional) coefficient of all variables in the objective function (default: 0.0)

        - ``names`` - optional list of names (default: ``None``)

        OUTPUT: The index of the variable created last.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")    # optional - Nonexistent_LP_solver
            sage: p.ncols()                                           # optional - Nonexistent_LP_solver
            0
            sage: p.add_variables(5)                                  # optional - Nonexistent_LP_solver
            4
            sage: p.ncols()                                           # optional - Nonexistent_LP_solver
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, names=['a','b']) # optional - Nonexistent_LP_solver
            6
        """
        raise NotImplementedError()

    @classmethod
    def _test_add_variables(cls, tester=None, **options):
        """
        Run tests on the method :meth:`.add_linear_constraints`.

        TEST::

            sage: from sage.numerical.backends.generic_backend import GenericBackend
            sage: p = GenericBackend()
            sage: p._test_add_variables()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        p = cls()                         # fresh instance of the backend
        if tester is None:
            tester = p._tester(**options)
        # Test from CVXOPT interface (part 1):
        ncols_added = 5
        ncols_before = p.ncols()
        add_variables_result = p.add_variables(ncols_added)
        ncols_after = p.ncols()
        tester.assertEqual(ncols_after, ncols_before+ncols_added, "Added the wrong number of columns")

    cpdef  set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1  Integer
            *  0  Binary
            *  -1  Continuous

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")   # optional - Nonexistent_LP_solver
            sage: p.ncols()                                        # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                  # optional - Nonexistent_LP_solver
            0
            sage: p.set_variable_type(0,1)                          # optional - Nonexistent_LP_solver
            sage: p.is_variable_integer(0)                          # optional - Nonexistent_LP_solver
            True
        """
        raise NotImplementedError()

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

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

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective
        function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            0
            sage: p.objective_coefficient(0)                         # optional - Nonexistent_LP_solver
            0.0
            sage: p.objective_coefficient(0,2)                       # optional - Nonexistent_LP_solver
            sage: p.objective_coefficient(0)                         # optional - Nonexistent_LP_solver
            2.0
        """
        raise NotImplementedError()

    cpdef objective_constant_term(self, d=None):
        """
        Set or get the constant term in the objective function

        INPUT:

        - ``d`` (double) -- its coefficient.  If `None` (default), return the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.objective_constant_term()                        # optional - Nonexistent_LP_solver
            0.0
            sage: p.objective_constant_term(42)                      # optional - Nonexistent_LP_solver
            sage: p.objective_constant_term()                        # optional - Nonexistent_LP_solver
            42.0
        """
        if d is None:
            return self.obj_constant_term
        else:
            self.obj_constant_term = d

    cpdef set_objective(self, list coeff, d = 0.0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose i-th element is the
          coefficient of the i-th variable in the objective function.

        - ``d`` (double) -- the constant term in the linear function (set to `0` by default)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")    # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                                 # optional - Nonexistent_LP_solver
            4
            sage: p.set_objective([1, 1, 2, 1, 3])                   # optional - Nonexistent_LP_solver
            sage: map(lambda x :p.objective_coefficient(x), range(5))  # optional - Nonexistent_LP_solver
            [1.0, 1.0, 2.0, 1.0, 3.0]

        Constants in the objective function are respected::

            sage: p = MixedIntegerLinearProgram(solver='Nonexistent_LP_solver') # optional - Nonexistent_LP_solver
            sage: x,y = p[0], p[1]                              # optional - Nonexistent_LP_solver
            sage: p.add_constraint(2*x + 3*y, max = 6)          # optional - Nonexistent_LP_solver
            sage: p.add_constraint(3*x + 2*y, max = 6)          # optional - Nonexistent_LP_solver
            sage: p.set_objective(x + y + 7)                    # optional - Nonexistent_LP_solver
            sage: p.set_integer(x); p.set_integer(y)            # optional - Nonexistent_LP_solver
            sage: p.solve()                                     # optional - Nonexistent_LP_solver
            9.0
        """
        raise NotImplementedError()

    cpdef set_verbosity(self, int level):
        """
        Set the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.set_verbosity(2)                                # optional - Nonexistent_LP_solver
        """
        raise NotImplementedError()

    cpdef remove_constraint(self, int i):
        r"""
        Remove a constraint.

        INPUT:

        - ``i`` -- index of the constraint to remove.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: v = p.new_variable(nonnegative=True)         # optional - Nonexistent_LP_solver
            sage: x,y = v[0], v[1]                             # optional - Nonexistent_LP_solver
            sage: p.add_constraint(2*x + 3*y, max = 6)         # optional - Nonexistent_LP_solver
            sage: p.add_constraint(3*x + 2*y, max = 6)         # optional - Nonexistent_LP_solver
            sage: p.set_objective(x + y + 7)                   # optional - Nonexistent_LP_solver
            sage: p.set_integer(x); p.set_integer(y)           # optional - Nonexistent_LP_solver
            sage: p.solve()                                    # optional - Nonexistent_LP_solver
            9.0
            sage: p.remove_constraint(0)                       # optional - Nonexistent_LP_solver
            sage: p.solve()                                    # optional - Nonexistent_LP_solver
            10.0
            sage: p.get_values([x,y])                          # optional - Nonexistent_LP_solver
            [0.0, 3.0]
        """
        raise NotImplementedError()

    cpdef remove_constraints(self, constraints):
        r"""
        Remove several constraints.

        INPUT:

        - ``constraints`` -- an iterable containing the indices of the rows to remove.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_constraint(p[0] + p[1], max = 10)           # optional - Nonexistent_LP_solver
            sage: p.remove_constraints([0])                         # optional - Nonexistent_LP_solver
        """
        if type(constraints) == int: self.remove_constraint(constraints)

        cdef int last = self.nrows() + 1

        for c in sorted(constraints, reverse=True):
            if c != last:
                self.remove_constraint(c)
                last = c

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` -- an iterable of pairs ``(i, v)``. In each
          pair, ``i`` is a variable index (integer) and ``v`` is a
          value (element of :meth:`base_ring`).

        - ``lower_bound`` -- element of :meth:`base_ring` or
          ``None``. The lower bound.

        - ``upper_bound`` -- element of :meth:`base_ring` or
          ``None``. The upper bound.

        - ``name`` -- string or ``None``. Optional name for this row.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")             # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                                           # optional - Nonexistent_LP_solver
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2.0, 2.0)  # optional - Nonexistent_LP_solver
            sage: p.row(0)                                                     # optional - Nonexistent_LP_solver
            ([0, 1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                                              # optional - Nonexistent_LP_solver
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo') # optional - Nonexistent_LP_solver
            sage: p.row_name(1)                                                           # optional - Nonexistent_LP_solver
            'foo'
        """
        raise NotImplementedError('add_linear_constraint')

    cpdef add_linear_constraint_vector(self, degree, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a vector-valued linear constraint.

        .. NOTE::

            This is the generic implementation, which will split the
            vector-valued constraint into components and add these
            individually. Backends are encouraged to replace it with
            their own optimized implementation.

        INPUT:

        - ``degree`` -- integer. The vector degree, that is, the
          number of new scalar constraints.

        - ``coefficients`` -- an iterable of pairs ``(i, v)``. In each
          pair, ``i`` is a variable index (integer) and ``v`` is a
          vector (real and of length ``degree``).

        - ``lower_bound`` -- either a vector or ``None``. The
          component-wise lower bound.

        - ``upper_bound`` -- either a vector or ``None``. The
          component-wise upper bound.

        - ``name`` -- string or ``None``. An optional name for all new
          rows.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: coeffs = ([0, vector([1, 2])], [1, vector([2, 3])])
            sage: upper = vector([5, 5])
            sage: lower = vector([0, 0])
            sage: p.add_variables(2)  # optional - Nonexistent_LP_solver
            1
            sage: p.add_linear_constraint_vector(2, coeffs, lower, upper, 'foo')  # optional - Nonexistent_LP_solver
        """
        for d in range(degree):
            coefficients_d = []
            for i, c in coefficients:
                coefficients_d.append((i, c[d]))
            lower_bound_d = None if lower_bound is None else lower_bound[d]
            upper_bound_d = None if upper_bound is None else upper_bound[d] 
            self.add_linear_constraint(coefficients_d, lower_bound_d, upper_bound_d, name=name)

    @classmethod
    def _test_add_linear_constraint_vector(cls, tester=None, **options):
        """
        Run tests on the method :meth:`.add_linear_constraint_vector`.

        TEST::

            sage: from sage.numerical.backends.generic_backend import GenericBackend
            sage: p = GenericBackend()
            sage: p._test_add_linear_constraint_vector()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        p = cls()                         # fresh instance of the backend
        if tester is None:
            tester = p._tester(**options)
        from sage.modules.all import vector
        # Ensure there are at least 2 variables
        p.add_variables(2)
        coeffs = ([0, vector([1, 2])], [1, vector([2, 3])])
        upper = vector([5, 5])
        p.add_linear_constraint_vector(2, coeffs, None, upper, 'foo')
        # FIXME: Tests here. Careful what we expect regarding ranged constraints with some solvers.

    cpdef add_col(self, list indices, list coeffs):
        """
        Add a column.

        INPUT:

        - ``indices`` (list of integers) -- this list contains the
          indices of the constraints in which the variable's
          coefficient is nonzero

        - ``coeffs`` (list of real values) -- associates a coefficient
          to the variable in each of the constraints in which it
          appears. Namely, the i-th entry of ``coeffs`` corresponds to
          the coefficient of the variable in the constraint
          represented by the i-th entry in ``indices``.

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
            sage: p.add_linear_constraints(5, 0, None)            # optional - Nonexistent_LP_solver
            sage: p.add_col(range(5), range(5))                   # optional - Nonexistent_LP_solver
            sage: p.nrows()                                       # optional - Nonexistent_LP_solver
            5
        """
        raise NotImplementedError()

    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=None):
        """
        Add constraints.

        INPUT:

        - ``number`` (integer) -- the number of constraints to add.

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``names`` - an optional list of names (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")   # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                                # optional - Nonexistent_LP_solver
            5
            sage: p.add_linear_constraints(5, None, 2)          # optional - Nonexistent_LP_solver
            sage: p.row(4)                                      # optional - Nonexistent_LP_solver
            ([], [])
            sage: p.row_bounds(4)                               # optional - Nonexistent_LP_solver
            (None, 2.0)
        """
        raise NotImplementedError()

    @classmethod
    def _test_add_linear_constraints(cls, tester=None, **options):
        """
        Run tests on the method :meth:`.add_linear_constraints`.

        TEST::

            sage: from sage.numerical.backends.generic_backend import GenericBackend
            sage: p = GenericBackend()
            sage: p._test_add_linear_constraints()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        p = cls()                         # fresh instance of the backend
        if tester is None:
            tester = p._tester(**options)
        nrows_before = p.nrows()
        nrows_added = 5
        p.add_linear_constraints(nrows_added, None, 2)
        nrows_after = p.nrows()
        # Test correct number of rows
        tester.assertEqual(nrows_after, nrows_before+nrows_added, "Added the wrong number of rows")
        # Test contents of the new rows are correct (sparse zero)
        for i in range(nrows_before, nrows_after):
            tester.assertEqual(p.row(i), ([], []))
            tester.assertEqual(p.row_bounds(i), (None, 2.0))
        # FIXME: Not sure if we should test that no new variables were added.
        # Perhaps some backend may need to introduce explicit slack variables?

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_linear_constraints(5, 0, None)             # optional - Nonexistent_LP_solver
            sage: p.add_col(range(5), range(5))                    # optional - Nonexistent_LP_solver
            sage: p.solve()                                        # optional - Nonexistent_LP_solver
            0
            sage: p.objective_coefficient(0,1)                 # optional - Nonexistent_LP_solver
            sage: p.solve()                                        # optional - Nonexistent_LP_solver
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """
        raise NotImplementedError()

    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            1
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3) # optional - Nonexistent_LP_solver
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

    cpdef best_known_objective_bound(self):
        r"""
        Return the value of the currently best known bound.

        This method returns the current best upper (resp. lower) bound on the
        optimal value of the objective function in a maximization
        (resp. minimization) problem. It is equal to the output of
        :meth:`get_objective_value` if the MILP found an optimal solution, but
        it can differ if it was interrupted manually or after a time limit (cf
        :meth:`solver_parameter`).

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: b = p.new_variable(binary=True)                      # optional - Nonexistent_LP_solver
            sage: for u,v in graphs.CycleGraph(5).edges(labels=False): # optional - Nonexistent_LP_solver
            ....:     p.add_constraint(b[u]+b[v]<=1)                   # optional - Nonexistent_LP_solver
            sage: p.set_objective(p.sum(b[x] for x in range(5)))       # optional - Nonexistent_LP_solver
            sage: p.solve()                                            # optional - Nonexistent_LP_solver
            2.0
            sage: pb = p.get_backend()                                 # optional - Nonexistent_LP_solver
            sage: pb.get_objective_value()                             # optional - Nonexistent_LP_solver
            2.0
            sage: pb.best_known_objective_bound()                      # optional - Nonexistent_LP_solver
            2.0
        """
        raise NotImplementedError()


    cpdef get_relative_objective_gap(self):
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

            sage: p = MixedIntegerLinearProgram(solver="Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: b = p.new_variable(binary=True)                      # optional - Nonexistent_LP_solver
            sage: for u,v in graphs.CycleGraph(5).edges(labels=False): # optional - Nonexistent_LP_solver
            ....:     p.add_constraint(b[u]+b[v]<=1)                   # optional - Nonexistent_LP_solver
            sage: p.set_objective(p.sum(b[x] for x in range(5)))       # optional - Nonexistent_LP_solver
            sage: p.solve()                                            # optional - Nonexistent_LP_solver
            2.0
            sage: pb = p.get_backend()                                 # optional - Nonexistent_LP_solver
            sage: pb.get_objective_value()                             # optional - Nonexistent_LP_solver
            2.0
            sage: pb.get_best_objective_value()                        # optional - Nonexistent_LP_solver
            2.0
            sage: pb.get_relative_objective_gap()                      # optional - Nonexistent_LP_solver
            0.0
        """
        raise NotImplementedError()


    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                              # optional - Nonexistent_LP_solver
            1
            sage: p.add_linear_constraint([(0,1), (1, 2)], None, 3) # optional - Nonexistent_LP_solver
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
        """
        Return the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            1
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            2
        """

        raise NotImplementedError()

    def _test_ncols_nonnegative(self, **options):
        tester = self._tester(**options)
        p = self
        tester.assertGreaterEqual(self.ncols(), 0)
    
    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: p.nrows()                                        # optional - Nonexistent_LP_solver
            0
            sage: p.add_linear_constraints(2, 2.0, None)         # optional - Nonexistent_LP_solver
            sage: p.nrows()                                      # optional - Nonexistent_LP_solver
            2
        """

        raise NotImplementedError()

    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

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
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")   # optional - Nonexistent_LP_solver
            sage: p.problem_name("There once was a french fry") # optional - Nonexistent_LP_solver
            sage: print p.problem_name()                        # optional - Nonexistent_LP_solver
            There once was a french fry
        """

        raise NotImplementedError()

    cpdef write_lp(self, char * name):
        """
        Write the problem to a ``.lp`` file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            2
            sage: p.add_linear_constraint([(0, 1], (1, 2)], None, 3) # optional - Nonexistent_LP_solver
            sage: p.set_objective([2, 5])                          # optional - Nonexistent_LP_solver
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))            # optional - Nonexistent_LP_solver
        """
        raise NotImplementedError()

    cpdef write_mps(self, char * name, int modern):
        """
        Write the problem to a ``.mps`` file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variables(2)                               # optional - Nonexistent_LP_solver
            2
            sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3) # optional - Nonexistent_LP_solver
            sage: p.set_objective([2, 5])                          # optional - Nonexistent_LP_solver
            sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))            # optional - Nonexistent_LP_solver
        """
        raise NotImplementedError()

    cpdef copy(self):
        """
        Returns a copy of self.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: b = p.new_variable() # optional - Nonexistent_LP_solver
            sage: p.add_constraint(b[1] + b[2] <= 6) # optional - Nonexistent_LP_solver
            sage: p.set_objective(b[1] + b[2]) # optional - Nonexistent_LP_solver
            sage: copy(p).solve() # optional - Nonexistent_LP_solver
            6.0
        """
        return self.__copy__()

    # Override this method in backends.
    cpdef __copy__(self):
        """
        Returns a copy of self.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: b = p.new_variable() # optional - Nonexistent_LP_solver
            sage: p.add_constraint(b[1] + b[2] <= 6) # optional - Nonexistent_LP_solver
            sage: p.set_objective(b[1] + b[2]) # optional - Nonexistent_LP_solver
            sage: cp = copy(p.get_backend()) # optional - Nonexistent_LP_solver
            sage: cp.solve() # optional - Nonexistent_LP_solver
            0
            sage: cp.get_objective_value() # optional - Nonexistent_LP_solver
            6.0
        """
        raise NotImplementedError()

    def __deepcopy__(self, memo={}):
        """
        Return a deep copy of ``self``.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: b = p.new_variable() # optional - Nonexistent_LP_solver
            sage: p.add_constraint(b[1] + b[2] <= 6) # optional - Nonexistent_LP_solver
            sage: p.set_objective(b[1] + b[2]) # optional - Nonexistent_LP_solver
            sage: cp = deepcopy(p.get_backend()) # optional - Nonexistent_LP_solver
            sage: cp.solve() # optional - Nonexistent_LP_solver
            0
            sage: cp.get_objective_value() # optional - Nonexistent_LP_solver
            6.0
        """
        return self.__copy__()

    cpdef row(self, int i):
        """
        Return a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variables(5)                               # optional - Nonexistent_LP_solver
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2) # optional - Nonexistent_LP_solver
            sage: p.row(0)                                     # optional - Nonexistent_LP_solver
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0]) ## FIXME: Why backwards?
            sage: p.row_bounds(0)                              # optional - Nonexistent_LP_solver
            (2.0, 2.0)
        """
        raise NotImplementedError()

    cpdef row_bounds(self, int index):
        """
        Return the bounds of a specific constraint.

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
            4
            sage: p.add_linear_constraint(range(5), range(5), 2, 2) # optional - Nonexistent_LP_solver
            sage: p.row(0)                                     # optional - Nonexistent_LP_solver
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0]) ## FIXME: Why backwards?
            sage: p.row_bounds(0)                              # optional - Nonexistent_LP_solver
            (2.0, 2.0)
        """
        raise NotImplementedError()

    cpdef col_bounds(self, int index):
        """
        Return the bounds of a specific variable.

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
            0
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                 # optional - Nonexistent_LP_solver
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, 5.0)
        """
        raise NotImplementedError()

    cpdef bint is_variable_binary(self, int index):
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            0
            sage: p.set_variable_type(0,0)                         # optional - Nonexistent_LP_solver
            sage: p.is_variable_binary(0)                          # optional - Nonexistent_LP_solver
            True

        """
        raise NotImplementedError()

    cpdef bint is_variable_integer(self, int index):
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            0
            sage: p.set_variable_type(0,1)                         # optional - Nonexistent_LP_solver
            sage: p.is_variable_integer(0)                         # optional - Nonexistent_LP_solver
            True
        """
        raise NotImplementedError()

    cpdef bint is_variable_continuous(self, int index):
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.ncols()                                       # optional - Nonexistent_LP_solver
            0
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            0
            sage: p.is_variable_continuous(0)                      # optional - Nonexistent_LP_solver
            True
            sage: p.set_variable_type(0,1)                         # optional - Nonexistent_LP_solver
            sage: p.is_variable_continuous(0)                      # optional - Nonexistent_LP_solver
            False

        """
        raise NotImplementedError()

    cpdef row_name(self, int index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])  # optional - Nonexistent_LP_solver
            sage: p.row_name(0)                                     # optional - Nonexistent_LP_solver
            'Empty constraint 1'

        """
        raise NotImplementedError()

    cpdef col_name(self, int index):
        """
        Return the ``index``-th column name

        INPUT:

        - ``index`` (integer) -- the column id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable(name="I am a variable")            # optional - Nonexistent_LP_solver
            1
            sage: p.col_name(0)                                     # optional - Nonexistent_LP_solver
            'I am a variable'
        """
        raise NotImplementedError()

    def _do_test_problem_data(self, tester, cp):
        """
        TESTS:

        Test, with an actual working backend, that comparing a problem with itself works::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver='GLPK')
            sage: tester = p._tester()
            sage: p._do_test_problem_data(tester, p)
        """
        tester.assertEqual(type(self), type(cp),
                           "Classes do not match")
        def assert_equal_problem_data(method):
            tester.assertEqual(getattr(self, method)(), getattr(cp, method)(),
                               "{} does not match".format(method))
        for method in ("ncols", "nrows", "objective_constant_term", "problem_name", "is_maximization"):
            assert_equal_problem_data(method)
        def assert_equal_col_data(method):
            for i in range(self.ncols()):
                tester.assertEqual(getattr(self, method)(i), getattr(cp, method)(i),
                                   "{}({}) does not match".format(method, i))
        for method in ("objective_coefficient", "is_variable_binary", "is_variable_binary", "is_variable_integer",
                       "is_variable_continuous", "col_bounds", "col_name"):
            # don't test variable_lower_bound, variable_upper_bound because we already test col_bounds.
            # TODO: Add a test elsewhere to ensure that variable_lower_bound, variable_upper_bound
            # are consistent with col_bounds.
            assert_equal_col_data(method)
        def assert_equal_row_data(method):
            for i in range(self.nrows()):
                tester.assertEqual(getattr(self, method)(i), getattr(cp, method)(i),
                                   "{}({}) does not match".format(method, i))
        for method in ("row_bounds", "row", "row_name"):
            assert_equal_row_data(method)
    
    def _test_copy(self, **options):
        """
        Test whether the backend can be copied
        and at least the problem data of the copy is equal to that of the original.
        Does not test whether solutions or solver parameters are copied.
        """
        tester = self._tester(**options)
        cp = copy(self)
        self._do_test_problem_data(tester, cp)

    def _test_copy_does_not_share_data(self, **options):
        """
        Test whether copy makes an independent copy of the backend.
        """
        tester = self._tester(**options)
        cp = copy(self)
        cpcp = copy(cp)
        del cp
        self._do_test_problem_data(tester, cpcp)

    # TODO: We should have a more systematic way of generating MIPs for testing.
    @classmethod
    def _test_copy_some_mips(cls, tester=None, **options):
        p = cls()                         # fresh instance of the backend
        if tester is None:
            tester = p._tester(**options)
        # From doctest of GenericBackend.solve:
        p.add_linear_constraints(5, 0, None)
        # p.add_col(range(5), range(5))     -- bad test because COIN sparsifies the 0s away on copy
        p.add_col(range(5), range(1, 6))
        # From doctest of GenericBackend.problem_name:
        p.problem_name("There once was a french fry")
        p._test_copy(**options)
        p._test_copy_does_not_share_data(**options)

    cpdef variable_upper_bound(self, int index, value = False):
        """
        Return or define the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``None``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            0
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)                 # optional - Nonexistent_LP_solver
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, 5.0)
        """
        raise NotImplementedError()

    cpdef variable_lower_bound(self, int index, value = False):
        """
        Return or define the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``None``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.add_variable()                                 # optional - Nonexistent_LP_solver
            0
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)                 # optional - Nonexistent_LP_solver
            sage: p.col_bounds(0)                              # optional - Nonexistent_LP_solver
            (5.0, None)
        """
        raise NotImplementedError()

    cpdef solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        .. NOTE::

           The list of available parameters is available at
           :meth:`~sage.numerical.mip.MixedIntegerLinearProgram.solver_parameter`.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "Nonexistent_LP_solver")  # optional - Nonexistent_LP_solver
            sage: p.solver_parameter("timelimit")                   # optional - Nonexistent_LP_solver
            sage: p.solver_parameter("timelimit", 60)               # optional - Nonexistent_LP_solver
            sage: p.solver_parameter("timelimit")                   # optional - Nonexistent_LP_solver
        """
        raise NotImplementedError()

    cpdef bint is_variable_basic(self, int index):
        """
        Test whether the given variable is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: x = p.new_variable(nonnegative=True)           # optional - Nonexistent_LP_solver
            sage: p.add_constraint(-x[0] + x[1] <= 2)            # optional - Nonexistent_LP_solver
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)    # optional - Nonexistent_LP_solver
            sage: p.set_objective(5.5 * x[0] - 3 * x[1])         # optional - Nonexistent_LP_solver
            sage: b = p.get_backend()                            # optional - Nonexistent_LP_solver
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()                                      # optional - Nonexistent_LP_solver
            0
            sage: b.is_variable_basic(0)                         # optional - Nonexistent_LP_solver
            True
            sage: b.is_variable_basic(1)                         # optional - Nonexistent_LP_solver
            False
        """
        raise NotImplementedError()

    cpdef bint is_variable_nonbasic_at_lower_bound(self, int index):
        """
        Test whether the given variable is nonbasic at lower bound.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: x = p.new_variable(nonnegative=True)           # optional - Nonexistent_LP_solver
            sage: p.add_constraint(-x[0] + x[1] <= 2)            # optional - Nonexistent_LP_solver
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)    # optional - Nonexistent_LP_solver
            sage: p.set_objective(5.5 * x[0] - 3 * x[1])         # optional - Nonexistent_LP_solver
            sage: b = p.get_backend()                            # optional - Nonexistent_LP_solver
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()                                      # optional - Nonexistent_LP_solver
            0
            sage: b.is_variable_nonbasic_at_lower_bound(0)       # optional - Nonexistent_LP_solver
            False
            sage: b.is_variable_nonbasic_at_lower_bound(1)       # optional - Nonexistent_LP_solver
            True
        """
        raise NotImplementedError()

    cpdef bint is_slack_variable_basic(self, int index):
        """
        Test whether the slack variable of the given row is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: x = p.new_variable(nonnegative=True)           # optional - Nonexistent_LP_solver
            sage: p.add_constraint(-x[0] + x[1] <= 2)            # optional - Nonexistent_LP_solver
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)    # optional - Nonexistent_LP_solver
            sage: p.set_objective(5.5 * x[0] - 3 * x[1])         # optional - Nonexistent_LP_solver
            sage: b = p.get_backend()                            # optional - Nonexistent_LP_solver
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()                                      # optional - Nonexistent_LP_solver
            0
            sage: b.is_slack_variable_basic(0)                   # optional - Nonexistent_LP_solver
            True
            sage: b.is_slack_variable_basic(1)                   # optional - Nonexistent_LP_solver
            False
        """
        raise NotImplementedError()

    cpdef bint is_slack_variable_nonbasic_at_lower_bound(self, int index):
        """
        Test whether the given variable is nonbasic at lower bound.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="Nonexistent_LP_solver") # optional - Nonexistent_LP_solver
            sage: x = p.new_variable(nonnegative=True)           # optional - Nonexistent_LP_solver
            sage: p.add_constraint(-x[0] + x[1] <= 2)            # optional - Nonexistent_LP_solver
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)    # optional - Nonexistent_LP_solver
            sage: p.set_objective(5.5 * x[0] - 3 * x[1])         # optional - Nonexistent_LP_solver
            sage: b = p.get_backend()                            # optional - Nonexistent_LP_solver
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()                                      # optional - Nonexistent_LP_solver
            0
            sage: b.is_slack_variable_nonbasic_at_lower_bound(0) # optional - Nonexistent_LP_solver
            False
            sage: b.is_slack_variable_nonbasic_at_lower_bound(1) # optional - Nonexistent_LP_solver
            True
        """
        raise NotImplementedError()

default_solver = None

def default_mip_solver(solver = None):
    """
    Returns/Sets the default MILP Solver used by Sage

    INPUT:

    - ``solver`` -- defines the solver to use:

        - GLPK (``solver="GLPK"``). See the `GLPK
          <http://www.gnu.org/software/glpk/>`_ web site.

        - GLPK's implementation of an exact rational simplex
          method (``solver="GLPK/exact"``).

        - COIN Branch and Cut (``solver="Coin"``). See the `COIN-OR
          <http://www.coin-or.org>`_ web site.

        - CPLEX (``solver="CPLEX"``). See the
          `CPLEX <http://www.ilog.com/products/cplex/>`_ web site.

        - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT
          <http://cvxopt.org/>`_ web site.

        - Gurobi (``solver="Gurobi"``). See the `Gurobi
          <http://www.gurobi.com/>`_ web site.

        - PPL (``solver="PPL"``). See the `PPL
          <http://bugseng.com/products/ppl/>`_ web site. This solver is
          an exact rational solver.

        - ``InteractiveLPProblem`` (``solver="InteractiveLP"``).  A didactical
          implementation of the revised simplex method in Sage.  It works over
          any exact ordered field, the default is ``QQ``.

        ``solver`` should then be equal to one of ``"GLPK"``,
        ``"Coin"``, ``"CPLEX"``,  ``"CVXOPT"``, ``"Gurobi"``, ``"PPL"`, or
        ``"InteractiveLP"``,

        - If ``solver=None`` (default), the current default solver's name is
          returned.

    OUTPUT:

    This function returns the current default solver's name if ``solver = None``
    (default). Otherwise, it sets the default solver to the one given. If this
    solver does not exist, or is not available, a ``ValueError`` exception is
    raised.

    EXAMPLE::

        sage: former_solver = default_mip_solver()
        sage: default_mip_solver("GLPK")
        sage: default_mip_solver()
        'Glpk'
        sage: default_mip_solver("PPL")
        sage: default_mip_solver()
        'Ppl'
        sage: default_mip_solver("GUROBI")
        Traceback (most recent call last):
        ...
        ValueError: Gurobi is not available. Please refer to the documentation to install it.
        sage: default_mip_solver("Yeahhhhhhhhhhh")
        Traceback (most recent call last):
        ...
        ValueError: 'solver' should be set to ...
        sage: default_mip_solver(former_solver)
    """
    global default_solver

    if solver is None:

        if default_solver is not None:
            return default_solver

        else:
            for s in ["Cplex", "Gurobi", "Coin", "Glpk"]:
                try:
                    default_mip_solver(s)
                    return s
                except ValueError:
                    pass

    solver = solver.capitalize()

    if solver == "Cplex":
        try:
            from sage.numerical.backends.cplex_backend import CPLEXBackend
            default_solver = solver
        except ImportError:
            raise ValueError("CPLEX is not available. Please refer to the documentation to install it.")

    elif solver == "Coin":
        try:
            from sage.numerical.backends.coin_backend import CoinBackend
            default_solver = solver
        except ImportError:
            raise ValueError("COIN is not available. Please refer to the documentation to install it.")

    elif solver == "Cvxopt":
        try:
            from sage.numerical.backends.cvxopt_backend import CVXOPTBackend
            default_solver = solver
        except ImportError:
            raise ValueError("CVXOPT is not available. Please refer to the documentation to install it.")

    elif solver == "Ppl":
        try:
            from sage.numerical.backends.ppl_backend import PPLBackend
            default_solver = solver
        except ImportError:
            raise ValueError("PPL is not available. Please refer to the documentation to install it.")

    elif solver == "Gurobi":
        try:
            from sage.numerical.backends.gurobi_backend import GurobiBackend
            default_solver = solver
        except ImportError:
            raise ValueError("Gurobi is not available. Please refer to the documentation to install it.")

    elif solver == "Glpk" or solver == "Glpk/exact":
        default_solver = solver

    elif solver == "Interactivelp":
        default_solver = solver

    else:
        raise ValueError("'solver' should be set to 'GLPK', 'Coin', 'CPLEX', 'CVXOPT', 'Gurobi', 'PPL', 'InteractiveLP', or None.")

cpdef GenericBackend get_solver(constraint_generation = False, solver = None, base_ring = None):
    """
    Return a solver according to the given preferences

    INPUT:

    - ``solver`` -- 7 solvers should be available through this class:

        - GLPK (``solver="GLPK"``). See the `GLPK
          <http://www.gnu.org/software/glpk/>`_ web site.

        - GLPK's implementation of an exact rational simplex
          method (``solver="GLPK/exact"``).

        - COIN Branch and Cut (``solver="Coin"``). See the `COIN-OR
          <http://www.coin-or.org>`_ web site.

        - CPLEX (``solver="CPLEX"``). See the
          `CPLEX <http://www.ilog.com/products/cplex/>`_ web site.

        - CVXOPT (``solver="CVXOPT"``). See the `CVXOPT
          <http://cvxopt.org/>`_ web site.

        - Gurobi (``solver="Gurobi"``). See the `Gurobi
          <http://www.gurobi.com/>`_ web site.

        - PPL (``solver="PPL"``). See the `PPL
          <http://bugseng.com/products/ppl/>`_ web site.  This solver is
          an exact rational solver.

        - ``InteractiveLPProblem`` (``solver="InteractiveLP"``).  A didactical
          implementation of the revised simplex method in Sage.  It works over
          any exact ordered field, the default is ``QQ``.

        ``solver`` should then be equal to one of the above strings,
        or ``None`` (default), in which case the default solver is used
        (see ``default_mip_solver`` method).

        ``solver`` can also be a callable, in which case it is called,
        and its result is returned.

    - ``base_ring`` -- If not ``None``, request a solver that works over this
        (ordered) field.  If ``base_ring`` is not a field, its fraction field
        is used.

        For example, is ``base_ring=ZZ`` is provided, the solver will work over
        the rational numbers.  This is unrelated to whether variables are
        constrained to be integers or not.

    - ``constraint_generation`` -- Only used when ``solver=None``.

      - When set to ``True``, after solving the ``MixedIntegerLinearProgram``,
        it is possible to add a constraint, and then solve it again.
        The effect is that solvers that do not support this feature will not be
        used.

      - Defaults to ``False``.

    .. SEEALSO::

        - :func:`default_mip_solver` -- Returns/Sets the default MIP solver.

    EXAMPLES::

        sage: from sage.numerical.backends.generic_backend import get_solver
        sage: p = get_solver()
        sage: p = get_solver(base_ring=RDF)
        sage: p.base_ring()
        Real Double Field
        sage: p = get_solver(base_ring=QQ); p
        <...sage.numerical.backends.ppl_backend.PPLBackend...>
        sage: p = get_solver(base_ring=ZZ); p
        <...sage.numerical.backends.ppl_backend.PPLBackend...>
        sage: p.base_ring()
        Rational Field
        sage: p = get_solver(base_ring=AA); p
        <...sage.numerical.backends.interactivelp_backend.InteractiveLPBackend...>
        sage: p.base_ring()
        Algebraic Real Field
        sage: d = polytopes.dodecahedron()
        sage: p = get_solver(base_ring=d.base_ring()); p
        <...sage.numerical.backends.interactivelp_backend.InteractiveLPBackend...>
        sage: p.base_ring()
        Number Field in sqrt5 with defining polynomial x^2 - 5
        sage: p = get_solver(solver='InteractiveLP', base_ring=QQ); p
        <...sage.numerical.backends.interactivelp_backend.InteractiveLPBackend...>
        sage: p.base_ring()
        Rational Field

    Passing a callable as the 'solver'::

        sage: from sage.numerical.backends.glpk_backend import GLPKBackend
        sage: p = get_solver(GLPKBackend); p
        <...sage.numerical.backends.glpk_backend.GLPKBackend...>

    Passing a callable that customizes a backend::

        sage: def glpk_exact_solver():
        ....:     from sage.numerical.backends.generic_backend import get_solver
        ....:     b = get_solver(solver="GLPK")
        ....:     b.solver_parameter("simplex_or_intopt", "exact_simplex_only")
        ....:     return b
        sage: delsarte_bound_additive_hamming_space(11,3,4,solver=glpk_exact_solver) # long time
        glp_exact...
        ...
        OPTIMAL SOLUTION FOUND
        8

    """
    if solver is None:

        solver = default_mip_solver()

        if base_ring is not None:
            base_ring = base_ring.fraction_field()
            from sage.rings.all import QQ, RDF
            if base_ring is QQ:
                solver = "Ppl"
            elif solver in ["Interactivelp", "Ppl"] and not base_ring.is_exact():
                solver = "Glpk"
            elif base_ring is not RDF:
                solver = "Interactivelp"

        # We do not want to use Coin for constraint_generation. It just does not
        # work
        if solver == "Coin" and constraint_generation:
            solver = "Glpk"

    elif callable(solver):
        kwds = {}
        if base_ring is not None:
            kwds['base_ring']=base_ring
        return solver(**kwds)

    else:
        solver = solver.capitalize()

    if solver == "Coin":
        from sage.numerical.backends.coin_backend import CoinBackend
        return CoinBackend()

    elif solver == "Glpk":
        from sage.numerical.backends.glpk_backend import GLPKBackend
        return GLPKBackend()

    elif solver == "Glpk/exact":
        from sage.numerical.backends.glpk_exact_backend import GLPKExactBackend
        return GLPKExactBackend()

    elif solver == "Cplex":
        from sage.numerical.backends.cplex_backend import CPLEXBackend
        return CPLEXBackend()

    elif solver == "Cvxopt":
        from sage.numerical.backends.cvxopt_backend import CVXOPTBackend
        return CVXOPTBackend()

    elif solver == "Gurobi":
        from sage.numerical.backends.gurobi_backend import GurobiBackend
        return GurobiBackend()

    elif solver == "Ppl":
        from sage.numerical.backends.ppl_backend import PPLBackend
        return PPLBackend(base_ring=base_ring)

    elif solver == "Interactivelp":
        from sage.numerical.backends.interactivelp_backend import InteractiveLPBackend
        return InteractiveLPBackend(base_ring=base_ring)

    else:
        raise ValueError("'solver' should be set to 'GLPK', 'GLPK/exact', 'Coin', 'CPLEX', 'CVXOPT', 'Gurobi', 'PPL', 'InteractiveLP', None (in which case the default one is used), or a callable.")
