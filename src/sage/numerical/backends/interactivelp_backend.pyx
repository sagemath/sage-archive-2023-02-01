r"""
COIN Backend

AUTHORS:

- Nathann Cohen (2010-10)      : generic_backend template
- Matthias Koeppe (2016-03)    : this backend

"""

#*****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#       Copyright (C) 2016 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.numerical.mip import MIPSolverException
from sage.numerical.interactive_simplex_method import InteractiveLPProblem, default_variable_name
from sage.modules.all import vector

cdef class InteractiveLPBackend:

    def __cinit__(self, maximization = True, base_ring = None):
        """
        Cython constructor

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")

        This backend can work with irrational algebraic numbers::

            sage: poly = polytopes.dodecahedron(base_ring=AA)
            sage: lp = poly.to_linear_program(solver='InteractiveLP')
            sage: b = lp.get_backend()
            sage: b.set_objective([1, 1, 1])
            sage: lp.solve()
            2.291796067500631?
        """

        if base_ring is None:
            from sage.rings.all import AA
            base_ring = AA

        self.lp = InteractiveLPProblem([], [], [], base_ring=base_ring)
        self.set_verbosity(0)

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

        self.obj_constant_term = 0

        self.row_names = []

    cpdef base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        A ring. The coefficients that the chosen solver supports.

        EXAMPLES::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.base_ring()
            Algebraic Real Field
        """
        return self.lp.base_ring()

    def _variable_type_from_bounds(self, lower_bound, upper_bound):
        """
        Return a string designating a variable type in `InteractiveLPProblem`.

        INPUT:

        - ``lower_bound`` - the lower bound of the variable

        - ``upper_bound`` - the upper bound of the variable

        OUTPUT:

        - a string, one of "", "<=", ">="

        The function raises an error if this pair of bounds cannot be
        represented by an `InteractiveLPProblem` variable type.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p._variable_type_from_bounds(0, None)
            '>='
            sage: p._variable_type_from_bounds(None, 0)
            '<='
            sage: p._variable_type_from_bounds(None, None)
            ''
            sage: p._variable_type_from_bounds(None, 5)
            Traceback (most recent call last):
            ...
            NotImplementedError: General variable bounds not supported
        """
        if lower_bound is None:
            if upper_bound is None:
                return ""
            elif upper_bound == 0:
                return "<="
            else:
                raise NotImplementedError("General variable bounds not supported")
        elif lower_bound == 0:
            if upper_bound is None:
                return ">="
            else:
                raise NotImplementedError("General variable bounds not supported")
        else:
            raise NotImplementedError("General variable bounds not supported")

    cpdef int add_variable(self, lower_bound=0, upper_bound=None,
                           binary=False, continuous=True, integer=False,
                           obj=None, name=None, coefficients=None) except -1:
        ## coefficients is an extension in this backend,
        ## and a proposed addition to the interface, to unify this with add_col.
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

        - ``obj`` - (optional) coefficient of this variable in the objective function (default: 0)

        - ``name`` - an optional name for the newly added variable (default: ``None``).

        - ``coefficients`` -- (optional) an iterable of pairs ``(i, v)``. In each
          pair, ``i`` is a variable index (integer) and ``v`` is a
          value (element of :meth:`base_ring`).

        OUTPUT: The index of the newly created variable

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable(continuous=True, integer=True)
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x',obj=1)
            1
            sage: p.col_name(1)
            'x'
            sage: p.objective_coefficient(1)
            1
        """
        A, b, c, x, constraint_types, variable_types, problem_type, ring = self._AbcxCVPR()
        cdef int vtype = int(binary) + int(continuous) + int(integer)
        if  vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")
        if not continuous:
            raise NotImplementedError("Integer variables are not supported")
        variable_types = variable_types + (self._variable_type_from_bounds(lower_bound, upper_bound),)
        col = vector(ring, self.nrows())
        if coefficients is not None:
            for (i, v) in coefficients:
                col[i] = v
        A = A.augment(col)
        if obj is None:
            obj = 0
        c = tuple(c) + (obj,)
        if name is None:
            var_names = default_variable_name("primal decision")
            name = "{}{:d}".format(var_names, self.ncols() + 1)
        x = tuple(x) + (name,)
        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring)
        return self.ncols() - 1

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

        - ``obj`` - (optional) coefficient of all variables in the objective function (default: 0)

        - ``names`` - optional list of names (default: ``None``)

        OUTPUT: The index of the variable created last.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variables(5)
            4
            sage: p.ncols()
            5
            sage: p.add_variables(2, names=['a','b'])
            6
        """
        for i in range(n):
            self.add_variable(lower_bound, upper_bound, binary, continuous, integer, obj,
                              None if names is None else names[i])
        return self.ncols() - 1

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,-1)
            sage: p.is_variable_continuous(0)
            True
        """
        if vtype == -1:
            pass
        else:
            raise NotImplementedError()

    def _AbcxCVPR(self):
        A, b, c, x = self.lp.Abcx()
        constraint_types = self.lp.constraint_types()
        variable_types = self.lp.variable_types()
        problem_type = self.lp.problem_type()
        base_ring = self.lp.base_ring()
        return A, b, c, x, constraint_types, variable_types, problem_type, base_ring

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        A, b, c, x, constraint_types, variable_types, problem_type, ring = self._AbcxCVPR()
        if sense == +1:
            problem_type = "max"
        else:
            problem_type = "min"
        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring)

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective
        function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2
        """
        if coeff is None:
            return self.lp.objective_coefficients()[variable]
        else:
            A, b, c, x, constraint_types, variable_types, problem_type, ring = self._AbcxCVPR()
            c = list(c)
            c[variable] = coeff
            self.lp = InteractiveLPProblem(A, b, c, x,
                                           constraint_types, variable_types,
                                           problem_type, ring)

    cpdef set_objective(self, list coeff, d = 0):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose i-th element is the
          coefficient of the i-th variable in the objective function.

        - ``d`` (real) -- the constant term in the linear function (set to `0` by default)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: map(lambda x :p.objective_coefficient(x), range(5))
            [1, 1, 2, 1, 3]

        Constants in the objective function are respected::

            sage: p = MixedIntegerLinearProgram(solver='InteractiveLP')
            sage: x,y = p[0], p[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.solve()
            47/5
        """
        A, b, c, x, constraint_types, variable_types, problem_type, ring = self._AbcxCVPR()
        c = coeff
        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring)
        self.obj_constant_term = d

    cpdef set_verbosity(self, int level):
        """
        Set the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.set_verbosity(2)
        """
        self.verbosity = level

    cpdef remove_constraint(self, int i):
        r"""
        Remove a constraint.

        INPUT:

        - ``i`` -- index of the constraint to remove.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="InteractiveLP")
            sage: v = p.new_variable(nonnegative=True)
            sage: x,y = v[0], v[1]
            sage: p.add_constraint(2*x + 3*y, max = 6)
            sage: p.add_constraint(3*x + 2*y, max = 6)
            sage: p.set_objective(x + y + 7)
            sage: p.solve()
            47/5
            sage: p.remove_constraint(0)
            sage: p.solve()
            10
            sage: p.get_values([x,y])
            [0, 3]
        """
        A, b, c, x, constraint_types, variable_types, problem_type, ring = self._AbcxCVPR()
        A = A.delete_rows((i,))
        b = list(b); del b[i]
        constraint_types=list(constraint_types); del constraint_types[i]
        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring)

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
            sage: p.row_bounds(0)
            (2, 2)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1, 1, name='foo')
            sage: p.row_name(1)
            'foo'
        """
        A, b, c, x, constraint_types, variable_types, problem_type, ring = self._AbcxCVPR()
        if lower_bound is None:
           if upper_bound is None:
               raise ValueError("At least one of lower_bound and upper_bound must be provided")
           else:
               constraint_types = constraint_types + ("<=",)
               b = tuple(b) + (upper_bound,)
        else:
            if upper_bound is None:
                constraint_types = constraint_types + (">=",)
                b = tuple(b) + (lower_bound,)
            elif lower_bound == upper_bound:
                constraint_types = constraint_types + ("==",)
                b = tuple(b) + (lower_bound,)
            else:
                raise NotImplementedError("Ranged constraints are not supported")

        row = vector(ring, self.ncols())
        for (i, v) in coefficients:
            row[i] = v
        A = A.stack(row)

        self.row_names.append(name)

        self.lp = InteractiveLPProblem(A, b, c, x,
                                       constraint_types, variable_types,
                                       problem_type, ring)


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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
            5
        """
        self.add_variable(coefficients = zip(indices, coeffs))

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.row(4)
            ([], [])
            sage: p.row_bounds(4)
            (0, None)
        """
        for i in range(number):
            self.add_linear_constraint(zip(range(self.ncols()),[0]*(self.ncols())),
                                       lower_bound, upper_bound,
                                       name=None if names is None else names[i])

    cpdef int solve(self) except -1:
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.solve()
            0
            sage: p.objective_coefficient(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """
        ## FIXME: standard_form should allow to pass slack names (which we would take from row_names).
        ## FIXME: Perhaps also pass the problem name as objective name
        lp_std_form = self.lp_std_form = self.lp.standard_form()
        output = lp_std_form.run_revised_simplex_method()
        ## FIXME: Display output as a side effect if verbosity is high enough
        d = self.final_dictionary = lp_std_form.final_revised_dictionary()
        if d.is_optimal():
            if lp_std_form.auxiliary_variable() in d.basic_variables():
                raise MIPSolverException("InteractiveLP: Problem has no feasible solution")
            else:
                return 0
        else:
            raise MIPSolverException("InteractiveLP: Problem is unbounded")

    cpdef get_objective_value(self):
        """
        Return the value of the objective function.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0,1), (1,2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            15/2
            sage: p.get_variable_value(0)
            0
            sage: p.get_variable_value(1)
            3/2
        """
        d = self.final_dictionary
        v = d.objective_value()
        if self.lp_std_form.is_negative():
            v = - v
        return self.obj_constant_term + v

    cpdef get_variable_value(self, int variable):
        """
        Return the value of a variable given by the solver.

        .. NOTE::

           Behavior is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([(0,1), (1, 2)], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            15/2
            sage: p.get_variable_value(0)
            0
            sage: p.get_variable_value(1)
            3/2
        """
        if str(self.lp.decision_variables()[variable]) != str(self.lp_std_form.decision_variables()[variable]):
            raise NotImplementedError("Undoing the standard-form transformation is not implemented")
        return self.final_dictionary.basic_solution()[variable]

    cpdef int ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            1
            sage: p.ncols()
            2
        """
        return self.lp.n_variables()

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(2, 0, None)
            sage: p.nrows()
            2
        """
        return self.lp.n_constraints()

    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        return self.lp.problem_type() == "max"

    cpdef problem_name(self, char * name = NULL):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.problem_name("There_once_was_a_french_fry")
            sage: print p.problem_name()
            There_once_was_a_french_fry
        """
        if name == NULL:
            if self.prob_name is not None:
                return self.prob_name
            else:
                return ""
        else:
            self.prob_name = str(name)

    ## cpdef write_lp(self, char * name):
    ##     """
    ##     Write the problem to a ``.lp`` file

    ##     INPUT:

    ##     - ``filename`` (string)

    ##     EXAMPLE::

    ##         sage: from sage.numerical.backends.generic_backend import get_solver
    ##         sage: p = get_solver(solver = "InteractiveLP")
    ##         sage: p.add_variables(2)
    ##         2
    ##         sage: p.add_linear_constraint([(0, 1], (1, 2)], None, 3)
    ##         sage: p.set_objective([2, 5])
    ##         sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))
    ##     """
    ##     raise NotImplementedError()

    ## cpdef write_mps(self, char * name, int modern):
    ##     """
    ##     Write the problem to a ``.mps`` file

    ##     INPUT:

    ##     - ``filename`` (string)

    ##     EXAMPLE::

    ##         sage: from sage.numerical.backends.generic_backend import get_solver
    ##         sage: p = get_solver(solver = "InteractiveLP")
    ##         sage: p.add_variables(2)
    ##         2
    ##         sage: p.add_linear_constraint([(0, 1), (1, 2)], None, 3)
    ##         sage: p.set_objective([2, 5])
    ##         sage: p.write_lp(os.path.join(SAGE_TMP, "lp_problem.lp"))
    ##     """
    ##     raise NotImplementedError()

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 0, None)
            sage: p.row(0)
            ([1, 2, 3, 4], [1, 2, 3, 4])
        """
        A, b, c, x = self.lp.Abcx()
        indices = []
        coeffs = []
        for j in range(self.ncols()):
            if A[i][j] != 0:
                indices.append(j)
                coeffs.append(A[i][j])
        return (indices, coeffs)

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row_bounds(0)
            (2, 2)
        """
        A, b, c, x = self.lp.Abcx()
        constraint_types = self.lp.constraint_types()
        if constraint_types[index] == '>=':
            return (b[index], None)
        elif constraint_types[index] == '<=':
            return (None, b[index])
        elif constraint_types[index] == '==':
            return (b[index], b[index])
        else:
            raise ValueError("Bad constraint_type")

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable(lower_bound=None)
            0
            sage: p.col_bounds(0)
            (None, None)
            sage: p.variable_lower_bound(0, 0)
            sage: p.col_bounds(0)
            (0, None)
        """
        t = self.lp.variable_types()[index]
        if t == ">=":
            return (0, None)
        elif t == "<=":
            return (None, 0)
        elif t == "":
            return (None, None)
        else:
            raise ValueError("Bad _variable_types")

    cpdef bint is_variable_binary(self, int index):
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_binary(0)
            False

        """
        return False

    cpdef bint is_variable_integer(self, int index):
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_integer(0)
            False
        """
        return False

    cpdef bint is_variable_continuous(self, int index):
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True

        """
        return True

    cpdef row_name(self, int index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])
            sage: p.row_name(0)
            'Empty constraint 1'

        """
        ## slack_variables = default_variable_name("primal slack")
        ## if slack_variables == default_variable_name("primal decision"):
        ##     index = index + self.ncols()
        ## return "{}{:d}".format(slack_variables, index + 1)
        return self.row_names[index]

    cpdef col_name(self, int index):
        """
        Return the ``index``-th column name

        INPUT:

        - ``index`` (integer) -- the column id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable(name="I_am_a_variable")
            0
            sage: p.col_name(0)
            'I_am_a_variable'
        """
        return str(self.lp.decision_variables()[index])

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
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable(lower_bound=None)
            0
            sage: p.col_bounds(0)
            (None, None)
            sage: p.variable_upper_bound(0) is None
            True
            sage: p.variable_upper_bound(0, 0)
            sage: p.col_bounds(0)
            (None, 0)
            sage: p.variable_upper_bound(0)
            0
            sage: p.variable_upper_bound(0, None)
            sage: p.variable_upper_bound(0) is None
            True
        """
        bounds = self.col_bounds(index)
        if value is False:
            return bounds[1]
        else:
            if value != bounds[1]:
                bounds = (bounds[0], value)
                A, b, c, x, constraint_types, variable_types, problem_type, ring = self._AbcxCVPR()
                variable_types = list(variable_types)
                variable_types[index] = self._variable_type_from_bounds(*bounds)
                self.lp = InteractiveLPProblem(A, b, c, x,
                                               constraint_types, variable_types,
                                               problem_type, ring)

    cpdef variable_lower_bound(self, int index, value = False):
        """
        Return or define the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has no lower bound. When set to ``None``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "InteractiveLP")
            sage: p.add_variable(lower_bound=None)
            0
            sage: p.col_bounds(0)
            (None, None)
            sage: p.variable_lower_bound(0) is None
            True
            sage: p.variable_lower_bound(0, 0)
            sage: p.col_bounds(0)
            (0, None)
            sage: p.variable_lower_bound(0)
            0
            sage: p.variable_lower_bound(0, None)
            sage: p.variable_lower_bound(0) is None
            True
        """
        bounds = self.col_bounds(index)
        if value is False:
            return bounds[0]
        else:
            if value != bounds[0]:
                bounds = (value, bounds[1])
                A, b, c, x, constraint_types, variable_types, problem_type, ring = self._AbcxCVPR()
                variable_types = list(variable_types)
                variable_types[index] = self._variable_type_from_bounds(*bounds)
                self.lp = InteractiveLPProblem(A, b, c, x,
                                               constraint_types, variable_types,
                                               problem_type, ring)

    ## cpdef solver_parameter(self, name, value = None):
    ##     """
    ##     Return or define a solver parameter

    ##     INPUT:

    ##     - ``name`` (string) -- the parameter

    ##     - ``value`` -- the parameter's value if it is to be defined,
    ##       or ``None`` (default) to obtain its current value.

    ##     .. NOTE::

    ##        The list of available parameters is available at
    ##        :meth:`~sage.numerical.mip.MixedIntegerLinearProgram.solver_parameter`.

    ##     EXAMPLE::

    ##         sage: from sage.numerical.backends.generic_backend import get_solver
    ##         sage: p = get_solver(solver = "InteractiveLP")
    ##         sage: p.solver_parameter("timelimit")
    ##         sage: p.solver_parameter("timelimit", 60)
    ##         sage: p.solver_parameter("timelimit")
    ##     """
    ##     raise NotImplementedError()

    cpdef bint is_variable_basic(self, int index):
        """
        Test whether the given variable is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="InteractiveLP")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_variable_basic(0)
            True
            sage: b.is_variable_basic(1)
            False
        """
        return self.lp_std_form.decision_variables()[index] in self.final_dictionary.basic_variables()

    cpdef bint is_variable_nonbasic_at_lower_bound(self, int index):
        """
        Test whether the given variable is nonbasic at lower bound.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="InteractiveLP")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_variable_nonbasic_at_lower_bound(0)
            False
            sage: b.is_variable_nonbasic_at_lower_bound(1)
            True
        """
        return self.lp_std_form.decision_variables()[index] in self.final_dictionary.nonbasic_variables()

    cpdef bint is_slack_variable_basic(self, int index):
        """
        Test whether the slack variable of the given row is basic.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="InteractiveLP")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_slack_variable_basic(0)
            True
            sage: b.is_slack_variable_basic(1)
            False
        """
        return self.lp_std_form.slack_variables()[index] in self.final_dictionary.basic_variables()

    cpdef bint is_slack_variable_nonbasic_at_lower_bound(self, int index):
        """
        Test whether the given variable is nonbasic at lower bound.

        This assumes that the problem has been solved with the simplex method
        and a basis is available.  Otherwise an exception will be raised.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(maximization=True,\
                                                solver="InteractiveLP")
            sage: x = p.new_variable(nonnegative=True)
            sage: p.add_constraint(-x[0] + x[1] <= 2)
            sage: p.add_constraint(8 * x[0] + 2 * x[1] <= 17)
            sage: p.set_objective(11/2 * x[0] - 3 * x[1])
            sage: b = p.get_backend()
            sage: # Backend-specific commands to instruct solver to use simplex method here
            sage: b.solve()
            0
            sage: b.is_slack_variable_nonbasic_at_lower_bound(0)
            False
            sage: b.is_slack_variable_nonbasic_at_lower_bound(1)
            True
        """
        return self.lp_std_form.slack_variables()[index] in self.final_dictionary.nonbasic_variables()

    cpdef dictionary(self):
        # Proposed addition to the general interface,
        # which would for other solvers return backend dictionaries (#18804)
        """
        Return a dictionary representing the current basis.
        """
        return self.final_dictionary

    cpdef interactive_lp_problem(self):

        return self.lp
