r"""
GLPK Backend
"""

from sage.numerical.mip import MIPSolverException

cdef class GLPKBackend(GenericBackend):

    def __cinit__(self):
        r"""
        Constructor

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
        """

        self.lp = glp_create_prob()
        self.iocp = new_c_glp_iocp()
        glp_init_iocp(self.iocp)
        self.iocp.presolve = GLP_ON
        glp_set_obj_dir(self.lp, GLP_MAX)
        self.set_log_level(0)
        #self.iocp.gmi_cuts = GLP_ON
        #self.iocp.fp_heur = GLP_ON
        #self.iocp.mir_cuts = GLP_ON

    cpdef int add_variable(self):
        r"""
        Adds a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.n_cols()
            0
            sage: p.add_variable()
            1
            sage: p.n_cols()
            1
        """

        glp_add_cols(self.lp, 1)
        cdef int n_var = glp_get_num_cols(self.lp)
        glp_set_col_bnds(self.lp, n_var, GLP_LO, 0, 0)
        glp_set_col_kind(self.lp, n_var, GLP_CV)

        return n_var

    cpdef int add_variables(self, int number):
        r"""
        Adds ``number`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.n_cols()
            0
            sage: p.add_variables(5)
            5
            sage: p.n_cols()
            5
        """

        glp_add_cols(self.lp, number)

        cdef int n_var
        n_var = glp_get_num_cols(self.lp)


        cdef int i

        for 0<= i < number:
            glp_set_col_bnds(self.lp, n_var-i, GLP_LO, 0, 0)
            glp_set_col_kind(self.lp, n_var-i, GLP_CV)

        return n_var

    cpdef set_variable_type(self, int variable, int vtype):
        r"""
        Sets the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1  Integer
            *  0  Binary
            * -1 Real

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.n_cols()
            0
            sage: p.add_variable()
            1
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True
        """

        if vtype==1:
            glp_set_col_kind(self.lp, variable+1, GLP_IV)

        elif vtype==0:
            glp_set_col_kind(self.lp, variable+1, GLP_BV)

        else:
            glp_set_col_kind(self.lp, variable+1, GLP_CV)

    cpdef set_direction(self, int sense):
        r"""
        Sets the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.is_maximization()
            True
            sage: p.set_direction(-1)
            sage: p.is_maximization()
            False
        """
        if sense == 1:
            glp_set_obj_dir(self.lp, GLP_MAX)
        else:
            glp_set_obj_dir(self.lp, GLP_MIN)

    cpdef set_objective_coeff(self, int variable, double coeff):
        r"""
        Sets the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variable()
            1
            sage: p.get_objective_coeff(0)
            0.0
            sage: p.set_objective_coeff(0,2)
            sage: p.get_objective_coeff(0)
            2.0
        """

        glp_set_obj_coef(self.lp, variable + 1, coeff)


    cpdef set_problem_name(self, char * name):
        r"""
        Sets the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.set_problem_name("There once was a french fry")
            sage: print p.get_problem_name()
            There once was a french fry
        """

        glp_set_prob_name(self.lp, name)

    cpdef  get_problem_name(self):
        r"""
        Returns the problem's name

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.set_problem_name("There once was a french fry")
            sage: print p.get_problem_name()
            There once was a french fry
        """

        cdef char * name = <char *> glp_get_prob_name(self.lp)
        if name == NULL:
            return ""
        else:
            return name

    cpdef set_objective(self, list coeff):
        r"""
        Sets the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variables(5)
            5
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: map(lambda x :p.get_objective_coeff(x), range(5))
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """

        cdef int i

        for i,v in enumerate(coeff):
            glp_set_obj_coef(self.lp, i+1, v)

    cpdef set_log_level(self, int level):
        r"""
        Sets the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.set_log_level(2)

        """
        if level == 0:
            self.iocp.msg_lev = GLP_MSG_OFF
        elif level == 1:
            self.iocp.msg_lev = GLP_MSG_ERR
        elif level == 2:
            self.iocp.msg_lev = GLP_MSG_ON
        else:
            self.iocp.msg_lev = GLP_MSG_ALL

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

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variables(5)
            5
            sage: p.add_constraints(5, +1, 2)
            sage: p.get_row(4)
            ([], [])
            sage: p.get_row_bounds(4)
            (None, 2.0)
        """

        glp_add_rows(self.lp, number)
        cdef int n = glp_get_num_rows(self.lp)

        if direction == +1:
            direction = GLP_UP
        elif direction == -1:
            direction = GLP_LO
        else:
            direction = GLP_FX

        cdef int i
        for 0<= i < number:
            glp_set_row_bnds(self.lp, n-i, direction, bound, bound)

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

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variables(5)
            5
            sage: p.add_constraint(range(5), range(5), 0, 2)
            sage: p.get_row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.get_row_bounds(0)
            (2.0, 2.0)
        """

        glp_add_rows(self.lp, 1)
        cdef int n = glp_get_num_rows(self.lp)

        if direction == +1:
            direction = GLP_UP
        elif direction == -1:
            direction = GLP_LO
        else:
            direction = GLP_FX

        cdef int * row_i
        cdef double * row_values

        row_i = <int *> sage_malloc((len(indices)+1) * sizeof(int))
        row_values = <double *> sage_malloc((len(indices)+1) * sizeof(double))

        for i,v in enumerate(indices):
            row_i[i+1] = v+1
        for i,v in enumerate(coeffs):
            row_values[i+1] = v

        glp_set_mat_row(self.lp, n, len(indices), row_i, row_values)
        glp_set_row_bnds(self.lp, n, direction, bound, bound)

    cpdef get_row(self, int index):
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

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variables(5)
            5
            sage: p.add_constraint(range(5), range(5), 0, 2)
            sage: p.get_row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.get_row_bounds(0)
            (2.0, 2.0)
        """
        cdef int n = glp_get_num_cols(self.lp)
        cdef int * c_indices = <int*> sage_malloc((n+1)*sizeof(int))
        cdef double * c_values = <double*> sage_malloc((n+1)*sizeof(double))
        cdef list indices = []
        cdef list values = []
        cdef int i,j

        i = glp_get_mat_row(self.lp, index + 1, c_indices, c_values)
        for 0 < j <= i:
            indices.append(c_indices[j]-1)
            values.append(c_values[j])

        sage_free(c_indices)
        sage_free(c_values)

        return (indices, values)

    cpdef get_row_bounds(self, int index):
        r"""
        Returns the bounds of a specific constraint.

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the constraint is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variables(5)
            5
            sage: p.add_constraint(range(5), range(5), 0, 2)
            sage: p.get_row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.get_row_bounds(0)
            (2.0, 2.0)
        """
        cdef double ub
        cdef double lb

        ub = glp_get_row_ub(self.lp, index + 1)
        lb = glp_get_row_lb(self.lp, index +1)

        return (
            (lb if lb != -DBL_MAX else None),
            (ub if ub != +DBL_MAX else None)
            )

    cpdef get_col_bounds(self, int index):
        r"""
        Returns the bounds of a specific variable.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        OUTPUT:

        A pair ``(lower_bound, upper_bound)``. Each of them can be set
        to ``None`` if the variable is not bounded in the
        corresponding direction, and is a real value otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variable()
            1
            sage: p.get_col_bounds(0)
            (0.0, None)
            sage: p.set_variable_max(0, 5)
            sage: p.get_col_bounds(0)
            (0.0, 5.0)
        """

        cdef double ub
        cdef double lb

        ub = glp_get_col_ub(self.lp, index +1)
        lb = glp_get_col_lb(self.lp, index +1)

        return (
            (lb if lb != -DBL_MAX else None),
            (ub if ub != +DBL_MAX else None)
            )

    cpdef double get_objective_coeff(self, int index):
        r"""
        Returns the coefficient of a variable in the objective
        function.

        INPUT:

        - ``index`` (integer) -- the variable's id.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variable()
            1
            sage: p.get_objective_coeff(0)
            0.0
            sage: p.set_objective_coeff(0,2)
            sage: p.get_objective_coeff(0)
            2.0
        """
        return glp_get_obj_coef(self.lp, index + 1)


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

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.n_cols()
            0
            sage: p.n_rows()
            0
            sage: p.add_constraints(5, -1, 0)
            sage: p.add_col(range(5), range(5))
            sage: p.n_rows()
            5
        """

        glp_add_cols(self.lp, 1)
        cdef int n = glp_get_num_cols(self.lp)

        cdef int * col_i
        cdef double * col_values

        col_i = <int *> sage_malloc((len(indices)+1) * sizeof(int))
        col_values = <double *> sage_malloc((len(indices)+1) * sizeof(double))

        for i,v in enumerate(indices):
            col_i[i+1] = v+1
        for i,v in enumerate(coeffs):
            col_values[i+1] = v

        glp_set_mat_col(self.lp, n, len(indices), col_i, col_values)
        glp_set_col_bnds(self.lp, n, GLP_LO, 0,0)


    cpdef int solve(self) except -1:
        r"""
        Solves the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_constraints(5, -1, 0)
            sage: p.add_col(range(5), range(5))
            sage: p.solve()
            0
            sage: p.set_objective_coeff(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """

        glp_intopt(self.lp, self.iocp)

        cdef int status = glp_mip_status(self.lp)
        if status == GLP_OPT:
            pass
        elif status == GLP_UNDEF:
            raise MIPSolverException("GLPK : Solution is undefined")
        elif status == GLP_FEAS:
            raise MIPSolverException("GLPK : Feasible solution found, while optimality has not been proven")
        elif status == GLP_NOFEAS:
            raise MIPSolverException("GLPK : There is no feasible integer solution to this Linear Program")

        return 0

    cpdef double get_objective_value(self):
        r"""
        Returns the value of the objective function.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variables(2)
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            7.5
            sage: p.get_variable_value(0)
            0.0
            sage: p.get_variable_value(1)
            1.5
        """
        return glp_mip_obj_val(self.lp)

    cpdef double get_variable_value(self, int variable):
        r"""
        Returns the value of a variable given by the solver.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variables(2)
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)
            sage: p.set_objective([2, 5])
            sage: p.solve()
            0
            sage: p.get_objective_value()
            7.5
            sage: p.get_variable_value(0)
            0.0
            sage: p.get_variable_value(1)
            1.5
        """
        return glp_mip_col_val(self.lp, variable+1)

    cpdef int n_cols(self):
        r"""
        Returns the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.n_cols()
            0
            sage: p.add_variables(2)
            2
            sage: p.n_cols()
            2
        """
        return glp_get_num_cols(self.lp)

    cpdef int n_rows(self):
        r"""
        Returns the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.n_rows()
            0
            sage: p.add_constraints(2, -1, 2)
            sage: p.n_rows()
            2
        """

        return glp_get_num_rows(self.lp)

    cpdef get_row_name(self, int index):
        r"""
        Returns the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variable()
            1
            sage: p.set_col_name(0, "I am a variable")
            sage: p.get_col_name(0)
            'I am a variable'
        """
        glp_create_index(self.lp)
        cdef char *  s = <char*> glp_get_row_name(self.lp, index + 1)

        if s != NULL:
            return s
        else:
            return ""


    cpdef set_col_name(self, int index, char * name):
        r"""
        Sets the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        - ``name`` (``char *``) -- its name

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variable()
            1
            sage: p.set_col_name(0, "I am a variable")
            sage: p.get_col_name(0)
            'I am a variable'
        """

        glp_set_col_name(self.lp, index + 1, name)

    cpdef get_col_name(self, int index):
        r"""
        Returns the ``index`` th variable name

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_constraints(1, -1, 2)
            sage: p.set_row_name(0, "Empty constraint 1")
            sage: p.get_row_name(0)
            'Empty constraint 1'
        """
        glp_create_index(self.lp)
        cdef char *  s = <char*> glp_get_col_name(self.lp, index + 1)

        if s != NULL:
            return s
        else:
            return ""

    cpdef set_row_name(self, int index, char * name):
        r"""
        Sets the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        - ``name`` (``char *``) -- its name

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_constraints(1, -1, 2)
            sage: p.set_row_name(0, "Empty constraint 1")
            sage: p.get_row_name(0)
            'Empty constraint 1'

        """

        glp_set_row_name(self.lp, index + 1, name)

    cpdef bint is_variable_binary(self, int index):
        r"""
        Tests whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.n_cols()
            0
            sage: p.add_variable()
            1
            sage: p.set_variable_type(0,0)
            sage: p.is_variable_binary(0)
            True

        """
        return glp_get_col_kind(self.lp, index + 1) == GLP_BV

    cpdef bint is_variable_integer(self, int index):
        r"""
        Tests whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.n_cols()
            0
            sage: p.add_variable()
            1
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True
        """
        return glp_get_col_kind(self.lp, index + 1) == GLP_IV

    cpdef bint is_variable_continuous(self, int index):
        r"""
        Tests whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.n_cols()
            0
            sage: p.add_variable()
            1
            sage: p.is_variable_continuous(0)
            True
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_continuous(0)
            False

        """
        return glp_get_col_kind(self.lp, index + 1) == GLP_CV

    cpdef bint is_maximization(self):
        r"""
        Tests whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.is_maximization()
            True
            sage: p.set_direction(-1)
            sage: p.is_maximization()
            False
        """

        return glp_get_obj_dir(self.lp) == GLP_MAX

    cpdef get_variable_max(self, int index):
        r"""
        Returns the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        OUTPUT:

        A real value if the variable has an upper bound, ``None``
        otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variable()
            1
            sage: p.get_variable_max(0) is None
            True
            sage: p.set_variable_max(0, 5)
            sage: p.get_variable_max(0)
            5.0

        """
        cdef double x = glp_get_col_ub(self.lp, index +1)
        if x == DBL_MAX:
            return None
        else:
            return x

    cpdef get_variable_min(self, int index):
        r"""
        Returns the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        OUTPUT:

        A real value if the variable has an lower bound, ``None``
        otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variable()
            1
            sage: p.get_variable_min(0)
            0.0
            sage: p.set_variable_min(0, 5)
            sage: p.get_variable_min(0)
            5.0
        """
        cdef double x = glp_get_col_lb(self.lp, index +1)
        if x == -DBL_MAX:
            return None
        else:
            return x

    cpdef set_variable_max(self, int index, value):
        r"""
        Sets the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variable()
            1
            sage: p.get_col_bounds(0)
            (0.0, None)
            sage: p.set_variable_max(0, 5)
            sage: p.get_col_bounds(0)
            (0.0, 5.0)
        """

        cdef double min = glp_get_col_lb(self.lp, index + 1)

        if value is None and min == -DBL_MAX:
            glp_set_col_bnds(self.lp, index + 1, GLP_FR, 0, 0)

        elif value is None:
            glp_set_col_bnds(self.lp, index + 1, GLP_LO, min, 0)

        elif min == -DBL_MAX:
            glp_set_col_bnds(self.lp, index + 1, GLP_UP, 0, value)

        elif min == value:
            glp_set_col_bnds(self.lp, index + 1, GLP_FX,  value, value)

        else:
            glp_set_col_bnds(self.lp, index + 1, GLP_DB, min, value)

    cpdef set_variable_min(self, int index, value):
        r"""
        Sets the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variable()
            1
            sage: p.get_col_bounds(0)
            (0.0, None)
            sage: p.set_variable_min(0, 5)
            sage: p.get_col_bounds(0)
            (5.0, None)
        """

        cdef double max = glp_get_col_ub(self.lp, index + 1)

        if value is None and max == DBL_MAX:
            glp_set_col_bnds(self.lp, index + 1, GLP_FR, 0.0, 0.0)

        elif value is None:
            glp_set_col_bnds(self.lp, index + 1, GLP_UP, 0.0, max)

        elif max == DBL_MAX:
            glp_set_col_bnds(self.lp, index + 1, GLP_LO, value, 0.0)

        elif max == value:
            glp_set_col_bnds(self.lp, index + 1, GLP_FX,  value, value)

        else:
            glp_set_col_bnds(self.lp, index + 1, GLP_DB, value, min)

    cpdef write_lp(self, char * filename):
        r"""
        Writes the problem to a .lp file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variables(2)
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)
            sage: p.set_objective([2, 5])
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp")
        """
        glp_write_lp(self.lp, NULL, filename)

    cpdef write_mps(self, char * filename, int modern):
        r"""
        Writes the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "GLPK")
            sage: p.add_variables(2)
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)
            sage: p.set_objective([2, 5])
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp")
        """
        glp_write_mps(self.lp, modern, NULL,  filename)

    def __dealloc__(self):
        r"""
        Destructor
        """
        glp_delete_prob(self.lp)
