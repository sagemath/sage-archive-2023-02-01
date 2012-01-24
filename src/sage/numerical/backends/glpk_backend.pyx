"""
GLPK Backend

AUTHORS:

- Nathann Cohen (2010-10): initial implementation
- John Perry (2012-01): glp_simplex preprocessing
"""

##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.numerical.mip import MIPSolverException

include "../../ext/interrupt.pxi"

cdef class GLPKBackend(GenericBackend):

    def __cinit__(self, maximization = True):
        """
        Constructor

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver="GLPK")
        """
        self.lp = glp_create_prob()
        self.preprocessing = 0
        self.smcp = <c_glp_smcp* > sage_malloc(sizeof(c_glp_smcp))
        glp_init_smcp(self.smcp)
        self.iocp = <c_glp_iocp* > sage_malloc(sizeof(c_glp_iocp))
        glp_init_iocp(self.iocp)
        self.iocp.presolve = GLP_ON
        self.set_verbosity(0)

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

    cpdef int add_variable(self, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, name=None) except -1:
        """
        Add a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive, real and the coefficient in the
        objective function is 0.0.

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.ncols()
            1
            sage: p.add_variable(binary=True)
            1
            sage: p.add_variable(lower_bound=-2.0, integer=True)
            2
            sage: p.add_variable(continuous=True, integer=True)
            Traceback (most recent call last):
            ...
            ValueError: ...
            sage: p.add_variable(name='x',obj=1.0)
            3
            sage: p.col_name(3)
            'x'
            sage: p.objective_coefficient(3)
            1.0
        """
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if  vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        glp_add_cols(self.lp, 1)
        cdef int n_var = glp_get_num_cols(self.lp)

        self.variable_lower_bound(n_var - 1, lower_bound)
        self.variable_upper_bound(n_var - 1, upper_bound)

        if continuous:
            glp_set_col_kind(self.lp, n_var, GLP_CV)
        elif binary:
            glp_set_col_kind(self.lp, n_var, GLP_BV)
        elif integer:
            glp_set_col_kind(self.lp, n_var, GLP_IV)

        if name is not None:
            glp_set_col_name(self.lp, n_var, name)

        if obj:
            self.objective_coefficient(n_var - 1, obj)

        return n_var - 1

    cpdef int add_variables(self, int number, lower_bound=0.0, upper_bound=None, binary=False, continuous=False, integer=False, obj=0.0, names=None) except -1:
        """
        Add ``number`` new variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive, real and theor coefficient in
        the objective function is 0.0.

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variables(5)
            4
            sage: p.ncols()
            5
            sage: p.add_variables(2, lower_bound=-2.0, integer=True, names=['a','b'])
            6
        """
        cdef int vtype = int(bool(binary)) + int(bool(continuous)) + int(bool(integer))
        if  vtype == 0:
            continuous = True
        elif vtype != 1:
            raise ValueError("Exactly one parameter of 'binary', 'integer' and 'continuous' must be 'True'.")

        glp_add_cols(self.lp, number)

        cdef int n_var
        n_var = glp_get_num_cols(self.lp)

        cdef int i

        for 0<= i < number:
            self.variable_lower_bound(n_var - i - 1, lower_bound)
            self.variable_upper_bound(n_var - i - 1, upper_bound)
            if continuous:
                glp_set_col_kind(self.lp, n_var - i, GLP_CV)
            elif binary:
                glp_set_col_kind(self.lp, n_var - i, GLP_BV)
            elif integer:
                glp_set_col_kind(self.lp, n_var - i, GLP_IV)

            if obj:
                self.objective_coefficient(n_var - i - 1, obj)

            if names is not None:
                glp_set_col_name(self.lp, n_var, names[number - i - 1])

        return n_var - 1

    cpdef set_variable_type(self, int variable, int vtype):
        """
        Set the type of a variable

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``vtype`` (integer) :

            *  1  Integer
            *  0  Binary
            * -1 Real

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
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

    cpdef set_sense(self, int sense):
        """
        Set the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """
        if sense == 1:
            glp_set_obj_dir(self.lp, GLP_MAX)
        else:
            glp_set_obj_dir(self.lp, GLP_MIN)

    cpdef objective_coefficient(self, int variable, coeff=None):
        """
        Set or get the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient or ``None`` for
          reading (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.objective_coefficient(0)
            0.0
            sage: p.objective_coefficient(0,2)
            sage: p.objective_coefficient(0)
            2.0
        """
        if coeff is None:
            return glp_get_obj_coef(self.lp, variable + 1)
        else:
            glp_set_obj_coef(self.lp, variable + 1, coeff)

    cpdef problem_name(self, char * name = NULL):
        """
        Return or define the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.problem_name("There once was a french fry")
            sage: print p.problem_name()
            There once was a french fry
        """
        cdef char * n

        if name == NULL:
            n =  <char *> glp_get_prob_name(self.lp)
            if n == NULL:
                return ""
            else:
                return n

        else:
            glp_set_prob_name(self.lp, name)

    cpdef set_objective(self, list coeff):
        """
        Set the objective function.

        INPUT:

        - ``coeff`` - a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.set_objective([1, 1, 2, 1, 3])
            sage: map(lambda x :p.objective_coefficient(x), range(5))
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """
        cdef int i

        for i,v in enumerate(coeff):
            glp_set_obj_coef(self.lp, i+1, v)

    cpdef set_verbosity(self, int level):
        """
        Set the verbosity level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.set_verbosity(2)

        """
        if level == 0:
            self.iocp.msg_lev = GLP_MSG_OFF
            self.smcp.msg_lev = GLP_MSG_OFF
        elif level == 1:
            self.iocp.msg_lev = GLP_MSG_ERR
            self.smcp.msg_lev = GLP_MSG_ERR
        elif level == 2:
            self.iocp.msg_lev = GLP_MSG_ON
            self.smcp.msg_lev = GLP_MSG_ON
        else:
            self.iocp.msg_lev = GLP_MSG_ALL
            self.smcp.msg_lev = GLP_MSG_ALL

    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=None):
        """
        Add a linear constraint.

        INPUT:

        - ``coefficients`` an iterable with ``(c,v)`` pairs where ``c``
          is a variable index (integer) and ``v`` is a value (real
          value).

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``name`` - an optional name for this row (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint( zip(range(5), range(5)), 2.0, 2.0)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
            (2.0, 2.0)
            sage: p.add_linear_constraint( zip(range(5), range(5)), 1.0, 1.0, name='foo')
            sage: p.row_name(1)
            'foo'
        """
        if lower_bound is None and upper_bound is None:
            raise ValueError("At least one of 'upper_bound' or 'lower_bound' must be set.")

        glp_add_rows(self.lp, 1)
        cdef int n = glp_get_num_rows(self.lp)

        cdef int * row_i
        cdef double * row_values

        row_i = <int *> sage_malloc((len(coefficients)+1) * sizeof(int))
        row_values = <double *> sage_malloc((len(coefficients)+1) * sizeof(double))

        for i,(c,v) in enumerate(coefficients):
            row_i[i+1] = c+1
            row_values[i+1] = v

        glp_set_mat_row(self.lp, n, len(coefficients), row_i, row_values)

        if upper_bound is not None and lower_bound is None:
            glp_set_row_bnds(self.lp, n, GLP_UP, upper_bound, upper_bound)
        elif lower_bound is not None and upper_bound is None:
            glp_set_row_bnds(self.lp, n, GLP_LO, lower_bound, lower_bound)
        elif upper_bound is not None and lower_bound is not None:
            if lower_bound == upper_bound:
                glp_set_row_bnds(self.lp, n, GLP_FX, lower_bound, upper_bound)
            else:
                glp_set_row_bnds(self.lp, n, GLP_DB, lower_bound, upper_bound)

        if name is not None:
            glp_set_row_name(self.lp, n, name)

        sage_free(row_i)
        sage_free(row_values)

    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=None):
        """
        Add ``'number`` linear constraints.

        INPUT:

        - ``number`` (integer) -- the number of constraints to add.

        - ``lower_bound`` - a lower bound, either a real value or ``None``

        - ``upper_bound`` - an upper bound, either a real value or ``None``

        - ``names`` - an optional list of names (default: ``None``)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraints(5, None, 2)
            sage: p.row(4)
            ([], [])
            sage: p.row_bounds(4)
            (None, 2.0)
            sage: p.add_linear_constraints(2, None, 2, names=['foo','bar'])
        """
        if lower_bound is None and upper_bound is None:
            raise ValueError("At least one of 'upper_bound' or 'lower_bound' must be set.")

        glp_add_rows(self.lp, number)
        cdef int n = glp_get_num_rows(self.lp)

        cdef int i
        for 0<= i < number:
            if upper_bound is not None and lower_bound is None:
                glp_set_row_bnds(self.lp, n-i, GLP_UP, upper_bound, upper_bound)
            elif lower_bound is not None and upper_bound is None:
                glp_set_row_bnds(self.lp, n-i, GLP_LO, lower_bound, lower_bound)
            elif upper_bound is not None and lower_bound is not None:
                if lower_bound == upper_bound:
                    glp_set_row_bnds(self.lp, n-i, GLP_FX, lower_bound, upper_bound)
                else:
                    glp_set_row_bnds(self.lp, n-i, GLP_DB, lower_bound, upper_bound)
            if names is not None:
                glp_set_row_name(self.lp, n-i, names[number-i-1])

    cpdef row(self, int index):
        r"""
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
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
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
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(5)
            4
            sage: p.add_linear_constraint(zip(range(5), range(5)), 2, 2)
            sage: p.row(0)
            ([4, 3, 2, 1], [4.0, 3.0, 2.0, 1.0])
            sage: p.row_bounds(0)
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
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
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

    cpdef add_col(self, list indices, list coeffs):
        """
        Add a column.

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
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.nrows()
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
        """
        Solve the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_linear_constraints(5, 0, None)
            sage: p.add_col(range(5), range(5))
            sage: p.solve()
            0
            sage: p.objective_coefficient(0,1)
            sage: p.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: ...

        .. WARNING::

            Sage uses GLPK's ``glp_intopt`` to find solutions.
            This routine sometimes FAILS CATASTROPHICALLY
            when given a system it cannot solve. (Ticket #12309.)
            Here, "catastrophic" can mean either "infinite loop" or
            segmentation fault. Upstream considers this behavior
            "essentially innate" to their design, and suggests
            preprocessing it with ``glpk_simplex`` first.
            Thus, if you suspect that your system is infeasible,
            set the ``preprocessing`` option first.

        EXAMPLE::

            sage: lp = MixedIntegerLinearProgram(solver = "GLPK")
            sage: lp.add_constraint( lp[1] +lp[2] -2.0 *lp[3] , max =-1.0)
            sage: lp.add_constraint( lp[0] -4.0/3 *lp[1] +1.0/3 *lp[2] , max =-1.0/3)
            sage: lp.add_constraint( lp[0] +0.5 *lp[1] -0.5 *lp[2] +0.25 *lp[3] , max =-0.25)
            sage: lp.solve()
            0.0
            sage: lp.add_constraint( lp[0] +4.0 *lp[1] -lp[2] +lp[3] , max =-1.0)
            sage: lp.solve()
            Traceback (most recent call last):
            ...
            RuntimeError: GLPK : Signal sent, try preprocessing option
            sage: lp.solver_parameter("preprocessing", True)
            sage: lp.solve()
            Traceback (most recent call last):
            ...
            MIPSolverException: 'GLPK : Simplex cannot find a feasible solution'
        """

        cdef int status
        if self.preprocessing:
            status = glp_simplex(self.lp, self.smcp)
            status = glp_get_prim_stat(self.lp)
            if status == GLP_OPT or status == GLP_FEAS:
                pass
            elif status == GLP_UNDEF or status == GLP_NOFEAS:
                raise MIPSolverException("GLPK : Simplex cannot find a feasible solution")

        sig_str('GLPK : Signal sent, try preprocessing option')
        glp_intopt(self.lp, self.iocp)
        sig_off()

        status = glp_mip_status(self.lp)
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
        """
        Returns the value of the objective function.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
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
        """
        Returns the value of a variable given by the solver.

        .. NOTE::

           Behaviour is undefined unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
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

    cpdef int ncols(self):
        """
        Return the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variables(2)
            1
            sage: p.ncols()
            2
        """
        return glp_get_num_cols(self.lp)

    cpdef int nrows(self):
        """
        Return the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.nrows()
            0
            sage: p.add_linear_constraints(2, 2, None)
            sage: p.nrows()
            2
        """

        return glp_get_num_rows(self.lp)

    cpdef col_name(self, int index):
        """
        Return the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable(name='I am a variable')
            0
            sage: p.col_name(0)
            'I am a variable'
        """
        cdef char * s

        glp_create_index(self.lp)
        s = <char*> glp_get_col_name(self.lp, index + 1)

        if s != NULL:
            return s
        else:
            return ""

    cpdef row_name(self, int index):
        """
        Return the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_linear_constraints(1, 2, None, names=['Empty constraint 1'])
            sage: p.row_name(0)
            'Empty constraint 1'
        """
        cdef char *  s

        glp_create_index(self.lp)
        s = <char*> glp_get_row_name(self.lp, index + 1)

        if s != NULL:
            return s
        else:
            return ""

    cpdef bint is_variable_binary(self, int index):
        """
        Test whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,0)
            sage: p.is_variable_binary(0)
            True

        """
        return glp_get_col_kind(self.lp, index + 1) == GLP_BV

    cpdef bint is_variable_integer(self, int index):
        """
        Test whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_integer(0)
            True
        """
        return glp_get_col_kind(self.lp, index + 1) == GLP_IV

    cpdef bint is_variable_continuous(self, int index):
        """
        Test whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.ncols()
            0
            sage: p.add_variable()
            0
            sage: p.is_variable_continuous(0)
            True
            sage: p.set_variable_type(0,1)
            sage: p.is_variable_continuous(0)
            False

        """
        return glp_get_col_kind(self.lp, index + 1) == GLP_CV

    cpdef bint is_maximization(self):
        """
        Test whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.is_maximization()
            True
            sage: p.set_sense(-1)
            sage: p.is_maximization()
            False
        """

        return glp_get_obj_dir(self.lp) == GLP_MAX

    cpdef variable_upper_bound(self, int index, value = False):
        """
        Return or define the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_upper_bound(0, 5)
            sage: p.col_bounds(0)
            (0.0, 5.0)
        """
        cdef double x
        cdef double min

        if value == False:
            x = glp_get_col_ub(self.lp, index +1)
            if x == DBL_MAX:
                return None
            else:
                return x

        else:
            min = glp_get_col_lb(self.lp, index + 1)

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


    cpdef variable_lower_bound(self, int index, value = False):
        """
        Return or define the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variable()
            0
            sage: p.col_bounds(0)
            (0.0, None)
            sage: p.variable_lower_bound(0, 5)
            sage: p.col_bounds(0)
            (5.0, None)
        """
        cdef double x
        cdef double max

        if value == False:
            x = glp_get_col_lb(self.lp, index +1)
            if x == -DBL_MAX:
                return None
            else:
                return x
        else:
            max = glp_get_col_ub(self.lp, index + 1)

            if value is None and max == DBL_MAX:
                glp_set_col_bnds(self.lp, index + 1, GLP_FR, 0.0, 0.0)

            elif value is None:
                glp_set_col_bnds(self.lp, index + 1, GLP_UP, 0.0, max)

            elif max == DBL_MAX:
                glp_set_col_bnds(self.lp, index + 1, GLP_LO, value, 0.0)

            elif max == value:
                glp_set_col_bnds(self.lp, index + 1, GLP_FX,  value, value)

            else:
                glp_set_col_bnds(self.lp, index + 1, GLP_DB, value, max)

    cpdef write_lp(self, char * filename):
        """
        Write the problem to a .lp file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp")
        """
        glp_write_lp(self.lp, NULL, filename)

    cpdef write_mps(self, char * filename, int modern):
        """
        Write the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.add_variables(2)
            1
            sage: p.add_linear_constraint([[0, 1], [1, 2]], None, 3)
            sage: p.set_objective([2, 5])
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp")
        """
        glp_write_mps(self.lp, modern, NULL,  filename)

    cpdef GLPKBackend copy(self):
        """
        Returns a copy of self.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = MixedIntegerLinearProgram(solver = "GLPK")
            sage: b = p.new_variable()
            sage: p.add_constraint(b[1] + b[2] <= 6)
            sage: p.set_objective(b[1] + b[2])
            sage: copy(p).solve()
            6.0
        """
        cdef GLPKBackend p = GLPKBackend(maximization = (1 if self.is_maximization() else -1))
        p.preprocessing = self.preprocessing
        p.iocp.tm_lim = self.iocp.tm_lim
        glp_copy_prob(p.lp, self.lp, 1)
        return p


    cpdef solver_parameter(self, name, value = None):
        """
        Return or define a solver parameter

        INPUT:

        - ``name`` (string) -- the parameter

        - ``value`` -- the parameter's value if it is to be defined,
          or ``None`` (default) to obtain its current value.

        .. NOTE::

           The list of available parameters is available at
           :meth:`sage.numerical.mip.MixedIntegerlinearProgram.solver_parameter`

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "GLPK")
            sage: p.solver_parameter("timelimit", 60)
            sage: p.solver_parameter("timelimit")
            60.0
        """
        if name == "timelimit":
            if value == None:
                return self.iocp.tm_lim/1000.0
            else:
                self.iocp.tm_lim = 1000 * value

        elif name == "preprocessing":
            if value == None:
                return self.preprocessing
            else:
                if value == True: value = 1
                elif value == False: value = 0
                self.preprocessing = value

        else:
            raise ValueError("This parameter is not available.")

    def __dealloc__(self):
        """
        Destructor
        """
        glp_delete_prob(self.lp)
        sage_free(self.iocp)
        sage_free(self.smcp)
