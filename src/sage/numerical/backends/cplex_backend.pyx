"""
CPLEX Backend

AUTHORS:

- Nathann Cohen (2010-10): initial implementation

"""

##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################


from sage.numerical.mip import MIPSolverException

cdef class CPLEXBackend(GenericBackend):

    def __cinit__(self, maximization = True):

        cdef int status
        self.env = CPXopenCPLEX (&status)
        check(status)

        cdef char * tmp = ""
        self.lp = CPXcreateprob (self.env, &status, tmp);
        check(status)

        if maximization:
            self.set_sense(+1)
        else:
            self.set_sense(-1)

    cpdef int add_variable(self):
        r"""
        Adds a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")                    # optional - CPLEX
            sage: p.ncols()                                         # optional - CPLEX
            0
            sage: p.add_variable()                                   # optional - CPLEX
            1
            sage: p.ncols()                                         # optional - CPLEX
            1
        """

        cdef int status
        status = CPXnewcols(self.env, self.lp, 1, NULL, NULL, NULL, NULL, NULL)
        check(status)
        return CPXgetnumcols(self.env, self.lp)

    cpdef int add_variables(self, int number):
        r"""
        Adds ``number`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")                    # optional - CPLEX
            sage: p.ncols()                                         # optional - CPLEX
            0
            sage: p.add_variables(5)                                 # optional - CPLEX
            5
            sage: p.ncols()                                         # optional - CPLEX
            5
        """

        cdef int status
        status = CPXnewcols(self.env, self.lp, number, NULL, NULL, NULL, NULL, NULL)
        check(status)
        return CPXgetnumcols(self.env, self.lp)

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

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")   # optional - CPLEX
            sage: p.ncols()                                        # optional - CPLEX
            0
            sage: p.add_variable()                                  # optional - CPLEX
            1
            sage: p.set_variable_type(0,1)                          # optional - CPLEX
            sage: p.is_variable_integer(0)                          # optional - CPLEX
            True
        """

        cdef int status

        cdef char type
        if vtype == 1:
            type = 'I'
        elif vtype == 0:
            type = 'B'
        else:
            type = 'C'

        status = CPXchgctype(self.env, self.lp, 1, &variable, &type)
        check(status)

    cpdef set_sense(self, int sense):
        r"""
        Sets the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.is_maximization()                              # optional - CPLEX
            True
            sage: p.set_sense(-1)                              # optional - CPLEX
            sage: p.is_maximization()                              # optional - CPLEX
            False
        """

        CPXchgobjsen(self.env, self.lp, -sense)

    cpdef set_objective_coefficient(self, int variable, double coeff):
        r"""
        Sets the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variable()                                 # optional - CPLEX
            1
            sage: p.get_objective_coefficient(0)                         # optional - CPLEX
            0.0
            sage: p.set_objective_coefficient(0,2)                       # optional - CPLEX
            sage: p.get_objective_coefficient(0)                         # optional - CPLEX
            2.0
        """

        cdef int status
        status = CPXchgobj(self.env, self.lp, 1, &variable, &coeff)
        check(status)

    cpdef problem_name(self, char * name = NULL):
        r"""
        Returns or defines the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name. When set to
          ``NULL`` (default), the method returns the problem's name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")   # optional - CPLEX
            sage: p.problem_name("There once was a french fry") # optional - CPLEX
            sage: print p.problem_name()                        # optional - CPLEX
            There once was a french fry
        """

        cdef int status
        cdef int zero
        cdef char * n
        if name == NULL:

            n = <char*> sage_malloc(500*sizeof(char))
            status = CPXgetprobname(self.env, self.lp, n, 500, &zero)
            check(status)
            s = str(n)
            sage_free(n)
            return s


        else:
            status = CPXchgprobname(self.env, self.lp, name)
            check(status)


    cpdef set_objective(self, list coeff):
        r"""
        Sets the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")    # optional - CPLEX
            sage: p.add_variables(5)                                 # optional - CPLEX
            5
            sage: p.set_objective([1, 1, 2, 1, 3])                   # optional - CPLEX
            sage: map(lambda x :p.get_objective_coefficient(x), range(5))  # optional - CPLEX
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """

        cdef int status
        cdef int n = self.ncols()
        cdef double * c_coeff = <double *> sage_malloc(n * sizeof(double))
        cdef int * c_indices = <int *> sage_malloc(n * sizeof(int))

        for i,v in enumerate(coeff):
            c_coeff[i] = v
            c_indices[i] = i

        status = CPXchgobj(self.env, self.lp, n, c_indices, c_coeff)
        check(status)


    cpdef set_verbosity(self, int level):
        r"""
        Sets the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")   # optional - CPLEX
            sage: p.set_verbosity(2)                                # optional - CPLEX

        """

        cdef int status
        if level == 0:
            status = CPXsetintparam (self.env, CPX_PARAM_SCRIND, 0)
            check(status)
        else:
            status = CPXsetintparam (self.env, CPX_PARAM_SCRIND, CPX_ON)
            check(status)

    cpdef add_linear_constraints(self, int number, int direction, double bound):
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
            sage: p = get_solver(solver = "CPLEX")   # optional - CPLEX
            sage: p.add_variables(5)                                # optional - CPLEX
            5
            sage: p.add_linear_constraints(5, +1, 2)                       # optional - CPLEX
            sage: p.row(4)                                      # optional - CPLEX
            ([], [])
            sage: p.row_bounds(4)                               # optional - CPLEX
            (None, 2.0)
        """

        cdef char * c_types = <char *> sage_malloc(number * sizeof(char))
        cdef double * c_bounds = <double *> sage_malloc(number * sizeof(double))
        cdef int i
        cdef int status

        c_bounds[0] = bound

        if direction == +1:
            c_types[0] = 'L'
        elif direction == -1:
            c_types[0] = 'G'
        else:
            c_types[0] = 'E'

        for 1<= i < number:
            c_types[i] = c_types[0]
            c_bounds[i] = c_bounds[0]


        status = CPXnewrows(self.env, self.lp, number, c_bounds, c_types, NULL, NULL)
        check(status)

    cpdef add_linear_constraint(self, list indices, list coeffs, int direction, double bound):
        r"""
        Adds a linear constraint.

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
            sage: p = get_solver(solver = "CPLEX") # optional - CPLEX
            sage: p.add_variables(5)                              # optional - CPLEX
            5
            sage: p.add_linear_constraint(range(5), range(5), 0, 2)      # optional - CPLEX
            sage: p.row(0)                                    # optional - CPLEX
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                             # optional - CPLEX
            (2.0, 2.0)
        """

        cdef int status
        cdef int i
        cdef int n = len(indices)
        cdef int nrows = self.nrows()
        cdef char sense

        if direction == 1:
            sense = 'L'
        elif direction == -1:
            sense = 'G'
        else:
            sense = 'E'

        status = CPXnewrows(self.env, self.lp, 1, &bound, &sense, NULL, NULL)
        check(status)

        cdef double * c_coeff = <double *> sage_malloc(n * sizeof(double))
        cdef int * c_indices = <int *> sage_malloc(n * sizeof(int))
        cdef int * c_row = <int *> sage_malloc(n * sizeof(int))

        for 0<= i < n:
            c_coeff[i] = coeffs[i]
            c_indices[i] = indices[i]
            c_row[i] = nrows


        status = CPXchgcoeflist(self.env, self.lp, n, c_row, c_indices, c_coeff)
        check(status)

    cpdef row(self, int index):
        r"""
        Returns a row

        INPUT:

        - ``index`` (integer) -- the constraint's id.

        OUTPUT:

        A pair ``(indices, coeffs)`` where ``indices`` lists the
        entries whose coefficient is nonzero, and to which ``coeffs``
        associates their coefficient on the model of the
        ``add_linear_constraint`` method.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(5)                               # optional - CPLEX
            5
            sage: p.add_linear_constraint(range(5), range(5), 0, 2)       # optional - CPLEX
            sage: p.row(0)                                     # optional - CPLEX
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                              # optional - CPLEX
            (2.0, 2.0)
        """

        cdef int status
        cdef int n,i
        cdef int zero
        cdef list indices = []
        cdef list values = []

        cdef double * c_coeff = <double *> sage_malloc((self.ncols()+10) * sizeof(double))
        cdef int * c_indices = <int *> sage_malloc((self.ncols()+10) * sizeof(int))

        status = CPXgetrows(self.env, self.lp, &n, &zero, c_indices, c_coeff, self.ncols()+3, &zero, index, index)

        check(status)

        for 0<= i<n:
            indices.append(c_indices[i])
            values.append(c_coeff[i])

        sage_free(c_coeff)
        sage_free(c_indices)

        return (indices, values)

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
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(5)                               # optional - CPLEX
            5
            sage: p.add_linear_constraint(range(5), range(5), 0, 2)       # optional - CPLEX
            sage: p.row(0)                                     # optional - CPLEX
            ([1, 2, 3, 4], [1.0, 2.0, 3.0, 4.0])
            sage: p.row_bounds(0)                              # optional - CPLEX
            (2.0, 2.0)
        """

        cdef int status
        cdef double value
        status = CPXgetrhs(self.env, self.lp, &value, index, index)
        check(status)

        cdef char direction
        status = CPXgetsense(self.env, self.lp, &direction, index, index)
        check(status)

        if direction == 'L':
            return (None, value)
        elif direction == 'G':
            return (value, None)
        else:
            return (value, value)

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
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variable()                                 # optional - CPLEX
            1
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, None)
            sage: p.variable_max(0, 5)                         # optional - CPLEX
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, 5.0)
        """

        cdef int status
        cdef double ub
        cdef double lb

        status = CPXgetub(self.env, self.lp, &ub, index, index)
        check(status)

        status = CPXgetlb(self.env, self.lp, &lb, index, index)
        check(status)

        return (lb if lb != -CPX_INFBOUND else None,
                ub if ub != +CPX_INFBOUND else None)

    cpdef double get_objective_coefficient(self, int index):
        r"""
        Sets the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variable()                                 # optional - CPLEX
            1
            sage: p.get_objective_coefficient(0)                         # optional - CPLEX
            0.0
            sage: p.set_objective_coefficient(0,2)                       # optional - CPLEX
            sage: p.get_objective_coefficient(0)                         # optional - CPLEX
            2.0
        """

        cdef int status
        cdef double value
        status = CPXgetobj(self.env, self.lp, &value, index, index)
        check(status)
        return value


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
            sage: p = get_solver(solver = "CPLEX")                  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.nrows()                                       # optional - CPLEX
            0
            sage: p.add_linear_constraints(5, -1, 0)                      # optional - CPLEX
            sage: p.add_col(range(5), range(5))                    # optional - CPLEX
            sage: p.nrows()                                       # optional - CPLEX
            5
        """

        cdef int status
        cdef int i
        cdef int n = len(indices)
        cdef int ncols = self.ncols()

        status = CPXnewcols(self.env, self.lp, 1, NULL, NULL, NULL, NULL, NULL)


        check(status)

        cdef double * c_coeff = <double *> sage_malloc(n * sizeof(double))
        cdef int * c_indices = <int *> sage_malloc(n * sizeof(int))
        cdef int * c_col = <int *> sage_malloc(n * sizeof(int))

        for 0<= i < n:
            c_coeff[i] = coeffs[i]
            c_indices[i] = indices[i]
            c_col[i] = ncols


        status = CPXchgcoeflist(self.env, self.lp, n, c_indices, c_col, c_coeff)
        check(status)

    cpdef int solve(self) except -1:
        r"""
        Solves the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX") # optional - CPLEX
            sage: p.add_linear_constraints(5, -1, 0)                     # optional - CPLEX
            sage: p.add_col(range(5), range(5))                   # optional - CPLEX
            sage: p.solve()                                       # optional - CPLEX
            0
            sage: p.set_objective_coefficient(0,1)                      # optional - CPLEX
            sage: p.solve()                                       # optional - CPLEX
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """

        cdef int status
        cdef int ptype
        cdef int solnmethod_p, solntype_p, pfeasind_p, dfeasind_p

        ptype = CPXgetprobtype(self.env, self.lp)

        if ptype == 1:
            status = CPXmipopt(self.env, self.lp)
        elif ptype == 0:
            status = CPXlpopt(self.env, self.lp)
        else:
            raise MIPSolverException("CPLEX: Unknown problem type")

        check(status)

        status = CPXsolninfo(self.env, self.lp, &solnmethod_p, &solntype_p, &pfeasind_p, &dfeasind_p)
        check(status)

        if not pfeasind_p:
            raise MIPSolverException("CPLEX: The primal has no feasible solution")
        if not dfeasind_p:
            raise MIPSolverException("CPLEX: The problem is unbounded")


        return 0


    cpdef double get_objective_value(self):
        r"""
        Returns the value of the objective function.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(2)                               # optional - CPLEX
            2
            sage: p.add_linear_constraint([0, 1], [1, 2], +1, 3)          # optional - CPLEX
            sage: p.set_objective([2, 5])                          # optional - CPLEX
            sage: p.solve()                                        # optional - CPLEX
            0
            sage: p.get_objective_value()                          # optional - CPLEX
            7.5
            sage: p.get_variable_value(0)                          # optional - CPLEX
            0.0
            sage: p.get_variable_value(1)                          # optional - CPLEX
            1.5
        """

        cdef int status
        cdef double value
        status = CPXgetobjval (self.env, self.lp, &value)
        check(status)

        return value


    cpdef double get_variable_value(self, int variable):
        r"""
        Returns the value of a variable given by the solver.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX") # optional - CPLEX
            sage: p.add_variables(2)                              # optional - CPLEX
            2
            sage: p.add_linear_constraint([0, 1], [1, 2], +1, 3)         # optional - CPLEX
            sage: p.set_objective([2, 5])                         # optional - CPLEX
            sage: p.solve()                                       # optional - CPLEX
            0
            sage: p.get_objective_value()                         # optional - CPLEX
            7.5
            sage: p.get_variable_value(0)                         # optional - CPLEX
            0.0
            sage: p.get_variable_value(1)                         # optional - CPLEX
            1.5
        """

        cdef int status
        cdef int zero
        cdef double value
        status = CPXgetx(self.env, self.lp, &value, variable, variable)
        check(status)

        return value


    cpdef int ncols(self):
        r"""
        Returns the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.add_variables(2)                               # optional - CPLEX
            2
            sage: p.ncols()                                       # optional - CPLEX
            2
        """

        return CPXgetnumcols(self.env, self.lp)

    cpdef int nrows(self):
        r"""
        Returns the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX") # optional - CPLEX
            sage: p.nrows()                                      # optional - CPLEX
            0
            sage: p.add_linear_constraints(2, -1, 2)                     # optional - CPLEX
            sage: p.nrows()                                      # optional - CPLEX
            2
        """

        return CPXgetnumrows(self.env, self.lp)

    cpdef row_name(self, int index, char * name = NULL):
        r"""
        Returns or defines the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_linear_constraints(1, -1, 2)                      # optional - CPLEX
            sage: p.row_name(0, "Empty constraint 1")          # optional - CPLEX
            sage: p.row_name(0)                                # optional - CPLEX
            'Empty constraint 1'

        """

        cdef int status
        cdef int zero
        cdef char * n

        if name == NULL:
            n = <char *>sage_malloc(500*sizeof(char))
            status = CPXgetrowname(self.env, self.lp, &n, n, 500, &zero, index, index)
            if status == 1219:
                sage_free(n)
                return ""
            check(status)

            s = str(n)
            sage_free(n)

            return s

            pass
        else:
            status = CPXchgrowname(self.env, self.lp, 1, &index, &name)
            check(status)

    cpdef col_name(self, int index, char * name = NULL):
        r"""
        Returns or defines the ``index`` th col name.

        INPUT:

        - ``index`` (integer) -- the col's id

        - ``name`` (``char *``) -- its name. When set to ``NULL``
          (default), the method returns the current name.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variable()                                 # optional - CPLEX
            1
            sage: p.col_name(0, "I am a variable")             # optional - CPLEX
            sage: p.col_name(0)                                # optional - CPLEX
            'I am a variable'
        """

        cdef int status
        cdef char * n
        cdef int zero

        if name == NULL:

            n = <char *>sage_malloc(500*sizeof(char))
            status = CPXgetcolname(self.env, self.lp, &n, n, 500, &zero, index, index)
            if status == 1219:
                sage_free(n)
                return ""
            check(status)

            s = str(n)
            sage_free(n)
            return s

        else:
            status = CPXchgcolname(self.env, self.lp, 1, &index, &name)
            check(status)



    cpdef bint is_variable_binary(self, int index):
        r"""
        Tests whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.add_variable()                                 # optional - CPLEX
            1
            sage: p.set_variable_type(0,0)                         # optional - CPLEX
            sage: p.is_variable_binary(0)                          # optional - CPLEX
            True

        """

        cdef int status
        cdef char ctype

        status = CPXgetctype(self.env, self.lp, &ctype, index, index)

        # status = 3003 when the problem is a LP and not a MILP
        if status == 3003:
            return False

        check(status)

        return ctype == 'B'


    cpdef bint is_variable_integer(self, int index):
        r"""
        Tests whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.add_variable()                                 # optional - CPLEX
            1
            sage: p.set_variable_type(0,1)                         # optional - CPLEX
            sage: p.is_variable_integer(0)                         # optional - CPLEX
            True
        """

        cdef int status
        cdef char ctype

        status = CPXgetctype(self.env, self.lp, &ctype, index, index)

        # status = 3003 when the problem is a LP and not a MILP
        if status == 3003:
            return False

        check(status)

        return ctype == 'I'


    cpdef bint is_variable_continuous(self, int index):
        r"""
        Tests whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.ncols()                                       # optional - CPLEX
            0
            sage: p.add_variable()                                 # optional - CPLEX
            1
            sage: p.is_variable_continuous(0)                      # optional - CPLEX
            True
            sage: p.set_variable_type(0,1)                         # optional - CPLEX
            sage: p.is_variable_continuous(0)                      # optional - CPLEX
            False

        """

        cdef int status
        cdef char ctype

        status = CPXgetctype(self.env, self.lp, &ctype, index, index)

        # status = 3003 when the problem is a LP and not a MILP
        if status == 3003:
            return True

        check(status)

        return ctype == 'C'


    cpdef bint is_maximization(self):
        r"""
        Tests whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX") # optional - CPLEX
            sage: p.is_maximization()                             # optional - CPLEX
            True
            sage: p.set_sense(-1)                             # optional - CPLEX
            sage: p.is_maximization()                             # optional - CPLEX
            False
        """

        return -1 == CPXgetobjsen(self.env, self.lp)

    cpdef variable_max(self, int index, value = False):
        r"""
        Returns or defines the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variable()                                 # optional - CPLEX
            1
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, None)
            sage: p.variable_max(0, 5)                         # optional - CPLEX
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, 5.0)
        """
        cdef int status
        cdef double ub
        cdef char x
        cdef double c_value

        if value == False:

            status = CPXgetub(self.env, self.lp, &ub, index, index)
            check(status)

            return ub if ub != CPX_INFBOUND else None

        else:

            x = 'U'
            c_value = value if value is not None else +CPX_INFBOUND
            status = CPXchgbds(self.env, self.lp, 1, &index, &x, &c_value)
            check(status)

    cpdef variable_min(self, int index, value = False):
        r"""
        Returns or defines the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound. When set to ``False``
          (default), the method returns the current value.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variable()                                 # optional - CPLEX
            1
            sage: p.col_bounds(0)                              # optional - CPLEX
            (0.0, None)
            sage: p.variable_min(0, 5)                         # optional - CPLEX
            sage: p.col_bounds(0)                              # optional - CPLEX
            (5.0, None)
        """
        cdef int status
        cdef double lb
        cdef char x
        cdef double c_value

        if value == False:
            status = CPXgetlb(self.env, self.lp, &lb, index, index)
            check(status)
            return lb if lb != -CPX_INFBOUND else None

        else:
            x = 'L'
            c_value = value if value is not None else -CPX_INFBOUND
            status = CPXchgbds(self.env, self.lp, 1, &index, &x, &c_value)
            check(status)

    cpdef write_lp(self, char * filename):
        r"""
        Writes the problem to a .lp file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(2)                               # optional - CPLEX
            2
            sage: p.add_linear_constraint([0, 1], [1, 2], +1, 3)          # optional - CPLEX
            sage: p.set_objective([2, 5])                          # optional - CPLEX
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp")            # optional - CPLEX
        """

        cdef int status
        cdef char * ext = "LP"
        status = CPXwriteprob(self.env, self.lp, filename, ext)
        check(status)

    cpdef write_mps(self, char * filename, int modern):
        r"""
        Writes the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import get_solver
            sage: p = get_solver(solver = "CPLEX")  # optional - CPLEX
            sage: p.add_variables(2)                               # optional - CPLEX
            2
            sage: p.add_linear_constraint([0, 1], [1, 2], +1, 3)          # optional - CPLEX
            sage: p.set_objective([2, 5])                          # optional - CPLEX
            sage: p.write_lp(SAGE_TMP+"/lp_problem.lp")            # optional - CPLEX
        """

        cdef int status
        cdef char * ext = "MPS"
        status = CPXwriteprob(self.env, self.lp, filename, ext)
        check(status)

    def __dealloc__(self):
        r"""
        Destructor for the class
        """
        CPXcloseCPLEX(self.env)

cdef check(number):
    r"""
    Given a number, raises the corresponding exception or does nothing
    if ``number == 0``.

    - ``number`` (integer) -- number corresponding to the error. If
      this number is unknown, the message contained in the raised
      exception will mention it.
    """

    # Below 1000 are 0 (no error), and some quality reports (but the
    # ERR* codes are above 1000)
    if number > 1000:
        from sage.numerical.mip import MIPSolverException
        default = "Error reported by the solver (unknown error number : "+str(number)+")"
        raise MIPSolverException("CPLEX: "+errors.get(number,default))

# Error codes
#
# Todo : when common error codes are returned, rewrite the message to
# be more meaningful

errors = {
    1001 : "CPXERR_NO_MEMORY",
    1002 : "CPXERR_NO_ENVIRONMENT",
    1003 : "CPXERR_BAD_ARGUMENT",
    1004 : "CPXERR_NULL_POINTER",
    1006 : "CPXERR_CALLBACK",
    1009 : "CPXERR_NO_PROBLEM",
    1012 : "CPXERR_LIMITS_TOO_BIG",
    1013 : "CPXERR_BAD_PARAM_NUM",
    1014 : "CPXERR_PARAM_TOO_SMALL",
    1015 : "CPXERR_PARAM_TOO_BIG",
    1016 : "CPXERR_RESTRICTED_VERSION",
    1017 : "CPXERR_NOT_FOR_MIP",
    1018 : "CPXERR_NOT_FOR_QP",
    1019 : "CPXERR_CHILD_OF_CHILD",
    1020 : "CPXERR_TOO_MANY_THREADS",
    1021 : "CPXERR_CANT_CLOSE_CHILD",
    1022 : "CPXERR_BAD_PROB_TYPE",
    1023 : "CPXERR_NOT_ONE_PROBLEM",
    1024 : "CPXERR_NOT_MILPCLASS",
    1026 : "CPXERR_STR_PARAM_TOO_LONG",
    1027 : "CPXERR_DECOMPRESSION",
    1028 : "CPXERR_BAD_PARAM_NAME",
    1029 : "CPXERR_NOT_MIQPCLASS",
    1031 : "CPXERR_NOT_FOR_QCP",
    1051 : "CPXERR_MSG_NO_CHANNEL",
    1052 : "CPXERR_MSG_NO_FILEPTR",
    1053 : "CPXERR_MSG_NO_FUNCTION",
    1101 : "CPXERR_PRESLV_INForUNBD",
    1103 : "CPXERR_PRESLV_NO_PROB",
    1106 : "CPXERR_PRESLV_ABORT",
    1107 : "CPXERR_PRESLV_BASIS_MEM",
    1108 : "CPXERR_PRESLV_COPYSOS",
    1109 : "CPXERR_PRESLV_COPYORDER",
    1110 : "CPXERR_PRESLV_SOLN_MIP",
    1111 : "CPXERR_PRESLV_SOLN_QP",
    1112 : "CPXERR_PRESLV_START_LP",
    1114 : "CPXERR_PRESLV_FAIL_BASIS",
    1115 : "CPXERR_PRESLV_NO_BASIS",
    1117 : "CPXERR_PRESLV_INF",
    1118 : "CPXERR_PRESLV_UNBD",
    1119 : "CPXERR_PRESLV_DUAL",
    1120 : "CPXERR_PRESLV_UNCRUSHFORM",
    1121 : "CPXERR_PRESLV_CRUSHFORM",
    1122 : "CPXERR_PRESLV_BAD_PARAM",
    1123 : "CPXERR_PRESLV_TIME_LIM",
    1200 : "CPXERR_INDEX_RANGE",
    1201 : "CPXERR_COL_INDEX_RANGE",
    1203 : "CPXERR_ROW_INDEX_RANGE",
    1205 : "CPXERR_INDEX_RANGE_LOW",
    1206 : "CPXERR_INDEX_RANGE_HIGH",
    1207 : "CPXERR_NEGATIVE_SURPLUS",
    1208 : "CPXERR_ARRAY_TOO_LONG",
    1209 : "CPXERR_NAME_CREATION",
    1210 : "CPXERR_NAME_NOT_FOUND",
    1211 : "CPXERR_NO_RHS_IN_OBJ",
    1215 : "CPXERR_BAD_SENSE",
    1216 : "CPXERR_NO_RNGVAL",
    1217 : "CPXERR_NO_SOLN",
    1219 : "CPXERR_NO_NAMES",
    1221 : "CPXERR_NOT_FIXED",
    1222 : "CPXERR_DUP_ENTRY",
    1223 : "CPXERR_NO_BARRIER_SOLN",
    1224 : "CPXERR_NULL_NAME",
    1225 : "CPXERR_NAN",
    1226 : "CPXERR_ARRAY_NOT_ASCENDING",
    1227 : "CPXERR_COUNT_RANGE",
    1228 : "CPXERR_COUNT_OVERLAP",
    1229 : "CPXERR_BAD_LUB",
    1230 : "CPXERR_NODE_INDEX_RANGE",
    1231 : "CPXERR_ARC_INDEX_RANGE",
    1232 : "CPXERR_NO_DUAL_SOLN",
    1233 : "CPXERR_DBL_MAX",
    1234 : "CPXERR_THREAD_FAILED",
    1251 : "CPXERR_INDEX_NOT_BASIC",
    1252 : "CPXERR_NEED_OPT_SOLN",
    1253 : "CPXERR_BAD_STATUS",
    1254 : "CPXERR_NOT_UNBOUNDED",
    1255 : "CPXERR_SBASE_INCOMPAT",
    1256 : "CPXERR_SINGULAR",
    1257 : "CPXERR_PRIIND",
    1258 : "CPXERR_NO_LU_FACTOR",
    1260 : "CPXERR_NO_SENSIT",
    1261 : "CPXERR_NO_BASIC_SOLN",
    1262 : "CPXERR_NO_BASIS",
    1263 : "CPXERR_ABORT_STRONGBRANCH",
    1264 : "CPXERR_NO_NORMS",
    1265 : "CPXERR_NOT_DUAL_UNBOUNDED",
    1266 : "CPXERR_TILIM_STRONGBRANCH",
    1267 : "CPXERR_BAD_PIVOT",
    1268 : "CPXERR_TILIM_CONDITION_NO",
    1292 : "CPXERR_BAD_METHOD",
    1421 : "CPXERR_NO_FILENAME",
    1422 : "CPXERR_FAIL_OPEN_WRITE",
    1423 : "CPXERR_FAIL_OPEN_READ",
    1424 : "CPXERR_BAD_FILETYPE",
    1425 : "CPXERR_XMLPARSE",
    1431 : "CPXERR_TOO_MANY_ROWS",
    1432 : "CPXERR_TOO_MANY_COLS",
    1433 : "CPXERR_TOO_MANY_COEFFS",
    1434 : "CPXERR_BAD_NUMBER",
    1435 : "CPXERR_BAD_EXPO_RANGE",
    1436 : "CPXERR_NO_OBJ_SENSE",
    1437 : "CPXERR_QCP_SENSE_FILE",
    1438 : "CPXERR_BAD_LAZY_UCUT",
    1439 : "CPXERR_BAD_INDCONSTR",
    1441 : "CPXERR_NO_NAME_SECTION",
    1442 : "CPXERR_BAD_SOS_TYPE",
    1443 : "CPXERR_COL_ROW_REPEATS",
    1444 : "CPXERR_RIM_ROW_REPEATS",
    1445 : "CPXERR_ROW_REPEATS",
    1446 : "CPXERR_COL_REPEATS",
    1447 : "CPXERR_RIM_REPEATS",
    1448 : "CPXERR_ROW_UNKNOWN",
    1449 : "CPXERR_COL_UNKNOWN",
    1453 : "CPXERR_NO_ROW_SENSE",
    1454 : "CPXERR_EXTRA_FX_BOUND",
    1455 : "CPXERR_EXTRA_FR_BOUND",
    1456 : "CPXERR_EXTRA_BV_BOUND",
    1457 : "CPXERR_BAD_BOUND_TYPE",
    1458 : "CPXERR_UP_BOUND_REPEATS",
    1459 : "CPXERR_LO_BOUND_REPEATS",
    1460 : "CPXERR_NO_BOUND_TYPE",
    1461 : "CPXERR_NO_QMATRIX_SECTION",
    1462 : "CPXERR_BAD_SECTION_ENDATA",
    1463 : "CPXERR_INT_TOO_BIG_INPUT",
    1464 : "CPXERR_NAME_TOO_LONG",
    1465 : "CPXERR_LINE_TOO_LONG",
    1471 : "CPXERR_NO_ROWS_SECTION",
    1472 : "CPXERR_NO_COLUMNS_SECTION",
    1473 : "CPXERR_BAD_SECTION_BOUNDS",
    1474 : "CPXERR_RANGE_SECTION_ORDER",
    1475 : "CPXERR_BAD_SECTION_QMATRIX",
    1476 : "CPXERR_NO_OBJECTIVE",
    1477 : "CPXERR_ROW_REPEAT_PRINT",
    1478 : "CPXERR_COL_REPEAT_PRINT",
    1479 : "CPXERR_RIMNZ_REPEATS",
    1480 : "CPXERR_EXTRA_INTORG",
    1481 : "CPXERR_EXTRA_INTEND",
    1482 : "CPXERR_EXTRA_SOSORG",
    1483 : "CPXERR_EXTRA_SOSEND",
    1484 : "CPXERR_TOO_MANY_RIMS",
    1485 : "CPXERR_TOO_MANY_RIMNZ",
    1486 : "CPXERR_NO_ROW_NAME",
    1487 : "CPXERR_BAD_OBJ_SENSE",
    1550 : "CPXERR_BAS_FILE_SHORT",
    1551 : "CPXERR_BAD_INDICATOR",
    1552 : "CPXERR_NO_ENDATA",
    1553 : "CPXERR_FILE_ENTRIES",
    1554 : "CPXERR_SBASE_ILLEGAL",
    1555 : "CPXERR_BAS_FILE_SIZE",
    1556 : "CPXERR_NO_VECTOR_SOLN",
    1560 : "CPXERR_NOT_SAV_FILE",
    1561 : "CPXERR_SAV_FILE_DATA",
    1562 : "CPXERR_SAV_FILE_WRITE",
    1563 : "CPXERR_FILE_FORMAT",
    1602 : "CPXERR_ADJ_SIGNS",
    1603 : "CPXERR_RHS_IN_OBJ",
    1604 : "CPXERR_ADJ_SIGN_SENSE",
    1605 : "CPXERR_QUAD_IN_ROW",
    1606 : "CPXERR_ADJ_SIGN_QUAD",
    1607 : "CPXERR_NO_OPERATOR",
    1608 : "CPXERR_NO_OP_OR_SENSE",
    1609 : "CPXERR_NO_ID_FIRST",
    1610 : "CPXERR_NO_RHS_COEFF",
    1611 : "CPXERR_NO_NUMBER_FIRST",
    1612 : "CPXERR_NO_QUAD_EXP",
    1613 : "CPXERR_QUAD_EXP_NOT_2",
    1614 : "CPXERR_NO_QP_OPERATOR",
    1615 : "CPXERR_NO_NUMBER",
    1616 : "CPXERR_NO_ID",
    1617 : "CPXERR_BAD_ID",
    1618 : "CPXERR_BAD_EXPONENT",
    1619 : "CPXERR_Q_DIVISOR",
    1621 : "CPXERR_NO_BOUND_SENSE",
    1622 : "CPXERR_BAD_BOUND_SENSE",
    1623 : "CPXERR_NO_NUMBER_BOUND",
    1627 : "CPXERR_NO_SOS_SEPARATOR",
    1650 : "CPXERR_INVALID_NUMBER",
    1660 : "CPXERR_PRM_DATA",
    1661 : "CPXERR_PRM_HEADER",
    1719 : "CPXERR_NO_CONFLICT",
    1720 : "CPXERR_CONFLICT_UNSTABLE",
    1801 : "CPXERR_WORK_FILE_OPEN",
    1802 : "CPXERR_WORK_FILE_READ",
    1803 : "CPXERR_WORK_FILE_WRITE",
    1804 : "CPXERR_IN_INFOCALLBACK",
    1805 : "CPXERR_MIPSEARCH_WITH_CALLBACKS",
    1806 : "CPXERR_LP_NOT_IN_ENVIRONMENT",
    1807 : "CPXERR_PARAM_INCOMPATIBLE",
    32000 : "CPXERR_LICENSE_MIN",
    32201 : "CPXERR_ILOG_LICENSE",
    32301 : "CPXERR_NO_MIP_LIC",
    32302 : "CPXERR_NO_BARRIER_LIC",
    32305 : "CPXERR_NO_MIQP_LIC",
    32018 : "CPXERR_BADLDWID",
    32023 : "CPXERR_BADPRODUCT",
    32024 : "CPXERR_ALGNOTLICENSED",
    32999 : "CPXERR_LICENSE_MAX",
    1701 : "CPXERR_IIS_NO_INFO",
    1702 : "CPXERR_IIS_NO_SOLN",
    1703 : "CPXERR_IIS_FEAS",
    1704 : "CPXERR_IIS_NOT_INFEAS",
    1705 : "CPXERR_IIS_OPT_INFEAS",
    1706 : "CPXERR_IIS_DEFAULT",
    1707 : "CPXERR_IIS_NO_BASIC",
    1709 : "CPXERR_IIS_NO_LOAD",
    1710 : "CPXERR_IIS_SUB_OBJ_LIM",
    1711 : "CPXERR_IIS_SUB_IT_LIM",
    1712 : "CPXERR_IIS_SUB_TIME_LIM",
    1713 : "CPXERR_IIS_NUM_BEST",
    1714 : "CPXERR_IIS_SUB_ABORT",
    3003 : "CPXERR_NOT_MIP",
    3006 : "CPXERR_BAD_PRIORITY",
    3007 : "CPXERR_ORDER_BAD_DIRECTION",
    3009 : "CPXERR_ARRAY_BAD_SOS_TYPE",
    3010 : "CPXERR_UNIQUE_WEIGHTS",
    3012 : "CPXERR_BAD_DIRECTION",
    3015 : "CPXERR_NO_SOS",
    3016 : "CPXERR_NO_ORDER",
    3018 : "CPXERR_INT_TOO_BIG",
    3019 : "CPXERR_SUBPROB_SOLVE",
    3020 : "CPXERR_NO_MIPSTART",
    3021 : "CPXERR_BAD_CTYPE",
    3023 : "CPXERR_NO_INT_X",
    3024 : "CPXERR_NO_SOLNPOOL",
    3301 : "CPXERR_MISS_SOS_TYPE",
    3412 : "CPXERR_NO_TREE",
    3413 : "CPXERR_TREE_MEMORY_LIMIT",
    3414 : "CPXERR_FILTER_VARIABLE_TYPE",
    3504 : "CPXERR_NODE_ON_DISK",
    3601 : "CPXERR_PTHREAD_MUTEX_INIT",
    3603 : "CPXERR_PTHREAD_CREATE",
    1212 : "CPXERR_UNSUPPORTED_CONSTRAINT_TYPE",
    1213 : "CPXERR_ILL_DEFINED_PWL",
    1530 : "CPXERR_NET_DATA",
    1531 : "CPXERR_NOT_MIN_COST_FLOW",
    1532 : "CPXERR_BAD_ROW_ID",
    1537 : "CPXERR_BAD_CHAR",
    1538 : "CPXERR_NET_FILE_SHORT",
    5002 : "CPXERR_Q_NOT_POS_DEF",
    5004 : "CPXERR_NOT_QP",
    5011 : "CPXERR_Q_DUP_ENTRY",
    5012 : "CPXERR_Q_NOT_SYMMETRIC",
    5014 : "CPXERR_Q_NOT_INDEF",
    6002 : "CPXERR_QCP_SENSE"
    }

