r"""
COIN Backend
"""

from sage.numerical.mip import MIPSolverException

cdef class CoinBackend(GenericBackend):

    def __cinit__(self, maximization = True):

        self.si = new_c_OsiCbcSolverInterface();
        # stuff related to the loglevel
        cdef c_CbcModel * model
        model = self.si.getModelPtr()

        import multiprocessing
        model.setNumberThreads(multiprocessing.cpu_count())

        self.set_log_level(0)

        if maximization:
            self.set_direction(+1)
        else:
            self.set_direction(-1)

    cpdef int add_variable(self):
        r"""
        Adds a variable.

        This amounts to adding a new column to the matrix. By default,
        the variable is both positive and real.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")    # optional - Coin
            sage: p.n_cols()                                         # optional - Coin
            0
            sage: p.add_variable()                                   # optional - Coin
            1
            sage: p.n_cols()                                         # optional - Coin
            1
        """

        self.si.addCol(0, NULL, NULL, 0, self.si.getInfinity(), 0)
        return self.si.getNumCols()

    cpdef int add_variables(self, int number):
        r"""
        Adds ``number`` variables.

        This amounts to adding new columns to the matrix. By default,
        the variables are both positive and real.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")    # optional - Coin
            sage: p.n_cols()                                         # optional - Coin
            0
            sage: p.add_variables(5)                                 # optional - Coin
            5
            sage: p.n_cols()                                         # optional - Coin
            5
        """

        cdef int i
        for 0<= i< number:
            self.si.addCol(0, NULL, NULL, 0, self.si.getInfinity(), 0)
        return self.si.getNumCols()

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
            sage: p = getSolver(solver = "Coin")   # optional - Coin
            sage: p.n_cols()                                        # optional - Coin
            0
            sage: p.add_variable()                                  # optional - Coin
            1
            sage: p.set_variable_type(0,1)                          # optional - Coin
            sage: p.is_variable_integer(0)                          # optional - Coin
            True
        """

        if vtype == 1:
            self.si.setInteger(variable)
        elif vtype == 0:
            self.si.setColLower(variable, 0)
            self.si.setInteger(variable)
            self.si.setColUpper(variable, 1)
        else:
            self.si.setContinuous(variable)

    cpdef set_direction(self, int sense):
        r"""
        Sets the direction (maximization/minimization).

        INPUT:

        - ``sense`` (integer) :

            * +1 => Maximization
            * -1 => Minimization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.is_maximization()                              # optional - Coin
            True
            sage: p.set_direction(-1)                              # optional - Coin
            sage: p.is_maximization()                              # optional - Coin
            False
        """
        self.si.setObjSense(-sense)

    cpdef set_objective_coeff(self, int variable, double coeff):
        r"""
        Sets the coefficient of a variable in the objective function

        INPUT:

        - ``variable`` (integer) -- the variable's id

        - ``coeff`` (double) -- its coefficient

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.get_objective_coeff(0)                         # optional - Coin
            0.0
            sage: p.set_objective_coeff(0,2)                       # optional - Coin
            sage: p.get_objective_coeff(0)                         # optional - Coin
            2.0
        """

        self.si.setObjCoeff(variable, coeff)

    cpdef set_objective(self, list coeff):
        r"""
        Sets the objective function.

        INPUT:

        - ``coeff`` -- a list of real values, whose ith element is the
          coefficient of the ith variable in the objective function.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")    # optional - Coin
            sage: p.add_variables(5)                                 # optional - Coin
            5
            sage: p.set_objective([1, 1, 2, 1, 3])                   # optional - Coin
            sage: map(lambda x :p.get_objective_coeff(x), range(5))  # optional - Coin
            [1.0, 1.0, 2.0, 1.0, 3.0]
        """

        cdef int i
        for i,v in enumerate(coeff):
            self.si.setObjCoeff(i, v)

    cpdef set_log_level(self, int level):
        r"""
        Sets the log (verbosity) level

        INPUT:

        - ``level`` (integer) -- From 0 (no verbosity) to 3.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")   # optional - Coin
            sage: p.set_log_level(2)                                # optional - Coin

        """

        cdef c_CbcModel * model
        cdef c_CoinMessageHandler * msg
        model = self.si.getModelPtr()
        msg = model.messageHandler()
        msg.setLogLevel(level)

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
            sage: p = getSolver(solver = "Coin")   # optional - Coin
            sage: p.add_variables(5)                                # optional - Coin
            5
            sage: p.add_constraints(5, +1, 2)                       # optional - Coin
            sage: p.get_row(4)                                      # optional - Coin
            ([], [])
            sage: p.get_row_bounds(4)                               # optional - Coin
            (None, 2.0)
        """

        cdef int i
        for 0<= i<number:
            self.add_constraint([],[],direction, bound)

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
            sage: p = getSolver(solver = "Coin") # optional - Coin
            sage: p.add_variables(5)                              # optional - Coin
            5
            sage: p.add_constraint(range(5), range(5), 0, 2)      # optional - Coin
            sage: p.get_row(0)                                    # optional - Coin
            ([0, 1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0, 4.0])
            sage: p.get_row_bounds(0)                             # optional - Coin
            (2.0, 2.0)
        """

        cdef int i
        cdef int n = len(indices)
        cdef c_CoinPackedVector* row
        row = new_c_CoinPackedVector();

        for 0<= i<n:
            row.insert(indices[i], coeffs[i])

        self.si.addRow (row[0],
                        bound if direction != +1 else -self.si.getInfinity(),
                        bound if direction != -1 else +self.si.getInfinity())


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
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variables(5)                               # optional - Coin
            5
            sage: p.add_constraint(range(5), range(5), 0, 2)       # optional - Coin
            sage: p.get_row(0)                                     # optional - Coin
            ([0, 1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0, 4.0])
            sage: p.get_row_bounds(0)                              # optional - Coin
            (2.0, 2.0)
        """

        cdef list indices = []
        cdef list values = []
        cdef int * c_indices
        cdef int i
        cdef double * c_values
        cdef c_CoinPackedMatrix * M = <c_CoinPackedMatrix *> self.si.getMatrixByRow()
        cdef c_CoinShallowPackedVector V = <c_CoinShallowPackedVector> M.getVector(index)
        cdef int n = V.getNumElements()

        c_indices = <int*> V.getIndices()
        c_values = <double*> V.getElements()

        for 0<= i < n:
            indices.append(c_indices[i])
            values.append(c_values[i])

        return (indices, values)

    cpdef get_row_bounds(self, int i):
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
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variables(5)                               # optional - Coin
            5
            sage: p.add_constraint(range(5), range(5), 0, 2)       # optional - Coin
            sage: p.get_row(0)                                     # optional - Coin
            ([0, 1, 2, 3, 4], [0.0, 1.0, 2.0, 3.0, 4.0])
            sage: p.get_row_bounds(0)                              # optional - Coin
            (2.0, 2.0)
        """

        cdef double * ub
        cdef double * lb

        ub = <double*> self.si.getRowUpper()
        lb = <double*> self.si.getRowLower()

        return (lb[i] if lb[i] != - self.si.getInfinity() else None,
                ub[i] if ub[i] != + self.si.getInfinity() else None)

    cpdef get_col_bounds(self, int i):
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
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.get_col_bounds(0)                              # optional - Coin
            (0.0, None)
            sage: p.set_variable_max(0, 5)                         # optional - Coin
            sage: p.get_col_bounds(0)                              # optional - Coin
            (0.0, 5.0)
        """

        cdef double * ub
        cdef double * lb

        ub = <double*> self.si.getColUpper()
        lb = <double*> self.si.getColLower()

        return (lb[i] if lb[i] != - self.si.getInfinity() else None,
                ub[i] if ub[i] != + self.si.getInfinity() else None)

    cpdef double get_objective_coeff(self, int index):
        return self.si.getObjCoefficients()[index]

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
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.n_cols()                                       # optional - Coin
            0
            sage: p.n_rows()                                       # optional - Coin
            0
            sage: p.add_constraints(5, -1, 0)                      # optional - Coin
            sage: p.add_col(range(5), range(5))                    # optional - Coin
            sage: p.n_rows()                                       # optional - Coin
            5
        """

        cdef int n = len(indices)
        cdef int * c_indices = <int*>sage_malloc(n*sizeof(int))
        cdef double * c_values  = <double*>sage_malloc(n*sizeof(double))
        cdef int i

        for 0<= i< n:
            c_indices[i] = indices[i]
            c_values[i] = coeffs[i]


        self.si.addCol (1, c_indices, c_values, 0, self.si.getInfinity(), 0)

    cpdef int solve(self) except -1:
        r"""
        Solves the problem.

        .. NOTE::

            This method raises ``MIPSolverException`` exceptions when
            the solution can not be computed for any reason (none
            exists, or the LP solver was not able to find it, etc...)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")    # optional - Coin
            sage: p.add_constraints(5, -1, 0)       # optional - Coin
            sage: p.add_col(range(5), [1,2,3,4,5])  # optional - Coin
            sage: p.solve()                         # optional - Coin
            0

        TESTS::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")    # optional - Coin
            sage: p.add_variable()                  # optional - Coin
            1
            sage: p.add_constraint([0], [1], +1, 4) # optional - Coin
            sage: p.add_constraint([0], [1], -1, 6) # optional - Coin
            sage: p.set_objective_coeff(0,1)        # optional - Coin
            sage: p.solve()                         # optional - Coin
            Traceback (most recent call last):
            ...
            MIPSolverException: ...
        """

        self.si.branchAndBound()

        cdef const_double_ptr solution

        if self.si.isAbandoned():
            raise MIPSolverException("CBC : The solver has abandoned!")

        elif self.si.isProvenPrimalInfeasible() or self.si.isProvenDualInfeasible():
            raise MIPSolverException("CBC : The problem or its dual has been proven infeasible!")

        elif (self.si.isPrimalObjectiveLimitReached() or self.si.isDualObjectiveLimitReached()):
            raise MIPSolverException("CBC : The objective limit has been reached for the problem or its dual!")

        elif self.si.isIterationLimitReached():
            raise MIPSolverException("CBC : The iteration limit has been reached!")

        elif not self.si.isProvenOptimal():
            raise MIPSolverException("CBC : Unknown error")

    cpdef double get_objective_value(self):
        r"""
        Returns the value of the objective function.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variables(2)                               # optional - Coin
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)          # optional - Coin
            sage: p.set_objective([2, 5])                          # optional - Coin
            sage: p.solve()                                        # optional - Coin
            0
            sage: p.get_objective_value()                          # optional - Coin
            7.5
            sage: p.get_variable_value(0)                          # optional - Coin
            0.0
            sage: p.get_variable_value(1)                          # optional - Coin
            1.5
        """

        return self.si.getObjValue()

    cpdef double get_variable_value(self, int variable):
        r"""
        Returns the value of a variable given by the solver.

        .. NOTE::

           Has no meaning unless ``solve`` has been called before.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin") # optional - Coin
            sage: p.add_variables(2)                              # optional - Coin
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)         # optional - Coin
            sage: p.set_objective([2, 5])                         # optional - Coin
            sage: p.solve()                                       # optional - Coin
            0
            sage: p.get_objective_value()                         # optional - Coin
            7.5
            sage: p.get_variable_value(0)                         # optional - Coin
            0.0
            sage: p.get_variable_value(1)                         # optional - Coin
            1.5
        """

        cdef double * solution
        solution = <double*> self.si.getColSolution()
        return solution[variable]

    cpdef int n_cols(self):
        r"""
        Returns the number of columns/variables.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.n_cols()                                       # optional - Coin
            0
            sage: p.add_variables(2)                               # optional - Coin
            2
            sage: p.n_cols()                                       # optional - Coin
            2
        """

        return self.si.getNumCols()

    cpdef int n_rows(self):
        r"""
        Returns the number of rows/constraints.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin") # optional - Coin
            sage: p.n_rows()                                      # optional - Coin
            0
            sage: p.add_constraints(2, -1, 2)                     # optional - Coin
            sage: p.n_rows()                                      # optional - Coin
            2
        """
        return self.si.getNumRows()


    cpdef bint is_variable_binary(self, int index):
        r"""
        Tests whether the given variable is of binary type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.n_cols()                                       # optional - Coin
            0
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.set_variable_type(0,0)                         # optional - Coin
            sage: p.is_variable_binary(0)                          # optional - Coin
            True

        """

        return (0 == self.si.isContinuous(index) and
                self.get_variable_min(index) == 0 and
                self.get_variable_max(index) == 1)

    cpdef bint is_variable_integer(self, int index):
        r"""
        Tests whether the given variable is of integer type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.n_cols()                                       # optional - Coin
            0
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.set_variable_type(0,1)                         # optional - Coin
            sage: p.is_variable_integer(0)                         # optional - Coin
            True
        """
        return (0 == self.si.isContinuous(index) and
                (self.get_variable_min(index) != 0 or
                self.get_variable_max(index) != 1))

    cpdef bint is_variable_continuous(self, int index):
        r"""
        Tests whether the given variable is of continuous/real type.

        INPUT:

        - ``index`` (integer) -- the variable's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.n_cols()                                       # optional - Coin
            0
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.is_variable_continuous(0)                      # optional - Coin
            True
            sage: p.set_variable_type(0,1)                         # optional - Coin
            sage: p.is_variable_continuous(0)                      # optional - Coin
            False

        """
        return 1 == self.si.isContinuous(index)


    cpdef bint is_maximization(self):
        r"""
        Tests whether the problem is a maximization

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin") # optional - Coin
            sage: p.is_maximization()                             # optional - Coin
            True
            sage: p.set_direction(-1)                             # optional - Coin
            sage: p.is_maximization()                             # optional - Coin
            False
        """

        return self.si.getObjSense() == -1


    cpdef get_variable_max(self, int i):
        r"""
        Returns the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        OUTPUT:

        A real value if the variable has an upper bound, ``None``
        otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.get_variable_max(0) is None                    # optional - Coin
            True
            sage: p.set_variable_max(0, 5)                         # optional - Coin
            sage: p.get_variable_max(0)                            # optional - Coin
            5.0

        """

        cdef double * ub
        ub = <double*> self.si.getColUpper()
        return ub[i] if ub[i] != + self.si.getInfinity() else None


    cpdef get_variable_min(self, int i):
        r"""
        Returns the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        OUTPUT:

        A real value if the variable has an lower bound, ``None``
        otherwise.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.get_variable_min(0)                            # optional - Coin
            0.0
            sage: p.set_variable_min(0, 5)                         # optional - Coin
            sage: p.get_variable_min(0)                            # optional - Coin
            5.0
        """

        cdef double * lb
        lb = <double*> self.si.getColLower()
        return lb[i] if lb[i] != - self.si.getInfinity() else None

    cpdef set_variable_max(self, int index, value):
        r"""
        Sets the upper bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not upper bound.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.get_col_bounds(0)                              # optional - Coin
            (0.0, None)
            sage: p.set_variable_max(0, 5)                         # optional - Coin
            sage: p.get_col_bounds(0)                              # optional - Coin
            (0.0, 5.0)
        """

        self.si.setColUpper(index, value if value is not None else +self.si.getInfinity())

    cpdef set_variable_min(self, int index, value):
        r"""
        Sets the lower bound on a variable

        INPUT:

        - ``index`` (integer) -- the variable's id

        - ``value`` -- real value, or ``None`` to mean that the
          variable has not lower bound.

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.get_col_bounds(0)                              # optional - Coin
            (0.0, None)
            sage: p.set_variable_min(0, 5)                         # optional - Coin
            sage: p.get_col_bounds(0)                              # optional - Coin
            (5.0, None)
        """

        self.si.setColLower(index, value if value is not None else -self.si.getInfinity())

    cpdef write_mps(self, char * filename, int modern):
        r"""
        Writes the problem to a .mps file

        INPUT:

        - ``filename`` (string)

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variables(2)                               # optional - Coin
            2
            sage: p.add_constraint([0, 1], [1, 2], +1, 3)          # optional - Coin
            sage: p.set_objective([2, 5])                          # optional - Coin
            sage: p.write_mps(SAGE_TMP+"/lp_problem.mps", 0)       # optional - Coin
        """

        cdef char * mps = "mps"
        self.si.writeMps(filename, mps, -1 if self.is_maximization() else 1)

    cpdef set_problem_name(self, char * name):
        r"""
        Sets the problem's name

        INPUT:

        - ``name`` (``char *``) -- the problem's name

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")   # optional - Coin
            sage: p.set_problem_name("There once was a french fry") # optional - Coin
            sage: print p.get_problem_name()                        # optional - Coin
            <BLANKLINE>
        """

        pass

    cpdef get_problem_name(self):
        r"""
        Returns the problem's name

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")   # optional - Coin
            sage: p.set_problem_name("There once was a french fry") # optional - Coin
            sage: print p.get_problem_name()                        # optional - Coin
            <BLANKLINE>
        """

        return ""

    cpdef get_row_name(self, int index):
        r"""
        Returns the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")           # optional - Coin
            sage: p.add_variable()                         # optional - Coin
            1
            sage: p.set_row_name(0, "I am a constraint")   # optional - Coin
            sage: print p.get_row_name(0)                  # optional - Coin
            <BLANKLINE>
        """

        return ""


    cpdef set_row_name(self, int index, char * name):
        r"""
        Sets the ``index`` th row name

        INPUT:

        - ``index`` (integer) -- the row's id

        - ``name`` (``char *``) -- its name

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_constraints(1, -1, 2)                      # optional - Coin
            sage: p.set_row_name(0, "Empty constraint 1")          # optional - Coin
            sage: print p.get_row_name(0)                          # optional - Coin
            <BLANKLINE>
        """

        pass

    cpdef set_objective_name(self, name):
        pass

    cpdef set_col_name(self, int index, char * name):
        r"""
        Sets the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        - ``name`` (``char *``) -- its name

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.set_col_name(0, "I am a variable")             # optional - Coin
            sage: print p.get_col_name(0)                          # optional - Coin
            <BLANKLINE>
        """
        pass

    cpdef get_col_name(self, int index):
        r"""
        Returns the ``index`` th col name

        INPUT:

        - ``index`` (integer) -- the col's id

        EXAMPLE::

            sage: from sage.numerical.backends.generic_backend import getSolver
            sage: p = getSolver(solver = "Coin")  # optional - Coin
            sage: p.add_variable()                                 # optional - Coin
            1
            sage: p.set_col_name(0, "I am a variable")             # optional - Coin
            sage: print p.get_col_name(0)                          # optional - Coin
            <BLANKLINE>
        """

        return ""
