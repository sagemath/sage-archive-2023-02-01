include '../../../../devel/sage/sage/ext/stdsage.pxi'
include '../ext/cdefs.pxi'

cdef int BINARY = 1
cdef int REAL = -1
cdef int INTEGER = 0

cdef class Osi_interface:
    cdef float osi_solve(self, LP, c_OsiSolverInterface * si,bool objective_only, bool is_cplex):
        from itertools import izip
        n_cols = len(LP._variables_type);

        # Definition of the objective function
        cdef double * objective = <double*> sage_malloc(n_cols*sizeof(double))

        # a hand-made memset
        memset(objective, 0, n_cols*sizeof(double))

        for (i,v) in izip(LP._objective_i,LP._objective_values):
            objective[i]=float(v)


        # Upper and lower bounds on the variables
        cdef double * col_lb = <double*> sage_malloc(n_cols*sizeof(double))
        cdef double * col_ub = <double*> sage_malloc(n_cols*sizeof(double))


        cdef int id
        id=0

        for (type,min,max) in izip(LP._variables_type,LP._variables_bounds_min,LP._variables_bounds_max):
            if type == BINARY:
                col_lb[id] = 0.0
                col_ub[id] = 1.0
            else:
                if min == None:
                    col_lb[id] = -1.0*si.getInfinity();
                else:
                    col_lb[id] = float(min)

                if max == None:
                    col_ub[id] = 1.0*si.getInfinity();
                else:
                    col_ub[id] = float(max)
            id+=1

        # Definition of the problem matrix
        n_rows = len(LP._constraints_bounds_min);
        cdef c_CoinPackedMatrix * matrix
        matrix = new_c_CoinPackedMatrix(False, 0, 0);
        matrix.setDimensions(0, n_cols);


        # Upper and lower bounds on the constraints
        cdef double * row_lb = <double*> sage_malloc(n_rows*sizeof(double))
        cdef double * row_ub = <double*> sage_malloc(n_rows*sizeof(double))

        # The constraint to be added
        cdef c_CoinPackedVector* row

        cdef int a
        a = 0

        for (min,max) in izip(LP._constraints_bounds_min,LP._constraints_bounds_max):
            if min == None:
                row_lb[a] = -1.0*si.getInfinity();
            else:
                row_lb[a] = float(min)

            if max == None:
                row_ub[a] = 1.0*si.getInfinity();
            else:
                row_ub[a] = float(max)
            a += 1

        cdef int last
        last = 0;
        row = new_c_CoinPackedVector()

        for (i,j,v) in izip (LP._constraints_matrix_i,LP._constraints_matrix_j,LP._constraints_matrix_values):
            if i != last:
                matrix.appendRow(row[0]);
                row = new_c_CoinPackedVector();
                last = i
            row.insert(j, float(v))

        matrix.appendRow(row[0]);

        # Direction of the problem : maximization, minimization...
        if LP._maximization == False:
            si.setObjSense(1.0)
        else:
            si.setObjSense(-1.0)

        si.assignProblem(matrix, col_lb, col_ub, objective, row_lb, row_ub);

        # Sets variables as integers when needed
        id = 0;
        for type in LP._variables_type:
            if type == REAL:
                si.setContinuous(id)
            else:
                si.setInteger(id)
            id += 1

        # Has there been any problem during the computations ?
        flag = True


        # Solves the MIP
        si.branchAndBound();

        cdef const_double_ptr solution

        if si.isAbandoned():
            exception_txt = "The solver has abandoned!"
        elif si.isProvenPrimalInfeasible() or si.isProvenDualInfeasible():
            exception_txt="The problem or its dual has been proven infeasible!"
        elif (not is_cplex) and (si.isPrimalObjectiveLimitReached() or si.isDualObjectiveLimitReached()):
            exception_txt = "The objective limit has been reached for the problem or its dual!"
        elif si.isIterationLimitReached():
            exception_txt = "The iteration limit has been reached!"
        elif si.isProvenOptimal() == True:
            flag = False
            return_value = si.getObjValue();
            if objective_only == False:
                solution = si.getColSolution();

                # Builds the returned dictionary of values
                for ((v,id),type) in izip(LP._variables.iteritems(),LP._variables_type):
                    LP._values[v] = solution[id] if type == REAL else round(solution[id])


        sage_free(objective)
        sage_free(col_lb)
        sage_free(col_ub)
        sage_free(row_lb)
        sage_free(row_ub)
        del_CoinPackedMatrix(matrix)

        if flag:
            from sage.numerical.mip import MIPSolverException
            raise MIPSolverException(exception_txt)
        return return_value
