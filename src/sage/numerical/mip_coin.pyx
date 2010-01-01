include 'mip_coin.pxd'
include '../../../../devel/sage/sage/ext/stdsage.pxi'

cdef int BINARY = 1
cdef int REAL = -1
cdef int INTEGER = 0


def solve_coin(self,log=False,objective_only=False):
    r"""
    Solves the MIP using COIN-OR CBC/CLP
    """

    from itertools import izip

    # Creates the solver interfaces
    cdef c_OsiCbcSolverInterface* si
    si = new_c_OsiCbcSolverInterface();

    n_cols = len(self._variables_type);

    # Definition of the objective function
    cdef double * objective = <double*> sage_malloc(n_cols*sizeof(double))

    # a hand-made memset
    for 0<= i < n_cols:
        objective[i] = 0

    for (i,v) in izip(self._objective_i,self._objective_values):
        objective[i]=float(v)


    # Upper and lower bounds on the variables
    cdef double * col_lb = <double*> sage_malloc(n_cols*sizeof(double))
    cdef double * col_ub = <double*> sage_malloc(n_cols*sizeof(double))


    cdef int id
    id=0

    for (type,min,max) in izip(self._variables_type,self._variables_bounds_min,self._variables_bounds_max):
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
    n_rows = len(self._constraints_bounds_min);
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

    for (min,max) in izip(self._constraints_bounds_min,self._constraints_bounds_max):
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

    for (i,j,v) in izip (self._constraints_matrix_i,self._constraints_matrix_j,self._constraints_matrix_values):
        if i != last:
            matrix.appendRow(row[0]);
            row = new_c_CoinPackedVector();
            last = i
        row.insert(j, float(v))

    matrix.appendRow(row[0]);

    # Direction of the problem : maximization, minimization...
    if self._maximization == False:
        si.setObjSense(1.0)
    else:
        si.setObjSense(-1.0)

    # stuff related to the loglevel
    cdef c_CoinMessageHandler * msg
    cdef c_CbcModel * model
    model = si.getModelPtr()
    msg = model.messageHandler()
    msg.setLogLevel(log if False != log else log)

    si.assignProblem(matrix, col_lb, col_ub, objective, row_lb, row_ub);

    # Sets variables as integers when needed
    id = 0;
    for type in self._variables_type:
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
        flag = True
        exception_txt = "Coin Branch and Cut: The solver has abandoned!"
    elif si.isProvenPrimalInfeasible() or si.isProvenDualInfeasible():
        flag = True
        exception_txt="Coin Branch and Cut: The problem or its dual has been proven infeasible!"
    elif si.isPrimalObjectiveLimitReached() or si.isDualObjectiveLimitReached():
        flag = True
        exception_txt = "Coin Branch and Cut: The objective limit has been reached for the problem or its dual!"
    elif si.isIterationLimitReached():
        flag = True
        exception_txt = "Coin Branch and Cut: The iteration limit has been reached!"
    elif si.isProvenOptimal() == True:
        flag = False
        return_value = si.getObjValue();
        if objective_only == False:
            solution = si.getColSolution();

            # Builds the returned dictionary of values
            for ((v,id),type) in izip(self._variables.iteritems(),self._variables_type):
                self._values[v] = solution[id] if type == REAL else round(solution[id])


    del_OsiCbcSolverInterface(si)
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
