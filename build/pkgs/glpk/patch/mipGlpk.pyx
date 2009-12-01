"""
AUTHOR:

- Nathan Cohen (2009): initial version
"""

include 'mipGlpk.pxd'
include '../../../../devel/sage/sage/ext/stdsage.pxi'

cdef int BINARY = 1
cdef int REAL = -1
cdef int INTEGER = 0

def solve_glpk(self, log=0, objective_only=False):
    """
    Solves the ``MixedIntegerLinearProgram`` using GLPK.
    Use ``solve()`` instead.

    INPUT:

    - ``log`` (integer) -- level of verbosity during the
      solving process. Set ``log=3`` for maximal verbosity and
      ``log=0`` for no verbosity at all.

    - ``objective_only`` (boolean) -- Indicates whether the
      values corresponding to an optimal assignment of the
      variables should be computed. Setting ``objective_only``
      to ``True`` saves some time on the computation when
      this information is not needed.

    OUTPUT:

    This function returns the optimal value of the objective
    function given by GLPK.
    """

    # Raises an exception if no objective has been defined
    if self._objective_i is None:
        raise ValueError("No objective function has been defined!")

    # Creates the structures
    cdef c_glp_prob * lp = glp_create_prob()
    cdef c_glp_iocp * iocp = new_c_glp_iocp()
    glp_init_iocp(iocp)

    # Builds the GLPK version of the problem
    build_glp_prob(lp, iocp, self,log, False)

    # Solves the problem
    glp_intopt(lp, iocp)

    cdef int status = glp_mip_status(lp)
    if status == GLP_UNDEF:
        glp_delete_prob(lp)
        from sage.numerical.mip import MIPSolverException
        raise MIPSolverException("GLPK : Solution is undefined")

    elif status == GLP_OPT:
        # Everything is fine !
        obj = glp_mip_obj_val(lp)

        # If the users is interested in the variables' value
        # fills self._values
        if objective_only == False:
            for (v, id) in self._variables.items():
                self._values[v] = glp_mip_col_val(lp, id+1) if self._variables_type[id]==REAL else round(glp_mip_col_val(lp, id+1))

        glp_delete_prob(lp)
        return obj

    elif status == GLP_FEAS:
        glp_delete_prob(lp)
        from sage.numerical.mip import MIPSolverException
        raise MIPSolverException("GLPK : Feasible solution found, while optimality has not been proven")

    elif status == GLP_NOFEAS:
        glp_delete_prob(lp)
        from sage.numerical.mip import MIPSolverException
        raise MIPSolverException("GLPK : There is no feasible integer solution to this Linear Program")

    else:
        glp_delete_prob(lp)
        raise ValueError("The value returned by the solver is unknown!")

def write_mps(self, filename,modern=True):
    """
    Write the linear program as a MPS file.

    INPUT:

    - ``filename`` -- The file in which you want the problem
      to be written.
    - ``modern`` -- Lets you choose between Fixed MPS and Free MPS
        - ``True`` -- Writes the problem in Free MPS
        - ``False`` -- Writes the problem in Fixed MPS

    OUTPUT:

    This function has no output.

    For information about the MPS file format :
    http://en.wikipedia.org/wiki/MPS_%28format%29
    """

    # Raises an exception if no objective has been defined
    if self._objective_i is None:
        raise ValueError("No objective function has been defined!")

    # Creates the structures
    cdef c_glp_prob * lp = glp_create_prob()
    cdef c_glp_iocp * iocp = new_c_glp_iocp()
    glp_init_iocp(iocp)

    # Builds the GLPK version of the problem
    build_glp_prob(lp, iocp, self, False, True)

    if modern:
        glp_write_mps(lp, GLP_MPS_FILE, NULL, filename)
    else:
        glp_write_mps(lp, GLP_MPS_DECK, NULL, filename)

def write_lp(self, filename):
    """
    Write the linear program as a MPS file.

    INPUT:

    - ``filename`` -- The file in which you want the problem
      to be written.

    OUTPUT:

    This function has no output.

    For more information about the LP file format :
    http://lpsolve.sourceforge.net/5.5/lp-format.htm
    """

    # Raises an exception if no objective has been defined
    if self._objective_i is None:
        raise ValueError("No objective function has been defined!")

    # Creates the structures
    cdef c_glp_prob * lp = glp_create_prob()
    cdef c_glp_iocp * iocp = new_c_glp_iocp()
    glp_init_iocp(iocp)

    # Builds the GLPK version of the problem
    build_glp_prob(lp, iocp, self,False,True)

    glp_write_lp(lp,NULL,filename)

cdef int build_glp_prob(c_glp_prob * lp, c_glp_iocp * iocp, LP, int log, bool names) except -1:
    """
    Builds the GLPK structure corresponding to the LP

    INPUT:

    - ``lp`` -- the ``c_glp_prob`` structure to be filled
    - ``iocp`` -- the ``c_glp_iocp`` structure to be filled
    - ``LP`` -- The data defining the LP
    - ``log`` -- integer indicating the level of verbosity
    - ``names`` -- should the names ( both variables and constraints ) be included
      in the problem ? ( only useful when the problem is to be written in a file )
    """
    glp_add_rows(lp, len(LP._constraints_bounds_max))
    glp_add_cols(lp, len(LP._variables_bounds_min))

    cdef int * c_matrix_i = <int*> sage_malloc((len(LP._constraints_matrix_i)+1)*sizeof(int))
    cdef int * c_matrix_j = <int*> sage_malloc((len(LP._constraints_matrix_j)+1)*sizeof(int))
    cdef double * c_matrix_values = <double*> sage_malloc((len(LP._constraints_matrix_values)+1)*sizeof(double))

    # iterator faster than "for 0<= i < len(LP.constraints.matrix.i)" ?

    cdef int c
    for c,v in enumerate(LP._constraints_matrix_i):
        c_matrix_i[c+1]=v+1

    for c,v in enumerate(LP._constraints_matrix_j):
        c_matrix_j[c+1]=v+1

    for c,v in enumerate(LP._constraints_matrix_values):
        c_matrix_values[c+1]=v

    glp_load_matrix(lp, len(LP._constraints_matrix_i),c_matrix_i,c_matrix_j,c_matrix_values)

    cdef int i
    for i in range(len(LP._objective_i)):
        glp_set_obj_coef(lp, LP._objective_i[i]+1, LP._objective_values[i])


    for i in range(len(LP._variables_bounds_min)):
        # Upper and lower bounds on the variables
        if LP._variables_bounds_min[i] is not None and LP._variables_bounds_max[i] is not None:
            glp_set_col_bnds(lp, i+1, GLP_DB, LP._variables_bounds_min[i], LP._variables_bounds_max[i])

        elif LP._variables_bounds_min[i] is None and LP._variables_bounds_max[i] is not None:
            glp_set_col_bnds(lp, i+1, GLP_UP, 0, LP._variables_bounds_max[i])

        elif LP._variables_bounds_min[i] is not None and LP._variables_bounds_max[i] is None:
            glp_set_col_bnds(lp, i+1, GLP_LO, LP._variables_bounds_min[i], 0)

        else:
            glp_set_col_bnds(lp, i+1, GLP_FR, 0, 0)

        # Types
        if LP._variables_type[i]==INTEGER:
            glp_set_col_kind(lp, i+1, GLP_IV)

        elif LP._variables_type[i]==BINARY:
            glp_set_col_kind(lp, i+1, GLP_BV)

        else:
            glp_set_col_kind(lp, i+1, GLP_CV)

    for i in range(len(LP._constraints_bounds_max)):
        if LP._constraints_bounds_min[i] is not None and LP._constraints_bounds_max[i] is not None:
            if LP._constraints_bounds_min[i] == LP._constraints_bounds_max[i]:
                glp_set_row_bnds(lp, i+1, GLP_FX, LP._constraints_bounds_min[i], LP._constraints_bounds_max[i])
            else:
                glp_set_row_bnds(lp, i+1, GLP_DB, LP._constraints_bounds_min[i], LP._constraints_bounds_max[i])

        elif LP._constraints_bounds_min[i] is None and LP._constraints_bounds_max[i] is not None:
            glp_set_row_bnds(lp, i+1, GLP_UP, 0, LP._constraints_bounds_max[i])

        elif LP._constraints_bounds_min[i] is not None and LP._constraints_bounds_max[i] is None:
            glp_set_row_bnds(lp, i+1, GLP_LO, LP._constraints_bounds_min[i], 0)

        else:
            glp_set_row_bnds(lp, i+1, GLP_FR, 0, 0)

    # Direction of the problem : maximization, minimization...
    if LP._maximization is False:
        glp_set_obj_dir(lp, GLP_MIN)
    else:
        glp_set_obj_dir(lp, GLP_MAX)

    # Writing the names if needed

    if names:
        if LP._name is not None:
            glp_set_prob_name(lp, LP._name)
        if LP._objective_name is not None:
            glp_set_obj_name(lp, LP._objective_name)

        for c,s in enumerate(LP._constraints_name):
            if s is not None:
                glp_set_row_name(lp,c+1,s)

        for c,s in enumerate(LP._variables_name):
            if s is not None:
                glp_set_col_name(lp,c+1,s)

    iocp.presolve = GLP_ON
    # stuff related to the loglevel
    if log == 0:
        iocp.msg_lev = GLP_MSG_OFF
    elif log == 1:
        iocp.msg_lev = GLP_MSG_ERR
    elif log == 2:
        iocp.msg_lev = GLP_MSG_ON
    else:
        iocp.msg_lev = GLP_MSG_ALL

    sage_free(c_matrix_i)
    sage_free(c_matrix_j)
    sage_free(c_matrix_values)
    return 0
