include 'mipGlpk.pxd'
include '../../../../devel/sage/sage/ext/stdsage.pxi'

cdef int BINARY = 1
cdef int REAL = -1
cdef int INTEGER = 0


def solve_glpk(self, log=False,objective_only=False):
    r"""
    Solves the MIP using GLPK. Use solve() instead.
    """

    # Raises an exception if no objective has been defined
    if self._objective_i == None:
        raise Exception("No objective function has been defined!")

    # Creates the structures
    cdef c_glp_prob * lp
    lp = glp_create_prob()
    cdef c_glp_iocp * iocp
    iocp = new_c_glp_iocp()
    glp_init_iocp(iocp)


    # Builds the GLPK version of the problem
    build_glp_prob(lp, iocp, self,log,False)

    # Solves the problem
    glp_intopt(lp, iocp)

    cdef int status = glp_mip_status(lp)
    if status == GLP_UNDEF:
        glp_delete_prob(lp)
        from sage.numerical.mip import MIPSolverException
        raise MIPSolverException, "GLPK : Solution is undefined"
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
        raise MIPSolverException, "GLPK : Feasible solution found, while optimality has not been proven"
    elif status == GLP_NOFEAS:
        glp_delete_prob(lp)
        from sage.numerical.mip import MIPSolverException
        raise MIPSolverException, "GLPK : There is no feasible integer solution to this Linear Program"
    else:
        glp_delete_prob(lp)
        raise ValueError, "The value returned by the solver is unknown !"



def write_mps(self, filename,modern=True):
    r"""
    Writes the MPS file corresponding to the LP
    """

    # Raises an exception if no objective has been defined
    if self._objective_i == None:
        raise Exception("No objective function has been defined!")

    # Creates the structures
    cdef c_glp_prob * lp
    lp = glp_create_prob()
    cdef c_glp_iocp * iocp
    iocp = new_c_glp_iocp()
    glp_init_iocp(iocp)

    # Builds the GLPK version of the problem
    build_glp_prob(lp, iocp, self,False,True)

    if modern:
        glp_write_mps(lp,GLP_MPS_FILE,NULL,filename)
    else:
        glp_write_mps(lp,GLP_MPS_DECK,NULL,filename)

def write_lp(self, filename):

    # Raises an exception if no objective has been defined
    if self._objective_i == None:
        raise Exception("No objective function has been defined!")

    # Creates the structures
    cdef c_glp_prob * lp
    lp = glp_create_prob()
    cdef c_glp_iocp * iocp
    iocp = new_c_glp_iocp()
    glp_init_iocp(iocp)

    # Builds the GLPK version of the problem
    build_glp_prob(lp, iocp, self,False,True)

    glp_write_lp(lp,NULL,filename)


cdef int build_glp_prob(c_glp_prob * lp, c_glp_iocp * iocp, LP, int log, bool names):
    r"""
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
    cdef int c = 1
    for v in LP._constraints_matrix_i:
        c_matrix_i[c]=v+1
        c=c+1

    c=1
    for v in LP._constraints_matrix_j:
        c_matrix_j[c]=v+1
        c=c+1

    c=1
    for v in LP._constraints_matrix_values:
        c_matrix_values[c]=v
        c=c+1


    glp_load_matrix(lp, len(LP._constraints_matrix_i),c_matrix_i,c_matrix_j,c_matrix_values)


    for 0<= i < len(LP._objective_i):
        glp_set_obj_coef(lp, LP._objective_i[i]+1, LP._objective_values[i])


    for 0<= i < len(LP._variables_bounds_min):
        # Upper and lower bounds on the variables
        if LP._variables_bounds_min[i] != None and LP._variables_bounds_max[i] != None:
            glp_set_col_bnds(lp, i+1, GLP_DB, LP._variables_bounds_min[i], LP._variables_bounds_max[i])
        elif LP._variables_bounds_min[i] == None and LP._variables_bounds_max[i] != None:
            glp_set_col_bnds(lp, i+1, GLP_UP, 0, LP._variables_bounds_max[i])
        elif LP._variables_bounds_min[i] != None and LP._variables_bounds_max[i] == None:
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



    for 0<= i < len(LP._constraints_bounds_max):
        if LP._constraints_bounds_min[i] != None and LP._constraints_bounds_max[i] != None:
            if LP._constraints_bounds_min[i] == LP._constraints_bounds_max[i]:
                glp_set_row_bnds(lp, i+1, GLP_FX, LP._constraints_bounds_min[i], LP._constraints_bounds_max[i])
            else:
                glp_set_row_bnds(lp, i+1, GLP_DB, LP._constraints_bounds_min[i], LP._constraints_bounds_max[i])
        elif LP._constraints_bounds_min[i] == None and LP._constraints_bounds_max[i] != None:
            glp_set_row_bnds(lp, i+1, GLP_UP, 0, LP._constraints_bounds_max[i])
        elif LP._constraints_bounds_min[i] != None and LP._constraints_bounds_max[i] == None:
            glp_set_row_bnds(lp, i+1, GLP_LO, LP._constraints_bounds_min[i], 0)
        else:
            glp_set_row_bnds(lp, i+1, GLP_FR, 0, 0)



    # Direction of the problem : maximization, minimization...
    if LP._maximization == False:
        glp_set_obj_dir(lp, GLP_MIN)
    else:
        glp_set_obj_dir(lp, GLP_MAX)

    # Writing the names if needed

    if names:
        if LP._name != None:
            glp_set_prob_name(lp, LP._name)
        if LP._objective_name != None:
            glp_set_obj_name(lp, LP._objective_name)


        c=1
        for s in LP._constraints_name:
            if s!=None:
                glp_set_row_name(lp,c,s)
            c+=1
        c=1
        for s in LP._variables_name:
            if s!=None:
                glp_set_col_name(lp,c,s)
            c+=1



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


