include '../../../../devel/sage/sage/ext/stdsage.pxi'

cdef int BINARY = 1
cdef int REAL = -1
cdef int INTEGER = 0



def solve_cplex(self,log=0,objective_only=False):
    r"""
    Solves the ``MixedIntegerLinearProgram`` using CPLEX
    *Use ``solve()`` instead*

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
    function given by CPLEX.

    EXAMPLE:

    Solving a simple Linear Program using CPLEX as a solver
    (Computation of a maximum stable set in Petersen's graph)::

        sage: from sage.numerical.mip_cplex import solve_cplex     # optional - requires Cplex
        sage: g = graphs.PetersenGraph()
        sage: p = MixedIntegerLinearProgram(maximization=True)
        sage: b = p.new_variable()
        sage: p.set_objective(sum([b[v] for v in g]))
        sage: for (u,v) in g.edges(labels=None):
        ...       p.add_constraint(b[u] + b[v], max=1)
        sage: p.set_binary(b)
        sage: solve_cplex(p,objective_only=True)     # optional - requires Cplex
        4.0
    """



    # Creates the solver interfaces
    cdef c_OsiCpxSolverInterface* si
    si = new_c_OsiCpxSolverInterface();

    from sage.numerical.mip import MIPSolverException
    cdef Osi_interface OSI

    OSI = Osi_interface()

    return_value = OSI.osi_solve(self, <c_OsiSolverInterface *> si, objective_only, True)
    del_OsiCpxSolverInterface(si)

    return return_value
