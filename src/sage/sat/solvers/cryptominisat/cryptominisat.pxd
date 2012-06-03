from decl cimport Solver, Var, Lit

cdef class CryptoMiniSat:
    cdef Solver *_solver
