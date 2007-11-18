cdef extern from "cremona/moddata.h":
    pass
cdef extern from "cremona/symb.h":
    pass
cdef extern from "cremona/cusp.h":
    pass

include "../../ext/interrupt.pxi"

cdef extern from "cremona/homspace.h":
    # From mat.h
    ctypedef struct mat "mat":
        pass

    # From homspace.h
    ctypedef struct homspace "homspace":
        # attributes
        long modulus
        int plusflag
        int cuspidal

        # member functions
        int (* h1dim)()
        mat (* heckeop)(long p, int dual, int display)

    # Constructors
    homspace *new_homspace "new homspace" (long n, int hp, int hcusp, int verbose)

    # General C++ stuff
    void delete_homspace "delete "(homspace* H)

cdef class ModularSymbols:
    cdef homspace* H
##     cpdef long level(self)
##     cpdef int dimension(self)
##     cpdef int sign(self)
##     cpdef bint is_cuspidal(self)


