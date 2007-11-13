cdef extern from "cremona/moddata.h":
    pass
cdef extern from "cremona/symb.h":
    pass
cdef extern from "cremona/cusp.h":
    pass

cdef extern from "cremona/homspace.h":
    ctypedef struct homspace "homspace":
        int (* h1dim)()

    homspace *new_homspace "new homspace" (long n, int hp, int hcusp, int verbose)

    void delete "delete "(void *o)


def dim(n, hp, hcusp, verbose):
    cdef homspace* H = new_homspace(n, hp, hcusp, verbose)
    cdef int d = H.h1dim()
    print "Dimension = ", d

