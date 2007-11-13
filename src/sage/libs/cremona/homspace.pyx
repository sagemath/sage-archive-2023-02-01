cdef extern from "cremona/moddata.h":
    pass
cdef extern from "cremona/symb.h":
    pass
cdef extern from "cremona/cusp.h":
    pass

include "../../ext/interrupt.pxi"

from sage.matrix.all import MatrixSpace
from sage.rings.all import GF
from sage.misc.all import cputime

characteristic = 134217689

cdef extern from "cremona/homspace.h":
    # From mat.h
    ctypedef int scalar "const int"
    ctypedef struct mat "mat":
        scalar* get_entries()   # TODO: possibly not int --

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
    void delete "delete "(void *o)

cdef class ModularSymbolsSpace:
    """
    Class of Cremona Modular Symbols
    """
    cdef homspace* H

    def __init__(self, long level, int sign=0, bint cusp_only=False, int verbose=0):
        if not (sign == 0 or sign==1):
            raise ValueError, "sign %s is not supported; use 0 or +1"%sign
        if verbose:
            print "WARNING: verbose loging currently broken!"
        _sig_on
        self.H = new_homspace(level, sign, cusp_only, verbose)
        _sig_off

    def __dealloc__(self):
        delete(<void*> self.H)

    def __repr__(self):
        return "Cremona %sModular Symbols space of dimension %s for Gamma_0(%s) of weight 2 with sign %s over GF(2^27-39)"%(
            'Cuspidal ' if self.is_cuspidal() else '',
            self.dimension(), self.level(), self.sign())

    cpdef long level(self):
        return self.H.modulus

    cpdef int dimension(self):
        return self.H.h1dim()

    cpdef int sign(self):
        return self.H.plusflag

    cpdef bint is_cuspidal(self):
        return self.H.cuspidal

    def hecke_matrix(self, long p, transpose=False, verbose=False, sparse=False, convert=True):
        cdef mat M
        cdef int i, j, k, n
        cdef scalar* v
        if sparse:
            raise NotImplementedError
        else:
            t = cputime()
            _sig_on
            M = self.H.heckeop(p, transpose, verbose)
            _sig_off
            print "Time to compute matrix: ", cputime(t)
            print "rest is overhead..."
            if not convert: return
            v = M.get_entries()
            n = self.dimension()
            T = MatrixSpace(GF(characteristic), n).zero_matrix()
            k = 0
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    T[i,j] = v[k]
                    k += 1
            return T





