from mat cimport MatrixFactory

cdef MatrixFactory MF = MatrixFactory()

cdef class ModularSymbols:
    """
    Class of Cremona Modular Symbols of given level, sign, etc.
    """
    def __init__(self, long level, int sign=0, bint cuspidal=False, int verbose=0):
        if not (sign == 0 or sign==1):
            raise ValueError, "sign %s is not supported; use 0 or +1"%sign
        if verbose:
            print "WARNING: verbose loging currently broken!"
        _sig_on
        self.H = new_homspace(level, sign, cuspidal, verbose)
        _sig_off

    def __dealloc__(self):
        delete_homspace(self.H)

    def __repr__(self):
        return "Cremona %sModular Symbols space of dimension %s for Gamma_0(%s) of weight 2 with sign %s"%(
            'Cuspidal ' if self.is_cuspidal() else '',
            self.dimension(), self.level(), self.sign())

    #cpdef long level(self):
    def level(self):
        return self.H.modulus

    #cpdef int dimension(self):
    def dimension(self):
        return self.H.h1dim()

    #cpdef int sign(self):
    def sign(self):
        return self.H.plusflag

    #cpdef bint is_cuspidal(self):
    def is_cuspidal(self):
        return self.H.cuspidal

    def hecke_matrix(self, long p, dual=False, verbose=False):
        """
        Return the matrix of the p-th Hecke operator acting on this space
        of modular symbols.

        INPUT:
            p -- a prime number, coprime to the level
            dual -- (default: False)
                    whether to compute the Hecke operator acting on the dual space, i.e., the
                    transpose of the Hecke operator
            verbose -- (default: False) print verbose output

        EXAMPLES:

        """
        if self.level() % p == 0:
            raise NotImplementedError, "computation of Hecke operator only implemented for p coprime to the level."
        _sig_on
        cdef mat M = self.H.heckeop(p, dual, verbose)
        _sig_off
        return MF.new_matrix(M)

