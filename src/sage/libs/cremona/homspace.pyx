cdef extern from "cremona/moddata.h":
    pass
cdef extern from "cremona/symb.h":
    pass
cdef extern from "cremona/cusp.h":
    pass

include "../../ext/interrupt.pxi"

from sage.matrix.all import MatrixSpace
from sage.rings.all import QQ
from sage.misc.all import cputime
from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse
from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense
from sage.rings.rational cimport Rational

characteristic = 134217689

cdef extern from "cremona/homspace.h":


    # From mat.h
    ctypedef int scalar   # TODO: int or long??

    ctypedef struct mat "mat":
        scalar* get_entries()   # TODO: possibly not int --

    long nrows(mat M)
    long ncols(mat M)

    # Constructors
    mat *new_mat "new mat" (mat m)

    # General C++ stuff
    void delete_mat "delete "(mat* m)

    #########################

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


cdef class Matrix

cdef class ModularSymbols:
    """
    Class of Cremona Modular Symbols of given level, sign, etc.
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
        delete_homspace(self.H)

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
        return new_Matrix(M)

cdef class Matrix:
    cdef mat* M

    def __init__(self):
        self.M = NULL

    def __repr__(self):
        return "Cremona matrix"

    cdef set(self, mat*  M):
        if self.M:
            raise RuntimeError, "self.M is already set."
        self.M = M

    def __dealloc__(self):
        if self.M:
            delete_mat(self.M)

    cpdef long nrows(self):
        return nrows(self.M[0])

    cpdef long ncols(self):
        return ncols(self.M[0])

    def sage_matrix_over_QQ(self, sparse=True):
        """
        Return corresponding Sage matrix over the rational numbers.

        INPUTS:
            sparse -- (default: True) whether the return matrix has a sparse representation
        """
        cdef long n = self.nrows()
        cdef long i, j, k
        cdef scalar* v = <scalar*> self.M.get_entries()   # coercion needed to deal with const

        cdef Matrix_rational_sparse T
        cdef Matrix_rational_dense Td

        # Ugly code -- and this is very slow -- it takes longer than computing
        # the matrix in the first splace is not so sparse...
        if sparse:
            T = MatrixSpace(QQ, n, sparse=sparse).zero_matrix()
            k = 0
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    if v[k]:
                        T.set_unsafe(i, j, Rational(v[k]))
                    k += 1
            return T
        else:
            Td = MatrixSpace(QQ, n, sparse=sparse).zero_matrix()
            k = 0
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    if v[k]:
                        Td.set_unsafe(i, j, Rational(v[k]))
                    k += 1
            return Td


cdef new_Matrix(mat M):
    cdef Matrix A = Matrix()
    A.set(new_mat(M))
    return A
