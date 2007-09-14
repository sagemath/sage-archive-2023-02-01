#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "../../ext/interrupt.pxi"
include "../../ext/stdsage.pxi"
include 'misc.pxi'
include 'decl.pxi'

from sage.libs.ntl.ntl_ZZ cimport ntl_ZZ
from sage.libs.ntl.ntl_ZZX cimport ntl_ZZX

from ntl_ZZ import unpickle_class_args

cdef make_ZZ(ZZ_c* x):
    cdef ntl_ZZ y
    y = ntl_ZZ()
    y.x = x[0]
    ZZ_delete(x)
    _sig_off
    return y

cdef make_mat_ZZ(mat_ZZ_c* x):
    cdef ntl_mat_ZZ y
    _sig_off
    y = ntl_mat_ZZ(_INIT)
    y.x = x[0]
    mat_ZZ_delete(x)
    y.__nrows = mat_ZZ_nrows(&y.x);
    y.__ncols = mat_ZZ_ncols(&y.x);
    return y

##############################################################################
#
# ntl_mat_ZZ: Matrices over the integers via NTL
#
##############################################################################

cdef class ntl_mat_ZZ:
    # see ntl.pxd for data members
    r"""
    The \class{mat_ZZ} class implements arithmetic with matrices over $\Z$.
    """
    def __init__(self, nrows=0,  ncols=0, v=None):
        if nrows == _INIT:
            return
        cdef unsigned long i, j
        cdef ntl_ZZ tmp
        if nrows == 0 and ncols == 0:
            return
        nrows = int(nrows)
        ncols = int(ncols)
        mat_ZZ_SetDims(&self.x, nrows, ncols)
        self.__nrows = nrows
        self.__ncols = ncols
        if v != None:
            for i from 0 <= i < nrows:
                for j from 0 <= j < ncols:
                    tmp = ntl_ZZ(v[i*ncols+j])
                    mat_ZZ_setitem(&self.x, i, j, &tmp.x)

    def __reduce__(self):
        """
        TEST:
            sage: m = ntl.mat_ZZ(3, 2, range(6)); m
            [[0 1]
            [2 3]
            [4 5]
            ]
            sage: loads(dumps(m))
            [[0 1]
            [2 3]
            [4 5]
            ]
        """
        return unpickle_class_args, (ntl_mat_ZZ, (self.__nrows, self.__ncols, self.list()))

    def __new__(self, nrows=0,  ncols=0, v=None):
        mat_ZZ_construct(&self.x)

    def __dealloc__(self):
        mat_ZZ_destruct(&self.x)

    def __repr__(self):
        _sig_on
        return string_delete(mat_ZZ_to_str(&self.x))

    def __mul__(ntl_mat_ZZ self, other):
        cdef ntl_mat_ZZ r = PY_NEW(ntl_mat_ZZ)
        if not PY_TYPE_CHECK(self, ntl_mat_ZZ):
            self = ntl_mat_ZZ(self)
        if not PY_TYPE_CHECK(other, ntl_mat_ZZ):
            other = ntl_mat_ZZ(other)
        _sig_on
        mat_ZZ_mul(r.x, (<ntl_mat_ZZ>self).x, (<ntl_mat_ZZ>other).x)
        _sig_off
        return r

    def __sub__(ntl_mat_ZZ self, other):
        cdef ntl_mat_ZZ r = PY_NEW(ntl_mat_ZZ)
        if not PY_TYPE_CHECK(self, ntl_mat_ZZ):
            self = ntl_mat_ZZ(self)
        if not PY_TYPE_CHECK(other, ntl_mat_ZZ):
            other = ntl_mat_ZZ(other)
        _sig_on
        mat_ZZ_sub(r.x, (<ntl_mat_ZZ>self).x, (<ntl_mat_ZZ>other).x)
        _sig_off
        return r

    def __add__(ntl_mat_ZZ self, other):
        cdef ntl_mat_ZZ r = PY_NEW(ntl_mat_ZZ)
        if not PY_TYPE_CHECK(self, ntl_mat_ZZ):
            self = ntl_mat_ZZ(self)
        if not PY_TYPE_CHECK(other, ntl_mat_ZZ):
            other = ntl_mat_ZZ(other)
        _sig_on
        mat_ZZ_add(r.x, (<ntl_mat_ZZ>self).x, (<ntl_mat_ZZ>other).x)
        _sig_off
        return r

    def __pow__(ntl_mat_ZZ self, long e, ignored):
        cdef ntl_mat_ZZ r = PY_NEW(ntl_mat_ZZ)
        _sig_on
        mat_ZZ_power(r.x, (<ntl_mat_ZZ>self).x, e)
        _sig_off
        return r

    def nrows(self):
        return self.__nrows

    def ncols(self):
        return self.__ncols

    def __setitem__(self, ij, x):
        cdef ntl_ZZ y
        cdef int i, j
        if not isinstance(x, ntl_ZZ):
            y = ntl_ZZ(x)
        else:
            y = x
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, 'ij must be a 2-tuple'
        i, j = int(ij[0]),int(ij[1])
        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError, "array index out of range"
        _sig_on
        mat_ZZ_setitem(&self.x, i, j, &y.x)
        _sig_off

    def __getitem__(self, ij):
        """
            sage: m = ntl.mat_ZZ(3, 2, range(6))
            sage: m[0,0]
            0
            sage: m[2,1]
            5
            sage: m[3,2] #  oops, 0 based
            Traceback (most recent call last):
            ...
            IndexError: array index out of range
        """
        cdef int i, j
        if not isinstance(ij, tuple) or len(ij) != 2:
            raise TypeError, 'ij must be a 2-tuple'
        i, j = ij
        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError, "array index out of range"
        return make_ZZ(mat_ZZ_getitem(&self.x, i+1, j+1))

    def list(self):
        """
        EXAMPLE:
            sage: m = ntl.mat_ZZ(3, 4, range(12)); m
            [[0 1 2 3]
            [4 5 6 7]
            [8 9 10 11]
            ]
            sage: m.list()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        """
        cdef int i, j
        return [make_ZZ(mat_ZZ_getitem(&self.x, i+1, j+1))
                    for i from 0 <= i < self.__nrows
                        for j from 0 <= j < self.__ncols]

    def determinant(self, deterministic=True):
        _sig_on
        return make_ZZ(mat_ZZ_determinant(&self.x, deterministic))

    def HNF(self, D=None):
        r"""
        The input matrix A=self is an n x m matrix of rank m (so n >=
        m), and D is a multiple of the determinant of the lattice L
        spanned by the rows of A.  The output W is the Hermite Normal
        Form of A; that is, W is the unique m x m matrix whose rows
        span L, such that

        - W is lower triangular,
        - the diagonal entries are positive,
        - any entry below the diagonal is a non-negative number
          strictly less than the diagonal entry in its column.

        This is implemented using the algorithm of [P. Domich,
        R. Kannan and L. Trotter, Math. Oper. Research 12:50-59,
        1987].

        TIMINGS:
        NTL isn't very good compared to MAGMA, unfortunately:

            sage.: import ntl
            sage.: a=MatrixSpace(Q,200).random_element()    # -2 to 2
            sage.: A=ntl.mat_ZZ(200,200)
            sage.: for i in xrange(a.nrows()):
               ....:     for j in xrange(a.ncols()):
               ....:         A[i,j] = a[i,j]
               ....:
            sage.: time d=A.determinant()
            Time.: 3.89 seconds
            sage.: time B=A.HNF(d)
            Time.: 27.59 seconds

        In comparison, MAGMA does this much more quickly:
        \begin{verbatim}
            > A := MatrixAlgebra(Z,200)![Random(-2,2) : i in [1..200^2]];
            > time d := Determinant(A);
            Time: 0.710
            > time H := HermiteForm(A);
            Time: 3.080
        \end{verbatim}

        Also, PARI is also faster than NTL if one uses the flag 1 to
        the mathnf routine.  The above takes 16 seconds in PARI.
        """
        cdef ntl_ZZ _D
        if D == None:
            _D = self.determinant()
        else:
            _D = ntl_ZZ(D)
        _sig_on
        return make_mat_ZZ(mat_ZZ_HNF(&self.x, &_D.x))

    def charpoly(self):
        """
        EXAMPLES:
            sage: ntl.mat_ZZ(2,2,[1,2,3,4]).determinant()
            -2
        """
        cdef ntl_ZZX r = ntl_ZZX()
        _sig_on
        mat_ZZ_CharPoly(r.x, self.x)
        _sig_off
        return r

    def LLL(self, a=3, b=4, return_U=False, verbose=False):
        r"""
        Performs LLL reduction of self (puts self in an LLL form).

        self is an m x n matrix, viewed as m rows of n-vectors.  m may
        be less than, equal to, or greater than n, and the rows need
        not be linearly independent. self is transformed into an
        LLL-reduced basis, and the return value is the rank r of self
        so as det2 (see below).  The first m-r rows of self are zero.

        More specifically, elementary row transformations are
        performed on self so that the non-zero rows of new-self form
        an LLL-reduced basis for the lattice spanned by the rows of
        old-self.  The default reduction parameter is $\delta=3/4$,
        which means that the squared length of the first non-zero
        basis vector is no more than $2^{r-1}$ times that of the
        shortest vector in the lattice.

        det2 is calculated as the \emph{square} of the determinant of
        the lattice---note that sqrt(det2) is in general an integer
        only when r = n.

        If return_U is True, a value U is returned which is the
        transformation matrix, so that U is a unimodular m x m matrix
        with U * old-self = new-self.  Note that the first m-r rows of
        U form a basis (as a lattice) for the kernel of old-B.

        The parameters a and b allow an arbitrary reduction parameter
        $\delta=a/b$, where $1/4 < a/b \leq 1$, where a and b are positive
        integers.  For a basis reduced with parameter delta, the
        squared length of the first non-zero basis vector is no more
        than $1/(delta-1/4)^{r-1}$ times that of the shortest vector in
        the lattice.

        The algorithm employed here is essentially the one in Cohen's
        book: [H. Cohen, A Course in Computational Algebraic Number
        Theory, Springer, 1993]

        INPUT:
           a        -- parameter a as described above (default: 3)
           b        -- parameter b as described above (default: 4)
           return_U -- return U as described above
           verbose  -- if True NTL will produce some verbatim messages on
                       what's going on internally (default: False)

        OUTPUT:
            (rank,det2,[U]) where rank,det2, and U are as described
            above and U is an optional return value if return_U is
            True.

        EXAMPLE:
            sage: M=ntl.mat_ZZ(3,3,[1,2,3,4,5,6,7,8,9])
            sage: M.LLL()
            (2, 54)
            sage: M
            [[0 0 0]
            [2 1 0]
            [-1 1 3]
            ]
            sage: M=ntl.mat_ZZ(4,4,[-6,9,-15,-18,4,-6,10,12,10,-16,18,35,-24,36,-46,-82]); M
            [[-6 9 -15 -18]
            [4 -6 10 12]
            [10 -16 18 35]
            [-24 36 -46 -82]
            ]
            sage: M.LLL()
            (3, 19140)
            sage: M
            [[0 0 0 0]
            [0 -2 0 0]
            [-2 1 -5 -6]
            [0 -1 -7 5]
            ]

        WARNING: This method modifies self. So after applying this method your matrix
        will be a vector of vectors.
        """
        cdef ZZ_c *det2
        cdef ntl_mat_ZZ U
        if return_U:
            U = PY_NEW(ntl_mat_ZZ)
            _sig_on
            rank = int(mat_ZZ_LLL_U(&det2, &self.x, &U.x, int(a), int(b), int(verbose)))
            _sig_off
            return rank, make_ZZ(det2), U
        else:
            _sig_on
            rank = int(mat_ZZ_LLL(&det2,&self.x,int(a),int(b),int(verbose)))
            _sig_off
            return rank,make_ZZ(det2)
