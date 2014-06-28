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

##############################################################################
#
# ntl_mat_GF2E: Matrices over the GF(2**x) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version (based on code by William Stein)
#
##############################################################################

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
include 'sage/ext/random.pxi'
include 'misc.pxi'
include 'decl.pxi'

from ntl_GF2E cimport ntl_GF2E
from ntl_GF2EContext import ntl_GF2EContext
from ntl_GF2EContext cimport ntl_GF2EContext_class
from sage.rings.integer cimport Integer

from sage.libs.ntl.ntl_ZZ import unpickle_class_args

cdef class ntl_mat_GF2E:
    r"""
    The \class{mat_GF2E} class implements arithmetic with matrices over $GF(2**x)$.
    """
    def __init__(self, modulus = None, nrows=0, ncols=0, v=None):
        """
        Constructs a matrix over ntl.GF2E.

        INPUT:
            modulus -- GF2E context
            nrows -- number of rows
            ncols -- number of columns
            v     -- either a list or a matrix over GF(2^x)

        EXAMPLES::

            sage: k.<a> = GF(2^4)
            sage: ctx = ntl.GF2EContext(k)
            sage: ntl.GF2XHexOutput(1)
            sage: ntl.mat_GF2E(ctx, 5,5, [0..24])
            [[0x0 0x1 0x2 0x3 0x4]
            [0x5 0x6 0x7 0x8 0x9]
            [0xa 0xb 0xc 0xd 0xe]
            [0xf 0x3 0x2 0x1 0x0]
            [0x7 0x6 0x5 0x4 0xb]
            ]
            sage: ntl.mat_GF2E(ctx, 5,5)
            [[0x0 0x0 0x0 0x0 0x0]
            [0x0 0x0 0x0 0x0 0x0]
            [0x0 0x0 0x0 0x0 0x0]
            [0x0 0x0 0x0 0x0 0x0]
            [0x0 0x0 0x0 0x0 0x0]
            ]
            sage: A= matrix(k,5,5,[k.fetch_int(_%(2^4)) for _ in range(25)])
            sage: ntl.mat_GF2E(ctx, A)
            [[0x0 0x1 0x2 0x3 0x4]
            [0x5 0x6 0x7 0x8 0x9]
            [0xa 0xb 0xc 0xd 0xe]
            [0xf 0x0 0x1 0x2 0x3]
            [0x4 0x5 0x6 0x7 0x8]
            ]
        """
        if modulus is None:
            raise ValueError, "You must specify a modulus when creating a GF2E."

        cdef unsigned long _nrows, _ncols
        cdef unsigned long i, j

        import sage.matrix.matrix
        if sage.matrix.matrix.is_Matrix(nrows):
            _nrows = nrows.nrows()
            _ncols = nrows.ncols()
            v     = nrows.list()
        else:
            _nrows = nrows
            _ncols = ncols

        self.x.SetDims(_nrows, _ncols)
        if v is not None:
            sig_on()
            for i from 0 <= i < _nrows:
                for j from 0 <= j < _ncols:
                    elem = v[i*_ncols+j]
                    if not PY_TYPE_CHECK(elem, ntl_GF2E):
                        elem = ntl_GF2E(elem, modulus)
                    mat_GF2E_setitem(&self.x, i, j, &(<ntl_GF2E>elem).x)
            sig_off()

    def __cinit__(self, modulus=None, nrows=0, ncols=0, v=None):
        #################### WARNING ###################
        ## Before creating a GF2E, you must create a  ##
        ## GF2EContext, and restore it.  In Python,   ##
        ## the error checking in __init__ will prevent##
        ## you from constructing an ntl_GF2E          ##
        ## inappropriately.  However, from Cython, you##
        ## could do r = PY_NEW(ntl_GF2E) without      ##
        ## first restoring a GF2EContext, which could ##
        ## have unfortunate consequences.  See _new  ##
        ## defined below for an example of the right  ##
        ## way to short-circuit __init__ (or just call##
        ## _new in your own code).                    ##
        ################################################
        if modulus is None:
            return
        if PY_TYPE_CHECK( modulus, ntl_GF2EContext_class ):
            self.c = <ntl_GF2EContext_class>modulus
            self.c.restore_c()
            mat_GF2E_construct(&self.x)
        else:
            self.c = <ntl_GF2EContext_class>ntl_GF2EContext(modulus)
            self.c.restore_c()
            mat_GF2E_construct(&self.x)

    cdef ntl_GF2E _new_element(self):
        cdef ntl_GF2E r
        self.c.restore_c()
        r = PY_NEW(ntl_GF2E)
        r.c = self.c
        return r

    cdef ntl_mat_GF2E _new(self):
        cdef ntl_mat_GF2E r
        self.c.restore_c()
        r = PY_NEW(ntl_mat_GF2E)
        r.x.SetDims(self.x.NumRows(),self.x.NumCols())
        r.c = self.c
        return r

    def modulus_context(self):
        """
        Returns the structure that holds the underlying NTL GF2E modulus.

        EXAMPLES::

            sage: ntl.GF2XHexOutput(0)
            sage: ctx = ntl.GF2EContext( ntl.GF2X([1,1,0,1,1,0,0,0,1]) )
            sage: a = ntl.GF2E(ntl.ZZ_pX([1,1,3],2), ctx)
            sage: A= ntl.mat_GF2E(ctx, 1, 1, [a])
            sage: cty = A.modulus_context(); cty
            NTL modulus [1 1 0 1 1 0 0 0 1]
            sage: ctx == cty
            True
        """
        return self.c

    def __dealloc__(self):
        # With NTL 6.0.0, mat_GF2E is a proper C++ class.
        # Therefore Cython automagically calls the class destructor.
        #mat_GF2E_destruct(&self.x)
        pass

    def __reduce__(self):
        """
        EXAMPLE::

            sage: k.<a> = GF(2^4)
            sage: ctx = ntl.GF2EContext(k)
            sage: A = ntl.mat_GF2E(ctx, 5,5, [0..24])
            sage: A == loads(dumps(A))
            True
        """
        return unpickle_class_args, (ntl_mat_GF2E, (self.modulus_context(), self.x.NumRows(), self.x.NumCols(), self.list()))

    def __repr__(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: ntl.GF2XHexOutput(1)
            sage: ntl.mat_GF2E(ctx, 2,2,range(4)).__repr__()
            '[[0x0 0x1]\n[0x0 0x1]\n]'
            sage: ntl.GF2XHexOutput(0)
            sage: ntl.mat_GF2E(ctx, 2,2,range(4)).__repr__()
            '[[[] [1]]\n[[] [1]]\n]'
        """
        self.c.restore_c()
        return mat_GF2E_to_PyString(&self.x)

    def __mul__(ntl_mat_GF2E self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: ntl.GF2XHexOutput(1)
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: n = ntl.mat_GF2E(ctx, 5,5,[3..27])
            sage: m*n ## indirect doctest
            [[0x87 0x04 0xc4 0xc7 0x87]
            [0x32 0x84 0x17 0x63 0x73]
            [0xa1 0x46 0x25 0xcd 0x2f]
            [0x1 0xcf 0xfb 0xd6 0x62]
            [0xcf 0x02 0x06 0xfd 0x79]
            ]
        """
        cdef ntl_mat_GF2E r = self._new()
        if not PY_TYPE_CHECK(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other, self.c)
        if not self.c is (<ntl_mat_GF2E>other).c:
            raise ValueError, "You can not perform arithmetic with matrices over different fields."
        sig_on()
        mat_GF2E_mul(r.x, self.x, (<ntl_mat_GF2E>other).x)
        sig_off()
        return r

    def __sub__(ntl_mat_GF2E self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: n = ntl.mat_GF2E(ctx, 5,5,[3..27])
            sage: ntl.GF2XHexOutput(0)
            sage: m-n ## indirect doctest
            [[[1 1] [1 0 1] [1 1 1] [1 0 1] [1 1]]
            [[1 0 1 1] [1 1 1 1] [1 0 1 1] [1 1] [1 0 1]]
            [[1 1 1] [1 0 1] [1 1] [1 0 1 1 1] [1 1 1 1 1]]
            [[1 0 1 1 1] [1 1] [1 0 1] [1 1 1] [1 0 1]]
            [[1 1] [1 0 1 1] [1 1 1 1] [1 0 1 1] [1 1]]
            ]
        """
        cdef ntl_mat_GF2E r = self._new()
        if not PY_TYPE_CHECK(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other, self.c)
        if not self.c is (<ntl_mat_GF2E>other).c:
            raise ValueError, "You can not perform arithmetic with matrices over different fields."
        sig_on()
        mat_GF2E_sub(r.x, self.x, (<ntl_mat_GF2E>other).x)
        sig_off()
        return r

    def __add__(ntl_mat_GF2E self, other):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: n = ntl.mat_GF2E(ctx, 5,5,[3..27])
            sage: m+n ## indirect doctest
            [[[1 1] [1 0 1] [1 1 1] [1 0 1] [1 1]]
            [[1 0 1 1] [1 1 1 1] [1 0 1 1] [1 1] [1 0 1]]
            [[1 1 1] [1 0 1] [1 1] [1 0 1 1 1] [1 1 1 1 1]]
            [[1 0 1 1 1] [1 1] [1 0 1] [1 1 1] [1 0 1]]
            [[1 1] [1 0 1 1] [1 1 1 1] [1 0 1 1] [1 1]]
            ]
        """
        cdef ntl_mat_GF2E r = self._new()
        if not PY_TYPE_CHECK(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other, self.c)
        if not self.c is (<ntl_mat_GF2E>other).c:
            raise ValueError, "You can not perform arithmetic with matrices over different fields."
        sig_on()
        mat_GF2E_add(r.x, self.x, (<ntl_mat_GF2E>other).x)
        sig_off()
        return r

    def __neg__(ntl_mat_GF2E self):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: -m == m ## indirect doctest
            True
        """
        cdef ntl_mat_GF2E r = self._new()
        sig_on()
        mat_GF2E_negate(r.x, self.x)
        sig_off()
        return r

    def __pow__(ntl_mat_GF2E self, long e, ignored):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: m**2 == m*m ## indirect doctest
            True
        """
        cdef ntl_mat_GF2E r = self._new()
        sig_on()
        mat_GF2E_power(r.x, self.x, e)
        sig_off()
        return r

    def __richcmp__(ntl_mat_GF2E self, other, op):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: n = ntl.mat_GF2E(ctx, 5,5,[3..27])
            sage: m == n ## indirect doctest
            False
            sage: m == m
            True
        """
        self.c.restore_c()

        if not PY_TYPE_CHECK(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other,self.c)

        if op != 2 and op != 3:
            raise TypeError, "Elements in GF2E are not ordered."

        cdef int t
        t = mat_GF2E_equal(self.x, (<ntl_mat_GF2E>other).x)
        if op == 2:
            return t == 1
        elif op == 3:
            return t == 0

    def NumRows(self):
        """
        Return the number of rows in self.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24]) ; m.NumRows()
            5
        """
        return int(self.x.NumRows())

    def NumCols(self):
        """
        Return the number of columns in self.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24]) ; m.NumCols()
            5
        """
        return int(self.x.NumCols())

    def __setitem__(self, ij, x):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: ntl.GF2XHexOutput(0)
            sage: m[0,0]
            []
            sage: m[0,0] = 1
            sage: m[0,0]
            [1]
        """
        cdef int i, j
        if not PY_TYPE_CHECK(x, ntl_GF2E):
            x = ntl_GF2E(x, self.c)

        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.x.NumCols()==1 and (PY_TYPE_CHECK(ij, Integer) or PY_TYPE_CHECK(ij,int)):
            i = ij
            j = 0
        elif self.x.NumRows()==1 and (PY_TYPE_CHECK(ij, Integer) or PY_TYPE_CHECK(ij,int)):
            i = 0
            j = ij
        else:
            raise TypeError, 'ij is not a matrix index'

        if i < 0 or i >= self.x.NumRows() or j < 0 or j >= self.x.NumCols():
            raise IndexError, "array index out of range"

        if not (<ntl_GF2E>x).c is self.c:
            raise ValueError, "You can not assign elements from different fields."

        self.c.restore_c()

        mat_GF2E_setitem(&self.x, i, j, &(<ntl_GF2E>x).x)

    def __getitem__(self, ij):
        """
        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: m[0,1]
            [1]
            sage: m[0,0] = 0
            sage: m[0,0]
            []
        """
        cdef int i, j
        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.x.NumCols() == 1 and (PY_TYPE_CHECK(ij, Integer) or PY_TYPE_CHECK(ij,int)):
            i = ij
            j = 0
        elif self.x.NumRows() == 1 and (PY_TYPE_CHECK(ij, Integer) or PY_TYPE_CHECK(ij,int)):
            i = 0
            j = ij
        else:
            raise TypeError, 'ij is not a matrix index'

        if i < 0 or i >= self.x.NumRows() or j < 0 or j >= self.x.NumCols():
            raise IndexError, "array index out of range"

        cdef ntl_GF2E e = self._new_element()
        e.x = self.x.get( i+1, j+1 )
        return e

    def determinant(self):
        """
        Returns the determinant.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: ntl.GF2XHexOutput(0)
            sage: ntl.mat_GF2E(ctx, 5,5,[0..24]).determinant()
            [0 1 0 1 1 1 1]
            sage: ntl.mat_GF2E(ctx, 5,5,[3..27]).determinant()
            [0 1 1 0 0 1]
        """
        cdef ntl_GF2E r = self._new_element()
        sig_on()
        r.x = mat_GF2E_determinant(self.x)
        sig_off()
        return r

    def gauss(self,ncols=-1):
        """
        Performs unitary row operations so as to bring this matrix
        into row echelon form.  If the optional argument \code{ncols}
        is supplied, stops when first ncols columns are in echelon
        form.  The return value is the rank (or the rank of the first
        ncols columns).

        INPUT:

        - ``ncols`` - number of columns to process (default: all)

        EXAMPLES::

            sage: m = ctx = ntl.GF2EContext(ntl.GF2X([1,1,0,1,1,0,0,0,1]))
            sage: ntl.mat_GF2E(ctx, 5,5,[3..27]).gauss()
            5
            sage: ntl.mat_GF2E(ctx, 5,5).gauss()
            0
            sage: ntl.mat_GF2E(ctx, 5,5,[3..27]).gauss(3)
            3
        """
        if ncols == -1:
            ncols = self.x.NumCols()
        return int(mat_GF2E_gauss(self.x, int(ncols)))

    def list(self):
        """
        Returns a list of the entries in this matrix

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 2,2,[ntl.GF2E_random(ctx) for x in xrange(2*2)])
            sage: ntl.GF2XHexOutput(0)
            sage: m.list()
            [[1 0 1 0 1 0 0 1], [1 0 1 1 1 0 0 1], [0 0 0 1 1 1 1], [1 1 1 1 1 1]]
        """
        return [self[i,j] for i in range(self.NumRows()) for j in range(self.x.NumCols())]

    def IsZero(self):
        """
        Return True if self is zero, and false otherwise.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: n = ntl.mat_GF2E(ctx, 5,5)
            sage: m.IsZero()
            False
            sage: n.IsZero()
            True
        """
        cdef long isZero
        sig_on()
        isZero = mat_GF2E_IsZero(self.x)
        sig_off()
        return bool(isZero)

    def _sage_(ntl_mat_GF2E self, k=None):
        """
        Returns a ``Matrix`` over a ``FiniteField`` representation
        of this element.

        INPUT:

        - ``k`` - optional GF(2**deg)

        OUTPUT:
            Matrix over k

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1])
            sage: m = ntl.mat_GF2E(ctx, 2,2,[3..6])
            sage: ntl.GF2XHexOutput(0)
            sage: m
            [[[1 1] [0 0 1]]
            [[1 0 1] [0 1 1]]
            ]
            sage: m._sage_()
            [  a + 1     a^2]
            [a^2 + 1 a^2 + a]
        """
        if k is None:
            from sage.rings.finite_rings.constructor import FiniteField
            f = self.c.m._sage_()
            e = GF2E_degree()
            k = FiniteField(2**e, name='a', modulus=f)

        l = [e._sage_(k) for e in self.list()] # we actually can do faster than this

        from sage.matrix.constructor import Matrix
        return Matrix(k,self.x.NumRows(),self.x.NumCols(),l)

    def transpose(ntl_mat_GF2E self):
        """
        Returns the transposed matrix of self.

        OUTPUT:
            transposed Matrix

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: n = m.transpose()
            sage: n == m
            False
            sage: n.transpose() == m
            True
        """
        cdef ntl_mat_GF2E r = self._new()
        sig_on()
        mat_GF2E_transpose(r.x, self.x)
        sig_off()
        return r

    def __invert__(self):
        """
        Return $X = A^{-1}$; an error is raised if A is singular.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: n = ~m
            sage: o = n*m
            sage: o.IsIdent()
            True
        """
        cdef ntl_mat_GF2E r = self._new()
        sig_on()
        mat_GF2E_inv(r.x, self.x)
        sig_off()
        return r

    def IsIdent(self, n = -1):
        """
        test if A is the n x n identity matrix

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 5,5,[0..24])
            sage: n = ~m
            sage: o = n*m
            sage: o.IsIdent()
            True
        """
        if n < 0:
            n = self.NumRows()
        return bool(mat_GF2E_IsIdent(self.x, n))

    def IsDiag(self, long n, ntl_GF2E d):
        """
        Test if X is an  n x n diagonal matrix with d on diagonal.

        EXAMPLES::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 3,3,[[0,1],0,0, 0,[0,1],0, 0,0,[0,1]])
            sage: m.IsDiag(2, ntl.GF2E([0,1],ctx))
            False
            sage: m.IsDiag(3, ntl.GF2E([0,1],ctx))
            True
        """
        return bool(mat_GF2E_IsDiag(self.x, n, d.x))

    def image(self):
        """
        The rows of X are computed as basis of A's row space.  X is is
        row echelon form.

        EXAMPLE::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 3,3,[0..24])
            sage: ntl.GF2XHexOutput(1)
            sage: m.image()
            [[0x3 0x4 0x5]
            [0x0 0x1 0x2]
            [0x0 0x0 0xc1]
            ]
        """
        cdef ntl_mat_GF2E X = self._new()
        sig_on()
        mat_GF2E_image(X.x, self.x)
        sig_off()
        return X

    def kernel(self):
        """
        Computes a basis for the kernel of the map ``x -> x*A``, where
        ``x`` is a row vector.

        EXAMPLE::

            sage: ctx = ntl.GF2EContext([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(ctx, 3,3,[0..24])
            sage: ntl.GF2XHexOutput(1)
            sage: m.kernel()
            []
        """
        cdef ntl_mat_GF2E X = self._new()
        sig_on()
        mat_GF2E_kernel(X.x, self.x)
        sig_off()
        return X

    def randomize(self, density=1, nonzero=False):
        """
        Randomize ``density`` proportion of the entries of this matrix,
        leaving the rest unchanged.

        INPUT:

        -  ``density`` - float; proportion (roughly) to be considered for
           changes
        -  ``nonzero`` - Bool (default: ``False``); whether the new entries
           are forced to be non-zero

        EXAMPLES::

            sage: k.<a> = GF(2^4)
            sage: ctx = ntl.GF2EContext(k)
            sage: ntl.GF2XHexOutput(1)
            sage: A = ntl.mat_GF2E(ctx, 100,100)
            sage: A.randomize()
            sage: len([e for e in A.list() if e!=0])
            9389

            sage: A = ntl.mat_GF2E(ctx, 100,100)
            sage: A.randomize(nonzero=True)
            sage: len([e for e in A.list() if e!=0])
            10000

            sage: A = ntl.mat_GF2E(ctx, 100,100)
            sage: A.randomize(nonzero=True, density=0.1)
            sage: len([e for e in A.list() if e!=0])
            994

        """
        cdef long i,j
        cdef GF2E_c tmp

        cdef float _density = density
        cdef randstate rstate = current_randstate()

        if _density <= 0:
            return
        if _density > 1:
            _density = 1.0

        if not nonzero:
            if _density == 1.0:
                for i in xrange(self.x.NumRows()):
                    for j in xrange(self.x.NumCols()):
                        tmp = GF2E_random()
                        mat_GF2E_setitem(&self.x, i, j, &tmp)
            else:
                for i in xrange(self.x.NumRows()):
                    for j in xrange(self.x.NumCols()):
                        if rstate.c_rand_double() <= _density:
                            tmp = GF2E_random()
                            mat_GF2E_setitem(&self.x, i, j, &tmp)
        else:
            if _density == 1.0:
                for i in xrange(self.x.NumRows()):
                    for j in xrange(self.x.NumCols()):
                        tmp = GF2E_random()
                        while GF2E_IsZero(tmp):
                            tmp = GF2E_random()
                        mat_GF2E_setitem(&self.x, i, j, &tmp)
            else:
                for i in xrange(self.x.NumRows()):
                    for j in xrange(self.x.NumCols()):
                        if rstate.c_rand_double() <= _density:
                            tmp = GF2E_random()
                            while GF2E_IsZero(tmp):
                                tmp = GF2E_random()
                            mat_GF2E_setitem(&self.x, i, j, &tmp)
