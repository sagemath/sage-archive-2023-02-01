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

from ntl_GF2E cimport ntl_GF2E

import ntl_GF2E as ntl_GF2E_module

from ntl_GF2E import ntl_GF2E_sage

##############################################################################
#
# ntl_mat_GF2E: Matrices over the GF(2**x) via NTL
#
# AUTHORS:
#  - Martin Albrecht <malb@informatik.uni-bremen.de>
#    2006-01: initial version (based on cody by William Stein)
#
##############################################################################


from ntl_ZZ import unpickle_class_args


cdef class ntl_mat_GF2E:
    r"""
    The \class{mat_GF2E} class implements arithmetic with matrices over $GF(2**x)$.
    """
    def __init__(self, nrows=0, ncols=0, v=None):
        """
        Constructs a matrix over ntl.GF2E.

        INPUT:
            nrows -- number of rows
            ncols -- nomber of columns
            v     -- either list or Matrix over GF(2**x)

        EXAMPLES:
            sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
            sage: m=ntl.mat_GF2E(10,10)
            sage: m=ntl.mat_GF2E(Matrix(GF(2**8, 'a'),10,10))
            sage: m=ntl.mat_GF2E(10,10,[ntl.GF2E_random() for x in xrange(10*10)])
        """

        if nrows is _INIT:
            return

        cdef unsigned long _nrows, _ncols

        import sage.matrix.matrix
        if sage.matrix.matrix.is_Matrix(nrows):
            _nrows = nrows.nrows()
            _ncols = nrows.ncols()
            v     = nrows.list()
        else:
            _nrows = nrows
            _ncols = ncols

        if not ntl_GF2E_module.__have_GF2E_modulus:
            raise "NoModulus"
        cdef unsigned long i, j
        cdef ntl_GF2E tmp

        mat_GF2E_SetDims(&self.x, _nrows, _ncols)
        self.__nrows = _nrows
        self.__ncols = _ncols
        if v != None:
            _sig_on
            for i from 0 <= i < _nrows:
                for j from 0 <= j < _ncols:
                    elem = v[i*_ncols+j]
                    if not isinstance(elem, ntl_GF2E):
                        tmp=ntl_GF2E(elem)
                    else:
                        tmp=elem
                    mat_GF2E_setitem(&self.x, i, j, &tmp.gf2e_x)
            _sig_off

    def __new__(self, nrows=0, ncols=0, v=None):
        mat_GF2E_construct(&self.x)

    def __dealloc__(self):
        mat_GF2E_destruct(&self.x)

    def __reduce__(self):
        """
        EXAMPLES:
            sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(5,5,range(25))
            sage: loads(dumps(m)) == m
            True
        """
        return unpickle_class_args, (ntl_mat_GF2E, (self.__nrows, self.__ncols, self.list()))

    def __repr__(self):
        _sig_on
        return string_delete(mat_GF2E_to_str(&self.x))

    def __mul__(ntl_mat_GF2E self, other):
        cdef ntl_mat_GF2E y
        cdef ntl_mat_GF2E r = ntl_mat_GF2E()
        if not isinstance(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other)
        y = other
        _sig_on
        mat_GF2E_mul(r.x, self.x, y.x)
        _sig_off
        return r

    def __sub__(ntl_mat_GF2E self, other):
        cdef ntl_mat_GF2E y
        cdef ntl_mat_GF2E r = ntl_mat_GF2E()
        if not isinstance(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other)
        y = other
        _sig_on
        mat_GF2E_sub(r.x, self.x, y.x)
        _sig_off
        return r

    def __add__(ntl_mat_GF2E self, other):
        cdef ntl_mat_GF2E y
        cdef ntl_mat_GF2E r = ntl_mat_GF2E()
        if not isinstance(other, ntl_mat_GF2E):
            other = ntl_mat_GF2E(other)
        y = other
        _sig_on
        mat_GF2E_add(r.x, self.x, y.x)
        _sig_off
        return r

    def __neg__(ntl_mat_GF2E self):
        cdef ntl_mat_GF2E r = ntl_mat_GF2E()
        _sig_on
        mat_GF2E_negate(r.x, self.x)
        _sig_off
        return r

    def __pow__(ntl_mat_GF2E self, long e, ignored):
        cdef ntl_mat_GF2E r = ntl_mat_GF2E()
        _sig_on
        mat_GF2E_power(r.x, self.x, e)
        _sig_off
        return r

    def __cmp__(self, other):
        c = cmp(type(self), type(other))
        if not c:
            if (self-other).is_zero():
                return 0
            else:
                return 1
        else:
            return c

    def nrows(self):
        """
        Number of rows
        """
        return self.__nrows

    def ncols(self):
        """
        Number of columns
        """
        return self.__ncols

    def __setitem__(self, ij, x):
        cdef ntl_GF2E y
        cdef int i, j
        if not isinstance(x, ntl_GF2E):
            x = ntl_GF2E(x)
        y = x

        from sage.rings.integer import Integer
        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.ncols()==1 and (isinstance(ij, Integer) or type(ij)==int):
            i, j = ij,0
        elif self.nrows()==1 and (isinstance(ij, Integer) or type(ij)==int):
            i, j = 0,ij
        else:
            raise TypeError, 'ij is not a matrix index'

        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError, "array index out of range"
        mat_GF2E_setitem(&self.x, i, j, &y.gf2e_x)

    def __getitem__(self, ij):
        cdef int i, j
        from sage.rings.integer import Integer
        if isinstance(ij, tuple) and len(ij) == 2:
            i, j = ij
        elif self.ncols()==1 and (isinstance(ij, Integer) or type(ij)==int):
            i, j = ij,0
        elif self.nrows()==1 and (isinstance(ij, Integer) or type(ij)==int):
            i, j = 0,ij
        else:
            raise TypeError, 'ij is not a matrix index'
        if i < 0 or i >= self.__nrows or j < 0 or j >= self.__ncols:
            raise IndexError, "array index out of range"
        cdef GF2E_c* elt = mat_GF2E_getitem(&self.x, i+1, j+1)
        cdef ntl_GF2E r = ntl_GF2E()
        r.gf2e_x = elt[0]
        GF2E_delete(elt)
        return r

    def determinant(self):
        """
        Returns the determinant.
        """
        cdef ntl_GF2E r = ntl_GF2E()
        _sig_on
        r.gf2e_x = mat_GF2E_determinant(self.x)
        _sig_off
        return r

    def echelon_form(self,ncols=0):
        """
        Performs unitary row operations so as to bring this matrix into row echelon
        form.  If the optional argument \code{ncols} is supplied, stops when first ncols
        columns are in echelon form.  The return value is the rank (or the
        rank of the first ncols columns).

        INPUT:
           ncols -- number of columns to process

        TODO: what is the output; does it return the rank?
        """
        cdef long w
        w = ncols
        _sig_on
        return mat_GF2E_gauss(&self.x, w)

    def list(self):
        """
        Returns a list of the entries in this matrix

        EXAMPLES:
            sage: ntl.GF2E_modulus([1,1,0,1,1,0,0,0,1])
            sage: m = ntl.mat_GF2E(2,2,[ntl.GF2E_random() for x in xrange(2*2)])
            sage: m.list()                   # random output
            [[1 0 1 0 1 0 1 1], [0 1 0 0 0 0 1], [0 0 0 1 0 0 1], [1 1 1 0 0 0 0 1]]
        """
        return [self[i,j] for i in range(self.__nrows) for j in range(self.__ncols)]

    def is_zero(self):
        cdef long isZero
        _sig_on
        isZero = mat_GF2E_is_zero(&self.x)
        _sig_off
        if isZero==0:
            return False
        else:
            return True

    def _sage_(ntl_mat_GF2E self, k=None, cache=None):
        """
        Returns a \class{Matrix} over a FiniteField representation
        of this element.

        INPUT:
            self  -- \class{mat_GF2E} element
            k     -- optional GF(2**deg)
            cache -- optional NTL to SAGE conversion dictionary

        OUTPUT:
            Matrix over k
        """
        if k==None:
            k = ntl_GF2E_sage()
        l = []
        for elem in self.list():
            l.append(elem._sage_(k, cache))

        from sage.matrix.constructor import Matrix
        return Matrix(k,self.nrows(),self.ncols(),l)

    def transpose(ntl_mat_GF2E self):
        """
        Returns the transposed matrix of self.

        OUTPUT:
            transposed Matrix
        """
        cdef ntl_mat_GF2E r = ntl_mat_GF2E()
        _sig_on
        mat_GF2E_transpose(r.x, self.x)
        _sig_off
        return r
