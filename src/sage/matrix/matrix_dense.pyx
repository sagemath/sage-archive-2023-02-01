r"""
Base class for dense matrices

TESTS:
    sage: R.<a,b> = QQ[]
    sage: m = matrix(R,2,[0,a,b,b^2])
    sage: loads(dumps(m)) == m
    True
"""

cimport matrix

from   sage.structure.element    cimport Element
import sage.matrix.matrix_space
import sage.structure.sequence

include '../ext/cdefs.pxi'
include '../ext/stdsage.pxi'

cdef class Matrix_dense(matrix.Matrix):
    cdef bint is_sparse_c(self):
        return 0

    cdef bint is_dense_c(self):
        return 1

    def __copy__(self):
        """
        Return a copy of this matrix.  Changing the entries of the
        copy will not change the entries of this matrix.
        """
        A = self.new_matrix(entries=self.list(), coerce=False, copy=False)
        if self.subdivisions is not None:
            A.subdivide(*self.get_subdivisions())
        return A

    def __hash__(self):
        """
        Return the hash of this matrix.

        Equal matrices should have equal hashes, even if one is sparse and
        the other is dense.

        EXAMPLES:
            sage: m = matrix(2, range(24), sparse=True)
            sage: m.set_immutable()
            sage: hash(m)
            976

            sage: d = m.dense_matrix()
            sage: d.set_immutable()
            sage: hash(d)
            976

            sage: hash(m) == hash(d)
            True
        """
        return self._hash()

    cdef long _hash(self) except -1:
        x = self.fetch('hash')
        if not x is None: return x

        if not self._mutability._is_immutable:
            raise TypeError, "mutable matrices are unhashable"

        v = self._list()
        cdef Py_ssize_t i
        cdef long h
        h = 0
        n = 1
        for i from 0 <= i < len(v):
            h = h ^ (i * hash(v[i]))

        self.cache('hash', h)
        return h

    def _multiply_classical(left, Matrix_dense right):
        """
        Multiply the matrices left and right using the classical $O(n^3)$
        algorithm.

        This method assumes that left and right have the same parent and
        compatable dimensions.
        """
        cdef Py_ssize_t i, j, k, l
        if left._ncols != right._nrows:
            raise IndexError, "Number of columns of left must equal number of rows of other."


        v = PyList_New(left._nrows * right._ncols)     # this is really sort of v = []..."
        zero = left.base_ring()(0)
        l = 0
        for i from 0 <= i < left._nrows:
            for j from 0 <= j < right._ncols:
                s = zero
                for k from 0 <= k < left._ncols:
                    s = s + left.get_unsafe(i,k) * right.get_unsafe(k,j)
                # This is really v.append(s)
                Py_INCREF(s); PyList_SET_ITEM(v, l, s)
                l = l + 1
        return left.new_matrix(left._nrows, right._ncols, entries = v, coerce=False, copy=False)




    def _pickle(self):
        version = -1
        data = self._list()  # linear list of all elements
        return data, version

    def _unpickle_generic(self, data, int version):
        cdef Py_ssize_t i, j, k
        if version == -1:
            # data is a *list* of the entries of the matrix.
            # TODO: Change the data[k] below to use the fast list access macros from the Python/C API
            k = 0
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    self.set_unsafe(i, j, data[k])
                    k = k + 1
        else:
            raise RuntimeError, "unknown matrix version (=%s)"%version

    cdef int _cmp_c_impl(self, Element right) except -2:
        """
        EXAMPLES:
            sage: P.<x> = QQ[]
            sage: m = matrix([[x,x+1],[1,x]])
            sage: n = matrix([[x+1,x],[1,x]])
            sage: o = matrix([[x,x],[1,x]])
            sage: m.__cmp__(n)
            -1
            sage: m.__cmp__(m)
            0
            sage: n.__cmp__(m)
            1
            sage: m.__cmp__(o)
            1
        """
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                res = cmp( self[i,j], right[i,j] )
                if res != 0:
                    return res
        return 0

    def transpose(self):
        """
        Returns the transpose of self, without changing self.

        EXAMPLES:
        We create a matrix, compute its transpose, and note that the
        original matrix is not changed.
            sage: M = MatrixSpace(QQ,  2)
            sage: A = M([1,2,3,4])
            sage: B = A.transpose()
            sage: print B
            [1 3]
            [2 4]
            sage: print A
            [1 2]
            [3 4]

            sage: A.subdivide(None, 1); A
            [1|2]
            [3|4]
            sage: A.transpose()
            [1 3]
            [---]
            [2 4]
        """
        f = []
        e = self.list()
        (nc, nr) = (self.ncols(), self.nrows())
        for j in xrange(nc):
            for i in xrange(nr):
                f.append(e[i*nc + j])
        trans = self.new_matrix(nrows = nc, ncols = nr,
                                entries = f, copy=False,
                                coerce=False)
        if self.subdivisions is not None:
            row_divs, col_divs = self.get_subdivisions()
            trans.subdivide(col_divs, row_divs)
        return trans


    def antitranspose(self):
        """
        Returns the anittranspose of self, without changing self.

        EXAMPLES:
            sage: A = matrix(2,3,range(6)); A
            [0 1 2]
            [3 4 5]
            sage: A.antitranspose()
            [5 2]
            [4 1]
            [3 0]

            sage: A.subdivide(1,2); A
            [0 1|2]
            [---+-]
            [3 4|5]
            sage: A.antitranspose()
            [5|2]
            [-+-]
            [4|1]
            [3|0]
        """
        f = []
        e = self.list()
        (nc, nr) = (self.ncols(), self.nrows())
        for j in reversed(xrange(nc)):
            for i in reversed(xrange(nr)):
                f.append(e[i*nc + j])
        trans = self.new_matrix(nrows = nc, ncols = nr,
                                entries = f, copy=False, coerce=False)
        if self.subdivisions is not None:
            row_divs, col_divs = self.get_subdivisions()
            trans.subdivide(list(reversed([nc - t for t in col_divs])),
                            list(reversed([nr - t for t in row_divs])))
        return trans

    def apply_morphism(self, phi):
        """
        Apply the morphism phi to the coefficients of this dense matrix.

        The resulting matrix is over the codomain of phi.

        INPUT:
            phi -- a morphism, so phi is callable and phi.domain()
                   and phi.codomain() are defined.  The codomain
                   must be a ring.

        OUTPUT:
            a matrix over the codomain of phi


        EXAMPLES:
            sage: m = matrix(ZZ, 3, range(9))
            sage: phi = ZZ.hom(GF(5))
            sage: m.apply_morphism(phi)
            [0 1 2]
            [3 4 0]
            [1 2 3]
            sage: parent(m.apply_morphism(phi))
            Full MatrixSpace of 3 by 3 dense matrices over Finite Field of size 5

        We apply a morphism to a matrix over a polynomial ring:
            sage: R.<x,y> = QQ[]
            sage: m = matrix(2, [x,x^2 + y, 2/3*y^2-x, x]); m
            [          x     x^2 + y]
            [2/3*y^2 - x           x]
            sage: phi = R.hom([y,x])
            sage: m.apply_morphism(phi)
            [          y     y^2 + x]
            [2/3*x^2 - y           y]

        """
        R = phi.codomain()
        M = sage.matrix.matrix_space.MatrixSpace(R, self._nrows,
                   self._ncols, sparse=False)
        image = M([phi(z) for z in self.list()])
        if self.subdivisions is not None:
            image.subdivide(*self.get_subdivisions())
        return image

    def apply_map(self, phi, R=None):
        """
        Apply the given map phi (an arbitrary Python function or
        callable object) to this dense matrix.  If R is not given,
        automatically determine the base ring of the resulting matrix.

        INPUT:
            phi -- arbitrary Python function or callable object
            R -- (optional) ring

        OUTPUT:
            a matrix over R

        EXAMPLES:
            sage: m = matrix(ZZ, 3, range(9))
            sage: k.<a> = GF(9)
            sage: f = lambda x: k(x)
            sage: n = m.apply_map(f); n
            [0 1 2]
            [0 1 2]
            [0 1 2]
            sage: n.parent()
            Full MatrixSpace of 3 by 3 dense matrices over Finite Field in a of size 3^2

        In this example, we explicitly specify
        the codomain.
            sage: s = GF(3)
            sage: f = lambda x: s(x)
            sage: n = m.apply_map(f, k); n
            [0 1 2]
            [0 1 2]
            [0 1 2]
            sage: n.parent()
            Full MatrixSpace of 3 by 3 dense matrices over Finite Field in a of size 3^2
        """
        v = [phi(z) for z in self.list()]
        if R is None:
            v = sage.structure.sequence.Sequence(v)
            R = v.universe()
        M = sage.matrix.matrix_space.MatrixSpace(R, self._nrows,
                   self._ncols, sparse=False)
        image = M(v)
        if self.subdivisions is not None:
            image.subdivide(*self.get_subdivisions())
        return image
