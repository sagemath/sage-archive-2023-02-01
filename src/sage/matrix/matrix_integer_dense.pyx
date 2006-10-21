"""
Dense matrices over the integers.
"""

######################################################################
#       Copyright (C) 2006 William Stein
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
######################################################################

from sage.misc.misc import verbose, get_verbose

include "../ext/gmp.pxi"

cimport sage.rings.integer
import sage.rings.integer

cimport matrix_integer
cimport matrix_generic

cdef class Matrix_integer_dense(matrix_integer.Matrix_integer):
    r"""
    Matrix over the integers.

    On a 32-bit machine, they can have at most $2^{32}-1$ rows or
    columns.  On a 64-bit machine, matrices can have at most
    $2^{64}-1$ rows or columns.

    EXAMPLES:
        sage: a = MatrixSpace(ZZ,3)(2); a
        [2 0 0]
        [0 2 0]
        [0 0 2]
        sage: a = matrix(ZZ,1,3, [1,2,-3]); a
        [ 1  2 -3]
        sage: a = MatrixSpace(ZZ,2,4)(2); a
        Traceback (most recent call last):
        ...
        TypeError: nonzero scalar matrix must be square
    """
    def __new__(self, parent, entries, coerce, copy):
        """
        INPUT:
            parent -- a matrix space
            entries -- list - create the matrix with those entries along the rows.
                       other -- a scalar; entries is coerced to an integer and the diagonal
                                entries of this matrix are set to that integer.
            coerce -- whether to coerce entries to the integers
            copy -- ignored (since integers are immutable)

        EXAMPLES:
        """
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        self._pivots = None
        matrix_generic.Matrix.__init__(self, parent)

        # Allocate an array where all the entries of the matrix are stored.
        self._entries = <mpz_t *>sage_malloc(sizeof(mpz_t) * (self._nrows * self._ncols))
        if self._entries == NULL:
            raise MemoryError, "out of memory allocating a matrix"

        # Allocate an array of pointers to the rows.
        self._matrix = <mpz_t **> sage_malloc(sizeof(mpz_t*)*self._nrows)
        if self._matrix == NULL:
            sage_free(self._entries)
            raise MemoryError, "out of memory allocating a matrix"

        # Set each of the pointers in the array self._matrix to point at the memory for
        # the corresponding row.
        cdef size_t i, k
        k = 0
        for i from 0 <= i < self._nrows:
            self._matrix[i] = self._entries + k
            k = k + self._ncols

    def __init__(self, parent, entries, coerce, copy):
        """
        INPUT:
            parent -- a matrix space
            entries -- list - create the matrix with those entries along the rows.
                       other -- a scalar; entries is coerced to an integer and the diagonal
                                entries of this matrix are set to that integer.
            coerce -- whether need to coerce entries to the integers (program may crash
                      if you get this wrong)
            copy -- ignored (since integers are immutable)

        EXAMPLES:

        """
        cdef size_t i, j
        cdef sage.rings.integer.Integer x

        if isinstance(entries, list):  # todo -- change to PyObject_TypeCheck???
            # Create the matrix whose entries are in the given entry list.
            if len(entries) != self._nrows * self._ncols:
                raise TypeError, "entries has the wrong length"
            if coerce:
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    x = sage.rings.integer.Integer(entries[i])   # todo -- see integer.pyx and the TODO there; perhaps this could be
                                     # sped up by creating a mpz_init_set_sage function.
                    mpz_init_set(self._entries[i], x.value)
            else:
                for i from 0 <= i < self._nrows * self._ncols:
                    # TODO: Should use an unsafe un-bounds-checked array access here.
                    x = entries[i]
                    mpz_init_set(self._entries[i], x.value)

        else:

            # Try to coerce entries to a scalar (an integer)
            x = sage.rings.integer.Integer(entries)

            # If x is zero, make the zero matrix and be done.
            if mpz_cmp_si(x.value, 0) == 0:
                self._zero_out_matrix()
                return

            # the matrix must be square:
            if self._nrows != self._ncols:
                sage_free(self._entries)
                sage_free(self._matrix)
                self._entries = NULL
                raise TypeError, "nonzero scalar matrix must be square"

            # Now we set all the diagonal entries to x and all other entries to 0.
            self._zero_out_matrix()
            j = 0
            for i from 0 <= i < self._nrows:
                mpz_init_set(self._entries[j], x.value)
                j = j + self._nrows + 1

    cdef void _zero_out_matrix(self):
        """
        Set this matrix to be the zero matrix.
        This is only for internal use.
        """
        # TODO: This is about 6-10 slower than MAGMA doing what seems to be the same thing.
        # Moreover, NTL can also do this quickly.  Why?
        _sig_on
        cdef size_t i
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init(self._entries[i])
        _sig_off


    def nrows(self):
        """
        Return the number of rows of this matrix.

        EXAMPLES:
            sage: a = MatrixSpace(ZZ,2,4)(range(8,16)); a
            [ 8  9 10 11]
            [12 13 14 15]
            sage: a.nrows()
            2
        """

        return sage.rings.integer.Integer(self._nrows)

    def ncols(self):
        """
        Return the number of columns of this matrix.

        EXAMPLES:
            sage: a = MatrixSpace(ZZ,2,7)(range(1,15)); a
            [ 1  2  3  4  5  6  7]
            [ 8  9 10 11 12 13 14]
            sage: a.nrows()
            2
            sage: a.ncols()
            7
        """
        return sage.rings.integer.Integer(self._ncols)

    def __reduce__(self):
        """
        Save this matrix to a binary stream.

        EXAMPLES:
            sage: a = matrix(ZZ, 2,2, [1,2,3,4])
            sage: loads(dumps(a)) == a
            True
        """
        # TODO: redo this to use the "mpz to bytes" stuff.
        import sage.matrix.reduce

        cdef int i, j, len_so_far, m, n
        cdef char *a
        cdef char *s, *t, *tmp

        if self._nrows == 0 or self._ncols == 0:
            entries = ''
        else:
            n = self._nrows*self._ncols*10
            s = <char*> sage_malloc(n * sizeof(char))
            t = s
            len_so_far = 0

            _sig_on
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    m = mpz_sizeinbase (self._matrix[i][j], 32)
                    if len_so_far + m + 1 >= n:
                        # copy to new string with double the size
                        n = 2*n + m + 1
                        tmp = <char*> sage_malloc(n * sizeof(char))
                        strcpy(tmp, s)
                        sage_free(s)
                        s = tmp
                        t = s + len_so_far
                    #endif
                    mpz_get_str(t, 32, self._matrix[i][j])
                    m = strlen(t)
                    len_so_far = len_so_far + m + 1
                    t = t + m
                    t[0] = <char>32
                    t[1] = <char>0
                    t = t + 1
            _sig_off
            entries = str(s)[:-1]
            free(s)

        return sage.matrix.reduce.make_Matrix_integer_dense, \
               (self.parent(), entries)


    def __cmp__(self, other):
        """
        Compare self to other.

        TODO -- this is dumb -- compares wrong if other isn't also a Matrix_integer_dense.

        EXAMPLES:
        """
        cdef size_t i
        cdef Matrix_integer_dense c_other
        # TODO: make this _cmp_siblings_ instead.  Definitely don't want to return -1 here -- ...
        if not isinstance(other, Matrix_integer_dense):
            return -1
        c_other = other
        if self._nrows != c_other._rows or self._ncols != c_other._ncols:
            return -1
        for i from 0 <= i < self._nrows * self._ncols:
            if mpz_cmp(self._entries[i], c_other._entries[i]):
                return -1
        return 0

    def __setitem__(self, ij, x):
        """
        Set position i,j of this matrix to x.

        EXAMPLES:
            sage: a = matrix(ZZ,2,3, range(6)); a
            [0 1 2]
            [3 4 5]

            sage: a[0,0] = 10
            ...
        """
        cdef sage.rings.integer.Integer z
        cdef size_t i, j

        if not PyTuple_Check(ij):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError, "index must be an integer or pair of integers"
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
                raise IndexError, "index out of bounds"
            # TODO: Use direct call to sage_to_mpz_t -- no need to make an Integer object.
            try:
                z = x
            except TypeError:
                z = sage.rings.integer.Integer(x)
            mpz_set(self._matrix[i][j], z.value)
        else:
            # If ij is not a tuple, coerce to an integer and set the row.
            i = ij
            for j from 0 <= j < self._ncols:
                # TODO: Use direct call to sage_to_mpz_t -- no need to make an Integer object.
                try:
                    z = x[j]
                except TypeError:
                    z = sage.rings.integer.Integer(x[j])
                mpz_set(self._matrix[i][j], z.value)

    def __getitem__(self, ij):
        """
        EXAMPLES:
            sage: a = MatrixSpace(ZZ,3)(range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a[1,2]
            5
            sage: a[0]
            (0, 1, 2)
            sage: a[4,7]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
            sage: a[-1,0]
            Traceback (most recent call last):
            ...
            IndexError: index out of bounds
        """
        cdef sage.rings.integer.Integer z
        cdef size_t i, j
        cdef object x

        if PyTuple_Check(ij):
            # ij is a tuple, so we get i and j efficiently, construct corresponding integer entry.
            if PyTuple_Size(ij) != 2:
                raise IndexError, "index must be an integer or pair of integers"
            i = <object> PyTuple_GET_ITEM(ij, 0)
            j = <object> PyTuple_GET_ITEM(ij, 1)
            if i < 0 or i >= self._nrows or j < 0 or j >= self._ncols:
                raise IndexError, "index out of bounds"
            z = sage.rings.integer.Integer()
            mpz_init_set(z.value, self._matrix[i][j])
            return z
        else:
            # If ij is not a tuple, coerce to an integer and get the row.
            i = ij
            return self.row(i)


    def  __dealloc__(self):
        """
        EXAMPLE:
        """
        if self._entries == NULL: return
        cdef size_t i
        for i from 0 <= i < (self._nrows * self._ncols):
            mpz_clear(self._entries[i])
        sage_free(self._entries)
        sage_free(self._matrix)


    def _mul_(Matrix_integer_dense self, Matrix_integer_dense other):
        """
        EXAMPLE:
            sage: n = 3
            sage: a = MatrixSpace(ZZ,n,n)(range(n^2))
            sage: a*a - a + a
            ...?
        """
        if self._ncols != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        cdef int i, j, k, nr, nc, snc
        cdef mpz_t *v
        cdef mpz_t s, z
        nr = self._nrows
        nc = other._ncols
        snc = self._ncols

        cdef Matrix_integer_dense M
        M = self.new_matrix(nr, nc, zero=False)

        cdef mpz_t **m
        m = M._matrix

        mpz_init(s); mpz_init(z)

        _sig_on
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                mpz_set_si(s,0)   # set s = 0
                v = self._matrix[i]
                for k from 0 <= k < snc:
                    mpz_mul(z, v[k], other._matrix[k][j])
                    mpz_add(s, s, z)
                mpz_init(m[i][j])
                mpz_set(m[i][j], s)
        _sig_off
        mpz_clear(s); mpz_clear(z)
        return M

    def _add_(Matrix_integer_dense self, Matrix_integer_dense other):
        if self._ncols != other._ncols:
            raise IndexError, "Number of columns of self must equal number of rows of other."
        if self._nrows != other._nrows:
            raise IndexError, "Number of columns of self must equal number of rows of other."

        cdef int i, j

        cdef Matrix_integer_dense M
        M = self.new_matrix(zero=False)

        cdef mpz_t **m
        m = M._matrix

        _sig_on
        for i from 0 <= i < self._ncols * self._nrows:
            mpz_init(m[i])
            mpz_add(m[i], self._entries[i], other._entries[i])
        _sig_off
        return M

    def transpose(self):
        """
        Returns the transpose of self.

        EXAMPLES:
            sage: a = matrix(ZZ,2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.transpose()
            [0 3]
            [1 4]
            [2 5]
            sage: a
            [0 1 2]
            [3 4 5]
            sage: a.transpose().transpose()
            [0 1 2]
            [3 4 5]
        """
        cdef size_t i, j
        cdef Matrix_integer_dense M

        P = self._parent.matrix_space(self._ncols, self._nrows)
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)
        _sig_on
        for i from 0 <= i < self._ncols:
            for j from 0 <= j < self._nrows:
                mpz_init_set(M._matrix[i][j], self._matrix[j][i])
        _sig_off
        return M

    def matrix_from_rows(self, rows):
        """
        Return the submatrix formed from the given rows.

        INPUT:
            rows -- list of int's
        OUTPUT:
            matrix created from the rows with given indexes

        EXAMPLES:

        """
        cdef int i, j, k, nc, nr
        cdef Matrix_integer_dense M

        if not isinstance(rows, list):
            raise TypeError, "rows (=%s) must be a list"%rows
        nr = len(rows)
        if nr == 0:
            return new_matrix(0, self._ncols)
        nc = self._ncols
        v = []
        for i in rows:
            v.append(int(i))
        rows = v
        if min(rows) < 0 or max(rows) >= self._nrows:
            raise IndexError, "invalid row indexes; rows don't exist"

        M = new_matrix(zero=False)
        cdef mpz_t **m
        m = M._matrix

        for i from 0 <= i < nr:
            k = rows[i]
            for j from 0 <= j < nc:
                mpz_init(m[i][j])
                mpz_set(m[i][j], self._matrix[k][j])

        return M



    def iterates(self, v, int n):
        """
        Let A be this matrix.   Return a matrix with rows
        $$
           v, Av, A^2v, ..., A^(n-1)v.
        $$
        """
        cdef int i, j, k, nr, nc
        cdef mpz_t s, z
        nr = n
        nc = self._ncols

        if self._nrows != self._ncols:
            raise ArithmeticError, "matrix must be square"
        if not isinstance(v, list):
            raise TypeError, "v must be a list"
        if len(v) != self._nrows:
            raise ArithmeticError, "incompatible matrix vector multiple"

        cdef Matrix_integer_dense M
        M = new_matrix(zero=False)

        cdef mpz_t **m
        m = M._matrix
        mpz_init(self.tmp)
        for j from 0 <= j < nc:
            string = str(v[j])
            r = mpz_set_str(self.tmp, string, 0)
            if r == -1:
                raise TypeError, "Invalid integer"
            mpz_init(m[0][j])
            mpz_set(m[0][j], self.tmp)

        mpz_init(s)
        mpz_init(z)
        for i from 1 <= i < nr:
            m[i] = <mpz_t *> sage_malloc(sizeof(mpz_t)*nc)
            if m[i] == <mpz_t*> 0:
                raise MemoryError, "Error allocating matrix"
            for j from 0 <= j < nc:
                mpz_set_si(s,0)  # set s = 0
                for k from 0 <= k < self._nrows:
                    mpz_mul(z, m[i-1][k], self._matrix[k][j])
                    mpz_add(s, s, z)
                mpz_init(m[i][j])
                mpz_set(m[i][j], s)
        mpz_clear(s); mpz_clear(z)
        return M


    def scalar_multiple(self, d):
        """
        Return the product self*d, as a new matrix.
        """
        cdef int i, j, nr, nc
        cdef sage.rings.integer.Integer z
        nr = self._nrows
        nc = self._ncols

        cdef mpz_t x
        try:
            z = d
        except (TypeError):
            z = sage.rings.integer.Integer(d)
        mpz_init(x)
        mpz_set(x, z.value)

        cdef Matrix_integer_dense M
        M = new_matrix(zero=False)

        cdef mpz_t *e
        e = M._entries

        _sig_on
        for i from 0 <= i < nr * nc:
            mpz_init(e[i])
            mpz_mul(e[i], self._entries[i], x)
        _sig_off
        mpz_clear(x)
        return M

    def copy(self):
        """
        EXAMPLE:
        """
        cdef int i, j, nr, nc
        nr = self._nrows; nc = self._ncols

        cdef Matrix_integer_dense M
        M = new_matrix(zero=False)

        for i from 0 <= i < nr * nc:
            mpz_set(M._entries[i], self._entries[i])

        for i from 0 <= i < nc:
            M._matrix[i] = self._matrix[i]

        return M

    def number_nonzero(self):
        """
        EXAMPLE:
        """
        cdef int i, j, n
        cdef mpz_t *v
        n = 0
        _sig_on
        for i from 0 <= i < self._nrows * self._ncols:
            if mpz_sgn(self._entries[i]):         # if nonzero
                n = n + 1
        _sig_off
        return n

    def list(self, int base=0):
        """
        EXAMPLE:
        """
        cdef int i, j
        cdef mpz_t* r
        cdef object v
        cdef sage.rings.integer.Integer x

        v = []
        _sig_on
        for i from 0 <= i < self._nrows:
            r = self._matrix[i]
            for j from 0 <= j < self._ncols:
                x = sage.rings.integer.Integer()
                x.set_from_mpz(r[j])
                v.append(x)
        _sig_off
        return v

    def rank(self):
        """
        Return the rank found during the last echelon operation on self.
        Of course if self is changed, and the echelon form of self is not
        recomputed, then the rank could be incorrect.
        """
        if self._pivots == None:
            raise ArithmeticError, "Echelon form has not yet been computed."
        return len(self._pivots)

    def pivots(self):
        """
        Return the pivots found during the last echelon operation on self.
        Of course if self is changed, and the echelon form of self is not
        recomputed, then the pivots could be incorrect.

        EXAMPLES:
        """
        if self._pivots == None:
            raise ArithmeticError, "Echelon form has not yet been computed."
        return self._pivots

    cdef int mpz_height(self, mpz_t height) except -1:
        """
        EXAMPLES:
        """
        cdef mpz_t x, h
        mpz_init(x)
        mpz_init_set_si(h, 0)
        cdef int i, j
        _sig_on
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                mpq_get_num(x,self._matrix[i][j])
                mpz_abs(x, x)
                if mpz_cmp(h,x) < 0:
                    mpz_set(h,x)
                mpq_get_den(x,self._matrix[i][j])
                mpz_abs(x, x)
                if mpz_cmp(h,x) < 0:
                    mpz_set(h,x)
        _sig_off
        mpz_set(height, h)
        mpz_clear(h)
        mpz_clear(x)
        return 0

    def height(self):
        """
        EXAMPLE:
        """
        cdef mpz_t h
        mpz_init(h)
        self.mpz_height(h)
        a = mpz_to_long(h)
        mpz_clear(h)
        return a


    def multiply_multi_modular(self, Matrix_integer_dense right):
        """
        Multiply this matrix by right using a multimodular algorithm
        and return the result.

        EXAMPLES:
        """
        raise NotImplementedError

cdef object mpz_to_long(mpz_t x):
    return long(mpz_to_str(x))

