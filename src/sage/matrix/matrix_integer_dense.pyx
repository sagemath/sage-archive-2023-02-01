"""
Dense matrices over the integers

TODO:
    * Implement a sensible __cmp__ method: this will be part of going through all the code in SAGE and
      creating a systematic __cmp__.

    * I implemented a fast matrix_from_rows, but not a fast matrix_from_columns.
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

include "../ext/stdsage.pxi"

cimport sage.rings.integer
import  sage.rings.integer

cimport matrix_generic
cimport matrix_pid_dense

def make_Matrix_integer_dense(parent, entries):
    cdef Matrix_integer_dense M
    # Create a new blank matrix with memory allocated but not yet initialized
    M = Matrix_integer_dense.__new__(Matrix_integer_dense, parent, 0,0,0)
    entries = entries.split()

    cdef size_t i, n
    n = M._nrows * M._ncols
    if len(entries) != n:
        raise RuntimeError, "invalid pickle data."
    for i from 0 <= i < n:
        s = entries[i]
        if mpz_init_set_str(M._entries[i], s, 32):
            raise RuntimeError, "invalid pickle data"

    return M

cdef class Matrix_integer_dense(matrix_pid_dense.Matrix_pid_dense):
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
        Create and allocate memory for the matrix.  Does not actually initialize
        any of the memory (e.g., mpz_init is not called).

        INPUT:
            parent, entries, coerce, copy -- as for __init__.

        EXAMPLES:
            sage: from sage.matrix.matrix_integer_dense import Matrix_integer_dense
            sage: a = Matrix_integer_dense.__new__(Matrix_integer_dense, Mat(ZZ,3), 0,0,0)
            sage: type(a)
            <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>

        WARNING: This is for internal use only, or if you really know what you're doing.
        """
        self._initialized = 0
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        self._pivots = None
        matrix.Matrix.__init__(self, parent)

        # Allocate an array where all the entries of the matrix are stored.
        self._entries = <mpz_t *>sage_malloc(sizeof(mpz_t) * (self._nrows * self._ncols))
        if self._entries == NULL:
            raise MemoryError, "out of memory allocating a matrix"

        # Allocate an array of pointers to the rows, which is useful for
        # certain algorithms.
        self._matrix = <mpz_t **> sage_malloc(sizeof(mpz_t*)*self._nrows)
        if self._matrix == NULL:
            sage_free(self._entries)
            self._entries = NULL
            raise MemoryError, "out of memory allocating a matrix"

        # Set each of the pointers in the array self._matrix to point
        # at the memory for the corresponding row.
        cdef size_t i, k
        k = 0
        for i from 0 <= i < self._nrows:
            self._matrix[i] = self._entries + k
            k = k + self._ncols

    def __init__(self, parent, entries, coerce, copy):
        r"""
        Initialize a dense matrix over the integers.

        INPUT:
            parent -- a matrix space
            entries -- list - create the matrix with those entries along the rows.
                       other -- a scalar; entries is coerced to an integer and the diagonal
                                entries of this matrix are set to that integer.
            coerce -- whether need to coerce entries to the integers (program may crash
                      if you get this wrong)
            copy -- ignored (since integers are immutable)

        EXAMPLES:

        The __init__ function is called implicitly in each of the
        examples below to actually fill in the values of the matrix.

        We create a $2 \times 2$ and a $1\times 4$ matrix:
            sage: matrix(ZZ,2,2,range(4))
            [0 1]
            [2 3]
            sage: Matrix(ZZ,1,4,range(4))
            [0 1 2 3]

        If the number of columns isn't given, it is determined from the number of
        elements in the list.
            sage: matrix(ZZ,2,range(4))
            [0 1]
            [2 3]
            sage: matrix(ZZ,2,range(6))
            [0 1 2]
            [3 4 5]

        Another way to make a matrix is to create the space of
        matrices and coerce lists into it.
            sage: A = Mat(ZZ,2); A
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: A(range(4))
            [0 1]
            [2 3]

        Actually it is only necessary that the input can be coerced to a list, so
        the following also works:
            sage: v = reversed(range(4)); type(v)
            <type 'listreverseiterator'>
            sage: A(v)
            [3 2]
            [1 0]

        Matrices can have many rows or columns (in fact, on a 64-bit machine they could
        have up to $2^64-1$ rows or columns):
            sage: v = Matrix(ZZ,1,10^5, 0)
            sage: v.parent()
            Full MatrixSpace of 1 by 100000 dense matrices over Integer Ring
        """
        cdef size_t i, j
        cdef int is_list
        cdef sage.rings.integer.Integer x

        if not isinstance(entries, list):  # todo -- change to PyObject_TypeCheck???
            try:
                entries = list(entries)
                is_list = 1
            except TypeError:
                try:
                    # Try to coerce entries to a scalar (an integer)
                    x = sage.rings.integer.Integer(entries)
                    is_list = 0
                except TypeError:
                    raise TypeError, "entries must be coercible to a list or integer"
        else:
            is_list = 1

        if is_list:

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
            self._initialized = 1

    cdef void _zero_out_matrix(self):
        """
        Set this matrix to be the zero matrix.
        This is only for internal use.
        """
        # TODO: This is about 6-10 slower than MAGMA doing what seems to be the same thing.
        # Moreover, NTL can also do this quickly.  Why?   I think both have specialized
        # small integer classes.
        _sig_on
        cdef size_t i
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init(self._entries[i])
        _sig_off
        self._initialized = 1

    cdef _new_unitialized_matrix(self, size_t nrows, size_t ncols):
        """
        Return a new matrix over the integers with the given number of rows and columns.
        All memory is allocated for this matrix, but its entries have not yet been
        filled in.
        """
        P = self._parent.matrix_space(nrows, ncols)
        return Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)

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
            sage: a = matrix(ZZ, 2,2, [-1,2,3,199])
            sage: loads(dumps(a)) == a
            True
        """
        # TODO: redo this to use mpz_import and mpz_export
        # from sec 5.14 of the GMP manual.
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
            for i from 0 <= i < self._nrows * self._ncols:
                m = mpz_sizeinbase (self._entries[i], 32)
                if len_so_far + m + 1 >= n:
                    # copy to new string with double the size
                    n = 2*n + m + 1
                    tmp = <char*> sage_malloc(n * sizeof(char))
                    strcpy(tmp, s)
                    sage_free(s)
                    s = tmp
                    t = s + len_so_far
                #endif
                mpz_get_str(t, 32, self._entries[i])
                m = strlen(t)
                len_so_far = len_so_far + m + 1
                t = t + m
                t[0] = <char>32
                t[1] = <char>0
                t = t + 1
            _sig_off
            entries = str(s)[:-1]
            free(s)

        return make_Matrix_integer_dense, (self.parent(), entries)

    def __cmp__(self, other):
        """
        Compare self to other.

        TODO -- this is dumb -- compares wrong if other isn't also a Matrix_integer_dense.

        EXAMPLES:
            sage: a = Mat(ZZ,2)([-1,1,10,3])
            sage: b = 2*a
            sage: a == b
            False
            sage: a < b
            False
            sage: b < a
            True
            sage: a + a == b
            True
        """
        # TODO: make this _cmp_siblings_ instead.  Definitely don't want to return -1 here -- ...

        cdef size_t i
        cdef int c
        cdef Matrix_integer_dense c_other

        try:
            c_other = other
        except TypeError:
            return -1
        if self._nrows != c_other._nrows or self._ncols != c_other._ncols:
            return -1
        for i from 0 <= i < self._nrows * self._ncols:
            c = mpz_cmp(self._entries[i], c_other._entries[i])
            if c < 0: return -1
            elif c > 0: return 1
        return 0

    def __setitem__(self, ij, x):
        """
        Set position i,j of this matrix to x.

        INPUT:
            ij -- tuple (i,j), where i is the row and j the column
        Alternatively, ij can be an integer, and the ij-th row is set.

        EXAMPLES:
            sage: a = matrix(ZZ,2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a[0,0] = 10
            sage: a
            [10  1  2]
            [ 3  4  5]
        """
        cdef sage.rings.integer.Integer z
        cdef size_t i, j

        if PyTuple_Check(ij):
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
            z = sage.rings.integer.Integer()  # todo: change to use new integer
            mpz_init_set(z.value, self._matrix[i][j])
            return z
        else:
            # If ij is not a tuple, coerce to an integer and get the row.
            i = ij
            return self.row(i)


    def  __dealloc__(self):
        """
        Frees all the memory allocated for this matrix.
        EXAMPLE:
            sage: a = Matrix(ZZ,2,[1,2,3,4])
            sage: del a
        """
        if self._entries == NULL: return
        cdef size_t i
        if self._initialized:
            for i from 0 <= i < (self._nrows * self._ncols):
                mpz_clear(self._entries[i])
        sage_free(self._entries)
        sage_free(self._matrix)


    def _mul_(Matrix_integer_dense self, Matrix_integer_dense other):
        """
        EXAMPLE:
            sage: n = 3
            sage: a = MatrixSpace(ZZ,n,n)(range(n^2))
            sage: a*a
            [ 15  18  21]
            [ 42  54  66]
            [ 69  90 111]
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
                mpz_init_set(m[i][j], s)
        _sig_off
        mpz_clear(s); mpz_clear(z)
        return M

    def _add_(Matrix_integer_dense self, Matrix_integer_dense other):
        """
        Add two dense matrices over ZZ.

        EXAMPLES:
            sage: a = MatrixSpace(ZZ,3)(range(9))
            sage: a+a
            [ 0  2  4]
            [ 6  8 10]
            [12 14 16]
        """
        if self._ncols != other._ncols:
            raise IndexError, "Number of columns of self must equal number of columns of other."
        if self._nrows != other._nrows:
            raise IndexError, "Number of rows of self must equal number of rows of other."

        cdef size_t i, j

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)

        _sig_on

        cdef mpz_t *entries
        entries = M._entries
        for i from 0 <= i < self._ncols * self._nrows:
            mpz_init(entries[i])
            mpz_add(entries[i], self._entries[i], other._entries[i])

        _sig_off
        return M

    def _sub_(Matrix_integer_dense self, Matrix_integer_dense other):
        """
        Subtract two dense matrices over ZZ.

        EXAMPLES:
            sage: M = Mat(ZZ,3)
            sage: a = M(range(9)); b = M(reversed(range(9)))
            sage: a - b
            [-8 -6 -4]
            [-2  0  2]
            [ 4  6  8]

        """
        if self._ncols != other._ncols:
            raise IndexError, "Number of columns of self must equal number of columns of other."
        if self._nrows != other._nrows:
            raise IndexError, "Number of rows of self must equal number of rows of other."

        cdef size_t i, j

        cdef Matrix_integer_dense M
        M = Matrix_integer_dense.__new__(Matrix_integer_dense, self._parent, None, None, None)

        _sig_on

        cdef mpz_t *entries
        entries = M._entries
        for i from 0 <= i < self._ncols * self._nrows:
            mpz_init(entries[i])
            mpz_sub(entries[i], self._entries[i], other._entries[i])

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

        M = self._new_unitialized_matrix(self._ncols, self._nrows)
        _sig_on
        for i from 0 <= i < self._ncols:
            for j from 0 <= j < self._nrows:
                mpz_init_set(M._matrix[i][j], self._matrix[j][i])
        _sig_off
        return M

    def matrix_from_rows(self, rows):
        """
        Return the matrix formed from the given rows.

        The output matrix need not be a submatrix of self.  See the
        examples below.

        INPUT:
            rows -- list of int's

        OUTPUT:
            A matrix created from the rows with given indexes

        EXAMPLES:
            sage: a = Mat(ZZ,2)([-1,1,10,3]); a
            [-1  1]
            [10  3]
            sage: a.matrix_from_rows([0])
            [-1  1]
            sage: a.matrix_from_rows([1,0])
            [10  3]
            [-1  1]
            sage: a.matrix_from_rows([1,1])
            [10  3]
            [10  3]
            sage: a.matrix_from_rows([])
            []
            sage: parent(a.matrix_from_rows([]))
            Full MatrixSpace of 0 by 2 dense matrices over Integer Ring
        """
        cdef int i, j, k, nc, nr
        cdef Matrix_integer_dense M

        if not isinstance(rows, (list, tuple)):
            raise TypeError, "rows (=%s) must be a list"%rows

        nr = len(rows)
        if nr == 0:
            return self.new_matrix(0, self._ncols)
        nc = self._ncols
        v = []
        for i in rows:
            v.append(int(i))
        rows = v
        if min(rows) < 0 or max(rows) >= self._nrows:
            raise IndexError, "invalid row indexes; rows don't exist"

        M = self._new_unitialized_matrix(nr, nc)

        cdef mpz_t **m
        m = M._matrix

        for i from 0 <= i < nr:
            k = rows[i]
            for j from 0 <= j < nc:
                mpz_init_set(m[i][j], self._matrix[k][j])

        return M

    def matrix_from_columns(self, columns):
        """
        Return the matrix formed from the given columns.

        The output matrix need not be a submatrix of self.  See the
        examples below.

        INPUT:
            columns -- list of int's

        OUTPUT:
            A matrix created from the columns with given indexes

        EXAMPLES:
            sage: a = Mat(ZZ,2)([-1,1,10,3]); a
            [-1  1]
            [10  3]
            sage: a.matrix_from_columns([0])
            [-1]
            [10]
            sage: a.matrix_from_columns([1,0])
            [ 1 -1]
            [ 3 10]
            sage: a.matrix_from_columns([1,0,1,1,0,1])
            [ 1 -1  1  1 -1  1]
            [ 3 10  3  3 10  3]
            sage: a.matrix_from_columns([1,1])
            [1 1]
            [3 3]
            sage: a.matrix_from_columns([])
            []
            sage: parent(a.matrix_from_columns([]))
            Full MatrixSpace of 2 by 0 dense matrices over Integer Ring
        """
        cdef int i, j, k, nc, nr
        cdef Matrix_integer_dense M

        columns = list(columns)
        if not isinstance(columns, (list, tuple)):
            raise TypeError, "columns (=%s) must be a list"%columns

        nc = len(columns)
        if nc == 0:
            return self.new_matrix(self._nrows, 0)
        nr = self._nrows
        v = []
        for i in columns:
            v.append(int(i))
        columns = v
        if min(columns) < 0 or max(columns) >= self._ncols:
            raise IndexError, "invalid column indexes; columns don't exist"

        M = self._new_unitialized_matrix(nr, nc)

        cdef mpz_t **m
        m = M._matrix

        for i from 0 <= i < nc:
            k = columns[i]
            for j from 0 <= j < nr:
                mpz_init_set(m[j][i], self._matrix[j][k])

        return M


    def iterates(self, v, int n):
        r"""
        Let A be this matrix.   Return a matrix with *rows* the iterates
        $$
          v, v A, v A^2, \ldots, v A^{n-1}.
        $$

        EXAMPLES:
            sage: a = matrix(ZZ,3, range(3^2))
            sage: a.iterates([1,0,0], 4)
            [  1   0   0]
            [  0   1   2]
            [ 15  18  21]
            [180 234 288]
            sage: v = (ZZ^3)([1,0,0])
            sage: v*(a^3)
            (180, 234, 288)
        """
        cdef int i, j, k, nr, nc
        cdef mpz_t s, z, tmp
        cdef sage.rings.integer.Integer x

        nr = n
        nc = self._ncols

        if self._nrows != self._ncols:
            raise ArithmeticError, "matrix must be square"

        v = list(v)
        if len(v) != self._nrows:
            raise ArithmeticError, "incompatible matrix vector multiple"

        cdef Matrix_integer_dense M
        M = self._new_unitialized_matrix(nr, nc)

        cdef mpz_t **m
        m = M._matrix

        for j from 0 <= j < nc:
            try:
                x = v[j]
            except TypeError:
                x = sage.rings.integer.Integer(v[j])
            mpz_init_set(m[0][j], x.value)

        mpz_init(s)
        mpz_init(z)
        for i from 1 <= i < nr:
            for j from 0 <= j < nc:
                mpz_set_si(s,0)  # set s = 0
                for k from 0 <= k < self._nrows:
                    mpz_mul(z, m[i-1][k], self._matrix[k][j])
                    mpz_add(s, s, z)
                mpz_init_set(m[i][j], s)
        mpz_clear(s)
        mpz_clear(z)

        return M


    def scalar_multiple(self, d):
        r"""
        Return the product self*d, as a new matrix.

        EXAMPLES:
            sage: a = Mat(ZZ,2,5)(range(10)); a
            [0 1 2 3 4]
            [5 6 7 8 9]
            sage: a.scalar_multiple(-2)
            [  0  -2  -4  -6  -8]
            [-10 -12 -14 -16 -18]
        """
        cdef size_t i, j
        cdef sage.rings.integer.Integer z
        cdef mpz_t x
        cdef Matrix_integer_dense M

        z = d
        mpz_init_set(x, z.value)
        M = self._new_unitialized_matrix(self._nrows, self._ncols)

        cdef mpz_t *e
        e = M._entries

        _sig_on
        for i from 0 <= i < self._nrows * self._ncols:
            mpz_init(e[i])
            mpz_mul(e[i], self._entries[i], x)
        _sig_off
        mpz_clear(x)
        return M

    def __copy__(self):
        """
        Return a copy of self.  Modifying entries of the copy does not affect self.

        EXAMPLE:
            sage: a = Mat(ZZ,2,5)(range(10)); a
            [0 1 2 3 4]
            [5 6 7 8 9]
            sage: copy(a)
            [0 1 2 3 4]
            [5 6 7 8 9]
            sage: copy(a) == a
            True
        """
        cdef size_t i, j, nr, nc
        nr = self._nrows; nc = self._ncols

        cdef Matrix_integer_dense M
        M = self._new_unitialized_matrix(self._nrows, self._ncols)
        _sig_on
        for i from 0 <= i < nr * nc:
            mpz_init_set(M._entries[i], self._entries[i])
        _sig_off

        return M

    def number_nonzero(self):
        """
        Return the number of nonzero entries in this matrix.

        EXAMPLE:
            sage: a = Mat(ZZ,2,5)(range(10)); a
            [0 1 2 3 4]
            [5 6 7 8 9]
            sage: a.number_nonzero()
            9
        """
        cdef size_t i, j, n
        n = 0
        _sig_on
        for i from 0 <= i < self._nrows * self._ncols:
            if mpz_sgn(self._entries[i]):         # if nonzero
                n = n + 1
        _sig_off
        return int(n)

    def list(self, int base=0):
        r"""
        Return a list of all the elements of self.

        This is the list form of the concatenation of the rows of self.

        EXAMPLE:
            sage: a = Mat(ZZ,2,5)(range(10)); a
            [0 1 2 3 4]
            [5 6 7 8 9]
            sage: a.list()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

        Note that {\tt a.list()} is different than {\tt list(a)} which returns the list of rows:
            sage: list(a)
            [(0, 1, 2, 3, 4), (5, 6, 7, 8, 9)]
        """
        cdef size_t i
        cdef sage.rings.integer.Integer x

        v = []
        _sig_on
        for i from 0 <= i < self._nrows * self._ncols:
            x = sage.rings.integer.Integer()   # todo: change to use new integer
            x.set_from_mpz(self._entries[i])
            v.append(x)
        _sig_off
        return v

    def rank(self):
        r"""
        Return the rank of this matrix.

        EXAMPLES:
             sage: ???
        """
        if self._pivots == None:
            X = self.echelon_form()
            return len(X._pivots)
        return len(self._pivots)

    def pivots(self):
        """
        Return the pivots found during the last echelon operation on self.
        Of course if self is changed, and the echelon form of self is not
        recomputed, then the pivots could be incorrect.

        EXAMPLES:
             sage: ???
        """
        if self._pivots == None:
            raise ArithmeticError, "Echelon form has not yet been computed."
        return self._pivots

    cdef int mpz_height(self, mpz_t height) except -1:
        """
        Used to compute the height of this matrix.

        INPUT:
             height -- a GMP mpz_t (that has not been initialized!)
        OUTPUT:
             sets the value of height to the height of this matrix, i.e., the max absolute
             value of the entries of the matrix.



        EXAMPLES:
             sage: ???
        """
        cdef mpz_t x, h
        cdef size_t i

        mpz_init_set_si(h, 0)
        mpz_init(x)

        _sig_on

        for i from 0 <= i < self._nrows * self._ncols:
            mpz_abs(x, self._entries[i])
            if mpz_cmp(h, x) < 0:
                mpz_set(h, x)

        _sig_off

        mpz_init_set(height, h)
        mpz_clear(h)
        mpz_clear(x)

        return 0   # no error occured.

    def height(self):
        """
        Return the height of this matrix, i.e., the max absolute value
        of the entries of the matrix.

        OUTPUT:
            A nonnegative integer.

        EXAMPLE:
            sage: a = Mat(ZZ,3)(range(9))
            sage: a.height()
            8
            sage: a = Mat(ZZ,2,3)([-17,3,-389,15,-1,0]); a
            [ -17    3 -389]
            [  15   -1    0]
            sage: a.height()
            389
        """
        cdef mpz_t h
        cdef sage.rings.integer.Integer x

        self.mpz_height(h)
        x = sage.rings.integer.Integer()   # todo: change to use new integer
        x.set_from_mpz(h)
        mpz_clear(h)

        return x

    def echelon_form(self, include_zero_rows=True):
        r"""
        Return the echelon form of this matrix over the integers.

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(),2)([1,2,3,4])
            sage: A.echelon_form()
            [1 0]
            [0 2]

            sage: A = MatrixSpace(IntegerRing(),5)(range(25))
            sage: A.echelon_form()
            [  5   0  -5 -10 -15]
            [  0   1   2   3   4]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
        """
        if self.nrows() == 0 or self.ncols() == 0:
            self.__rank = 0
            return self
        # The following complicated sequence of column reversals
        # and transposes is needed since PARI's Hermite Normal Form
        # does column operations instead of row operations.
        n = self.ncols()
        r = []
        for i from 0 <= i < n:
            r.append(n-i)
        v = self._pari_()
        v = v.vecextract(r) # this reverses the order of columns
        v = v.mattranspose()
        w = v.mathnf(1)

        cdef Matrix_integer_dense H_m
        H = convert_parimatrix(w[0])
        if self.ncols() == 1:
            H = [H]

        # We can do a 'fast' change of the above into a list of ints,
        # since we know all entries are ints:
        (nr,nc) = (self.nrows(), self.ncols())
        num_missing_rows = (nr*nc - len(H)) / nc
        rank = nr - num_missing_rows
        if include_zero_rows:
            H = H + ['0']*(num_missing_rows*nc)
            H_m = self.new_matrix(nrows=nr, ncols=nc, entries=H, coerce=True)
        else:
            H_m = self.new_matrix(nrows=rank, ncols=nc, entries=H, coerce=True)
        H_m.__rank = rank
        H_m.set_immutable()
        return H_m

    def elementary_divisors(self):
        """
        Return the elementary divisors of self, in order.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of this matrix.  They are
        ordered in reverse by divisibility.

        INPUT:
            matrix
        OUTPUT:
            list of int's

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
            sage: A.elementary_divisors()
            [0, 3, 1]
            sage: C = MatrixSpace(ZZ,4)([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: C.elementary_divisors()
            [687, 1, 1, 1]

        SEE ALSO: smith_form
        """
        try:
            return self.__elementary_divisors
        except AttributeError:
            if self.nrows() == 0 or self.ncols() == 0:
                return []
            d = self._pari_().matsnf(0).python()
            if self.is_immutable():
                self.__elementary_divisors = d
            return d

    def smith_form(self, transformation=False):
        """
        Returns matrices S, U, and V such that S = U*self*V, and S
        is in Smith normal form.  Thus S is diagonal with diagonal
        entries the ordered elementary divisors of S.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of this matrix.  They are
        ordered in reverse by divisibility.

        EXAMPLES:
            sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
            sage: D, U, V = A.smith_form()
            sage: D
            [0 0 0]
            [0 3 0]
            [0 0 1]
            sage: U
            [-1  2 -1]
            [ 0 -1  1]
            [ 0  1  0]
            sage: V
            [ 1  4 -1]
            [-2 -3  1]
            [ 1  0  0]
            sage: U*A*V
            [0 0 0]
            [0 3 0]
            [0 0 1]

        It also makes sense for nonsquare matrices:

            sage: A = Matrix(ZZ,3,2,range(6))
            sage: D, U, V = A.smith_form()
            sage: D
            [0 0]
            [2 0]
            [0 1]
            sage: U
            [-1  2 -1]
            [ 0 -1  1]
            [ 0  1  0]
            sage: V
            [ 3 -1]
            [-2  1]
            sage: U * A * V
            [0 0]
            [2 0]
            [0 1]

        SEE ALSO: elementary_divisors
        """
        v = self._pari_().matsnf(1).python()
        print "v =", v
        D = self.matrix_space()(v[2])
        U = self.matrix_space(ncols = self.nrows())(v[0])
        V = self.matrix_space(nrows = self.ncols())(v[1])
        return D, U, V

    def frobenius(self,flag=0):
        """
        Return the Frobenius form (rational canonical form) of this matrix.

        If flag is 1, return only the elementary divisors.  If flag is
        2, return a two-components vector [F,B] where F is the
        Frobenius form and B is the basis change so that $M=B^{-1}FB$.

        INPUT:
           flag -- 0 (default), 1 or 2 as described above

        ALGORITHM: uses pari's matfrobenius()

        EXAMPLE:
           sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
           sage: A.frobenius(0)
           [ 0  0  0]
           [ 1  0 18]
           [ 0  1 12]
           sage: A.frobenius(1)
           [x^3 - 12*x^2 - 18*x]
           sage: A.frobenius(2)
           ([ 0  0  0]
           [ 1  0 18]
           [ 0  1 12],
           [    -1      2     -1]
           [     0  23/15 -14/15]
           [     0  -2/15   1/15])

        AUTHOR:
           -- 2006-04-02: Martin Albrecht

        TODO:
           -- move this to work for more general matrices than just over Z.
              This will require fixing how PARI polynomials are coerced
              to SAGE polynomials.
        """
        if self.nrows()!=self.ncols():
            raise ArithmeticError, "frobenius matrix of non-square matrix not defined."

        v = self._pari_().matfrobenius(flag)
        if flag==0:
            return self.matrix_space()(v.python())
        elif flag==1:
            r = polynomial_ring.PolynomialRing(self.base_ring())
            #TODO: this should be handled in PolynomialRing not here
            retr = []
            for x in v:
                retr.append(eval(str(x).replace("^","**"), {}, r.gens_dict()))
            return retr
        elif flag==2:
            F = matrix_space.MatrixSpace(rational_field.RationalField(),
                                         self.nrows())(v[0].python())
            B = matrix_space.MatrixSpace(rational_field.RationalField(),
                                         self.nrows())(v[1].python())
            return F,B

    def kernel(self, LLL=False):
        r"""
        Return the kernel of this matrix, as a module over the integers.

        INPUT:
           LLL -- bool (default: False); if True the basis is an LLL
                  reduced basis; otherwise, it is an echelon basis.

        EXAMPLES:
            sage: M = MatrixSpace(IntegerRing(),4,2)(range(8))
            sage: M.kernel()
            Free module of degree 4 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -3  2]
            [ 0  1 -2  1]
        """

        Z = self.base_ring()

        if self.nrows() == 0:    # from a 0 space
            M = sage.modules.free_module.FreeModule(Z, self.nrows())
            return M.zero_subspace()

        elif self.ncols() == 0:  # to a 0 space
            return sage.modules.free_module.FreeModule(Z, self.nrows())

        A = self._pari_().mattranspose()
        B = A.matkerint()


        n = self.nrows()
        Z = sage.rings.integer_ring.IntegerRing()
        M = sage.modules.free_module.FreeModule(Z, n)

        if B.ncols() == 0:
            return M.zero_submodule()

        # Now B is a basis for the LLL-reduced integer kernel as a
        # PARI object.  The basis vectors or B[0], ..., B[n-1],
        # where n is the dimension of the kernel.
        B = []
        for b in B:
          tmp = []
          for x in b:
              tmp.append(Z(x))
          B.append(M(tmp))
#        B = [M([Z(x) for x in b]) for b in B]
        if LLL:
            return M.span_of_basis(B)
        else:
            return M.span(B)

    def _adjoint(self):
        """assumes self is a square matrix (checked in adjoint)"""
        return self.parent()(self._pari_().matadjoint().python())

    def _lllgram(self):
        """assumes self is a square matrix (checked in lllgram)"""
        Z = sage.rings.integer_ring.IntegerRing()
        n = self.nrows()
        # pari does not like negative definite forms
        if n > 0 and self[0,0] < 0:
            self = -self
        # maybe should be /unimodular/ matrices ?
        MS = matrix_space.MatrixSpace(Z,n,n)
        try:
            U = MS(self._pari_().lllgramint().python())
        except (RuntimeError, ArithmeticError):
            raise ValueError, "not a definite matrix"
        # Fix last column so that det = +1
        if U.det() == -1:
            for i in range(n):
                U[i,n-1] = - U[i,n-1]
        return U

    def _ntl_(self):
        """
        ntl.mat_ZZ representation of self.

        \note{NTL only knows dense matrices, so if you provide a
        sparse matrix NTL will allocate memory for every zero entry.}
        """
        return mat_ZZ(self.nrows(),self.ncols(),self.list())



    def multiply_multi_modular(self, Matrix_integer_dense right):
        """
        Multiply this matrix by right using a multimodular algorithm
        and return the result.

        EXAMPLES:
            sage: ???
        """
        raise NotImplementedError


cdef object mpz_to_long(mpz_t x):
    return long(mpz_to_str(x))

def _parimatrix_to_strlist(A):
    s = str(A)
    s = s.replace('Mat(','').replace(')','')
    s = s.replace(';',',').replace(' ','')
    s = s.replace(",", "','")
    s = s.replace("[", "['")
    s = s.replace("]", "']")
    return eval(s)

def _parimatrix_to_reversed_strlist(A):
    s = str(A)
    if s.find('Mat') != -1:
        return _parimatrix_to_strlist(A)
    s = s.replace('[','').replace(']','').replace(' ','')
    v = s.split(';')
    v.reverse()
    s = "['" + (','.join(v)) + "']"
    s = s.replace(",", "','")
    return eval(s)

def convert_parimatrix(z):
    n = z.ncols();
    r = []
    for i from 0 <= i < n:
        r.append(n-i)
    z = z.vecextract(r)
    z = z.mattranspose()
    n = z.ncols();
    r = []
    for i from 0 <= i < n:
        r.append(n-i)
    z = z.vecextract(r)
    return _parimatrix_to_strlist(z)

