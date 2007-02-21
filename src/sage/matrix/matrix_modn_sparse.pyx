r"""
Sparse matrices over $\Z/n\Z$ for $n$ small.

This is a compiled implementation of sparse matrices over $\Z/n\Z$
for $n$ small.

AUTHOR:
   -- William Stein (2004): first version
   -- William Stein (2006-02-12): added set_row_to_multiple_of_row
   -- William Stein (2006-03-04): added multimodular echelon, __reduce__, etc.
   -- William Wtein (2006-08-09): Pyrexification

TODO:
   -- move vectors into a pyrex vector class
   -- add _add_ and _mul_ methods.

EXAMPLES:
    sage: a = matrix(Integers(37),3,3,range(9),sparse=True); a
    [0 1 2]
    [3 4 5]
    [6 7 8]
    sage: type(a)
    <type 'sage.matrix.matrix_modn_sparse.Matrix_modn_sparse'>
    sage: parent(a)
    Full MatrixSpace of 3 by 3 sparse matrices over Ring of integers modulo 37
    sage: a^2
    [15 18 21]
    [ 5 17 29]
    [32 16  0]
    sage: a+a
    [ 0  2  4]
    [ 6  8 10]
    [12 14 16]
    sage: b = a.new_matrix(2,3,range(6)); b
    [0 1 2]
    [3 4 5]
    sage: a*b
    Traceback (most recent call last):
    ...
    TypeError: incompatible dimensions
    sage: b*a
    [15 18 21]
    [ 5 17 29]

    sage: a == loads(dumps(a))
    True
    sage: b == loads(dumps(b))
    True

    sage: a.echelonize(); a
    [ 1  0 36]
    [ 0  1  2]
    [ 0  0  0]
    sage: b.echelonize(); b
    [ 1  0 36]
    [ 0  1  2]
    sage: a.pivots()
    [0, 1]
    sage: b.pivots()
    [0, 1]
    sage: a.rank()
    2
    sage: b.rank()
    2
    sage: a[2,2] = 5
    sage: a.rank()
    3
"""

#############################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

include "../ext/cdefs.pxi"
include '../ext/interrupt.pxi'
include '../ext/stdsage.pxi'

cimport matrix
cimport matrix_sparse
from sage.rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract

from sage.misc.misc import verbose, get_verbose

################
# TODO: change this to use extern cdef's methods.
from sage.ext.arith cimport arith_int
cdef arith_int ai
ai = arith_int()
################

import sage.ext.multi_modular
MAX_MODULUS = sage.ext.multi_modular.MAX_MODULUS

cdef class Matrix_modn_sparse(matrix_sparse.Matrix_sparse):

    ########################################################################
    # LEVEL 1 functionality
    # x * __new__
    # x * __dealloc__
    # x * __init__
    # x * set_unsafe
    # x * get_unsafe
    # x * __richcmp__    -- always the same
    ########################################################################
    def __new__(self, parent, entries, copy, coerce):
        # allocate memory
        cdef Py_ssize_t i, nr, nc
        cdef int p

        nr = parent.nrows()
        nc = parent.ncols()
        p = parent.base_ring().order()

        _sig_on
        self.rows = <c_vector_modint*> sage_malloc(nr*sizeof(c_vector_modint))
        _sig_off
        if not self.rows:
            raise MemoryError, "error allocating memory for sparse matrix"

        for i from 0 <= i < nr:
            init_c_vector_modint(&self.rows[i], p, nc, 0)


    def __dealloc__(self):
        cdef int i
        for i from 0 <= i < self._nrows:
            clear_c_vector_modint(&self.rows[i])


    def __init__(self, parent, entries, copy, coerce):
        """
        Create a sparse matrix modulo n.

        INPUT:
            parent -- a matrix space
            entries -- * a Python list of triples (i,j,x), where 0 <= i < nrows,
                         0 <= j < ncols, and x is coercible to an int.  The i,j
                         entry of self is set to x.  The x's can be 0.
                       * Alternatively, entries can be a list of *all* the entries
                         of the sparse matrix (so they would be mostly 0).
            copy -- ignored
            coerce -- ignored
        """
        cdef int s, y, z, p
        cdef Py_ssize_t i, j, k

        cdef object seq
        cdef void** X

        matrix.Matrix.__init__(self, parent)

        self.p = parent.base_ring().order()
        p = self.p

        if isinstance(entries, dict):
            # Sparse input format.
            R = self._base_ring
            for ij, x in entries.iteritems():
                y = x
                z = R(y)
                if z != 0:
                    i, j = ij  # nothing better to do since this is user input, which may be bogus.
                    if i < 0 or j < 0 or i >= self._nrows or j >= self._ncols:
                        raise IndexError, "invalid entries list"
                    set_entry(&self.rows[i], j, z)
        elif isinstance(entries, list):
            # Dense input format
            if len(entries) != self._nrows * self._ncols:
                raise TypeError, "list of entries must be a dictionary of (i,j):x or a dense list of n * m elements"
            seq = PySequence_Fast(entries,"expected a list")
            X = PySequence_Fast_ITEMS(seq)
            k = 0
            R = self._base_ring
            # Get fast access to the entries list.
            for i from 0 <= i < self._nrows:
                for  j from 0 <= j < self._ncols:
                    set_entry(&self.rows[i], j, R(<object>X[k]))
                    k = k + 1
        else:
            # scalar?
            s = int(self._base_ring(entries))
            if s == 0:
                return
            if self._nrows != self._ncols:
                raise TypeError, "matrix must be square to initialize with a scalar."
            for i from 0 <= i < self._nrows:
                set_entry(&self.rows[i], i, s)


    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        set_entry(&self.rows[i], j, (<IntegerMod_int> value).ivalue)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef IntegerMod_int n
        n =  IntegerMod_int.__new__(IntegerMod_int)
        IntegerMod_abstract.__init__(n, self._base_ring)
        n.ivalue = get_entry(&self.rows[i], j)
        return n

    def __richcmp__(matrix.Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        return self._hash()


    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    #   * cdef _add_c_impl
    #   * cdef _mul_c_impl
    #   * cdef _cmp_c_impl
    #   * __neg__
    #   * __invert__
    #   * __copy__
    #   * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################
    # def _pickle(self):
    # def _unpickle(self, data, int version):   # use version >= 0
    # cdef ModuleElement _add_c_impl(self, ModuleElement right):
    # cdef _mul_c_impl(self, Matrix right):
    # cdef int _cmp_c_impl(self, Matrix right) except -2:
    # def __neg__(self):
    # def __invert__(self):
    # def __copy__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):
    # def _dict(self):


    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_c_impl
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    # x      - echelon form in place
    ########################################################################
    def swap_rows(self, r1, r2):
        self.check_bounds_and_mutability(r1,0)
        self.check_bounds_and_mutability(r2,0)
        self.swap_rows_c(r1, r2)

    cdef swap_rows_c(self, Py_ssize_t n1, Py_ssize_t n2):
        """
        Swap the rows in positions n1 and n2.  No bounds checking.
        """
        cdef c_vector_modint tmp
        tmp = self.rows[n1]
        self.rows[n1] = self.rows[n2]
        self.rows[n2] = tmp

    def _echelon_in_place_classical(self):
        """
        Replace self by its reduction to reduced row echelon form.

        ALGORITHM: We use Gauss elimination, in a slightly intelligent
        way, in that we clear each column using a row with the minimum
        number of nonzero entries.

        TODO: Implement switching to a dense method when the matrix
        gets dense.
        """
        x = self.fetch('in_echelon_form')
        if not x is None and x: return  # already known to be in echelon form
        self.check_mutability()

        cdef Py_ssize_t i, r, c, min, min_row,  start_row
        cdef int a0, a_inverse, b, do_verb
        cdef c_vector_modint tmp
        start_row = 0
        pivots = []
        fifth = self._ncols / 10 + 1
        tm = verbose(caller_name = 'sparse_matrix_pyx matrix_modint echelon')
        do_verb = (get_verbose() >= 2)
        for c from 0 <= c < self._ncols:
            if do_verb and (c % fifth == 0 and c>0):
                tm = verbose('on column %s of %s'%(c, self._ncols),
                             level = 2,
                             caller_name = 'matrix_modn_sparse echelon')
            #end if
            min = self._ncols + 1
            min_row = -1
            for r from start_row <= r < self._nrows:
                if self.rows[r].num_nonzero > 0 and self.rows[r].num_nonzero < min:
                    # Since there is at least one nonzero entry, the first entry
                    # of the positions list is defined.  It is the first position
                    # of a nonzero entry, and it equals c precisely if row r
                    # is a row we could use to clear column c.
                    if self.rows[r].positions[0] == c:
                        min_row = r
                        min = self.rows[r].num_nonzero
                    #endif
                #endif
            #endfor

            if min_row != -1:
                r = min_row
                #print "min number of entries in a pivoting row = ", min
                pivots.append(c)
                _sig_on
                # Since we can use row r to clear column c, the
                # entry in position c in row r must be the first nonzero entry.
                a = self.rows[r].entries[0]
                if a != 1:
                    a_inverse = ai.c_inverse_mod_int(a, self.p)
                    scale_c_vector_modint(&self.rows[r], a_inverse)
                self.swap_rows_c(r, start_row)
                for i from 0 <= i < self._nrows:
                    if i != start_row:
                        b = get_entry(&self.rows[i], c)
                        if b != 0:
                            add_c_vector_modint_init(&tmp, &self.rows[i],
                                                     &self.rows[start_row], self.p - b)
                            clear_c_vector_modint(&self.rows[i])
                            self.rows[i] = tmp
                start_row = start_row + 1
                _sig_off

        self.cache('pivots',pivots)
        self.cache('in_echelon_form',True)





###########################################################################################
# A mini C-level library just for sparse vectors modulo an integer.
#
#  IMPLEMENTATION NOTES:
#  A *sparse* c_vector_modint is stored using the following data
#  structure.  The entries and positions arrays are of the same length.
#  The entries array contains the nonzero elements of the vector, and
#  positions contains the SORTED list of the positions of those elements.
## The num_nonzero integer is the number of nonzero entries of the
## vector, which is the length of the arries entries and positions.  The
## degree is the dimension of the ambient vector space, and p is the
## prime (so it's a vector modulo p).
##
###########################################################################################
cdef int allocate_c_vector_modint(c_vector_modint* v, Py_ssize_t num_nonzero) except -1:
    """
    Allocate memory for a c_vector_modint -- most user code won't call this.
    """
    v.entries = <int*>sage_malloc(num_nonzero*sizeof(int))
    if v.entries == NULL:
        raise MemoryError, "Error allocating memory"
    v.positions = <Py_ssize_t*>sage_malloc(num_nonzero*sizeof(Py_ssize_t))
    if v.positions == NULL:
        raise MemoryError, "Error allocating memory"
    return 0

cdef int init_c_vector_modint(c_vector_modint* v, int p, Py_ssize_t degree,
                              Py_ssize_t num_nonzero) except -1:
    """
    Initialize a c_vector_modint.
    """
    if (allocate_c_vector_modint(v, num_nonzero) == -1):
        raise MemoryError, "Error allocating memory for sparse vector."
    if p > 46340:
        raise OverflowError, "The prime must be <= 46340."
    v.num_nonzero = num_nonzero
    v.degree = degree
    v.p = p
    return 0

cdef void clear_c_vector_modint(c_vector_modint* v):
    sage_free(v.entries)
    sage_free(v.positions)

cdef Py_ssize_t binary_search0(Py_ssize_t* v, Py_ssize_t n, int x):
    """
    Find the position of the int x in the array v, which has length n.
    Returns -1 if x is not in the array v.
    """
    if n == 0:
        return -1

    cdef Py_ssize_t i, j, k
    i = 0
    j = n-1
    while i<=j:
        if i == j:
            if v[i] == x:
                return i
            return -1
        k = (i+j)/2
        if v[k] > x:
            j = k-1
        elif v[k] < x:
            i = k+1
        else:   # only possibility is that v[k] == x
            return k
    return -1

cdef Py_ssize_t binary_search(Py_ssize_t* v, Py_ssize_t n, int x, Py_ssize_t* ins):
    """
    Find the position of the integer x in the array v, which has length n.
    Returns -1 if x is not in the array v, and in this case ins is
    set equal to the position where x should be inserted in order to
    obtain an ordered array.
    """
    if n == 0:
        ins[0] = 0
        return -1

    cdef Py_ssize_t i, j, k
    i = 0
    j = n-1
    while i<=j:
        if i == j:
            if v[i] == x:
                ins[0] = i
                return i
            if v[i] < x:
                ins[0] = i + 1
            else:
                ins[0] = i
            return -1
        k = (i+j)/2
        if v[k] > x:
            j = k-1
        elif v[k] < x:
            i = k+1
        else:   # only possibility is that v[k] == x
            ins[0] = k
            return k
    # end while
    ins[0] = j+1
    return -1

cdef int get_entry(c_vector_modint* v, Py_ssize_t n) except -1:
    """
    Returns the n-th entry of the sparse vector v.  This
    would be v[n] in Python syntax.
    """
    if n >= v.degree or n < 0:
        raise IndexError, "Index must be between 0 and the degree minus 1."
        return -1
    cdef Py_ssize_t m
    m = binary_search0(v.positions, v.num_nonzero, n)
    if m == -1:
        return 0
    return v.entries[m]

cdef object to_list(c_vector_modint* v):
    """
    Returns a Python list of 2-tuples (i,x), where x=v[i] runs
    through the nonzero elements of x, in order.
    """
    cdef object X
    cdef Py_ssize_t i
    X = []
    for i from 0 <= i < v.num_nonzero:
        X.append( (v.positions[i], v.entries[i]) )
    return X

cdef int set_entry(c_vector_modint* v, Py_ssize_t n, int x) except -1:
    """
    Set the n-th component of the sparse vector v equal to x.
    This would be v[n] = x in Python syntax.
    """
    if n < 0 or n >= v.degree:
        raise IndexError, "Index (=%s) must be between 0 and %s."%(n, v.degree-1)
        return -1
    cdef Py_ssize_t i, m, ins
    cdef Py_ssize_t m2, ins2
    cdef Py_ssize_t *pos
    cdef int *e

    x = x % v.p
    if x<0: x = x + v.p
    m = binary_search(v.positions, v.num_nonzero, n, &ins)

    if m != -1:
        # The position n was found in the array of positions.
        # Now there are two cases:
        #   1. x =/= 0, which is easy, and
        #   2. x = 0, in which case we have to recopy
        #      positions and entries, without the m-th
        #      element, and change num_nonzero.
        if x != 0:   # case 1
            v.entries[m] = x
        else:        # case 2
            e = v.entries
            pos = v.positions
            allocate_c_vector_modint(v, v.num_nonzero - 1)
            for i from 0 <= i < m:
                v.entries[i] = e[i]
                v.positions[i] = pos[i]
            for i from m < i < v.num_nonzero:
                v.entries[i-1] = e[i]
                v.positions[i-1] = pos[i]
            sage_free(e)
            sage_free(pos)
            v.num_nonzero = v.num_nonzero - 1
    else:
        # Allocate new memory and copy over elements from the
        # old array.  This is similar to case 2 above,
        # except we are inserting a new entry rather than
        # deleting an old one.  The new entry should be inserted
        # at position ins, which was computed using binary search.
        #
        # There is one exception -- if the new entry is 0, we
        # do nothing and return.
        if x == 0:
            return 0
        v.num_nonzero = v.num_nonzero + 1
        e = v.entries
        pos = v.positions
        allocate_c_vector_modint(v, v.num_nonzero)
        for i from 0 <= i < ins:
            v.entries[i] = e[i]
            v.positions[i] = pos[i]
        v.entries[ins] = x
        v.positions[ins] = n
        for i from ins < i < v.num_nonzero:
            v.entries[i] = e[i-1]
            v.positions[i] = pos[i-1]
        sage_free(e)
        sage_free(pos)

cdef int add_c_vector_modint_init(c_vector_modint* sum, c_vector_modint* v,
                                  c_vector_modint* w, int multiple) except -1:
    """
    Set sum = v + multiple*w.
    """
    if v.p != w.p:
        raise ArithmeticError, "The vectors must be modulo the same prime."
        return -1
    if v.degree != w.degree:
        raise ArithmeticError, "The vectors must have the same degree."
        return -1

    cdef int s
    cdef Py_ssize_t nz, i, j, k
    cdef c_vector_modint* z

    multiple = multiple % v.p    # need this to avoid overflow.
    if multiple < 0:
        multiple = multiple + v.p

    z = sum
    # ALGORITHM:
    # 1. Allocate enough memory to hold the union of the two
    #    lists of positions.  We allocate the sum of the number
    #    of positions of both (up to the degree), to avoid
    #    having to make two passes.  This might be slightly wasteful of
    #    memory, but is faster.
    # 2. Move along the entries of v and w, copying them into the
    #    new position / entry array.  When position are the same,
    #    add modulo p.
    # 3. Set num_nonzero and return success code.

    # 1. Allocate memory:
    nz = v.num_nonzero + w.num_nonzero
    if nz > v.degree: nz = v.degree
    init_c_vector_modint(z, v.p, v.degree, nz)
    # 2. Merge entries
    i = 0  # index into entries of v
    j = 0  # index into entries of w
    k = 0  # index into z (the vector we are creating)
    while i < v.num_nonzero or j < w.num_nonzero:
        if i >= v.num_nonzero:   # just copy w in
            z.positions[k] = w.positions[j]
            z.entries[k] = (multiple*w.entries[j])%v.p
            j = j + 1
            k = k + 1
        elif j >= w.num_nonzero:  # just copy v in
            z.positions[k] = v.positions[i]
            z.entries[k] = v.entries[i]
            i = i + 1
            k = k + 1
        elif v.positions[i] < w.positions[j]:  # copy entry from v in
            z.positions[k] = v.positions[i]
            z.entries[k] = v.entries[i]
            i = i + 1
            k = k + 1
        elif v.positions[i] > w.positions[j]: # copy entry from w in
            s = (multiple*w.entries[j])%v.p
            if s != 0:
                z.positions[k] = w.positions[j]
                z.entries[k] = s
                k = k + 1
            j = j + 1
        else:                                 # equal, so add and copy
            s = (v.entries[i] + multiple*w.entries[j]) % v.p
            if s != 0:
                z.positions[k] = v.positions[i]
                z.entries[k] = s
                k = k + 1     # only increment if sum is nonzero!
            i = i + 1
            j = j + 1
        #end if
    # end while
    z.num_nonzero = k
    return 0

cdef int scale_c_vector_modint(c_vector_modint* v, int scalar) except -1:
    scalar = scalar % v.p
    if scalar == 0:
        clear_c_vector_modint(v)
        init_c_vector_modint(v, v.p, v.degree, 0)
        return 0
    if scalar < 0:
        scalar = scalar + v.p
    cdef Py_ssize_t i
    for i from 0 <= i < v.num_nonzero:
        v.entries[i] = (v.entries[i] * scalar) % v.p
    return 0

