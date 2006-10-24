r"""
Sparse matrices over $\F_p$ and $\Q$

This is a compiled implementation of sparse matrix algebra over small
prime finite fields and the rational numbers, which is used mainly
internally by other classes.

IMPLEMENTATION NOTES:
A *sparse* c_vector_modint is stored using the following data
structure.  The entries and positions arrays are of the same length.
The entries array contains the nonzero elements of the vector, and
positions contains the SORTED list of the positions of those elements.
The num_nonzero integer is the number of nonzero entries of the
vector, which is the length of the arries entries and positions.  The
degree is the dimension of the ambient vector space, and p is the
prime (so it's a vector modulo p).

AUTHOR:
   -- William Stein (2004): first version
   -- William Stein (2006-02-12): added set_row_to_multiple_of_row
   -- William Stein (2006-03-04): added multimodular echelon, __reduce__, etc.
   -- William Wtein (2006-08-09): Pyrexification

TODO:
   -- move vectors into a pyrex vector class
   -- add _add_ and _mul_ methods.
"""

#############################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#############################################################################


import random

import sage.misc.all
import sage.rings.finite_field
import sage.rings.arith
import sage.rings.integer

cimport matrix_generic
import matrix_generic

cimport sage.ext.arith
import sage.ext.arith
cdef sage.ext.arith.arith_int ai
ai = sage.ext.arith.arith_int()

include "../ext/cdefs.pxi"
include '../ext/interrupt.pxi'

cdef int allocate_c_vector_modint(c_vector_modint* v, int num_nonzero) except -1:
    """
    Allocate memory for a c_vector_modint -- most user code won't call this.
    """
    v.entries = <int*>PyMem_Malloc(num_nonzero*sizeof(int))
    if v.entries == <int*> 0:
        raise MemoryError, "Error allocating memory"
    v.positions = <int*>PyMem_Malloc(num_nonzero*sizeof(int))
    if v.positions == <int*> 0:
        raise MemoryError, "Error allocating memory"
    return 0

cdef int init_c_vector_modint(c_vector_modint* v, int p, int degree,
                              int num_nonzero) except -1:
    """
    Initialize a c_vector_modint -- most user code *will* call this.
    """
    if (allocate_c_vector_modint(v, num_nonzero) == -1):
        return -1
    if p > 46340:
        raise OverflowError, "The prime must be <= 46340."
    v.num_nonzero = num_nonzero
    v.degree = degree
    v.p = p

cdef void clear_c_vector_modint(c_vector_modint* v):
    PyMem_Free(v.entries)
    PyMem_Free(v.positions)

cdef int linear_search0(int* v, int n, int x):
    """
    Find the position of the integer x in the array v, which has length n.
    Returns -1 if x is not in the array v.

    (This is a drop-in replacement for binary_search0, which we do
    not use in practice.  I just wrote it to make sure binary_search0
    was working correctly.)
    """
    if n == 0:
        return -1
    cdef int k
    for k from 0 <= k < n:
        if v[k] == x:
            return k
    return -1

cdef int binary_search0(int* v, int n, int x):
    """
    Find the position of the integer x in the array v, which has length n.
    Returns -1 if x is not in the array v.
    """
    if n == 0:
        return -1

    cdef int i, j, k
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


cdef int linear_search(int* v, int n, int x, int* ins):
    """
    Find the position of the integer x in the array v, which has length n.
    Returns -1 if x is not in the array v, and in this case ins is
    set equal to the position where x should be inserted in order to
    get the array ordered.

    (This is a drop-in replacement for binary_search, which we do
    not use in practice.  I just wrote it to make sure binary_search
    was working.)
    """
    if n == 0:
        ins[0] = 0
        return -1

    cdef int k
    for k from 0 <= k < n:
        if v[k] >= x:
            ins[0] = k
            if v[k] == x:
                return k
            else:
                return -1
    ins[0] = n
    return -1

cdef int binary_search(int* v, int n, int x, int* ins):
    """
    Find the position of the integer x in the array v, which has length n.
    Returns -1 if x is not in the array v, and in this case ins is
    set equal to the position where x should be inserted in order to
    obtain an ordered array.
    """
    if n == 0:
        ins[0] = 0
        return -1

    cdef int i, j, k
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

def bs(v, x):
    cdef int* w
    w = <int*>PyMem_Malloc(sizeof(int)*len(v))
    for i from 0 <= i < len(v):
        w[i] = v[i]
    cdef int ins, b
    b = binary_search(w, len(v), x, &ins)
    s1 = (b, ins)
    b = linear_search(w, len(v), x, &ins)
    s2 = (b, ins)
    PyMem_Free(w)
    return s1, s2

cdef int get_entry(c_vector_modint* v, int n) except -1:
    """
    Returns the n-th entry of the sparse vector v.  This
    would be v[n] in Python syntax.
    """
    if n >= v.degree or n < 0:
        raise IndexError, "Index must be between 0 and the degree minus 1."
        return -1
    cdef int m
    m = binary_search0(v.positions, v.num_nonzero, n)
    #m = linear_search0(v.positions, v.num_nonzero, n)
    if m == -1:
        return 0
    return v.entries[m]

cdef object to_list(c_vector_modint* v):
    """
    Returns a Python list of 2-tuples (i,x), where x=v[i] runs
    through the nonzero elements of x, in order.
    """
    cdef object X
    cdef int i
    X = []
    for i from 0 <= i < v.num_nonzero:
        X.append( (v.positions[i], v.entries[i]) )
    return X

cdef int set_entry(c_vector_modint* v, int n, int x) except -1:
    """
    Set the n-th component of the sparse vector v equal to x.
    This would be v[n] = x in Python syntax.
    """
    if n >= v.degree:
        raise IndexError, "Index must be between 0 and the degree minus 1."
        return -1
    cdef int i, m, ins
    cdef int m2, ins2
    cdef int *e, *pos

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
            PyMem_Free(e)
            PyMem_Free(pos)
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
        PyMem_Free(e)
        PyMem_Free(pos)

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

    cdef int nz, i, j, k, s
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
    cdef int i
    for i from 0 <= i < v.num_nonzero:
        v.entries[i] = (v.entries[i] * scalar) % v.p
    return 0


cdef class Vector_modint:
    """
    Vector_modint -- a sparse vector modulo p.  This is a Python
    extension type that wraps the C implementation of sparse vectors
    modulo a small prime.
    """
    cdef c_vector_modint v

    def __init__(self, int p, int degree, int num_nonzero=0):
        init_c_vector_modint(&self.v, p, degree, num_nonzero)

    def __dealloc__(self):
        clear_c_vector_modint(&self.v)

    def __getitem__(self, int n):
        return get_entry(&self.v, n)

    def __setitem__(self, int n, int x):
        set_entry(&self.v, n, x)

    def __repr__(self):
        return str(list(self))

    def p(self):
        return self.v.p

    def degree(self):
        return self.v.degree

    def num_nonzero(self):
        return self.v.num_nonzero

    def __list__(self):
        return to_list(&self.v)
#         cdef int i
#         v = []
#         for i from 0 <= i < v.num_nonzero:
#             v.append((self.v.positions[i], self.v.entries[i]))
#         return v

    def rescale(self, int x):
        scale_c_vector_modint(&self.v, x)

    def __add__(Vector_modint self, Vector_modint other):
        cdef c_vector_modint z1, *z2
        cdef Vector_modint w

        add_c_vector_modint_init(&z1, &self.v, &other.v, 1)
        w = Vector_modint(self.v.p, self.v.degree)
        z2 = &(w.v)
        clear_c_vector_modint(z2)   # free memory wasted on allocated w
        z2.entries = z1.entries
        z2.positions = z1.positions
        z2.num_nonzero = z1.num_nonzero
        # at this point we do *not* free z1, since it is referenced by w.
        return w

    def __sub__(Vector_modint self, Vector_modint other):
        return self + other*(-1)

    def copy(self):
        cdef int i
        cdef Vector_modint w
        w = Vector_modint(self.v.p, self.v.degree, self.v.num_nonzero)
        for i from 0 <= i < self.v.num_nonzero:
            w.v.entries[i] = self.v.entries[i]
            w.v.positions[i] = self.v.positions[i]
        return w

    def __mul__(x, y):
        if isinstance(x, Vector_modint) and isinstance(y,int):
            self = x
            other = y
        elif isinstance(y, Vector_modint) and isinstance(x, int):
            self = y
            other = x
        else:
            raise TypeError, "Invalid types."
        z = self.copy()
        z.rescale(other)
        return z

    def randomize(self, int sparcity, bound=3):
        """
        randomize(self, int sparcity, exact=False):

        The sparcity is a bound on the number of nonzeros per row.
        """
        cdef int i
        for i from 0 <= i < sparcity:
            self[random.randrange(self.v.degree)] = random.randrange(1,bound)

cdef class Matrix_modn_sparse(matrix_sparse.Matrix_sparse):

    def __new__(self, parent, int p, int nrows, int ncols, object entries=[]):
        # allocate memory
        cdef int i
        self.rows = <c_vector_modint*> PyMem_Malloc(nrows*sizeof(c_vector_modint))
        for i from 0 <= i < nrows:
            init_c_vector_modint(&self.rows[i], p, ncols, 0)

    def __dealloc__(self):
        cdef int i
        for i from 0 <= i < self.nr:
            clear_c_vector_modint(&self.rows[i])


    def __init__(self, parent, int p, int nrows, int ncols, object entries=[]):
        """
        p should be a prime < 46340 and entries should be a list of triples
        (i,j,x), where 0 <= i < nrows, 0 <= j < ncols, and x is an integer
        which need not be reduced modulo p.  Then the i,j entry of the matrix
        is set to x.  It is OK for some x to be zero.
        """
        cdef int y, z

        matrix_generic.Matrix.__init__(self, parent)

        self.p = p
        self.nr = nrows
        self.nc = ncols
        self.__pivots = []

        if len(entries) == 0:
            return

        _sig_on
        for i, j, x in entries:
            y = x
            z = y % p
            if z < 0:
                z = z + p
            if z:
                self[i,j] = z
        _sig_off

    def pivots(self):
        return self.__pivots

    def lift_row_to_dict(self, int i):
        """
        Return an associative arrow of pairs
               n:x
        where the keys n run through the nonzero positions of the row,
        and the x are nonzero and of type Integer.
        """
        cdef int j, n
        cdef object entries
        if i < 0 or i >= self.nr: raise IndexError
        X = {}
        for j from 0 <= j < self.rows[i].num_nonzero:
            n = self.rows[i].positions[j]
            x = sage.rings.integer.Integer()
            x.set_si(self.rows[i].entries[j])
            X[n] = x
        return X

    def lift_to_dict(self):
        """
        Return an associative arrow of pairs
               (i,j):x
        where the keys (i,j) run through the nonzero positions of the matrix
        and the x are nonzero and of type Integer.
        """
        cdef int i, j, n
        cdef object entries
        X = {}
        for i from 0 <= i < self.nr:
            for j from 0 <= j < self.rows[i].num_nonzero:
                n = self.rows[i].positions[j]
                x = sage.rings.integer.Integer()
                x.set_si(self.rows[i].entries[j])
                X[(i,n)] = x
        return X

    def randomize(self, int sparcity, exact=False):
        """
        randomize(self, int sparcity, exact=False):

        INPUT:
            sparcity -- integer
            exact -- bool (default: False)

        The sparcity is (an upper bound) on the number of nonzeros per row.
        """
        cdef int i, j, k, r
        for i from 0 <= i < self.nr:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            if exact:
                r = sparcity
            else:
                r = random.randrange(sparcity)
            for j from 0 <= j <= r:
                self[i, random.randrange(0,self.nc)] = \
                        random.randrange(1,self.p)

    def __repr__(self):
        cdef int i, j
        s = "[\n"
        for i from 0 <= i < self.nr:
            if PyErr_CheckSignals(): raise KeyboardInterrupt
            for j from 0 <= j < self.nc:
                s = s + str(get_entry(&self.rows[i],j)) + ", "
            s = s + "\n"
        s = s[:-3] + "\n]"
        return s

    def list(self):
        """
        Return list of entries of self.
        """
        R = self.base_ring()
        X = [R(0)]*(self.nr*self.nc)
        cdef int i
        for i from 0 <= i < self.nr:
            for j, x in to_list(&self.rows[i]):
                X[i*self.nc + j] = x
        return X

    def __getitem__(self, t):
        if not isinstance(t, tuple) or len(t) != 2:
            raise IndexError, "Index of matrix item must be a row and a column."
        cdef int i, j
        i, j = t
        return self.base_ring()(get_entry(&self.rows[i], j))   # make this quickly

    def __setitem__(self, t, x):
        if not isinstance(t, tuple) or len(t) != 2:
            raise IndexError, "Index for setting matrix item must be a row and a column."
        cdef int i, j
        i, j = t
        if i<0 or i >= self.nr or j<0 or j >= self.nc:
            raise IndexError, "Array index out of bounds."
        set_entry(&self.rows[i], j, x)

    def nrows(self):
        return self.nr

    def ncols(self):
        return self.nc

    def swap_rows(self, int n1, int n2):
        """
        Swap the rows in positions n1 and n2
        """
        if n1 < 0 or n1 >= self.nr or n2 < 0 or n2 >= self.nr:
            raise IndexError, "Invalid row number."
        if n1 == n2:
            return
        cdef c_vector_modint tmp
        tmp = self.rows[n1]
        self.rows[n1] = self.rows[n2]
        self.rows[n2] = tmp

    def echelon(self):
        """
        Replace self by its reduction to reduced row echelon form.

        ALGORITHM: We use Gauss elimination, which is slightly
        intelligent, in these sense that we clear each column using a
        row with the minimum number of nonzero entries.

        TODO: Implement switching to a dense method at some point.
        """
        cdef int i, r, c, a0, a_inverse, b, min, min_row, start_row, tenth, do_verb
        cdef c_vector_modint tmp
        start_row = 0
        self.__pivots = []
        fifth = self.nc / 10 + 1
        tm = sage.misc.all.verbose(caller_name = 'sparse_matrix_pyx matrix_modint echelon')
        do_verb = (sage.misc.all.get_verbose() >= 2)
        for c from 0 <= c < self.nc:
            if c % fifth == 0 and c>0:
                if do_verb:
                    tm = sage.misc.all.verbose('on column %s of %s'%(c, self.nc),
                             level = 2,
                             caller_name = 'matrix_modn_sparse echelon')
            min = self.nc + 1
            min_row = -1
            for r from start_row <= r < self.nr:
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
            #if min_row != -1:
            #    sage.misc.all.verbose('c=%s; using row %s with %s nonzero entries'%(c, min_row, min), level=2)

            if min_row != -1:
                r = min_row
                #print "min number of entries in a pivoting row = ", min
                self.__pivots.append(c)
                _sig_on
                # Since we can use row r to clear column c, the
                # entry in position c in row r must be the first nonzero entry.
                a = self.rows[r].entries[0]
                if a != 1:
                    a_inverse = ai.c_inverse_mod_int(a, self.p)
                    scale_c_vector_modint(&self.rows[r], a_inverse)
                self.swap_rows(r, start_row)
                for i from 0 <= i < self.nr:
                    if i != start_row:
                        b = get_entry(&self.rows[i], c)
                        if b != 0:
                            add_c_vector_modint_init(&tmp, &self.rows[i],
                                                     &self.rows[start_row], self.p - b)
                            clear_c_vector_modint(&self.rows[i])
                            self.rows[i] = tmp
                start_row = start_row + 1
                _sig_off

    def rank(self):
        """
        Return the rank found during the last echelon operation on self.
        Of course if self is changed, then the rank could be incorrect.
        """
        if self.__pivots is None:
            raise ArithmeticError, "Echelon form has not yet been computed."
        return len(self.__pivots)

    def linear_combination_of_rows(self, Vector_modint v):
        if self.nr != v.degree():
            raise ArithmeticError, "Incompatible vector * matrix multiply."
        cdef c_vector_modint w, sum, sum2
        cdef int i, r, nr, p
        cdef Vector_modint ans
        nr = self.nr
        p = self.p
        w = v.v
        init_c_vector_modint(&sum, p, nr, 0)
        _sig_on
        for i from 0 <= i < w.num_nonzero:
            r = w.positions[i]
            add_c_vector_modint_init(&sum2, &sum, &self.rows[r], w.entries[i])
            # Now sum2 is initialized and equals sum + w[i]*self.rows[i]
            # We want sum to equal this.
            clear_c_vector_modint(&sum)
            sum = sum2
        _sig_off
        # Now sum is a sparse C-vector that gives the linear combination of rows.
        # Convert to a Vector_modint and return.
        ans = Vector_modint(p, nr)
        clear_c_vector_modint(&ans.v)
        ans.v = sum
        return ans

