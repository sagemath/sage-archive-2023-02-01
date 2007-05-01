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
include '../modules/vector_modn_sparse_c.pxi'

cimport matrix
cimport matrix_sparse
from sage.rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract

from sage.misc.misc import verbose, get_verbose, graphics_filename

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
        matrix.Matrix.__init__(self, parent)

        # allocate memory
        cdef Py_ssize_t i, nr, nc
        cdef int p

        nr = parent.nrows()
        nc = parent.ncols()
        p = parent.base_ring().order()
        self.p = p


        self.rows = <c_vector_modint*> sage_malloc(nr*sizeof(c_vector_modint))
        if self.rows == NULL:
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
        cdef int s, z, p
        cdef Py_ssize_t i, j, k

        cdef void** X

        if isinstance(entries, dict):
            # Sparse input format.
            R = self._base_ring
            for ij, x in entries.iteritems():
                z = R(x)
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
                    z = R(<object>X[k])
                    if z != 0:
                        set_entry(&self.rows[i], j, z)
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
    #   * x _dict -- sparse dictionary of underlying elements (need not be a copy)
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

    def _dict(self):
        """
        Unsafe version of the dict method, mainly for internal use.
        This may return the dict of elements, but as an *unsafe*
        reference to the underlying dict of the object.  It is might
        be dangerous if you change entries of the returned dict.
        """
        d = self.fetch('dict')
        if not d is None:
            return d

        cdef Py_ssize_t i, j, k
        d = {}
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self.rows[i].num_nonzero:
                d[(int(i),int(self.rows[i].positions[j]))] = self.rows[i].entries[j]
        self.cache('dict', d)
        return d


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
                # Since we can use row r to clear column c, the
                # entry in position c in row r must be the first nonzero entry.
                a = self.rows[r].entries[0]
                if a != 1:
                    a_inverse = ai.c_inverse_mod_int(a, self.p)
                    scale_c_vector_modint(&self.rows[r], a_inverse)
                self.swap_rows_c(r, start_row)
                _sig_on
                for i from 0 <= i < self._nrows:
                    if i != start_row:
                        b = get_entry(&self.rows[i], c)
                        if b != 0:
                            add_c_vector_modint_init(&tmp, &self.rows[i],
                                                     &self.rows[start_row], self.p - b)
                            clear_c_vector_modint(&self.rows[i])
                            self.rows[i] = tmp
                _sig_off
                start_row = start_row + 1


        self.cache('pivots',pivots)
        self.cache('in_echelon_form',True)

    def visualize_structure(self, filename=None, maxsize=512):
        """
        Write a PNG image to 'filename' which visualizes self by putting
        black pixels in those positions which have nonzero entries.

        White pixels are put at positions with zero entries. If 'maxsize'
        is given, then the maximal dimension in either x or y direction is
        set to 'maxsize' depending on which is bigger. If the image is
        scaled, the darkness of the pixel reflects how many of the
        represented entries are nonzero. So if e.g. one image pixel
        actually represents a 2x2 submatrix, the dot is darker the more of
        the four values are nonzero.

        INPUT:
            filename -- either a path or None in which case a filename in
                        the current directory is chosen automatically
                        (default:None)
            maxsize -- maximal dimension in either x or y direction of the resulting
                       image (default: 512)

        """
        import gd
        import os

        cdef Py_ssize_t i, j, k
        cdef float _maxsize = float(maxsize)
        cdef float blk,invblk
        cdef int delta
        cdef int x,y,r,g,b

        mr, mc = self.nrows(), self.ncols()
        if max(mr,mc) > _maxsize:
            ir = int(mc * _maxsize/max(mr,mc))
            ic = int(mr * _maxsize/max(mr,mc))
            blk = max(mr,mc)/_maxsize
            invblk = _maxsize/max(mr,mc)
        else:
            ir = mc
            ic = mr
            blk = 1.0
            invblk = 1.0

        delta = 255.0 / blk**2

        im = gd.image((ir,ic),1)
        white = im.colorExact((255,255,255))
        im.fill((0,0),white)

        colorComponents = im.colorComponents
        getPixel = im.getPixel
        setPixel = im.setPixel
        colorExact = im.colorExact

        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self.rows[i].num_nonzero:
                x = <int>(invblk * self.rows[i].positions[j])
                y = <int>(invblk * i)
                r,g,b = colorComponents( getPixel((x,y)))
                setPixel( (x,y), colorExact((r-delta,g-delta,b-delta)) )

        if filename is None:
            filename = sage.misc.misc.graphics_filename()

        im.writePng(filename)

    def density(self):
        """
        Return the density of self.

        By density we understand the ration of the number of nonzero
        positions and the self.nrows() * self.ncols(), i.e. the number
        of possible nonzero positions.


        EXAMPLE:

            First, note that the density parameter does not ensure
            the density of a matrix, it is only an upper bound.

            sage: A = random_matrix(GF(127),200,200,density=0.3, sparse=True)
            sage: A.density() # somewhat random
            643/2500

            sage: A = matrix(QQ,3,3,[0,1,2,3,0,0,6,7,8],sparse=True)
            sage: A.density()
            2/3
        """
        from sage.rings.rational_field import QQ

        cdef Py_ssize_t i, j, k
        k = 0
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self.rows[i].num_nonzero:
                k+=1
        return QQ(k)/QQ(self.nrows()*self.ncols())
