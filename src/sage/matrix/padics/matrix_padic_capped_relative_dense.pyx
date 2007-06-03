"""
Dense matrices over either pAdicRingCappedRelative or pAdicFieldCappedRelative.

AUTHOR:
David Roe (roed@math.harvard.edu) -- 3/30/07
"""

##############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include "../../ext/stdsage.pxi"
include "../../ext/interrupt.pxi"
include "../../ext/cdefs.pxi"
include "../../ext/gmp.pxi"

ctypedef unsigned long ulong

from sage.rings.padics.misc import min

cimport sage.matrix.matrix_dense
cimport sage.rings.integer
cimport sage.matrix.matrix_integer_dense
cimport sage.structure.element
cimport sage.matrix.matrix

from sage.structure.element cimport RingElement
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.rings.integer cimport Integer
from sage.matrix.matrix cimport Matrix
from sage.matrix.matrix_dense cimport Matrix_dense

from sage.rings.padics.padic_capped_relative_element import pAdicCappedRelativeElement
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ

cdef class Matrix_padic_capped_relative_dense(Matrix_dense):

    ########################################################################
    # LEVEL 0 functionality (overridden: y = for speed purposes, c = change of functionality)
    # y * __setitem__
    # y * swap_columns_c
    # y * swap_rows_c
    # y * add_multiple_of_row_c
    # y * add_multiple_of_column_c
    # y * _coerce_element
    # y * rescale_row_c
    # y * rescale_column_c
    # y * set_row_to_multiple_of_row
    # y * linear_combination_of_rows
    # y * linear_combination_of_columns
    # y * _nonzero_positions_by_row
    # y * _nonzero_positions_by_column
    # y * _nonzero_positions_in_row
    # y * _nonzero_positions_in_column
    # y * _vector_times_matrix_c_impl
    # y * _matrix_times_vector_c_impl
    # y * __mod__
    # y * mod
    # y * _lmul_c_impl
    # y * __nonzero__
    # y * denominator
    # c * lift
    ## Lots of other functions that use get_unsafe and set_unsafe
    ## because get_unsafe has to construct a p-adic element and then
    ## set_unsafe just discards it
    ########################################################################

    ########################################################################
    # LEVEL 1 functionality
    # x * __new__
    # x * __dealloc__
    # x * __init__
    # x * set_unsafe
    # x * get_unsafe
    # x * __richcmp__    -- always the same
    # x * __hash__       -- alway simple
    ########################################################################
    def __new__(self, parent, entries, copy, coerce):
        Matrix_dense.__init__(self, parent)
        self._nrows = parent.nrows()
        self._ncols = parent.ncols()
        # Allocate an array for the precisions of the elements
        _sig_on
        self._relprecs = <ulong *>sage_malloc(size_of(ulong) * (self._nrows * self._ncols))
        _sig_off
        if self._relprecs == NULL:
            raise MemoryError, "out of matrix allocating precisions"

        # Allocate an array of pointers to the rows of _relprecs,
        # which is useful for certain algorithms.
        ##################################
        # *IMPORTANT*: FOR pAdicCappedRelative MATRICES, WE ALWAYS ASSUME THAT
        # THIS ARRAY IS *not* PERMUTED.  This should be OK, since all
        # algorithms operate on the underlying integer matrix anyway.
        ##################################
        self._relprec_mat = <ulong **> sage_malloc(sizeof(ulong*)*self._nrows)
        # Should eventually switch to unsigned long ints
        if self._relprec_mat == NULL:
            sage_free(self._relprecs)
            self._relprecs = NULL
            raise MemoryError, "out of memory allocating precisions"

        # Set each of the pointers in the array self._matrix to point
        # at the memory for the corresponding row.
        cdef Py_ssize_t i, k
        k = 0
        for i from 0 <= i < self._nrows:
            self._relprec_mat[i] = self._relprecs + k
            k = k + self._ncols

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._initialized:
            sage_free(self._relprecs) #shouldn't this part always be done, even if self is not initialized?
            sage_free(self._relprec_mat)

    def __richcmp__(Matrix self, right, int op):  # always need for mysterious reasons.
        return self._richcmp(right, op)
    def __hash__(self):
        x = self.fetch('hash')
        if not x is None: return x

        if not self._mutability._is_immutable:
            raise TypeError, "mutable matrices are unhashable"

        v = self._relprecs
        cdef Py_ssize_t i
        cdef long h
        h = 0
        cdef PyObject** w
        w = FAST_SEQ_UNSAFE(v)
        for i from 0 <= i < len(v):
            h = h ^ (i * PyObject_Hash( <object> w[i] ))
        h = h ^ self._value_matrix._hash()
        self.cache('hash', h)
        return h

    cdef _comp_valaddeds(self):
        cdef Py_ssize_t i, j
        cdef int n
        cdef Integer p
        n = 0
        p = self._base_ring.prime()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self._valaddeds[n] = self._value_matrix.get_unsafe(i, j)._valuation(p)
                n += 1


    def __init__(self, parent, entries, copy, coerce):
        cdef Py_ssize_t i, j, Icur
        cdef int is_list
        cdef object coerced_list

        self._padic_values = {}
        self._valaddeds = {}

        if entries is None:
            x = self.base_ring()(0)
            is_list = 0
        elif isinstance(entries, (int,long)) or is_Element(entries):
            try:
                x = self.base_ring()(entries)
            except TypeError:
                self._initialized = False
                raise TypeError, "unable to coerce entry to a p-adic"
            is_list = 0
        elif isinstance(entries, tuple) and len(entries) == 4 and entries[3] is None: # The None entry is a flag to use this constructor
            # Format is (int_matrix, valbase, prec_info, valaddeds,
            #First we fill the _value_matrix with values from entries[0]
            self._value_matrix = Matrix_integer_dense(MatrixSpace(ZZ, self._nrows, self._ncols), entries[0], copy=copy, coerce = coerce)
            #Next we set the common valuation.
            self._valbase = entries[2]
            #Now we fill in precision information
            if entries[1] is None:
                #Use precision_cap of base_ring
                self._comp_valaddeds()
                self._adjust_prec_info_global(infinity, self._base_ring.precision_cap())
                self._initialized = True
                return
            elif isinstance(entries[1], (int, long)) or is_Element(entries[1]):
                try:
                    a = ZZ(entries[1])
                except TypeError:
                    self._initialized = False
                    raise TypeError, "unable to create an integer out of desired precision"
                self._adjust_prec_info_global(a, self._base_ring.precision_cap())
            elif isinstance(entries[1], tuple) and len(entries[1]) == 2:
                alist = False
                blist = False
                if entries[1][0] is None or entries[1][0] is infinity:
                    a = infinity
                elif isinstance(entries[1][0], (list, tuple)) and len(entries[1][0]) == self._nrows * self._ncols:
                    a = entries[1][0]
                    alist = True
                else:
                    try:
                        a = ZZ(entries[1][0])
                    except TypeError:
                        self._initialized = False
                        raise TypeError, "unable to create an integer out of desired precision"
                if entries[1][1] is None or entries[1][1] is infinity:
                    b = infinity
                elif isinstance(entries[1][1], (list, tuple)) and len(entries[1][1]) == self._nrows * self._ncols:
                    a = entries[1][1]
                    blist = True
                else:
                    try:
                        b = ZZ(entries[1][1])
                    except TypeError:
                        self._initialized = False
                        raise TypeError, "unable to create an integer out of desired precision"
                if alist:
                    if blist:
                        self._adjust_prec_info_local_local(a, b, coerce)
                    else:
                        self._adjust_prec_info_local_global(a, min(b, self._base_ring.precision_cap()), coerce)
                else:
                    if blist:
                        self._adjust_prec_info_global_local(a, b, coerce)
                    else:
                        self._adjust_prec_info_global(a, min(b, self._base_ring.precision_cap()))
                return
            elif isinstance(entries[1], list) and len(entries[1]) == self._nrows*self._ncols:
                self._adjust_prec_info_local(entries[1], coerce)
                return
            else:
                raise TypeError, "Unrecognizable precision information in creating p-adic matrix"
        elif not isinstance(entries, (list, tuple)):
            entries = list(entries) # does this actually handle dicts correctly?
            is_list = 1
        else:
            is_list = 1

        if is_list:

            # Create the matrix whose entries are in the given entry list.
            if len(entries) != self._nrows * self._ncols:
                sage_free(self._relprecs)
                sage_free(self._relprec_mat)
                self._relprecs = NULL
                # Do I also need to set self._relprec_mat to NULL?
                raise TypeError, "entries has the wrong length"
            if coerce:
                # How do I declare a list in Pyrex?
                coerced_list = [self.base_ring()(e) for e in entries]
                self._valbase = min([e.valuation() for e in coerced_list]) #should eventually not make this a list before passing to min
                self._value_matrix = Matrix_integer_dense(MatrixSpace(ZZ, self._nrows, self._ncols), 0, copy = False, coerce = True)
                Icur = 0
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < self._ncols:
                        self._value_matrix.set_unsafe(i, j, coerced_list[Icur + j].__rshift__(self._valbase)._integer_())
                        self._relprec_mat[i][j] = e.precision_absolute() - self._valbase
                    Icur += self._nrows
                self._initialized = True
            else:
                Icur = 0
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < self._ncols:
                        self._value_matrix.set_unsafe(i, j, entries[Icur + j].__rshift__(self._valbase)._integer_())
                        self._relprec_mat[i][j] = e.precision_absolute() - self._valbase
                    Icur += self._nrows
                self._initialized = True
        else:

            # If x is zero, make the zero matrix and be done.
            if x == 0:
                self._zero_out_matrix()
                self._initialized = True #?
                return

            # the matrix must be square:
            if self._nrows != self._ncols:
                sage_free(self._relprecs)
                sage_free(self._relprec_mat)
                self._relprecs = NULL
                raise TypeError, "nonzero scalar matrix must be square"

            # Now we set all the diagonal entries to x and all other entries to 0.
            self._zero_out_matrix(x.valuation(), precs = False)
            xint = x.unit_part()
            xrprec = x.precision_relative()
            for i from 0 <= i < self._nrows:
                self._value_matrix.set_unsafe(i, i, xint)
                self._relprecs_mat[i][i] = xrprec
            self._initialized = True

    cdef void _comp_valaddeds(self):
        cdef Py_ssize_t i, j
        cdef int n
        cdef Integer p
        n = 0
        p = self._base_ring.prime()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                self._valaddeds[n] = self._value_matrix.get_unsafe(i, j)._valuation(p)
                n += 1

    cdef void _adjust_prec_info_global(self, RingElement absolute, RingElement relative):
        pass

    cdef void _adjust_prec_info_global_local(self, RingElement absolute, object relative, coerce):
        pass

    cdef void _adjust_prec_info_local_global(self, object absolute, RingElement relative, coerce):
        pass

    cdef void _adjust_prec_info_local_local(self, object absolute, object relative, coerce):
        pass

    def _zero_out_matrix(self, val = None, precs = True):
        self._value_matrix._zero_out_matrix() #this may have been done already...
        if val is None:
            self._valbase = Integer(0)
        else:
            self._valbase = val
        if precs:
            self._relprecs = [infinity] * (self._nrows * self._ncols)
            self._valaddeds = [infinity] * (self._nrows * self._ncols)
        self._padic_values = {}

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object value):
        r"""
        Sets self[i][j] to value.  No error checking is done.

        INPUT:
            i -- the row
            j -- the column
            value -- either (1 is faster)
               1) a tuple (a, c) where the desired value is p^(self._valbase) * a
                  and c + self._valbase is the desired absolute precision
               2) an element of the base ring.

        NOTE:

        Using the second option for value can result in having to
        rescale the whole matrix to have a lower valuation.  For
        example, if one creates a matrix filled with
        pAdicFieldCappedRelativeElements that all have positive
        valuation, then set_unsafe one of the entries to have negative
        valuation, this will be slow.  One should particularly avoid
        doing this multiple times.

        You have three options.  Either
          1) set the matrix to originally have base valuation that is the minimum you need.
          2) call self.set_base_valuation(min_valuation) appropriately before doing a lot of set_unsafeing.
          3) set the element with minimum valuation first.
        """
        cdef Integer a
        cdef Integer b, c # should eventually be an unsigned long int
        cdef Integer ppow
        cdef Py_ssize_t i, j
        if PyTupleCheck(value):
            a = <Integer> PyTuple_GET_ITEM(value, 0) # integer value above self._valbase
            c = <Integer> PyTuple_GET_ITEM(value, 1) # precision above self._valbase
            self._value_matrix.set_unsafe(i, j, a)
            self._value_matrix.reduce_entry_unsafe(i, j, self._base_ring.prime_pow(c))
            self._relprec_mat[i][j] =  mpz_get_ui(c.value)
        else:
            #d = self._base_ring(value)
            a = <Integer> d._unit_part().lift()
            b = <Integer> d.valuation() - self._valbase
            c = <Integer> d.precision_absolute() - self._valbase
            if b < 0:
                # We are trying to set an entry to a value that would make
                # the corresponding entry of self._value_matrix not be an
                # integer.  We must therefore rescale every entry of the
                # matrix.  To avoid doing this multiple times, call
                # self.set_base_valuation(min_valuation) with the minimum
                # valuation of any element before setting lots of entries

                # We rescale c back up, because this entry will have
                # _valadded 0, and we're decreasing _valbase by -b.  Note that since c >= b, now c >= 0.
                c -= b
                self._valbase += b
                # We now rescale every entry.  This should be made more efficient
                ppow = <Integer> self._base_ring.prime_pow(-b)
                self._value_matrix = self._value_matrx * ppow
                for i from 0 <= i < self._nrows * self._ncols:
                    self._relprecs[i] -= b

                self._relprec_mat[i][j] = c
                self._value_matrix.set_unsafe(i, j, a)
            elif b == 0:
                self._relprec_mat[i][j] = c
                self._value_matrix.set_unsafe(i, j, a)
            else: # b > 0
                self._relprec_mat[i][j] = c
                self._value_matrix.set_unsafe(i, j, a * self._base_ring.prime_pow(b))

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef int n
        cdef Integer m, relprec
        n = i * self._ncols + j
        if self._padic_values.has_key(n):
            return self._padic_values[n]
        if self._cache.has_key('list'):
            return self.fetch('list')[n]
        if not self._valaddeds.has_key(n):
            m = self._value_matrix.get_unsafe(i, j)
            self._valaddeds[n] = m.valuation()
        relprec = self._relprecs[n] - self._valaddeds[n]
        return pAdicCappedRelativeElement(self._base_ring, (self._valbase + self._valaddeds[n], Mod(self._value_matrix.get_unsafe(i, j) // self._base_ring.prime_pow(self._valaddeds[n]), self._base_ring.prime_pow(relprec)), relprec), construct = True)

    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    #   * cdef _add_c_impl
    #   * cdef _sub_c_impl
    #   * cdef _mul_c_impl
    #   * cdef _cmp_c_impl
    #   * __neg__
    #   * __invert__
    #   * __copy__
    #   * _multiply_classical
    #   * _multiply_split_prec
    #   * _matrix_times_matrix_c_impl
    #   * _list -- list of underlying elements (need not be a copy)
    # y * _dict -- sparse dictionary of underlying elements (need not be a copy)
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
    #   * cdef _sub_c_impl
    #   * __deepcopy__
    #   * __invert__
    #   * _determinant_split_prec
    #   * _charpoly_split_prec
    #   * kernel
    #   * echelonize
    #   * _adjoint
    ########################################################################
