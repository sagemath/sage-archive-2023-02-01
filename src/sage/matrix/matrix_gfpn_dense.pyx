r"""
Dense Matrices over `\mathbb F_q`, with `q<255` odd and not prime

This module is a wrapper for version 2.4.24 of the Aachen
`C-MeatAxe <http://www.math.rwth-aachen.de/homes/MTX/download.html>`_,
improved by an implementation of the Winograd-Strassen multiplication
algorithm. It provides matrices over the finite field `\mathbb F_q`,
where `q\le 255` is odd and not prime.

AUTHORS:

- Simon King (2015-09-18): initial version

"""

#*****************************************************************************
#       Copyright (C) 2015 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


## Define an environment variable that enables MeatAxe to find
## its multiplication tables.

from sage.env import DOT_SAGE
import os
cdef extern from "Python.h":
    object PyString_FromStringAndSize(char *s, Py_ssize_t len)
    char* PyString_AsString(object string)
MtxLibDir = PyString_AsString(os.path.join(DOT_SAGE,'meataxe'))

####################
#
# import sage types
#
####################

from sage.rings.integer import Integer
from sage.rings.finite_rings.constructor import GF
from sage.rings.finite_rings.integer_mod import IntegerMod_int
from sage.matrix.constructor import random_matrix
from sage.rings.arith import is_prime_power, factor
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.randstate import current_randstate
from sage.misc.cachefunc import cached_method, cached_function
from sage.structure.element cimport Element, ModuleElement, RingElement, Matrix

include 'sage/ext/stdsage.pxi'

####################
#
# auxiliary functions
#
####################
import sys
from libc.string cimport memcpy

cdef inline int setfield(long n) except -1:
    # This is a wrapper around FfSetField, but
    # we guard it against MTX_Error, which would immediately
    # crash the Sage session.
    if n == FfOrder:
        return 0
    if not (0 < n < 255 and is_prime_power(n)):
        raise ValueError("Only finite fields of order at most 255 are supported")
    return FfSetField(n)

# Fast conversion from field to int and int to field
cdef class FieldConverter_class:
    """
    An auxiliary class, used to convert between <int> and finite field element

    This class is for non-prime fields only. The method
    :meth:`int_to_field` exists for speed. The method
    :meth:`field_to_int` exists in order to have a common interface
    for elements of prime and non-prime fields; see
    :class:`PrimeFieldConverter_class`.

    EXAMPLE::

        sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class
        sage: F.<y> = GF(125)
        sage: C = FieldConverter_class(F)
        sage: C.int_to_field(15)
        3*y
        sage: F.fetch_int(15)
        3*y
        sage: %timeit C.int_to_field(15)    #not tested
        625 loops, best of 3: 1.04 µs per loop
        sage: %timeit F.fetch_int(15)       #not tested
        625 loops, best of 3: 3.97 µs per loop
        sage: C.field_to_int(y)
        5
        sage: y.integer_representation()
        5

    """
    def __init__(self, field):
        """
        INPUT:

        A finite *non-prime* field. This assumption is not tested.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class
            sage: F.<y> = GF(125)
            sage: C = FieldConverter_class(F)
            sage: C.int_to_field(15)
            3*y
            sage: F.fetch_int(15)
            3*y
            sage: C.field_to_int(y)
            5
            sage: y.integer_representation()
            5

        """
        self.field = field._cache.fetch_int
    cdef object int_to_field(self, int x):
        """
        Fetch a python int into the field.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class
            sage: F.<y> = GF(125)
            sage: C = FieldConverter_class(F)
            sage: C.int_to_field(15)
            3*y
            sage: F.fetch_int(15)
            3*y

        """
        return self.field(x)
    cdef int field_to_int(self, x):
        """
        Represent a field element by a python int.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import FieldConverter_class
            sage: F.<y> = GF(125)
            sage: C = FieldConverter_class(F)
            sage: C.field_to_int(y)
            5
            sage: y.integer_representation()
            5

        """
        return x.integer_representation()

cdef class PrimeFieldConverter_class(FieldConverter_class):
    """
    An auxiliary class, used to convert between <int> and finite field element

    This class is for prime fields only. The methods
    :meth:`int_to_field` and :meth:`field_to_int` exist in order to
    have a common interface for elements of prime and non-prime fields;
    see :class:`FieldConverter_class`.

    EXAMPLE::

        sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class
        sage: F = GF(5)
        sage: C = PrimeFieldConverter_class(F)
        sage: C.int_to_field(int(2))
        2
        sage: F(2)
        2
        sage: C.field_to_int(F(2))
        2
        sage: int(F(2))
        2

    """
    def __init__(self, field):
        """
        INPUT:

        A finite *prime* field. This assumption is not tested.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class
            sage: F = GF(5)
            sage: C = PrimeFieldConverter_class(F)
            sage: C.int_to_field(int(2))
            2
            sage: F(2)
            2
            sage: C.field_to_int(F(2))
            2
            sage: int(F(2))
            2

        """
        self.field = field
    cdef object int_to_field(self, int x):
        """
        Fetch a python int into the field.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class
            sage: F = GF(5)
            sage: C = PrimeFieldConverter_class(F)
            sage: C.int_to_field(int(2))
            2
            sage: F(2)
            2

        """
        return IntegerMod_int(self.field, x)
    cdef int field_to_int(self, x):
        """
        Represent a field element by a python int.

        EXAMPLE::

            sage: from sage.matrix.matrix_gfpn_dense import PrimeFieldConverter_class
            sage: F = GF(5)
            sage: C = PrimeFieldConverter_class(F)
            sage: C.field_to_int(F(2))
            2
            sage: int(F(2))
            2

        """
        return int(x)

cdef dict _converter_cache = {}
cdef FieldConverter_class FieldConverter(field):
    """
    Return a :class:`FieldConverter_class` or :class:`PrimeFieldConverter_class` instance,
    depending whether the field is prime or not.

    EXAMPLE::

        sage: MS = MatrixSpace(GF(5^3,'y'),2)
        sage: A = MS.random_element()
        sage: A*2 == A+A    # indirect doctest
        True
        sage: A = MS.random_element()
        sage: A*2 == A+A
        True

    """
    try:
        return _converter_cache[field]
    except KeyError:
        if field.is_prime_field():
            return _converter_cache.setdefault(field, PrimeFieldConverter_class(field))
        return _converter_cache.setdefault(field, FieldConverter_class(field))



cdef class Matrix_gfpn_dense(Matrix_dense):
    r"""
    Dense matrices over `\mathbb F_q`, `q<255` odd and not prime.

    NOTE:

    This class uses a major modification of the Aachen C-MeatAxe
    as backend. In principle, it would also work for prime fields
    and in characteristic two. However, other matrices in Sage,
    relying on linbox, m4ri or m4rie, are more efficient in these
    cases.

    EXAMPLES::

        sage: M = MatrixSpace(GF(25,'z'),2,3)([1,2,3,4,5,6])
        sage: print M
        [1 2 3]
        [4 0 1]
        sage: type(M)
        <type 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>

    The documentation of the ``__init__`` methods shows further
    ways of creating a :class:`Matrix_gfpn_dense` instance.
    However, these should only be of internal use.

    """
##################
## Init, Dealloc, Copy
    def __cinit__(self, parent=None, entries=None, *args, **kwds):
        """
        TESTS::

            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense
            sage: Matrix_gfpn_dense.__new__(Matrix_gfpn_dense)   # indirect doctest
            []
            sage: Matrix_gfpn_dense(MatrixSpace(GF(64,'z'),4), None)
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

        """
        if parent is None:  # this makes Matrix_gfpn_dense.__new__(Matrix_gfpn_dense) work,
                            # returning a non-initialised matrix
            return
        if isinstance(parent, basestring): # this allows to provide a file when initialising a matrix
            return
        cdef int f = parent.base_ring().order()
        cdef int nrows = parent.nrows()
        cdef int ncols = parent.ncols()
        self.Data = MatAlloc(f, nrows, ncols)

    def __init__(self, parent, data=None, mutable=True, copy=False, coerce=False):
        """
        Matrix extension class using libmeataxe as backend

        INPUT:

        Instances of this class can be created by providing one of
        the following input data, where ``q<255`` is a prime power,
        ``m,n`` are non-negative integers, and `a_{11},...,a_{mn}`
        can be coerced into ``GF(q)``. Note that a user should
        create these instances via the matrix constructors; what
        we explain here is for internal use only!

        - None => empty matrix over an unspecified field (used for unpickling)
        - a string ``f`` ==> load matrix from the file named ``f``
        - A matrix space of `m\\times n` matrices over GF(q) and either

          - a list `[a_{11},a_{12},...,a_{1n},a_{21},...,a_{m1},...,a_{mn}]`,
            which results in a matrix with the given marks
          - ``None``, which is the fastest way to creata a zero matrix.
          - an element of GF(q), which results in a diagonal matrix with the
            given element on the diagonal.

        If the optional parameter ``mutable`` is ``False`` (by default,
        it is ``True``), the resulting matrix can not be changed, and
        it can be used as key in a Python dictionary.

        The arguments ``copy`` and ``coerce`` are ignored, they are only
        here for a common interface with :class:`~sage.matrix.matrix.Matrix`.

        EXAMPLES::

            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense

        1. Creating an empty matrix::

            sage: Matrix_gfpn_dense(None)
            []

        2. Creating a zero (3x2)-matrix::

            sage: Matrix_gfpn_dense(MatrixSpace(GF(4,'z'),3,2))
            [0 0]
            [0 0]
            [0 0]

        3. Creating a matrix from a list or list of lists::

            sage: Matrix_gfpn_dense(MatrixSpace(GF(5),2,3),[1,2,3,4,5,6])
            [1 2 3]
            [4 0 1]
            sage: Matrix_gfpn_dense(MatrixSpace(GF(5),2,3),[[1,2,3],[4,5,6]])  # indirect doctest
            [1 2 3]
            [4 0 1]

        4. Creating a diagonal matrix::

            sage: M = Matrix_gfpn_dense(MatrixSpace(GF(7),5),2); M
            [2 0 0 0 0]
            [0 2 0 0 0]
            [0 0 2 0 0]
            [0 0 0 2 0]
            [0 0 0 0 2]

        5. Creating a matrix from a file in MeatAxe format.

           First, we have to create that file; we use a temporary file,
           that will be removed when leaving Sage. Note that the method
           :meth:`msave` must be used, which does not use Python pickling
           but relies on the intrinsic C--MeatAxe way of saving.
           ::

            sage: f = tmp_filename()
            sage: M.msave(f)
            sage: Matrix_gfpn_dense(f)
            [2 0 0 0 0]
            [0 2 0 0 0]
            [0 0 2 0 0]
            [0 0 0 2 0]
            [0 0 0 0 2]

        TESTS::

            sage: MS = MatrixSpace(GF(125,'y'),2)
            sage: A = MS(0)
            sage: A.left_kernel()
            Vector space of degree 2 and dimension 2 over Finite Field in y of size 5^3
            Basis matrix:
            [1 0]
            [0 1]
            sage: A.right_kernel()
            Vector space of degree 2 and dimension 2 over Finite Field in y of size 5^3
            Basis matrix:
            [1 0]
            [0 1]

        """
        if parent is None:
            self._is_immutable = False
            self._ncols = 0
            self._nrows = 0
            self._cache = {}
            return
        if isinstance(parent, basestring): # load from file
            FILE = os.path.realpath(parent)
            try:
                fsock = open(FILE,"rb",0)
                fsock.close()
            except (OSError,IOError):
                return
            self.Data = MatLoad(FILE)
            if FfSetField(self.Data.Field):
                raise ValueError("Invalid data in file {}".format(FILE))
            B = GF(self.Data.Field, 'z')
            parent = MatrixSpace(B, self.Data.Nor, self.Data.Noc)
            self._is_immutable = False
            self._parent = parent
            self._base_ring = B
            self._converter = FieldConverter(B)
            self._ncols = self.Data.Noc
            self._nrows = self.Data.Nor
            self._cache = {}
            return

        if not self.Data: # should have been initialised by __cinit__
            raise MemoryError, "Error allocating memory for MeatAxe matrix"
        Matrix_dense.__init__(self, parent)
        self._is_immutable = not mutable
        B = self._base_ring
        self._converter = FieldConverter(B)
        if data is None:
            return

        cdef int i,j
        cdef FEL f
        cdef PTR x
        if not isinstance(data,list):
            if not data:
                return
            if self._nrows != self._ncols:
                raise ValueError("Cannot initialise non-square matrix from {}".format(data))
            f = FfFromInt(self._converter.field_to_int(self._coerce_element(data)))
            x = self.Data.Data
            for j from 0 <= j < self.Data.Noc:
                FfInsert(x,j,f)
                FfStepPtr(&x)
            return

        x = self.Data.Data
        cdef int nr = self.Data.Nor
        cdef int nc = self.Data.Noc
        assert self._ncols == nc
        assert self._nrows == nr
        if nr==0 or nc==0:
            return
        if len(data)<nr:
            raise ValueError, "Expected a list of size at least the number of rows"
        cdef list dt, dt_i
        FfSetField(self.Data.Field)
        FfSetNoc(nc)
        if isinstance(data[0],list):
            # The matrix is given by a list of rows
            dt = data
            for i from 0 <= i < nr:
                idx = 0
                dt_i = dt[i]
                for j from 0 <= j < nc:
                    FfInsert(x, j, FfFromInt(self._converter.field_to_int(self._coerce_element(dt_i[j]))))
                FfStepPtr(&(x))
        else:
            # It is supposed to be a flat list of all entries, sorted by rows
            dtnext = data.__iter__().next
            for i from 0 <= i < nr:
                for j from 0 <= j < nc:
                    bla = self._converter.field_to_int(self._coerce_element(dtnext()))
                    FfInsert(x, j, FfFromInt(bla))
                FfStepPtr(&(x))

    def rowsize(self):
        return self.Data.RowSize

    def __dealloc__(self):
        if self.Data != NULL:
            MatFree(self.Data)
            self.Data = NULL

    def __copy__(self):
        """
        Return a copy of this matrix.

        EXAMPLES::

            sage: M=MatrixSpace(GF(25,'x')([20*[0],20*[0],[1]+19*[0]])
            sage: N=copy(M)
            sage: print N
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            sage: N==M
            True
            sage: N is M
            False
            sage: from sage.matrix.matrix_gfpn_dense import Matrix_gfpn_dense
            sage: M=Matrix_gfpn_dense('')
            sage: N=copy(M)
            sage: N
            Empty MTX matrix
            sage: N==M
            True
            sage: N is M
            False
        """
        cdef Matrix_gfpn_dense retval = type(self).__new__(type(self))
        # Do the initialisation "manually"
        retval._is_immutable = False  # a copy of a matrix is mutable!
        retval._parent = self._parent
        retval._base_ring = self._base_ring
        retval._converter = self._converter
        retval._ncols = self._ncols
        retval._nrows = self._nrows
        retval._cache = dict(self._cache.iteritems()) if self._cache is not None else {}
        if self.Data:
            retval.Data = MatDup(self.Data)
            if not retval.Data:
                raise MemoryError, "Error copying a %s instance"%repr(type(self))
        else:
            retval.Data = NULL
        return retval

    ##########################
    ## Saving should be done via pickling
    ## However, we keep a method that relies on MeatAxe matsave:
    def msave(self,f):
        """
        M.msave('filename') ==> save matrix into file <filename>

        It can be reloaded with ``Matrix_gfpn_dense('filename')``.
        """
        MatSave(self.Data,f)

    ## Pickling and string representation is taken care of by implementing get_unsafe
    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get an element without checking.

        TEST::

            sage: F.<z> = GF(9)
            sage: M = MatrixSpace(F,3)(list(F))
            sage: type(M)
            <type 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
            sage: M    # indirect doctest
            [      0     2*z   z + 1]
            [  z + 2       2       z]
            [2*z + 2 2*z + 1       1]

        """
        if self.Data == NULL:
            raise IndexError, "Matrix is empty"
        return self._converter.int_to_field(FfToInt(FfExtract(MatGetPtr(self.Data,i), j)))

    cdef inline int get_unsafe_int(self, Py_ssize_t i, Py_ssize_t j):
        # NOTE:
        # It is essential that you call FfSetField and FfSetNoc YOURSELF
        # and that you assert that the matrix is not empty!
        # This method is here for speed!
        return FfToInt(FfExtract(FfGetPtr(self.Data.Data,i) ,j))

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        # ASSUMPTION: value's parent is the base ring
        if self.Data == NULL:
            raise IndexError, "Matrix is empty"
        FfInsert(MatGetPtr(self.Data,i), j, FfFromInt(self._converter.field_to_int(value)))

    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value):
        # NOTE:
        # It is essential that you call FfSetField and FfSetNoc YOURSELF
        # and that you assert that the matrix is not empty!
        # This method is here for speed!
        FfInsert(FfGetPtr(self.Data.Data,i), j, FfFromInt(value))

    def randomize(self, density=None, nonzero=False, *args, **kwds):
        """
        Fill the matrix with random values.

        INPUT:

        - ``density`` (optional real number between zero and one) --
          the expected density of the resulting matrix
        - ``nonzero`` (optional bool, default ``False``) --
          If true, all inserted marks are non-zero.

        EXAMPLE::

            sage: MS = MatrixSpace(GF(27,'z'),6,6)
            sage: M = MS.random_element(); M    # indirect doctest
            [              1           z + 1     z^2 + z + 1             z^2       2*z^2 + z           z + 1]
            [2*z^2 + 2*z + 2   2*z^2 + z + 2         z^2 + 1 2*z^2 + 2*z + 2         z^2 + z   2*z^2 + z + 1]
            [        2*z + 2     z^2 + z + 2           z + 2 2*z^2 + 2*z + 2           2*z^2           2*z^2]
            [  2*z^2 + z + 2             z^2           z + 2         z^2 + z       2*z^2 + 2         z^2 + 2]
            [      2*z^2 + z             2*z 2*z^2 + 2*z + 1       2*z^2 + 1 2*z^2 + 2*z + 1       2*z^2 + z]
            [        2*z + 1         z^2 + z             z^2             z^2     2*z^2 + 2*z           z + 1]
            sage: type(M)
            <type 'sage.matrix.matrix_gfpn_dense.Matrix_gfpn_dense'>
            sage: MS.random_element(nonzero=True)
            [            2*z               1   z^2 + 2*z + 1   2*z^2 + z + 1             z^2     z^2 + z + 1]
            [    2*z^2 + 2*z   2*z^2 + z + 2         2*z + 1       z^2 + 2*z     2*z^2 + 2*z             z^2]
            [        z^2 + z     z^2 + z + 2 2*z^2 + 2*z + 1         z^2 + 2               1           2*z^2]
            [              z     2*z^2 + 2*z           2*z^2         2*z + 1           z + 2           z + 2]
            [        z^2 + z             z^2           z + 2     2*z^2 + 2*z         2*z + 1         z^2 + z]
            [    z^2 + z + 2       2*z^2 + z             z^2           z + 1     2*z^2 + 2*z   z^2 + 2*z + 1]
            sage: MS.random_element(density=0.5)
            [        z^2 + 2               0   z^2 + 2*z + 2       2*z^2 + z               0     z^2 + z + 2]
            [              0               1               0               0               0               0]
            [  2*z^2 + z + 1   2*z^2 + z + 2               0     z^2 + z + 2               0     z^2 + z + 1]
            [              0               0               0               0               0               0]
            [2*z^2 + 2*z + 2               0               0   2*z^2 + z + 2               0         2*z + 1]
            [              0       2*z^2 + z               0               1               0   2*z^2 + z + 1]

        """
        self.check_mutability()
        cdef int fl = self.Data.Field
        density = float(density)
        if density <= 0:
            return
        if density > 1:
            density = float(1)

        self.clear_cache()

        cdef PTR x
        cdef unsigned char *y
        x = self.Data.Data
        cdef int nr = self.Data.Nor
        cdef int nc = self.Data.Noc
        cdef int i, j, k

        FfSetField(fl)
        FfSetNoc(nc)
        cdef int O, MPB, tmp
        randint = current_randstate().c_random
        randdouble = current_randstate().c_rand_double

        if not nonzero:
            if density == 1:
                MPB = 0
                tmp = fl
                while tmp <= 256:
                    MPB += 1
                    tmp *= fl
                O = (fl**MPB)
                sig_on()
                if nc%MPB:
                    for i from 0 <= i < nr:
                        y = <unsigned char*>x
                        for j from 0 <= j < FfCurrentRowSizeIo-1:
                            y[j] = randint()%O
                        y[FfCurrentRowSizeIo-1] = randint()%(fl**(nc%MPB))
                        FfStepPtr(&(x))
                else:
                    for i from 0 <= i < nr:
                        y = <unsigned char*>x
                        for j from 0 <= j < FfCurrentRowSizeIo:
                            y[j] = randint()%O
                        FfStepPtr(&(x))
                sig_off()
            else:
                sig_on()
                for i from 0 <= i < nr:
                    for j from 0 <= j < nc:
                        if randdouble() < density:
                            FfInsert(x, j, FfFromInt( (randint()%fl) ))
                    FfStepPtr(&(x))
                sig_off()
        else:
            if density == 1:
                fl -= 1
                sig_on()
                for i from 0 <= i < nr:
                    for j from 0 <= j < nc:
                        FfInsert(x, j, FfFromInt( (randint()%fl)+1 ))
                    FfStepPtr(&(x))
                sig_off()
            else:
                fl -= 1
                sig_on()
                for i from 0 <= i < nr:
                    for j from 0 <= j < nc:
                        if randdouble() < density:
                            FfInsert(x, j, FfFromInt( (randint()%fl)+1 ))
                    FfStepPtr(&(x))
                sig_off()

    def show_contents(self, r=None):
        FfSetField(self.Data.Field)
        FfSetNoc(self.Data.Noc)
        cdef PTR p
        cdef size_t i, j
        if r is not None:
            r_min = r
            r_max = r+1
        else:
            r_min = 0
            r_max = self.Data.Nor
        for i in range(r_min, r_max):
            p = FfGetPtr(self.Data.Data, i)
            for j from 0<=j<self.Data.RowSize:
                print "%3.3d"%p[j],
            print

##################
## comparison
    cpdef int _cmp_(left, Element right) except -2:
        """
        Compare two Matrix_gfpn_dense matrices

        Of course, '<' and '>' doesn't make much sense for matrices.

        EXAMPLES::

            sage: M = MatrixSpace(GF(125,'x'),[20*[0],20*[0],[1]+19*[0]])
            sage: N = copy(M)
            sage: M == N
            True
            sage: M != N
            False
            sage: print M < N
            None
            sage: N[2,19] = 1
            sage: M == N
            False
            sage: M != N
            True
        """
        cdef Matrix_gfpn_dense self = left
        cdef Matrix_gfpn_dense N = right
        cdef char* d1
        cdef char* d2
        if self.Data.Field != N.Data.Field:
            if self.Data.Field > N.Data.Field:
                return 1
            return -1
        if self.Data.Noc != N.Data.Noc:
            if self.Data.Noc > N.Data.Noc:
                return 1
            return -1
        if self.Data.Nor != N.Data.Nor:
            if self.Data.Nor > N.Data.Nor:
                return 1
            return -1
        d1 = <char*>(self.Data.Data)
        d2 = <char*>(N.Data.Data)
        cdef str s1 = PyString_FromStringAndSize(d1,self.Data.RowSize * self.Data.Nor)
        cdef str s2 = PyString_FromStringAndSize(d2,N.Data.RowSize * N.Data.Nor)
        if s1 != s2:
            if s1 > s2:
                return 1
            return -1
        return 0

    def _rowlist_(self, i, j=-1):
        "M._rowlist_(i): Return row <i> as a list of python ints"
        cdef int k
        if self.Data:
            FfSetField(self.Data.Field)
        else:
            raise ValueError("Matrix is empty")
        if (i<0) or (i>=self.Data.Nor):
            raise IndexError("Index {} out of range 0..{}",format(i,self.Data.Nor-1))
        cdef PTR p
        p = MatGetPtr(self.Data,i)
        L = [FfToInt(FfExtract(p,k)) for k in range(self.Data.Noc)]
        if j!=-1:
            if not(isinstance(j,int) or isinstance(j,Integer)):
                raise TypeError, "Second index must be an integer"
            if j >= self.Data.Nor:
                raise IndexError, "Index out of range"
            for k from i < k <= j:
                FfStepPtr(&(p)) # This is only called after MatGetPtr, hence, after FfSetNoc.
                L.extend([FfToInt(FfExtract(p,l)) for l in range(self.Data.Noc)])
        return L

    def _matlist_(self):
        "M._matlist_(): Return M as a list of lists of python ints"
        cdef int i
        if self.Data:
            FfSetField(self.Data.Field)
            FfSetNoc(self.Data.Noc)
        else:
            raise IndexError, "Matrix is empty"
        cdef PTR p
        p = self.Data.Data
        l_out=[]
        for i from 1<=i<self.Data.Nor:
            l_out.append([FfToInt(FfExtract(p,j)) for j in range(self.Data.Noc)])
            FfStepPtr(&(p))
        l_out.append([FfToInt(FfExtract(p,j)) for j in range(self.Data.Noc)])
        return l_out

#########################
## Arithmetics
    cdef rescale_row_c(self, Py_ssize_t i, s, Py_ssize_t start_col):
        if start_col != 0 or self.Data == NULL:
            raise ValueError
        cdef PTR = MatGetPtr(self.Data, i)
        FfMulRow(PTR, FfFromInt(self._converter.field_to_int(s)))

    cpdef ModuleElement _add_(self, ModuleElement right):
        "add two MTX matrices of equal size"
        cdef Matrix_gfpn_dense Self = self
        cdef Matrix_gfpn_dense Right = right
        assert Self is not None
        assert Right is not None
        if Self.Data == NULL or Right.Data == NULL:
            raise NotImplementedError, "The matrices must not be empty"
        cdef Matrix_gfpn_dense Left = Self.__copy__()
        if MatAdd(Left.Data, Right.Data) != NULL:
            return Left
        else:
            raise ArithmeticError("Matrix sizes or fields not compatible")

    cpdef ModuleElement _sub_(self, ModuleElement right):
        "subtract two MTX matrices of equal size"
        cdef Matrix_gfpn_dense Self = self
        cdef Matrix_gfpn_dense Right = right
        assert Self is not None
        assert Right is not None
        if Self.Data == NULL or Right.Data == NULL:
            raise NotImplementedError, "The matrices must not be empty"
        cdef Matrix_gfpn_dense Left = Self.__copy__()
        Left._is_immutable = False
        if MatAddMul(Left.Data, Right.Data, mtx_taddinv[1]) != NULL:
            return Left
        else:
            raise ArithmeticError, "Matrix sizes or fields not compatible"

    def __neg__(self):
        "return negation of a MTX matrix: -M == M.__neg__()"
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        return self._rmul_(self._base_ring(-1))

    cpdef ModuleElement _rmul_(self, RingElement left):
        "Scalar multiplication"
        if self.Data == NULL:
            return self.__copy__()
        FfSetField(self.Data.Field)
        cdef Matrix_gfpn_dense OUT = self.__copy__()
        if MatMulScalar(OUT.Data, FfFromInt(self._converter.field_to_int(left))) != NULL:
            return OUT
        raise ArithmeticError("Matrix sizes or fields not compatible")

    cpdef ModuleElement _lmul_(self, RingElement right):
        "Scalar multiplication"
        if self.Data == NULL:
            return self.__copy__()
        FfSetField(self.Data.Field)
        cdef Matrix_gfpn_dense OUT = self.__copy__()
        if MatMulScalar(OUT.Data, FfFromInt(self._converter.field_to_int(right))) != NULL:
            return OUT
        raise ArithmeticError("Matrix sizes or fields not compatible")

    cdef Matrix _matrix_times_matrix_(self, Matrix right):
        # Surprisingly, Winograd-Strassen can compete with school book
        # multiplication for smallish matrices, and of course it is
        # asymptotically faster. So, we used it by default.
        return self._multiply_strassen(right)

    cpdef Matrix_gfpn_dense _multiply_classical(Matrix_gfpn_dense self, Matrix_gfpn_dense right):
        "multiply two meataxe matrices by the school book algorithm"
        if self.Data == NULL or right.Data == NULL:
            raise ValueError("The matrices must not be empty")
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")
        MS = self.matrix_space(self._nrows, right._ncols, False)
        cdef Matrix_gfpn_dense OUT = Matrix_gfpn_dense.__new__(Matrix_gfpn_dense)
        sig_on()
        OUT.Data = MatDup(self.Data)
        if OUT.Data == NULL:
            sig_off()
            raise MemoryError
        if not MatMul(OUT.Data,right.Data):
            sig_off()
            raise ArithmeticError("Matrix sizes or fields not compatible")
        sig_off()
        OUT._nrows = OUT.Data.Nor
        OUT._ncols = OUT.Data.Noc
        OUT._is_immutable = False
        OUT._parent = MS
        OUT._base_ring = self._base_ring
        OUT._converter = self._converter
        OUT._cache = {}
        return OUT

    cpdef Matrix_gfpn_dense _multiply_strassen(Matrix_gfpn_dense self, Matrix_gfpn_dense right, cutoff=0):
        """
        cutoff is NOT the number of rows/columns, but the rowsize expressed in bytes.
        If `cutoff==0` then the default ``sizeof(long)^2/2`` is chosen.
        """
        if self.Data == NULL or right.Data == NULL:
            raise ValueError("The matrices must not be empty")
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")
        MS = self.matrix_space(self._nrows, right._ncols, False)
        cdef Matrix_gfpn_dense OUT = Matrix_gfpn_dense(MS, None)
        # Now, OUT.Data is initialised, which is neede for MatrixMulStrassen to work.
        cutoff = cutoff//sizeof(long)
        StrassenSetCutoff(cutoff)
        sig_on()
        if MatMulStrassen(OUT.Data, self.Data, right.Data) == NULL:
            raise ArithmeticError("Error multiplying matrices by Strassen-Winograd algorithm")
        sig_off()
        return OUT

    cdef ModuleElement _mul_long(self, long n):
        "multiply an MTX matrix with a field element represented by an integer"
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        cdef Matrix_gfpn_dense left
        cdef FEL r
        if n < 0:
            r = mtx_taddinv[FfFromInt(-n)]
        else:
            r = FfFromInt(n)
        left = self.__copy__()
        if MatMulScalar(left.Data, r) != NULL:
            return left
        raise ArithmeticError("Matrix sizes or fields not compatible")

    def __div__(Matrix_gfpn_dense self, p):
        "divide an MTX matrix by a field element represented by an integer"
        if self.Data == NULL:
            return self.__copy__()
        if not p:
            raise ZeroDivisionError
        if p not in self._base_ring:
            raise ValueError("{} is not a scalar".format(p))
        FfSetField(self.Data.Field)
        cdef Matrix_gfpn_dense OUT = self.__copy__()
        cdef FEL r = mtx_tmultinv[FfFromInt(self._converter.field_to_int(p))]
        if MatMulScalar(OUT.Data, r) != NULL:
            return OUT
        raise ArithmeticError("Matrix sizes or fields not compatible")

    def __pow__(Matrix_gfpn_dense self, n, ignored):
        "M.__pow__(n): return M^n"
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        if ignored is not None:
            raise RuntimeError("__pow__ third argument not used")
        cdef Matrix_gfpn_dense OUT
        cdef Matrix_gfpn_dense SELFINV
        OUT = type(self).__new__(type(self))
        OUT._nrows = self._nrows
        OUT._ncols = self._ncols
        OUT._is_immutable = False
        OUT._parent = self._parent
        OUT._base_ring = self._base_ring
        OUT._converter = self._converter
        OUT._cache = {}
        if n>=0:
            OUT.Data = MatPower(self.Data,n)
        else:
            SELFINV = self.__invert__()
            OUT.Data = MatPower(SELFINV.Data,-n)
        if OUT.Data != NULL:
            return OUT
        raise ArithmeticError("Failure in exponentiating a matrix")

    def __invert__(Matrix_gfpn_dense self):
        "M__invert__(): return M^(-1)"
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        cdef Matrix_gfpn_dense OUT = type(self).__new__(type(self))
        OUT._nrows = self._nrows
        OUT._ncols = self._ncols
        OUT._is_immutable = False
        OUT._parent = self._parent
        OUT._base_ring = self._base_ring
        OUT._converter = self._converter
        OUT._cache = {}
        OUT.Data = MatInverse(self.Data)
        if OUT.Data != NULL:
            return OUT
        raise ArithmeticError("This matrix is not invertible")

    def transpose(Matrix_gfpn_dense self):
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        cdef Matrix_gfpn_dense OUT = type(self).__new__(type(self))
        OUT._nrows = self._ncols
        OUT._ncols = self._nrows
        OUT._is_immutable = False
        OUT._parent = self.matrix_space(self._ncols, self._nrows, False)
        OUT._base_ring = self._base_ring
        OUT._converter = self._converter
        OUT._cache = {}
        OUT.Data = MatTransposed(self.Data)
        return OUT

    def order(self):
        "M.order(): return multiplicative order of M"
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        if (self.Data.Nor <> self.Data.Noc):
            raise ValueError("only defined for square matrices")
        o = MatOrder(self.Data)
        if o==-1:
            raise ArithmeticError("order too large")
        else:
            return o

###################
## Gauss algorithm

    def nullity(self):
        "M.nullity(): return the nullity of M"
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        return MatNullity(self.Data)

    def left_kernel_matrix(self):
        """M.left_kernel_matrix(): return the null space of M

        M.left_kernel_matrix()*M is a null matrix
        """
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        cdef Matrix_gfpn_dense OUT = type(self).__new__(type(self))
        OUT.Data = MatNullSpace(self.Data)
        if OUT.Data == NULL:
            return OUT
        OUT._nrows = OUT.Data.Nor
        OUT._ncols = OUT.Data.Noc
        OUT._is_immutable = False
        OUT._parent = self.matrix_space(OUT._nrows, OUT._ncols, False)
        OUT._base_ring = self._base_ring
        OUT._converter = self._converter
        OUT._cache = {}
        return OUT

    def lead(self):
        """
(f,i) = M.lead() <=> f=M[0,i] is the first non-zero coefficient in the first row of M

If the first row of M has no non-zero entry then f==0
        """
        cdef int i
        cdef int fe
        if self.Data == NULL:
            raise ValueError("The matrix must not be empty")
        FfSetField(self.Data.Field)
        for i from 0 <= i < self.Data.Noc:
            fe = FfToInt(FfExtract(self.Data.Data,i))
            if fe:
                return fe, i
        return 0, self.Data.Noc

########################
### String representations
#    def __repr__(self):
#        "return a short description of an MTX matrix"
#        if self.Data == NULL:
#            return 'Empty MTX matrix'
#        return '(%s x %s) MTX matrix over GF(%s)'%(self.Data.Nor, self.Data.Noc, self.Data.Field)
#
#    def __str__(self):
#        "return a string showing the contents of an MTX matrix"
#        # cdef long i,j
#        if self.Data == NULL:
#            return '[]'
#        nc = self.Data.Noc
#        nr = self.Data.Nor
#        setfield(self.Data.Field)
#        fln = len(str(FfOrder))
#        matL = self._matlist_()
#        return "\n".join(["["+" ".join([str(el).rjust(fln) for el in matL[i]])+"]" \
#                                   for i in range(nr)])

###############################################################################
# Further features may be added later
###############################################################################

