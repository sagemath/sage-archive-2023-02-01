r"""
Dense vectors using a numpy backend.  This serves as a base class for
dense vectors over Real Double Field and Complex Double Field

EXAMPLES:
    sage: v = vector(CDF,[(1,-1), (2,pi), (3,5)])
    sage: v
    (1.0 - 1.0*I, 2.0 + 3.14159265359*I, 3.0 + 5.0*I)
    sage: type(v)
    <type 'sage.modules.vector_complex_double_dense.Vector_complex_double_dense'>
    sage: parent(v)
    Vector space of dimension 3 over Complex Double Field
    sage: v[0] = 5
    sage: v
    (5.0, 2.0 + 3.14159265359*I, 3.0 + 5.0*I)
    sage: loads(dumps(v)) == v
    True
    sage: v = vector(RDF, [1,2,3,4]); v
    (1.0, 2.0, 3.0, 4.0)
    sage: loads(dumps(v)) == v
    True

AUTHORS:
    -- Jason Grout, Oct 2008: switch to numpy backend, factored out
       Vector_double_dense class
    -- Josh Kantor
    -- William Stein
"""

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###############################################################################

cimport free_module_element
import  free_module_element

from sage.structure.element cimport Element, ModuleElement, RingElement, Vector

from sage.rings.complex_double import CDF
from sage.rings.complex_double cimport ComplexDoubleElement, new_ComplexDoubleElement

cimport sage.ext.numpy as cnumpy

numpy = None
scipy = None

# This is for the Numpy C API to work
cnumpy.import_array()

cdef class Vector_double_dense(free_module_element.FreeModuleElement):
    """
    Base class for vectors over the Real Double Field and the Complex
    Double Field.  These are supposed to be fast vector operations
    using C doubles. Most operations are implemented using numpy which
    will call the underlying BLAS, if needed, on the system.

    This class cannot be instantiated on its own.  The numpy vector
    creation depends on several variables that are set in the
    subclasses.

    EXAMPLES:
        sage: v = vector(RDF, [1,2,3,4]); v
        (1.0, 2.0, 3.0, 4.0)
        sage: v*v
        30.0
    """
    def __new__(self, parent, entries, coerce=True, copy=True):
        """
        Set up a new vector

        EXAMPLE:
        sage: v = vector(RDF, range(3))
        sage: v.__new__(v.__class__, v.parent(), [1,1,1]) # random output
        (3.00713073107e-261, 3.06320700422e-322, 2.86558074588e-322)
        sage: v = vector(CDF, range(3))
        sage: v.__new__(v.__class__, v.parent(), [1,1,1]) # random output
        (-2.26770549592e-39, 5.1698223615e-252*I, -5.9147262602e-62 + 4.63145528786e-258*I)
        """
        self._is_mutable = 1
        self._degree = parent.degree()
        self._parent = parent

    cdef Vector_double_dense _new(self, cnumpy.ndarray vector_numpy):
        """
        Return a new vector with same parent as self.
        """
        cdef Vector_double_dense v
        v = self.__class__.__new__(self.__class__,self._parent,None,None,None)
        v._is_mutable = 1
        v._parent = self._parent
        v._degree = self._parent.degree()

        v._vector_numpy = vector_numpy
        return v

    def __create_vector__(self):
        """
        Create a new uninitialized numpy array to hold the data for the class.

        This function assumes that self._numpy_dtypeint and
        self._nrows and self._ncols have already been initialized.

        EXAMPLE:
        In this example, we throw away the current array and make a
        new unitialized array representing the data for the class.
            sage: a=vector(RDF, range(9))
            sage: a.__create_vector__()
        """
        cdef cnumpy.npy_intp dims[1]
        dims[0] = self._degree
        self._vector_numpy = cnumpy.PyArray_SimpleNew(1, dims, self._numpy_dtypeint)
        return

    cdef bint is_dense_c(self):
        """
        Return True (i.e., 1) if self is dense.
        """
        return 1

    cdef bint is_sparse_c(self):
        """
        Return True (i.e., 1) if self is sparse.
        """
        return 0


    def __copy__(self, copy=True):
        """
        Return a copy of the vector

        EXAMPLE:
            sage: a = vector(RDF, range(9))
            sage: a == a.copy()
            True
        """
        if self._degree == 0:
            return self
        return self._new(self._vector_numpy.copy())

    def __init__(self, parent, entries, coerce = True, copy = True):
        """
        Fill the vector with entries.

        The numpy array must have already been allocated.

        EXAMPLES:
            sage: vector(RDF, range(9))
            (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0)
            sage: vector(CDF, 5)
            (0, 0, 0, 0, 0)

        TESTS:
            sage: vector(CDF, 0)
            ()
            sage: vector(RDF, 0)
            ()
            sage: vector(CDF, 4)
            (0, 0, 0, 0)
            sage: vector(RDF, 4)
            (0.0, 0.0, 0.0, 0.0)
            sage: vector(CDF, [CDF(1+I)*j for j in range(4)])
            (0, 1.0 + 1.0*I, 2.0 + 2.0*I, 3.0 + 3.0*I)
            sage: vector(RDF, 4, range(4))
            (0.0, 1.0, 2.0, 3.0)
        """
        cdef Py_ssize_t i,j
        if isinstance(entries,(tuple, list)):
            if len(entries)!=self._degree:
                    raise TypeError, "entries has wrong length"

            if coerce:
                for i from 0<=i<self._degree:
                    self.set_unsafe(i,self._python_dtype(entries[i]))
            else:
                for i from 0<=i<self._degree:
                    self.set_unsafe(i,entries[i])

        elif isinstance(entries, cnumpy.ndarray):
            self._replace_self_with_numpy(entries)
        else:
            cnumpy.PyArray_FILLWBYTE(self._vector_numpy, 0)
            if entries is None:
                print "entries is None"
                return
            else:
                try:
                    z = self._python_dtype(entries)
                except TypeError:
                    raise TypeError, "cannot coerce entry to type %s"%self._python_dtype
                if z != 0:
                    for i from 0<=i<self._degree:
                        self.set_unsafe(i,z)

    def __len__(self):
        """
        Return the length of the vector.

        EXAMPLE:
            sage: v = vector(RDF, 5); v
            (0.0, 0.0, 0.0, 0.0, 0.0)
            sage: len(v)
            5
        """
        return self._degree

    def __setitem__(self, Py_ssize_t i, object value):
        """
        Set the ith entry of self.

        EXAMPLES:
            sage: v = vector(CDF, [1,CDF(3,2), -1]); v
            (1.0, 3.0 + 2.0*I, -1.0)
            sage: v[1] = 2
            sage: v[-1] = I
            sage: v[5] = 2
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: v
            (1.0, 2.0, 1.0*I)
        """
        if not self._is_mutable:
            raise ValueError, "vector is immutable; please change a copy instead"
        if i < 0:
            i += self._degree

        if i < 0 or i >= self._degree:
            raise IndexError, 'index out of range'

        self.set_unsafe(i, value)

    def __getitem__(self, Py_ssize_t i):
        """
        Return the ith entry of self.

        EXAMPLES:
            sage: v = vector(CDF, [1,CDF(3,2), -1]); v
            (1.0, 3.0 + 2.0*I, -1.0)
            sage: v[1]
            3.0 + 2.0*I
            sage: v[-1]
            -1.0
            sage: v[5]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
            sage: v[-5]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        if i < 0:
            i += self._degree
        if i < 0 or i >= self._degree:
            raise IndexError, 'index out of range'
        return self.get_unsafe(i)

    cdef set_unsafe(self, Py_ssize_t i, object value):
        """
        Set the ith entry to value without any bounds checking,
        mutability checking, etc.
        """
        # We assume that Py_ssize_t is the same as npy_intp

        # We must patch the ndarrayobject.h file so that the SETITEM
        # macro does not have a semicolon at the end for this to work.
        # Cython wraps the macro in a function that converts the
        # returned int to a python object, which leads to compilation
        # errors because after preprocessing you get something that
        # looks like "););".  This is bug
        # http://scipy.org/scipy/numpy/ticket/918

        # We call the self._python_dtype function on the value since
        # numpy does not know how to deal with complex numbers other
        # than the builtin complex number type.
        cdef int status
        status = cnumpy.PyArray_SETITEM(self._vector_numpy,
                        cnumpy.PyArray_GETPTR1(self._vector_numpy, i),
                        self._python_dtype(value))
        #TODO: Throw an error if status == -1


    cdef get_unsafe(self, Py_ssize_t i):
        """
        Get the (i,j) entry without any bounds checking, etc.
        """
        # We assume that Py_ssize_t is the same as npy_intp
        return self._sage_dtype(cnumpy.PyArray_GETITEM(self._vector_numpy,
                                                cnumpy.PyArray_GETPTR1(self._vector_numpy, i)))


    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add two vectors together.

        EXAMPLES:
            sage: A = vector(RDF, range(3))
            sage: A+A
            (0.0, 2.0, 4.0)
        """
        if self._degree == 0:
            return self.copy()

        cdef Vector_double_dense _right, _left
        _right = right
        _left = self

        return self._new(_left._vector_numpy + _right._vector_numpy)

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Return self - right

        EXAMPLES:
            sage: A = vector(RDF, range(3))
            sage: (A-A).is_zero()
            True
        """
        if self._degree == 0:
            return self.copy()

        cdef Vector_double_dense _right, _left
        _right = right
        _left = self

        return self._new(_left._vector_numpy - _right._vector_numpy)

    cpdef Element _dot_product_(self, Vector right):
        """
        Dot product of self and right.

        EXAMPLES:
            sage: v = vector(RDF, [1,2,3]); w = vector(RDF, [2, 4, -3])
            sage: v*w
            1.0
            sage: w*v
            1.0
        """
        if not right.parent() == self.parent():
            right = self.parent().ambient_module()(right)
        if self._degree == 0:
            return self.copy()

        cdef Vector_double_dense _right, _left
        _right = right
        _left = self

        global scipy
        if scipy is None:
            import scipy
        return self._sage_dtype(scipy.dot(_left._vector_numpy, _right._vector_numpy))

    cpdef Vector _pairwise_product_(self, Vector right):
        """
        Return the component-wise product of self and right.

        EXAMPLES:
            sage: v = vector(CDF, [1,2,3]); w = vector(CDF, [2, 4, -3])
            sage: v.pairwise_product(w)
            (2.0, 8.0, -9.0)
        """
        if not right.parent() == self.parent():
            right = self.parent().ambient_module()(right)

        if self._degree == 0:
            return self.copy()

        cdef Vector_double_dense _right, _left
        _right = right
        _left = self

        return self._new(_left._vector_numpy * _right._vector_numpy)

    cpdef ModuleElement _rmul_(self, RingElement left):
        """
        Multiply a scalar and vector

        EXAMPLE:
            sage: v = vector(CDF, range(3))
            sage: 3*v
            (0, 3.0, 6.0)
        """
        if self._degree == 0:
            return self.copy()

        return self._new(self._python_dtype(left)*self._vector_numpy)


    cpdef ModuleElement _lmul_(self, RingElement right):
        """
        Multiply a scalar and vector

        EXAMPLE:
            sage: v = vector(CDF, range(3))
            sage: v*3
            (0, 3.0, 6.0)
        """
        if self._degree == 0:
            return self.copy()

        return self._new(self._vector_numpy*self._python_dtype(right))


    def inv_fft(self,algorithm="radix2", inplace=False):
        """
        This performs the inverse fast fourier transform on the vector.

        The fourier transform can be done in place using the keyword
        inplace=True

        This will be fastest if the vector's length is a power of 2.

        EXAMPLES:
            sage: v = vector(CDF,[1,2,3,4])
            sage: w = v.fft()
            sage: max(v - w.inv_fft()) < 1e-12
            True
        """
        return self.fft(direction="backward",algorithm=algorithm,inplace=inplace)

    def fft(self, direction = "forward", algorithm = "radix2", inplace=False):
        """
        This performs a fast fourier transform on the vector.

        INPUT:
           direction -- 'forward' (default) or 'backward'

           The algorithm and inplace arguments are ignored.

        This function is fastest if the vector's length is a power of 2.

        EXAMPLES:
            sage: v = vector(CDF,[1+2*I,2,3*I,4])
            sage: v.fft()
            (7.0 + 5.0*I, 1.0 + 1.0*I, -5.0 + 5.0*I, 1.0 - 3.0*I)
            sage: v.fft(direction='backward')
            (1.75 + 1.25*I, 0.25 - 0.75*I, -1.25 + 1.25*I, 0.25 + 0.25*I)
            sage: v.fft().fft(direction='backward')
            (1.0 + 2.0*I, 2.0, 3.0*I, 4.0)
            sage: v.fft().parent()
            Vector space of dimension 4 over Complex Double Field
            sage: v.fft(inplace=True)
            sage: v
            (7.0 + 5.0*I, 1.0 + 1.0*I, -5.0 + 5.0*I, 1.0 - 3.0*I)

            sage: v = vector(RDF,4,range(4)); v
            (0.0, 1.0, 2.0, 3.0)
            sage: v.fft()
            (6.0, -2.0 + 2.0*I, -2.0, -2.0 - 2.0*I)
            sage: v.fft(direction='backward')
            (1.5, -0.5 - 0.5*I, -0.5, -0.5 + 0.5*I)
            sage: v.fft().fft(direction='backward')
            (0, 1.0, 2.0, 3.0)
            sage: v.fft().parent()
            Vector space of dimension 4 over Complex Double Field
            sage: v.fft(inplace=True)
            Traceback (most recent call last):
            ...
            ValueError: inplace can only be True for CDF vectors
        """
        if direction not in ('forward', 'backward'):
            raise ValueError, "direction must be either 'forward' or 'backward'"

        if self._degree == 0:
            return self
        global scipy
        if scipy is None:
            import scipy
        import scipy.fftpack

        if inplace:
            if self._sage_dtype is not CDF:
                raise ValueError, "inplace can only be True for CDF vectors"
            if direction == 'forward':
                self._vector_numpy = scipy.fftpack.fft(self._vector_numpy, overwrite_x = True)
            else:
                self._vector_numpy = scipy.fftpack.ifft(self._vector_numpy, overwrite_x = True)
        else:
            V = CDF**self._degree
            from vector_complex_double_dense import Vector_complex_double_dense
            if direction == 'forward':
                return Vector_complex_double_dense(V, scipy.fft(self._vector_numpy))
            else:
                return Vector_complex_double_dense(V, scipy.ifft(self._vector_numpy))


    def numpy(self):
        """
        Return numpy array corresponding to this vector.

        EXAMPLES:
            sage: v = vector(CDF,4,range(4))
            sage: v.numpy()
            array([ 0.+0.j,  1.+0.j,  2.+0.j,  3.+0.j])
            sage: v = vector(CDF,0)
            sage: v.numpy()
            array([], dtype=complex128)
            sage: v = vector(RDF,4,range(4))
            sage: v.numpy()
            array([ 0.,  1.,  2.,  3.])
            sage: v = vector(RDF,0)
            sage: v.numpy()
            array([], dtype=float64)
        """
        return self._vector_numpy.copy()

    cdef _replace_self_with_numpy(self,cnumpy.ndarray numpy_array):
        """
        Replace the underlying numpy array with numpy_array.
        """
        if self._degree == 0:
            return
        shape = numpy_array.shape
        if len(shape)!=1 or len(self._vector_numpy) != shape[0]:
            raise ValueError, "vector lengths are not the same"

        self._vector_numpy = numpy_array.astype(self._numpy_dtype)

    def n(self, *args, **kwargs):
        """
        Returns a numerical approximation of self by calling the n()
        method on all of its entries.

        EXAMPLES:
            sage: v = vector(CDF, [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Complex Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Complex Field with 75 bits of precision
            sage: v = vector(RDF, [1,2,3])
            sage: v.n()
            (1.00000000000000, 2.00000000000000, 3.00000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 53 bits of precision
            sage: v.n(prec=75)
            (1.000000000000000000000, 2.000000000000000000000, 3.000000000000000000000)
            sage: _.parent()
            Vector space of dimension 3 over Real Field with 75 bits of precision
        """
        return free_module_element.vector( [e.n(*args, **kwargs) for e in self] )


    def complex_vector(self):
        """
        Return the associated complex vector, i.e., this vector but with
        coefficients viewed as complex numbers.

        EXAMPLES:
            sage: v = vector(RDF,4,range(4)); v
            (0.0, 1.0, 2.0, 3.0)
            sage: v.complex_vector()
            (0, 1.0, 2.0, 3.0)
            sage: v = vector(RDF,0)
            sage: v.complex_vector()
            ()
        """
        return self.change_ring(CDF)


    #############################
    # statistics
    #############################
    def mean(self):
        """
        Calculate the arithmetic mean of the vector.

        EXAMPLE:
            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.mean()
            4.0
            sage: w.mean()
            4.0 + 5.0*I
        """
        global scipy
        if scipy is None:
            import scipy
        import scipy.stats
        return self._sage_dtype(scipy.stats.mean(self._vector_numpy))

    def variance(self, population=True):
        """
        Calculate the variance of entries of the vector.

        INPUT:
            population -- If False, calculate the sample variance.

        EXAMPLE:
            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.variance()
            7.5
            sage: v.variance(population=False)
            6.66666666667
            sage: w.variance()
            -15.0*I
            sage: w.variance(population=False)
            -13.3333333333*I
        """
        global scipy
        if scipy is None:
            import scipy
        import scipy.stats

        if population is False:
            return self._sage_dtype(scipy.stats.samplevar(self._vector_numpy))
        else:
            return self._sage_dtype(scipy.stats.var(self._vector_numpy))

    def standard_deviation(self, population=True):
        """
        Calculate the standard deviation of entries of the vector.

        INPUT:
            population -- If False, calculate the sample standard deviation.

        EXAMPLES:
            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.standard_deviation()
            2.73861278753
            sage: v.standard_deviation(population=False)
            2.58198889747
            sage: w.standard_deviation()
            2.73861278753 - 2.73861278753*I
            sage: w.standard_deviation(population=False)
            2.58198889747 - 2.58198889747*I
        """
        global scipy
        if scipy is None:
            import scipy
        import scipy.stats

        if population is False:
            return self._sage_dtype(scipy.stats.samplestd(self._vector_numpy))
        else:
            return self._sage_dtype(scipy.stats.std(self._vector_numpy))


    def stats_kurtosis(self):
        """
        Compute the kurtosis of a dataset.

        Kurtosis is the fourth central moment divided by the square of
        the variance. Since we use Fisher's definition, 3.0 is
        subtracted from the result to give 0.0 for a normal
        distribution. (Paragraph from the scipy.stats docstring.)

        EXAMPLE:
            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.stats_kurtosis()
            -1.23
            sage: w.stats_kurtosis()
            -1.23
        """
        global scipy
        if scipy is None:
            import scipy
        import scipy.stats
        return self._sage_dtype(scipy.stats.kurtosis(self._vector_numpy))

    def prod(self):
        """
        Return the product of the entries of self.

        EXAMPLES:
            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.prod()
            0.0
            sage: w.prod()
            57204225.0*I
        """
        return self._sage_dtype(self._vector_numpy.prod())

    def sum(self):
        """
        Return the sum of the entries of self.

        EXAMPLES:
            sage: v = vector(RDF, range(9))
            sage: w = vector(CDF, [k+(9-k)*I for k in range(9)])
            sage: v.sum()
            36.0
            sage: w.sum()
            36.0 + 45.0*I
        """
        return self._sage_dtype(self._vector_numpy.sum())

