"""
Time Series

This is a module for working with discrete floating point time series.
It is designed so that every operation is very fast, typically much
faster than with other generic code, e.g., Python lists of doubles or
even numpy arrays.  The semantics of time series is more similar to
Python lists of doubles than Sage real double vectors or numpy 1d
arrays.   In particular, time series are not endowed with much
algebraic structure and are always mutable.

NOTES: Numpy arrays are faster at slicing, since slices return
references, and numpy arrays have strides.  However, this speed at
slicing makes numpy slower at certain other operations.

EXAMPLES:
    sage: set_random_seed(1)
    sage: t = finance.TimeSeries([random()-0.5 for _ in xrange(10)]); t
    [0.3294, 0.0959, -0.0706, -0.4646, 0.4311, 0.2275, -0.3840, -0.3528, -0.4119, -0.2933]
    sage: t.sums()
    [0.3294, 0.4253, 0.3547, -0.1099, 0.3212, 0.5487, 0.1647, -0.1882, -0.6001, -0.8933]
    sage: t.exponential_moving_average(0.7)
    [0.0000, 0.3294, 0.1660, 0.0003, -0.3251, 0.2042, 0.2205, -0.2027, -0.3078, -0.3807]
    sage: t.standard_deviation()
    0.33729638212891383
    sage: t.mean()
    -0.089334255069294391
    sage: t.variance()
    0.11376884939725425

AUTHOR:
    -- William Stein
"""

include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"
include "../ext/python_string.pxi"
include "../ext/random.pxi"

cdef extern from "math.h":
    double exp(double)
    double floor(double)
    double log(double)
    double pow(double, double)
    double sqrt(double)

cdef extern from "string.h":
    void* memcpy(void* dst, void* src, size_t len)

cdef extern from "arrayobject.h":
    cdef enum:
        NPY_OWNDATA = 0x0004 #bit mask so numpy does not free array contents when its destroyed
    object PyArray_FromDimsAndData(int,int*,int,double *)
    ctypedef int intp
    void import_array()
    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef int nd
        cdef int flags
        cdef intp *dimensions
        cdef char* data

from sage.rings.integer import Integer
from sage.rings.real_double import RDF
from sage.modules.real_double_vector cimport RealDoubleVectorSpaceElement

max_print = 10
digits = 4

cdef class TimeSeries:
    def __new__(self, values=None):
        """
        Create new empty uninitialized time series.

        EXAMPLES:
        This implicitly calls new.
            sage: finance.TimeSeries([1,3,-4,5])
            [1.0000, 3.0000, -4.0000, 5.0000]
        """
        self._values = NULL

    def __init__(self, values):
        """
        Initialize new time series.

        INPUT:
            values -- integer (number of values) or an iterable of floats

        EXAMPLES:
        This implicity calls init.
            sage: finance.TimeSeries([pi, 3, 18.2])
            [3.1416, 3.0000, 18.2000]

        Conversion from a numpy 1d array, which is very fast.
            sage: v = finance.TimeSeries([1..5])
            sage: w = v.numpy()
            sage: finance.TimeSeries(w)
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000]

        Conversion from an n-dimensional numpy array also works:
            sage: import numpy
            sage: v = numpy.array([[1,2], [3,4]], dtype=float); v
            array([[ 1.,  2.],
                   [ 3.,  4.]])
            sage: finance.TimeSeries(v)
            [1.0000, 2.0000, 3.0000, 4.0000]
        """
        cdef RealDoubleVectorSpaceElement z
        cdef ndarray np
        if isinstance(values, (int, long, Integer)):
            values = [0.0]*values
        elif PY_TYPE_CHECK(values, ndarray):
            np = values
            if np.nd != 1:
                np = values.reshape([values.size])
            self._length = np.dimensions[0]
            self._values = <double*> sage_malloc(sizeof(double) * self._length)
            if self._values == NULL:
                raise MemoryError
            memcpy(self._values, np.data, sizeof(double)*self._length)
            return
        elif PY_TYPE_CHECK(values, RealDoubleVectorSpaceElement):
            # Fast constructor from real double vector.
            z = values
            self._length = z.v.size
            self._values = <double*> sage_malloc(sizeof(double) * self._length)
            if self._values == NULL:
                raise MemoryError
            memcpy(self._values, z.v.data, sizeof(double)*self._length)
            return
        else:
            values = [float(x) for x in values]
        self._length = len(values)
        self._values = <double*> sage_malloc(sizeof(double) * self._length)
        if self._values == NULL:
            raise MemoryError
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            self._values[i] = values[i]

    def __reduce__(self):
        """
        Used in pickling time series.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-3.5])
            sage: v.__reduce__()
            (<built-in function unpickle_time_series_v1>,
             ('\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x0c\xc0', 2))
            sage: loads(dumps(v)) == v
            True

        Note that dumping and loading with compress False is much faster, though
        dumping with compress True can save a lot of space.
            sage: v = finance.TimeSeries([1..10^5])
            sage: loads(dumps(v, compress=False),compress=False) == v
            True
        """
        buf = PyString_FromStringAndSize(<char*>self._values, self._length*sizeof(double)/sizeof(char))
        return unpickle_time_series_v1, (buf, self._length)

    def __cmp__(self, _other):
        """
        Compare self and other.  This has the same semantics
        as list comparison.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,2,3]); w = finance.TimeSeries([1,2])
            sage: v < w
            False
            sage: w < v
            True
            sage: v == v
            True
            sage: w == w
            True
        """
        cdef TimeSeries other
        cdef Py_ssize_t c, i
        cdef double d
        if not isinstance(_other, TimeSeries):
            _other = TimeSeries(_other)
        other = _other
        for i from 0 <= i < min(self._length, other._length):
            d = self._values[i] - other._values[i]
            if d:
                return -1 if d < 0 else 1
        c = self._length - other._length
        if c < 0:
            return -1
        elif c > 0:
            return 1
        return 0

    def  __dealloc__(self):
        """
        Free up memory used by a time series.

        EXAMPLES:
        This tests __dealloc__ implicitly:
            sage: v = finance.TimeSeries([1,3,-4,5])
            sage: del v
        """
        if self._values:
            sage_free(self._values)

    def vector(self):
        """
        Return real double vector whose entries are the values of this
        time series.  This is useful since vectors have standard
        algebraic structure and play well with matrices.

        OUTPUT:
            a real double vector

        EXAMPLES:
            sage: v = finance.TimeSeries([1..10])
            sage: v.vector()
            (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0)
        """
        V = RDF**self._length
        cdef RealDoubleVectorSpaceElement x = RealDoubleVectorSpaceElement(V, 0)
        memcpy(x.v.data, self._values, sizeof(double)*self._length)
        return x

    def __repr__(self):
        """
        Return string representation of self.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,3.1908439,-4,5.93932])
            sage: v.__repr__()
            '[1.0000, 3.1908, -4.0000, 5.9393]'

        By default 4 digits after the decimal point are displayed.  To
        change this change self.finance.time_series.digits.
            sage: sage.finance.time_series.digits = 2
            sage: v.__repr__()
            '[1.00, 3.19, -4.00, 5.94]'
            sage: v
            [1.00, 3.19, -4.00, 5.94]
            sage: sage.finance.time_series.digits = 4
            sage: v
            [1.0000, 3.1908, -4.0000, 5.9393]
        """
        return self._repr()

    def _repr(self, prec=None):
        """
        Print representation of a time series.

        INPUT:
            prec -- number of digits of precision or None; if None
                    use the default sage.finance.time_series.digits
        OUTPUT:
             a string

        EXAMPLES:
            sage: v = finance.TimeSeries([1,3.1908439,-4,5.93932])
            sage: v._repr()
            '[1.0000, 3.1908, -4.0000, 5.9393]'
            sage: v._repr(10)
            '[1.0000000000, 3.1908439000, -4.0000000000, 5.9393200000]'
            sage: v._repr(2)
            '[1.00, 3.19, -4.00, 5.94]'
        """
        if prec is None: prec = digits
        format = '%.' + str(prec) + 'f'
        v = self.list()
        if len(v) > max_print:
            v0 = v[:max_print//2]
            v1 = v[-max_print//2:]
            return '[' + ', '.join([format%x for x in v0]) + ' ... ' + \
                         ', '.join([format%x for x in v1]) + ']'
        else:
            return '[' + ', '.join([format%x for x in v]) + ']'

    def __len__(self):
        """
        Return the number of entries in this time series.

        OUTPUT:
            Python integer

        EXAMPLES:
            sage: v = finance.TimeSeries([1,3.1908439,-4,5.93932])
            sage: v.__len__()
            4
            sage: len(v)
            4
        """
        return self._length

    def __getitem__(self, i):
        """
        Return i-th entry or slice of self.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-4,3,-2.5,-4,3])
            sage: v[2]
            3.0
            sage: v[-1]
            3.0
            sage: v[-10]
            Traceback (most recent call last):
            ...
            IndexError: TimeSeries index out of range
            sage: v[5]
            3.0
            sage: v[6]
            Traceback (most recent call last):
            ...
            IndexError: TimeSeries index out of range

        Some slice examples:
            sage: v[-3:]
            [-2.5000, -4.0000, 3.0000]
            sage: v[-3:-1]
            [-2.5000, -4.0000]
            sage: v[::2]
            [1.0000, 3.0000, -4.0000]
            sage: v[3:20]
            [-2.5000, -4.0000, 3.0000]
            sage: v[3:2]
            []

        Make a copy:
            sage: v[:]
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]

        Reverse the time series:
            sage: v[::-1]
            [3.0000, -4.0000, -2.5000, 3.0000, -4.0000, 1.0000]
        """
        cdef Py_ssize_t start, stop, step, j
        cdef TimeSeries t
        if PySlice_Check(i):
            start = 0 if (i.start is None) else i.start
            stop = self._length if (i.stop is None) else i.stop
            step = 1 if (i.step is None) else i.step
            if start < 0:
                start += self._length
                if start < 0: start = 0
            elif start >= self._length:
                start = self._length - 1
            if stop < 0:
                stop += self._length
                if stop < 0: stop = 0
            elif stop > self._length:
                stop = self._length
            if start >= stop:
                return new_time_series(0)
            if step < 0:
                step = -step
                t = new_time_series((stop-start)/step)
                for j from 0 <= j < (stop-start)/step:
                    t._values[j] = self._values[stop-1 - j*step]
            elif step > 1:
                t = new_time_series((stop-start)/step)
                for j from 0 <= j < (stop-start)/step:
                    t._values[j] = self._values[j*step+start]
            else:
                t = new_time_series(stop-start)
                memcpy(t._values, self._values + start, sizeof(double)*t._length)
            return t
        else:
            if i < 0:
                i += self._length
                if i < 0:
                    raise IndexError, "TimeSeries index out of range"
            elif i >= self._length:
                raise IndexError, "TimeSeries index out of range"
            return self._values[i]

    def __setitem__(self, Py_ssize_t i, double x):
        """
        Set the i-th entry of self to x.

        INPUT:
            i -- a nonnegative integer
            x -- a float

        EXAMPLES:
            sage: v = finance.TimeSeries([1,3,-4,5.93932]); v
            [1.0000, 3.0000, -4.0000, 5.9393]
            sage: v[0] = -5.5; v
            [-5.5000, 3.0000, -4.0000, 5.9393]
            sage: v[-1] = 3.2; v
            [-5.5000, 3.0000, -4.0000, 3.2000]
            sage: v[10]
            Traceback (most recent call last):
            ...
            IndexError: TimeSeries index out of range
            sage: v[-5]
            Traceback (most recent call last):
            ...
            IndexError: TimeSeries index out of range
        """
        if i < 0:
            i += self._length
            if i < 0:
                raise IndexError, "TimeSeries index out of range"
        elif i >= self._length:
            raise IndexError, "TimeSeries index out of range"
        self._values[i] = x

    def __copy__(self):
        """
        Return a copy of self.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-4,3,-2.5,-4,3])
            sage: v.__copy__()
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]
            sage: copy(v)
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]
            sage: copy(v) is v
            False
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        memcpy(t._values, self._values , sizeof(double)*self._length)
        return t

    def __add__(left, right):
        """
        Concatenate the time series self and right.

        NOTE: To add a single number to the entries of a time series,
        use the add_scalar method, and to add componentwise use
        the add_entries method.

        INPUT:
            right -- a time series
        OUTPUT:
            a time series

        EXAMPLES:
            sage: v = finance.TimeSeries([1,2,3]); w = finance.TimeSeries([1,2])
            sage: v + w
            [1.0000, 2.0000, 3.0000, 1.0000, 2.0000]
            sage: v = finance.TimeSeries([1,2,-5]); v
            [1.0000, 2.0000, -5.0000]

        Note that both summands must be a time series:
            sage: v + xrange(4)
            Traceback (most recent call last):
            ...
            TypeError: right operand must be a time series
            sage: [1,5] + v
            Traceback (most recent call last):
            ...
            TypeError: left operand must be a time series
        """
        if not isinstance(right, TimeSeries):
            raise TypeError, "right operand must be a time series"
        if not isinstance(left, TimeSeries):
            raise TypeError, "left operand must be a time series"
        cdef TimeSeries R = right
        cdef TimeSeries L = left
        cdef TimeSeries t = new_time_series(L._length + R._length)
        memcpy(t._values, L._values, sizeof(double)*L._length)
        memcpy(t._values + L._length, R._values, sizeof(double)*R._length)
        return t

    def __mul__(left, right):
        """
        Multiply a time series by an integer n, which (like for lists)
        results in the time series concatenated with itself n times.

        NOTE: To multiply all the entries of a time series by a single
        scalar, use the scale method.

        INPUT:
            left, right -- an integer and a time series
        OUTPUT:
            a time series

        EXAMPLES:
            sage: v = finance.TimeSeries([1,2,-5]); v
            [1.0000, 2.0000, -5.0000]
            sage: v*3
            [1.0000, 2.0000, -5.0000, 1.0000, 2.0000, -5.0000, 1.0000, 2.0000, -5.0000]
            sage: 3*v
            [1.0000, 2.0000, -5.0000, 1.0000, 2.0000, -5.0000, 1.0000, 2.0000, -5.0000]
            sage: v*v
            Traceback (most recent call last):
            ...
            TypeError: 'sage.finance.time_series.TimeSeries' object cannot be interpreted as an index
        """
        cdef Py_ssize_t n, i
        cdef TimeSeries T
        if isinstance(left, TimeSeries):
            T = left
            n = right
        else:
            T = right
            n = left
        # Make n copies of T concatenated together
        cdef TimeSeries v = new_time_series(T._length * n)
        for i from 0 <= i < n:
            memcpy(v._values + i*T._length, T._values, sizeof(double)*T._length)
        return v

    def linear_filter(self, M):
        """
        Return a linear filter using the first M autocovariance values
        of self.  This is useful in forcasting by taking a weighted
        average of previous values of a time series.

        WARNING: The input sequence is assumed to be stationary, which
        means that the autocovariance $\langle y_j y_k \rangle$ depends
        only on the difference $|j-k|$.

        INPUT:
            M -- integer

        OUTPUT:
            TimeSeries -- the weights in the linear filter.

        EXAMPLES:
            sage: set_random_seed(0)
            sage: v = finance.TimeSeries(10^4).randomize('normal').sums()
            sage: F = v.linear_filter(100)
            sage: v
            [0.6767, 0.2756, 0.6332, 0.0469, -0.8897 ... 87.6759, 87.6825, 87.4120, 87.6639, 86.3194]
            sage: v.linear_forecast(F)
            86.017728504280015
            sage: F
            [1.0148, -0.0029, -0.0105, 0.0067, -0.0232 ... -0.0106, -0.0068, 0.0085, -0.0131, 0.0092]
        """
        acvs = [self.autocovariance(i) for i in range(M+1)]
        return linear_filter(acvs)

    def linear_forecast(self, filter):
        """
        Given a linear filter as output by the linear_filter command,
        compute the forecast for the next term in the series.

        INPUT:
            filter -- a time series output by the linear filter command.

        EXAMPLES:
            sage: set_random_seed(0)
            sage: v = finance.TimeSeries(100).randomize('normal').sums()
            sage: F = v[:-1].linear_filter(5); F
            [1.0019, -0.0524, -0.0643, 0.1323, -0.0539]
            sage: v.linear_forecast(F)
            11.782029861181114
            sage: v
            [0.6767, 0.2756, 0.6332, 0.0469, -0.8897 ... 9.2447, 9.6709, 10.4037, 10.4836, 12.1960]
        """
        cdef TimeSeries filt
        if isinstance(filter, TimeSeries):
            filt = filter
        else:
            filt = TimeSeries(filter)

        cdef double f = 0
        cdef Py_ssize_t i
        for i from 0 <= i < min(self._length, filt._length):
            f += self._values[self._length - i - 1] * filt._values[i]
        return f

    def reversed(self):
        """
        Return new time series obtain from this time series by
        reversing the order of the entries in this time series.

        OUTPUT:
            time series

        EXAMPLES:
            sage: v = finance.TimeSeries([1..5])
            sage: v.reversed()
            [5.0000, 4.0000, 3.0000, 2.0000, 1.0000]
        """
        cdef Py_ssize_t i, n = self._length-1
        cdef TimeSeries t = new_time_series(self._length)

        for i from 0 <= i < self._length:
            t._values[i] = self._values[n - i]
        return t

    def extend(self, right):
        """
        Extend this time series by appending elements from the iterable right.

        INPUT:
            right -- iterable that can be converted to a time series

        EXAMPLES:
            sage: v = finance.TimeSeries([1,2,-5]); v
            [1.0000, 2.0000, -5.0000]
            sage: v.extend([-3.5, 2])
            sage: v
            [1.0000, 2.0000, -5.0000, -3.5000, 2.0000]
        """
        if not isinstance(right, TimeSeries):
            right = TimeSeries(right)
        if len(right) == 0:
            return
        cdef TimeSeries T = right
        cdef double* z = <double*> sage_malloc(sizeof(double)*(self._length + T._length))
        if z == NULL:
            raise MemoryError
        memcpy(z, self._values, sizeof(double)*self._length)
        memcpy(z + self._length, T._values, sizeof(double)*T._length)
        sage_free(self._values)
        self._values = z
        self._length = self._length + T._length

    def list(self):
        """
        Return list of elements of self.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-4,3,-2.5,-4,3])
            sage: v.list()
            [1.0, -4.0, 3.0, -2.5, -4.0, 3.0]
        """
        v = [0.0]*self._length
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            v[i] = self._values[i]
        return v

    def log(self):
        """
        Return new time series got by taking the logarithms of all the
        terms in the time series.

        OUTPUT:
            a new time series.

        EXAMPLES:
        We exponentiate then log a time seris and get back
        the original series.
            sage: v = finance.TimeSeries([1,-4,3,-2.5,-4,3]); v
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]
            sage: v.exp()
            [2.7183, 0.0183, 20.0855, 0.0821, 0.0183, 20.0855]
            sage: v.exp().log()
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]
        """
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            if self._values[i] <= 0:
                raise ValueError, "every entry of self must be positive."

        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = log(self._values[i])
        return t

    def exp(self):
        """
        Return new time series got by applying the exponential map to
        all the terms in the time series.

        OUTPUT:
            a new time series.

        EXAMPLES:
            sage: v = finance.TimeSeries([1..5]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000]
            sage: v.exp()
            [2.7183, 7.3891, 20.0855, 54.5982, 148.4132]
            sage: v.exp().log()
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000]
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = exp(self._values[i])
        return t

    def abs(self):
        """
        Return new time series got by replacing all entries
        of self by their absolute value.

        OUTPUT:
            a new time series

        EXAMPLES:
            sage: v = finance.TimeSeries([1,3.1908439,-4,5.93932])
            sage: v
            [1.0000, 3.1908, -4.0000, 5.9393]
            sage: v.abs()
            [1.0000, 3.1908, 4.0000, 5.9393]
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = self._values[i] if self._values[i] >= 0 else -self._values[i]
        return t

    def diffs(self, Py_ssize_t k = 1):
        """
        Return the new time series got by taking the differences of
        successive terms in the time series.  So if self is the time
        series $X_0, X_1, X_2, ...$, then this function outputs the
        series $X_1 - X_0, X_2 - X_1, ...$.  The output series has one
        less term than the input series.  If the optional parameter
        $k$ is given, return $X_k - X_0, X_{2k} - X_k, ...$.

        INPUT:
            k -- positive integer (default: 1)

        OUTPUT:
            a new time series.

        EXAMPLES:
            sage: v = finance.TimeSeries([5,4,1.3,2,8]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000]
            sage: v.diffs()
            [-1.0000, -2.7000, 0.7000, 6.0000]
        """
        if k != 1:
            return self.scale_time(k).diffs()
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length - 1)
        for i from 1 <= i < self._length:
            t._values[i-1] = self._values[i] - self._values[i-1]
        return t

    def scale_time(self, Py_ssize_t k):
        """
        Return the new time series at scale k.  If the input
        time series is $X_0, X_1, X_2, ...$, then this function
        returns the shorter time series $X_0, X_k, X_{2k}, ...$.

        INPUT:
            k -- a positive integer

        OUTPUT:
            a new time series.

        EXAMPLES:
            sage: v = finance.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.scale_time(1)
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.scale_time(2)
            [5.0000, 1.3000, 8.0000, 3.0000]
            sage: v.scale_time(3)
            [5.0000, 2.0000]
            sage: v.scale_time(10)
            []

        A series of odd length:
            sage: v = finance.TimeSeries([1..5]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000]
            sage: v.scale_time(2)
            [1.0000, 3.0000, 5.0000]

        TESTS:
            sage: v.scale_time(0)
            Traceback (most recent call last):
            ...
            ValueError: k must be positive
            sage: v.scale_time(-1)
            Traceback (most recent call last):
            ...
            ValueError: k must be positive
        """
        if k <= 0:
            raise ValueError, "k must be positive"

        cdef Py_ssize_t i, n
        n = self._length/k
        if self._length % 2:
            n += 1
        cdef TimeSeries t = new_time_series(n)  # in C / is floor division.
        for i from 0 <= i < n:
            t._values[i] = self._values[i*k]
        return t

    def scale(self, double s):
        """
        Return new time series obtained by multiplying every value in the series by s.

        INPUT:
            s -- float
        OUTPUT:
            a new time series with all values multiplied by s.

        EXAMPLES:
            sage: v = finance.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.scale(0.5)
            [2.5000, 2.0000, 0.6500, 1.0000, 4.0000, 5.0000, 1.5000, -2.5000]
        """
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = self._values[i] * s
        return t

    def add_scalar(self, double s):
        """
        Return new time series obtained by adding a scalar to every
        value in the series.

        NOTE: To add componentwise use the add_entries method.

        INPUT:
            s -- float
        OUTPUT:
            a new time series with s added to all values.

        EXAMPLES:
            sage: v = finance.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.add_scalar(0.5)
            [5.5000, 4.5000, 1.8000, 2.5000, 8.5000, 10.5000, 3.5000, -4.5000]
        """
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = self._values[i] + s
        return t

    def add_entries(self, t):
        """
        Add corresponding entries of self and t together,
        extending either self or t by 0's if they do
        not have the same length.

        NOTE: To add a single number to the entries of a time series,
        use the add_scalar method.

        INPUT:
            t -- a time seris
        OUTPUT:
            a time series with length the maxima of the lengths of
            self and t.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,2,-5]); v
            [1.0000, 2.0000, -5.0000]
            sage: v.add_entries([3,4])
            [4.0000, 6.0000, -5.0000]
            sage: v.add_entries(v)
            [2.0000, 4.0000, -10.0000]
            sage: v.add_entries([3,4,7,18.5])
            [4.0000, 6.0000, 2.0000, 18.5000]
        """
        if not PY_TYPE_CHECK(t, TimeSeries):
            t = TimeSeries(t)
        cdef Py_ssize_t i, n
        cdef TimeSeries T = t, shorter, longer
        if self._length > T._length:
            shorter = T
            longer = self
        else:
            shorter = self
            longer = T

        n = longer._length
        cdef TimeSeries v = new_time_series(n)

        # Iterate up to the length of the shorter one
        for i from 0 <= i < shorter._length:
            v._values[i] = shorter._values[i] + longer._values[i]

        # Copy over the rest of the values
        if shorter._length != n:
            memcpy(v._values + shorter._length,
                   longer._values + shorter._length,
                   sizeof(double)*(v._length - shorter._length))

        return v

    def show(self, *args, **kwds):
        """
        Calls plot and passes all arguments onto the plot function.  This is thus
        just an alias for plot.

        EXAMPLES:
        Draw a plot of a time series:
            sage: finance.TimeSeries([1..10]).show()
        """
        return self.plot(*args, **kwds)

    def plot(self, Py_ssize_t plot_points=1000, points=False, **kwds):
        """
        Return a plot of this time series as a line or points through
        (i,T(i)), where i ranges over nonnegative integers up to the
        length of self.

        INPUT:
            plot_points -- (default: 1000) 0 or positive integer; only
                           plot the given number of equally spaced
                           points in the time series; if 0, plot all points
            points -- bool (default: False) -- if True, return just the
                      points of the time series
            **kwds -- passed to the line or point command

        EXAMPLES:
            sage: v = finance.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.plot()
            sage: v.plot(points=True)
            sage: v.plot() + v.plot(points=True, rgbcolor='red')
            sage: v.plot() + v.plot(points=True, rgbcolor='red',pointsize=50)
        """
        from sage.plot.all import line, point
        cdef Py_ssize_t s

        if self._length < plot_points:
            plot_points = 0

        if plot_points > 0:
            s = self._length/plot_points
            if plot_points <= 0:
                raise ValueError, "plot_points must be a positive integer"
            v = self.scale_time(s).list()[:plot_points]
        else:
            s = 1
            v = self.list()
        w = [(i * s, y) for i,y in enumerate(v)]
        if points:
            L = point(w, **kwds)
        else:
            L = line(w, **kwds)
        L.ymin(min(v))
        L.ymax(max(v))
        L.xmin(0)
        L.xmax(len(v)*s)
        return L

    def simple_moving_average(self, Py_ssize_t k):
        """
        Return the moving average time series over the last k time units.
        Assumes the input time series was constant with its starting value
        for negative time.  The t-th step of the output is the sum of
        the previous k-1 steps of self and the kth step divided by k.
        Thus k values are avaraged at each point.

        INPUT:
            k -- positive integer

        OUTPUT:
            a time series with the same number of steps as self.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.simple_moving_average(0)
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.simple_moving_average(1)
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.simple_moving_average(2)
            [1.0000, 1.0000, 1.0000, 1.5000, 2.5000]
            sage: v.simple_moving_average(3)
            [1.0000, 1.0000, 1.0000, 1.3333, 2.0000]
        """
        if k == 0 or k == 1:
            return self.__copy__()
        if k <= 0:
            raise ValueError, "k must be positive"
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        if self._length == 0:
            return t

        cdef double s = self._values[0] * k
        for i from 0 <= i < self._length:
            if i >= k:
                s -= self._values[i-k]
            else:
                s -= self._values[0]
            # I have a slight concern about accumulated rounding error given how
            # this algorithm adds and subtracts.
            s += self._values[i]
            t._values[i] = s/k
        return t

    def exponential_moving_average(self, double alpha):
        """
        Return the exponential moving average time series.  Assumes
        the input time series was constant with its starting value for
        negative time.  The t-th step of the output is the sum of the
        previous k-1 steps of self and the kth step divided by k.

        The 0th term is formally undefined, so we define it to be 0,
        and we define the first term to be self[0].

        INPUT:
            alpha -- float; a smoothing factor with 0 <= alpha <= 1

        OUTPUT:
            a time series with the same number of steps as self.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.exponential_moving_average(0)
            [0.0000, 1.0000, 1.0000, 1.0000, 1.0000]
            sage: v.exponential_moving_average(1)
            [0.0000, 1.0000, 1.0000, 1.0000, 2.0000]
            sage: v.exponential_moving_average(0.5)
            [0.0000, 1.0000, 1.0000, 1.0000, 1.5000]

        Some more examples:
            sage: v = finance.TimeSeries([1,2,3,4,5])
            sage: v.exponential_moving_average(1)
            [0.0000, 1.0000, 2.0000, 3.0000, 4.0000]
            sage: v.exponential_moving_average(0)
            [0.0000, 1.0000, 1.0000, 1.0000, 1.0000]
        """
        if alpha < 0 or alpha > 1:
            raise ValueError, "alpha must be between 0 and 1"
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        if self._length == 0:
            return t
        t._values[0] = 0
        if self._length == 1:
            return t
        t._values[1] = self._values[0]
        for i from 2 <= i < self._length:
            t._values[i] = alpha * self._values[i-1] + (1-alpha) *t._values[i-1]
        return t

    def sums(self, double s=0):
        """
        Return the new time series got by taking the running partial
        sums of the terms of this time series.

        INPUT:
            s -- starting value for partial sums
        OUTPUT:
            TimeSeries

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.sums()
            [1.0000, 2.0000, 3.0000, 5.0000, 8.0000]
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            s += self._values[i]
            t._values[i] = s
        return t

    def sum(self):
        """
        Return the sum of all the entries of self.  If self has
        length 0, returns 0.

        OUTPUT:
            double

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.sum()
            8.0
        """
        cdef double s = 0
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            s += self._values[i]
        return s

    def prod(self):
        """
        Return the prod of all the entries of self.  If self has
        length 0, returns 1.

        OUTPUT:
            double

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.prod()
            6.0
        """
        cdef double s = 1
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            s *= self._values[i]
        return s


    def mean(self):
        """
        Return the mean (average) of the elements of self.

        OUTPUT:
            double

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.mean()
            1.6000000000000001
        """
        return self.sum() / self._length

    def pow(self, double k):
        """
        Return new time series with every elements of self raised to the
        kth power.

        INPUT:
            k -- float
        OUTPUT:
            time series

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.pow(2)
            [1.0000, 1.0000, 1.0000, 4.0000, 9.0000]
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = pow(self._values[i], k)
        return t

    def moment(self, int k):
        """
        Return the k-th moment of self, which is just the
        mean of the k-th powers of the elements of self.

        INPUT:
            k -- a positive integer

        OUTPUT:
            double

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.moment(1)
            1.6000000000000001
            sage: v.moment(2)
            3.2000000000000002
        """
        if k <= 0:
            raise ValueError, "k must be positive"
        if k == 1:
            return self.mean()
        cdef double s = 0
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            s += pow(self._values[i], k)
        return s / self._length

    def central_moment(self, int k):
        """
        Return the k-th central moment of self, which is just the mean
        of the k-th powers of the differences self[i]-mu, where mu is
        the mean of self.

        INPUT:
            k -- a positive integer
        OUTPUT:
            double

        EXAMPLES:
            sage: v = finance.TimeSeries([1,2,3])
            sage: v.central_moment(2)
            0.66666666666666663

        Note that the central moment is different than the moment
        here, since the mean is not $0$:
            sage: v.moment(2)     # doesn't subtract off mean
            4.666666666666667

        We compute the central moment directly:
            sage: mu = v.mean(); mu
            2.0
            sage: ((1-mu)^2 + (2-mu)^2 + (3-mu)^2) / 3
            0.66666666666666663
        """
        if k == 1:
            return float(0)
        mu = self.mean()
        # We could make this slightly faster by doing the add scalar
        # and moment calculation together.  But that would be nasty.
        return self.add_scalar(-mu).moment(k)

    def covariance(self, TimeSeries other):
        r"""
        Return the covariance of the time series self and other.

        INPUT:
            self, other -- time series

        Whichever time series has more terms is truncated.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-2,3]); w = finance.TimeSeries([4,5,-10])
            sage: v.covariance(w)
            -11.777777777777779
        """
        cdef double mu = self.mean(), mu2 = other.mean()
        cdef double s = 0
        cdef Py_ssize_t i
        cdef Py_ssize_t n = min(self._length, other._length)
        for i from 0 <= i < n:
            s += (self._values[i] - mu)*(other._values[i] - mu2)
        return s / n
        # NOTE: There is also a formula
        #       self._values[i]*other._values[i] - mu*mu2
        # but that seems less numerically stable (?), and when tested
        # was not noticeably faster.

    def autocovariance(self, Py_ssize_t k=0):
        r"""
        Return the k-th autocovariance function $\gamma(k)$ of self.
        This is the covariance of self with self shifted by $k$, i.e.,
        $$
        ( \sum_{t=0}^{n-k-1} (x_t - \mu)(x_{t + k} - \mu) ) / n,
        $$
        where $n$ is the length of self.

        Note the denominator of $n$, which gives a "better" sample
        estimatator.

        INPUT:
            k -- a nonnegative integer (default: 0)

        OUTPUT:
            float

        EXAMPLES:
            sage: v = finance.TimeSeries([13,8,15,4,4,12,11,7,14,12])
            sage: v.autocovariance(0)
            14.4
            sage: mu = v.mean(); sum([(a-mu)^2 for a in v])/len(v)
            14.4
            sage: v.autocovariance(1)
            -2.7000000000000002
            sage: mu = v.mean(); sum([(v[i]-mu)*(v[i+1]-mu) for i in range(len(v)-1)])/len(v)
            -2.7000000000000002
            sage: v.autocovariance(1)
            -2.7000000000000002

        We illustrate with a random sample that an independently and
        identically distributed distribution with zero mean and
        variance $\sigma^2$ has autocovariance function $\gamma(h)$
        with $\gamma(0) = \sigma^2$ and $\gamma(h) = 0$ for $h\neq 0$.
            sage: set_random_seed(0)
            sage: v = finance.TimeSeries(10^6)
            sage: v.randomize('normal', 0, 5)
            [3.3835, -2.0055, 1.7882, -2.9319, -4.6827 ... -5.1868, 9.2613, 0.9274, -6.2282, -8.7652]
            sage: v.autocovariance(0)
            24.954106897195892
            sage: v.autocovariance(1)
            -0.0050839047886276651
            sage: v.autocovariance(2)
            0.022056325329509487
            sage: v.autocovariance(3)
            -0.019020003743134766
        """
        cdef double mu = self.mean()
        cdef double s = 0
        cdef Py_ssize_t i
        cdef Py_ssize_t n = self._length - k
        for i from 0 <= i < n:
            s += (self._values[i] - mu)*(self._values[i+k] - mu)
        return s / self._length

    def correlation(self, TimeSeries other):
        """
        Return the correlation of self and other, which is the
        covariance of self and other divided by the product of their
        standard deviation.

        INPUT:
            self, other -- time series

        Whichever time series has more terms is truncated.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-2,3]); w = finance.TimeSeries([4,5,-10])
            sage: v.correlation(w)
            -0.55804160922502144
            sage: v.covariance(w)/(v.standard_deviation() * w.standard_deviation())
            -0.55804160922502144
        """
        return self.covariance(other) / (self.standard_deviation() * other.standard_deviation())

    def autocorrelation(self, Py_ssize_t k=1):
        r"""
        Return the $k$th sample autocorrelation of this time series
        $x_i$.

        Let $\mu$ be the sample mean.  Then the sample autocorrelation
        function is
           $$
              \frac{\sum_{t=0}^{n-k-1} (x_t - \mu)(x_{t+k} - \mu) }
                   {\sum_{t=0}^{n-1}   (x_t - \mu)^2}
           $$

        Note that the variance must be nonzero or you will get a
        ZeroDivisionError.

        INPUT:
            k -- a nonnegative integer (default: 1)

        OUTPUT:
            Time series

        EXAMPLE:
            sage: v = finance.TimeSeries([13,8,15,4,4,12,11,7,14,12])
            sage: v.autocorrelation()
            -0.1875
            sage: v.autocorrelation(1)
            -0.1875
            sage: v.autocorrelation(0)
            1.0
            sage: v.autocorrelation(2)
            -0.20138888888888887
            sage: v.autocorrelation(3)
            0.18055555555555555

            sage: finance.TimeSeries([1..1000]).autocorrelation()
            0.997
        """
        return self.autocovariance(k) / self.variance(bias=True)

    def variance(self, bias=False):
        """
        Return the variance of the elements of self, which is the mean
        of the squares of the differences from the mean.

        INPUT:
            bias -- bool (default: False); if False, divide by
                    self.length() - 1 instead of self.length()
                    to give a less biased estimator for the variance
        OUTPUT:
            double

        EXAMPLE:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.variance()
            0.80000000000000004
            sage: v.variance(bias=True)
            0.64000000000000001
        """
        cdef double mu = self.mean()
        cdef double s = 0
        cdef double a
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            a = self._values[i] - mu
            s += a * a
        if bias:
            return s / self._length
        else:
            return s / (self._length - 1)

    def standard_deviation(self, bias=False):
        """
        Return the standard deviation of the entries of self.

        INPUT:
            bias -- bool (default: False); if False, divide by
                    self.length() - 1 instead of self.length()
                    to give a less biased estimator for the variance
        OUTPUT:
            double

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.standard_deviation()
            0.89442719099991586
            sage: v.standard_deviation(bias=True)
            0.80000000000000004
        """
        return sqrt(self.variance(bias=bias))

    def range_statistic(self):
        """
        Return the range statistic of self, which is the difference
        between the maximum and minimum values of this time series,
        divided by the standard deviation of the series of differences.

        OUTPUT:
            float

        EXAMPLES:
        Notice that if we make a Brownian motion random walk, there
        is no difference if we change the standard deviation.
            sage: set_random_seed(0); finance.TimeSeries(10^6).randomize('normal').sums().range_statistic()
            1873.9206979719115
            sage: set_random_seed(0); finance.TimeSeries(10^6).randomize('normal',0,100).sums().range_statistic()
            1873.920697971955
        """
        return (self.max() - self.min())/self.diffs().standard_deviation()

    def hurst_exponent(self):
        """
        Returns a very simple naive estimate of the Hurst exponent of
        this time series.  There are many possible ways to estimate
        this number, and this is perhaps the most naive.  The estimate
        we take here is the log of the range statistic divided by the
        length of this time series.

        We define the Hurst exponent of a constant time series to be 1.

        NOTE: This is only a very rough estimator.  There are supposed
        to be better ones that use wavelets.

        EXAMPLES:
        The Hurst exponent of Brownian motion is 1/2.  We approximate
        it with some specific samples:
            sage: set_random_seed(0)
            sage: bm = finance.TimeSeries(10^5).randomize('normal').sums(); bm
            [0.6767, 0.2756, 0.6332, 0.0469, -0.8897 ... 152.2437, 151.5327, 152.7629, 152.9169, 152.9084]
            sage: bm.hurst_exponent()
            0.5174890556918027

        We compute the Hurst exponent of a simulated fractional Brownian
        motion with Hurst parameter 0.7.  This function estimates the
        Hurst exponent as 0.6678...
            sage: set_random_seed(0)
            sage: fbm = finance.fractional_brownian_motion_simulation(0.7,0.1,10^5,1)[0].sums()
            sage: fbm.hurst_exponent()
            0.66787027921443409

        Another example with small Hurst exponent (notice how bad the prediction is...):
            sage: fbm = finance.fractional_brownian_motion_simulation(0.2,0.1,10^5,1)[0].sums()
            sage: fbm.hurst_exponent()
            0.30450273560706259

        The above example illustrate that this is not a very good
        estimate of the Hurst exponent.
        """
        cdef double r = self.range_statistic()
        if r == 0:
            return float(1)
        return log(r)/log(self._length)

    def min(self, bint index=False):
        """
        Return the smallest value in this time series. If this series
        has length 0 we raise a ValueError.

        INPUT:
            index -- bool (default: False); if True, also return index of
                     minimal entry.
        OUTPUT:
            float -- smallest value
            integer -- index of smallest value; only returned if index=True

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-4,3,-2.5,-4])
            sage: v.min()
            -4.0
            sage: v.min(index=True)
            (-4.0, 1)
        """
        if self._length == 0:
            raise ValueError, "min() arg is an empty sequence"
            return 0.0, -1
        cdef Py_ssize_t i, j
        cdef double s = self._values[0]
        j = 0
        for i from 1 <= i < self._length:
            if self._values[i] < s:
                s = self._values[i]
                j = i
        if index:
            return s, j
        else:
            return s

    def max(self, bint index=False):
        """
        Return the largest value in this time series. If this series
        has length 0 we raise a ValueError

        INPUT:
            index -- bool (default: False); if True, also return index of
                     maximum entry.
        OUTPUT:
            float -- largest value
            integer -- index of largest value; only returned if index=True

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-4,3,-2.5,-4,3])
            sage: v.max()
            3.0
            sage: v.max(index=True)
            (3.0, 2)
        """
        if self._length == 0:
            raise ValueError, "max() arg is an empty sequence"
        cdef Py_ssize_t i, j
        cdef double s = self._values[0]
        for i from 1 <= i < self._length:
            if self._values[i] > s:
                s = self._values[i]
                j = i
        if index:
            return s, j
        else:
            return s

    def clip_remove(self, min=None, max=None):
        """
        Return new time series obtained from self by removing all
        values that are less than or equal to a certain miminum value
        or greater than or equal to a certain maximum.

        INPUT:
            min -- None or double
            max -- None or double

        OUTPUT:
            time series

        EXAMPLES:
            sage: v = finance.TimeSeries([1..10])
            sage: v.clip_remove(3,7)
            [3.0000, 4.0000, 5.0000, 6.0000, 7.0000]
            sage: v.clip_remove(3)
            [3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
            sage: v.clip_remove(max=7)
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000]
        """
        cdef Py_ssize_t i, j
        cdef TimeSeries t
        cdef double mn, mx
        cdef double x
        if min is None and max is None:
            return self.copy()
        # This code is ugly but I think as fast as possible.
        if min is None and max is not None:
            # Return everything <= max
            j = 0
            mx = max
            for i from 0 <= i < self._length:
                if self._values[i] <= mx:
                    j += 1

            t = TimeSeries(j)
            j = 0
            for i from 0 <= i < self._length:
                if self._values[i] <= mx:
                    t._values[j] = self._values[i]
                    j += 1
        elif max is None and min is not None:
            # Return everything >= min
            j = 0
            mn = min
            for i from 0 <= i < self._length:
                if self._values[i] >= mn:
                    j += 1
            t = TimeSeries(j)
            j = 0
            for i from 0 <= i < self._length:
                if self._values[i] >= mn:
                    t._values[j] = self._values[i]
                    j += 1
        else:
            # Return everything between min and max
            j = 0
            mn = min; mx = max
            for i from 0 <= i < self._length:
                x = self._values[i]
                if x >= mn and x <= mx:
                    j += 1
            t = TimeSeries(j)
            j = 0
            for i from 0 <= i < self._length:
                x = self._values[i]
                if x >= mn and x <= mx:
                    t._values[j] = x
                    j += 1
        return t

    def histogram(self, Py_ssize_t bins=50, bint normalize=False):
        """
        Return the frequency histogram of the values in
        this time series divided up into the given
        numberof bins.

        INPUT:
            bins -- a positive integer (default: 50)

        OUTPUT:
            counts -- list of counts of numbers of elements in
                      each bin
            endpoints -- list of 2-tuples (a,b) that give the
                      endpoints of the bins

        EXAMPLES:
            sage: v = finance.TimeSeries([5,4,1.3,2,8,10,3,-5]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000, 10.0000, 3.0000, -5.0000]
            sage: v.histogram(3)
            ([1, 4, 3], [(-5.0, 0.0), (0.0, 5.0), (5.0, 10.0)])
        """
        if bins <= 0:
            raise ValueError, "bins must be positive"

        cdef double mn = self.min(), mx = self.max()
        cdef double r = mx - mn, step = r/bins

        if r == 0:
            raise ValueError, "bins have 0 width"

        v = [(mn + j*step, mn + (j+1)*step) for j in range(bins)]
        if self._length == 0:
            return [], v

        if step == 0:
            counts = [0]*bins
            counts[0] = [self._length]
            return counts, v

        cdef Py_ssize_t i
        cdef Py_ssize_t* cnts = <Py_ssize_t*>sage_malloc(sizeof(Py_ssize_t)*bins)
        if cnts == NULL:
            raise MemoryError
        for i from 0 <= i < bins:
            cnts[i] = 0

        cdef Py_ssize_t j
        for i from 0 <= i < self._length:
            j = int((self._values[i] - mn)/step)
            if j >= bins:
                j = bins-1
            cnts[j] += 1

        b = 1.0/(self._length * step)
        if normalize:
            counts = [cnts[i]*b for i in range(bins)]
        else:
            counts = [cnts[i] for i in range(bins)]
        sage_free(cnts)

        return counts, v

    def plot_histogram(self, bins=50, normalize=True, **kwds):
        """
        Return histogram plot of this time series with given number of bins.

        INPUT:
            bins -- positive integer (default: 50)
            normalize -- (default: True) whether to normalize so the total
                         area in the bars of the histogram is 1.
            **kwds -- passed to the bar_chart function
        OUTPUT:
            a histogram plot

        EXAMPLES:
            sage: v = finance.TimeSeries([1..50])
            sage: v.plot_histogram(bins=10)
        """
        from sage.plot.all import bar_chart, polygon
        counts, intervals = self.histogram(bins, normalize=normalize)
        s = 0
        for i, (x0,x1) in enumerate(intervals):
            s += polygon([(x0,0), (x0,counts[i]), (x1,counts[i]), (x1,0)], **kwds)
        if len(intervals) > 0:
            s.xmin(intervals[0][0])
            s.xmax(intervals[-1][1])
            s.ymin(0)
            s.ymax(max(counts))
        return s

    def numpy(self, copy=True):
        """
        Return version of this time series in numpy.

        NOTE: If copy is False, return a numpy 1d array reference to
        exactly the same block of memory as this time series.  This is
        very very fast, and makes it easy to quickly use all
        numpy/scipy functionality on self.  However, it is dangerous
        because when this time series goes out of scope and is garbage
        collected, the corresponding numpy reference object will point
        to garbage.

        INPUT:
            copy -- bool (default: True)

        OUTPUT:
            a numpy 1-d array

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-3,4.5,-2])
            sage: w = v.numpy(copy=False); w
            array([ 1. , -3. ,  4.5, -2. ])
            sage: type(w)
            <type 'numpy.ndarray'>
            sage: w.shape
            (4,)

        Notice that changing w also changes v too!
            sage: w[0] = 20
            sage: w
            array([ 20. ,  -3. ,   4.5,  -2. ])
            sage: v
            [20.0000, -3.0000, 4.5000, -2.0000]

        If you want a separate copy do not give the copy=False option.
            sage: z = v.numpy(); z
            array([ 20. ,  -3. ,   4.5,  -2. ])
            sage: z[0] = -10
            sage: v
            [20.0000, -3.0000, 4.5000, -2.0000]
        """
        import_array() #This must be called before using the numpy C/api or you will get segfault
        cdef int dims[1]
        dims[0] = self._length
        cdef ndarray n = PyArray_FromDimsAndData(1, dims, 12, self._values)
        if copy:
            return n.copy()
        else:
            return n

    def randomize(self, distribution='uniform', loc=0, scale=1, **kwds):
        """
        Randomize the entries in this time series, and return a reference
        to self.  Thus this function both changes self in place, and returns
        a copy of self, for convenience.

        INPUT:
            distribution -- 'uniform':    from loc to loc + scale
                            'normal':     mean loc and standard deviation scale
                            'semicircle': with center at loc (scale is ignored)
            loc   -- float (default: 0)
            scale -- float (default: 1)

        NOTE: All random numbers are generated using algorithms that
        build on the high quality GMP random number function
        gmp_urandomb_ui.  Thus this function fully respects the Sage
        set_random_state command.  It's not quite as fast at the C
        library random number generator, but is of course much better
        quality, and is platform independent.

        EXAMPLES:
        We generate 5 uniform random numbers in the interval [0,1]:
            sage: set_random_seed(0)
            sage: finance.TimeSeries(5).randomize()
            [0.8685, 0.2816, 0.0229, 0.1456, 0.7314]

        We generate 5 uniform random numbers from 5 to 5+2=7:
            sage: set_random_seed(0)
            sage: finance.TimeSeries(10).randomize('uniform', 5, 2)
            [6.7369, 5.5632, 5.0459, 5.2911, 6.4628, 5.2412, 5.2010, 5.2761, 5.5813, 5.5439]

        We generate 5 normal random values with mean 0 and variance 1.
            sage: set_random_seed(0)
            sage: finance.TimeSeries(5).randomize('normal')
            [0.6767, -0.4011, 0.3576, -0.5864, -0.9365]

        We generate 10 normal random values with mean 5 and variance 2.
            sage: set_random_seed(0)
            sage: finance.TimeSeries(10).randomize('normal', 5, 2)
            [6.3534, 4.1978, 5.7153, 3.8273, 3.1269, 2.9598, 3.7484, 6.7472, 3.8986, 4.6271]

        We generate 5 values using the semicircle distribution.
            sage: set_random_seed(0)
            sage: finance.TimeSeries(5).randomize('semicircle')
            [0.7369, -0.9541, 0.4628, -0.7990, -0.4187]

        We generate 1 million normal random values and create a frequency histogram.
            sage: set_random_seed(0)
            sage: a = finance.TimeSeries(10^6).randomize('normal')
            sage: a.histogram(10)[0]
            [36, 1148, 19604, 130699, 340054, 347870, 137953, 21290, 1311, 35]

        We take the above values, and compute the proportion that lie within
        1, 2, 3, 5, and 6 standard deviations of the mean (0):
            sage: s = a.standard_deviation()
            sage: len(a.clip_remove(-s,s))/float(len(a))
            0.68309399999999998
            sage: len(a.clip_remove(-2*s,2*s))/float(len(a))
            0.95455900000000005
            sage: len(a.clip_remove(-3*s,3*s))/float(len(a))
            0.997228
            sage: len(a.clip_remove(-5*s,5*s))/float(len(a))
            0.99999800000000005

        There were no "six sigma events":
            sage: len(a.clip_remove(-6*s,6*s))/float(len(a))
            1.0
        """
        if distribution == 'uniform':
            self._randomize_uniform(loc, loc + scale)
        elif distribution == 'normal':
            self._randomize_normal(loc, scale)
        elif distribution == 'semicircle':
            self._randomize_semicircle(loc)
        elif distribution == 'lognormal':
            self._randomize_lognormal(loc, scale)
        else:
            raise NotImplementedError
        return self

    def _randomize_uniform(self, double left, double right):
        """
        Generates a uniform random distribution of doubles between left and
        right and stores values in place.

        INPUT:
            left -- left bound on random distribution
            right -- right bound on random distribution

        EXAMPLES:
        We generate 5 values distributed with respect to the uniform
        distribution over the interval [0,1].
            sage: v = finance.TimeSeries(5)
            sage: set_random_seed(0)
            sage: v.randomize('uniform')
            [0.8685, 0.2816, 0.0229, 0.1456, 0.7314]

        We now test that the mean is indeed 0.5.
            sage: v = finance.TimeSeries(10^6)
            sage: set_random_seed(0)
            sage: v.randomize('uniform').mean()
            0.50069085504319877
        """
        if left >= right:
            raise ValueError, "left must be less than right"

        cdef randstate rstate = current_randstate()
        cdef Py_ssize_t k
        cdef double d = right - left
        for k from 0 <= k < self._length:
            self._values[k] = rstate.c_rand_double() * d + left

    def _randomize_normal(self, double m, double s):
        """
        Generates a normal random distribution of doubles with mean m and
        standard deviation s and stores values in place. Uses the
        Box-Muller algorithm.

        INPUT:
            m -- mean
            s -- standard deviation

        EXAMPLES:
        We generate 5 values distributed with respect to the normal
        distribution with mean 0 and standard deviation 1.
            sage: set_random_seed(0)
            sage: v = finance.TimeSeries(5)
            sage: v.randomize('normal')
            [0.6767, -0.4011, 0.3576, -0.5864, -0.9365]

        We now test that the mean is indeed 0.
            sage: set_random_seed(0)
            sage: v = finance.TimeSeries(10^6)
            sage: v.randomize('normal').mean()
            6.2705472723385207e-05

        The same test with mean equal to 2 and standard deviation equal
        to 5.
            sage: set_random_seed(0)
            sage: v = finance.TimeSeries(10^6)
            sage: v.randomize('normal', 2, 5).mean()
            2.0003135273636117
        """
        # Ported from http://users.tkk.fi/~nbeijar/soft/terrain/source_o2/boxmuller.c
        # This the box muller algorithm.
        cdef randstate rstate = current_randstate()

        cdef double x1, x2, w, y1, y2
        cdef Py_ssize_t k
        for k from 0 <= k < self._length:
            while 1:
                x1 = 2*rstate.c_rand_double() - 1
                x2 = 2*rstate.c_rand_double() - 1
                w = x1*x1 + x2*x2
                if w < 1: break
            w = sqrt( (-2*log(w))/w )
            y1 = x1 * w
            y2 = x2 * w
            self._values[k] = m + y1*s
            k += 1
            if k < self._length:
                self._values[k] = m + y2*s

    def _randomize_semicircle(self, double center):
        """
        Generates a semicircle random distribution of doubles about center
        and stores values in place. Uses the acceptance-rejection method.

        INPUT:
            center -- the center of the semicircle distribution

        EXAMPLES:
        We generate 5 values distributed with respect to the semicircle
        distribution located at center.
            sage: v = finance.TimeSeries(5)
            sage: set_random_seed(0)
            sage: v.randomize('semicircle')
            [0.7369, -0.9541, 0.4628, -0.7990, -0.4187]

        We now test that the mean is indeed the center.
            sage: v = finance.TimeSeries(10^6)
            sage: set_random_seed(0)
            sage: v.randomize('semicircle').mean()
            0.00072074971804614557

        The same test with center equal to 2.
            sage: v = finance.TimeSeries(10^6)
            sage: set_random_seed(0)
            sage: v.randomize('semicircle', 2).mean()
            2.0007207497179227
        """
        cdef Py_ssize_t k
        cdef double x, y, s, d = 2, left = center - 1, z
        cdef randstate rstate = current_randstate()
        z = d*d
        s = 1.5707963267948966192  # pi/2
        for k from 0 <= k < self._length:
            while 1:
                x = rstate.c_rand_double() * d - 1
                y = rstate.c_rand_double() * s
                if y*y + x*x < 1:
                    break
            self._values[k] = x + center

    def _randomize_lognormal(self, double m, double s):
        """
        Generates a log-normal random distribution of doubles with mean m
        and standard deviation s. Usues Box-Muller algorithm and the identity:
        if Y is a random variable with normal distribution then X = exp(Y)
        is a random variable with log-normal distribution.

        INPUT:
            m -- mean
            s -- standard deviation

        EXAMPLES:
        We generate 5 values distributed with respect to the lognormal
        distribution with mean 0 and standard deviation 1.
            sage: set_random_seed(0)
            sage: v = finance.TimeSeries(5)
            sage: v.randomize('lognormal')
            [1.9674, 0.6696, 1.4299, 0.5563, 0.3920]

        We now test that the mean is indeed sqrt(e).
            sage: set_random_seed(0)
            sage: v = finance.TimeSeries(10^6)
            sage: v.randomize('lognormal').mean()
            1.6473519736548801
            sage: e^0.5
            1.648721270700128

        A log-normal distribution can be simply thought of as the logarithm
        of a normally distributed dataset. We test that here by generating
        5 values distributed with respect to the normal distribution with mean
        0 and standard deviation 1.
            sage: set_random_seed(0)
            sage: w = finance.TimeSeries(5)
            sage: w.randomize('normal')
            [0.6767, -0.4011, 0.3576, -0.5864, -0.9365]
            sage: exp(w)
            [1.9674, 0.6696, 1.4299, 0.5563, 0.3920]
        """
        # Ported from http://users.tkk.fi/~nbeijar/soft/terrain/source_o2/boxmuller.c
        # This the box muller algorithm.
        cdef randstate rstate = current_randstate()
        cdef double x1, x2, w, y1, y2
        cdef Py_ssize_t k
        for k from 0 <= k < self._length:
            while 1:
                x1 = 2*rstate.c_rand_double() - 1
                x2 = 2*rstate.c_rand_double() - 1
                w = x1*x1 + x2*x2
                if w < 1: break
            w = sqrt( (-2*log(w))/w )
            y1 = x1 * w
            y2 = x2 * w
            self._values[k] = exp(m + y1*s)
            k += 1
            if k < self._length:
                self._values[k] = exp(m + y2*s)

    def fft(self, bint overwrite=False):
        """
        Return the real discrete fast fourier transform of self, as a
        real time series:

          [y(0),Re(y(1)),Im(y(1)),...,Re(y(n/2))]              if n is even

          [y(0),Re(y(1)),Im(y(1)),...,Re(y(n/2)),Im(y(n/2))]   if n is odd

        where

          y(j) = sum[k=0..n-1] x[k] * exp(-sqrt(-1)*j*k* 2*pi/n)

        for j = 0..n-1.  Note that y(-j) = y(n-j).

        EXAMPLES:
            sage: v = finance.TimeSeries([1..9]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000]
            sage: w = v.fft(); w
            [45.0000, -4.5000, 12.3636, -4.5000, 5.3629, -4.5000, 2.5981, -4.5000, 0.7935]

        We get just the series of real parts of
            sage: finance.TimeSeries([w[0]]) + w[1:].scale_time(2)
            [45.0000, -4.5000, -4.5000, -4.5000, -4.5000]

        An example with an even number of terms.
            sage: v = finance.TimeSeries([1..10]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
            sage: w = v.fft(); w
            [55.0000, -5.0000, 15.3884, -5.0000, 6.8819, -5.0000, 3.6327, -5.0000, 1.6246, -5.0000]
            sage: v.fft().ifft()
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
        """
        import scipy.fftpack
        if overwrite:
            y = self
        else:
            y = self.__copy__()
        w = y.numpy(copy=False)
        scipy.fftpack.rfft(w, overwrite_x=1)
        return y

    def ifft(self, bint overwrite=False):
        """
        Return the real discrete inverse fast fourier transform of
        self, which is also a real time series.

        This is the inverse of fft().

        The returned real array contains

                         [y(0),y(1),...,y(n-1)]

        where for n is even

          y(j) = 1/n (sum[k=1..n/2-1] (x[2*k-1]+sqrt(-1)*x[2*k])
                                       * exp(sqrt(-1)*j*k* 2*pi/n)
                      + c.c. + x[0] + (-1)**(j) x[n-1])

        and for n is odd

          y(j) = 1/n (sum[k=1..(n-1)/2] (x[2*k-1]+sqrt(-1)*x[2*k])
                                       * exp(sqrt(-1)*j*k* 2*pi/n)
                      + c.c. + x[0])

        c.c. denotes complex conjugate of preceeding expression.

        EXAMPLES:
            sage: v = finance.TimeSeries([1..10]); v
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
            sage: v.ifft()
            [5.1000, -5.6876, 1.4764, -1.0774, 0.4249, -0.1000, -0.2249, 0.6663, -1.2764, 1.6988]
            sage: v.ifft().fft()
            [1.0000, 2.0000, 3.0000, 4.0000, 5.0000, 6.0000, 7.0000, 8.0000, 9.0000, 10.0000]
        """
        import scipy.fftpack
        if overwrite:
            y = self
        else:
            y = self.__copy__()
        w = y.numpy(copy=False)
        scipy.fftpack.irfft(w, overwrite_x=1)
        return y


cdef new_time_series(Py_ssize_t length):
    """
    Return a new uninitialized time series of the given length.
    The entries of the time series are garbage.

    INPUT:
        length -- integer
    OUTPUT:
        TimeSeries

    EXAMPLES:
    This uses new_time_series implicitly:
        sage: v = finance.TimeSeries([1,-3,4.5,-2])
        sage: v.__copy__()
        [1.0000, -3.0000, 4.5000, -2.0000]
    """
    if length < 0:
        raise ValueError, "length must be nonnegative"
    cdef TimeSeries t = PY_NEW(TimeSeries)
    t._length = length
    t._values = <double*> sage_malloc(sizeof(double)*length)
    return t

def unpickle_time_series_v1(v, Py_ssize_t n):
    """
    Version 0 unpickle method.

    INPUT:
        v -- a raw char buffer

    EXAMPLES:
        sage: v = finance.TimeSeries([1,2,3])
        sage: s = v.__reduce__()[1][0]; s
        '\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00@\x00\x00\x00\x00\x00\x00\x08@'
        sage: type(s)
        <type 'str'>
        sage: sage.finance.time_series.unpickle_time_series_v1(s,3)
        [1.0000, 2.0000, 3.0000]
        sage: sage.finance.time_series.unpickle_time_series_v1(s+s,6)
        [1.0000, 2.0000, 3.0000, 1.0000, 2.0000, 3.0000]
        sage: sage.finance.time_series.unpickle_time_series_v1('',0)
        []
    """
    cdef TimeSeries t = new_time_series(n)
    memcpy(t._values, PyString_AsString(v), n*sizeof(double))
    return t




def linear_filter(acvs):
    """
    Create a linear filter with given autocovariance sequence.

    EXAMPLES:

    In this example we consider the multifractal cascade random walk
    of length 1000, and use simultations to estimate the
    expected first few autocovariance parameters for this model, then
    use them to construct a linear filter that works vastly better
    than a linear filter constructed from the same data but not using
    this model. The Monte-Carlo method illustrated below should work for
    predicting one "time step" into the future for any
    model that can be simultated.  To predict k time steps into the
    future would require using a similar technique but would require
    scaling time by k.

    We create 100 simultations ofa multifractal random walk.  This
    models the logarithms of a stock price sequence.
        sage: set_random_seed(0)
        sage: y = finance.multifractal_cascade_random_walk_simulation(3700,0.02,0.01,0.01,1000,100)

    For each walk below we replace the walk by the walk but where each
    step size is replaced by its absolute value -- this is what we
    expect to be able to predict given the model, which is only a
    model for predicting volatility.  We compute the first 200
    autocovariance values for every random walk:
        sage: c = [[a.diffs().abs().sums().autocovariance(i) for a in y] for i in range(200)]

    We make a time series out of the expected values of the
    autocovariances:
        sage: ac = finance.TimeSeries([finance.TimeSeries(z).mean() for z in c])
        sage: ac
        [3.9962, 3.9842, 3.9722, 3.9601, 3.9481 ... 1.7144, 1.7033, 1.6922, 1.6812, 1.6701]

    Note: ac looks like a line -- one could best fit it to yield a lot
    more approximate autocovariances.

    We compute the linear filter given by the above autocovariances:
        sage: F = finance.linear_filter(ac); F
        [0.9982, -0.0002, -0.0002, 0.0003, 0.0001 ... 0.0002, -0.0002, -0.0000, -0.0002, -0.0014]

    Note that the sum is close to 1.
        sage: sum(F)
        0.99593284089454...

    Now we make up an 'out of sample' sequence:
        sage: y2 = finance.multifractal_cascade_random_walk_simulation(3700,0.02,0.01,0.01,1000,1)[0].diffs().abs().sums()
        sage: y2
        [0.0013, 0.0059, 0.0066, 0.0068, 0.0184 ... 6.8004, 6.8009, 6.8063, 6.8090, 6.8339]

    And we forecast the very last value using our linear filter; the forecast is close:
        sage: y2[:-1].linear_forecast(F)
        6.7836741372407...

    In fact it is closer than we would get by forecasting using a
    linear filter made from all the autocovariances of our sequence:
        sage: y2[:-1].linear_forecast(y2[:-1].linear_filter(len(y2)))
        6.7701687056683...

    We record the last 20 forecasts, always using all correct values up to the
    one we are forecasting:
        sage: s1 = sum([(y2[:-i].linear_forecast(F)-y2[-i])^2 for i in range(1,20)])

    We do the same, but using the autocovariances of the sample sequence:
        sage: F2 = y2[:-100].linear_filter(len(F))
        sage: s2 = sum([(y2[:-i].linear_forecast(F2)-y2[-i])^2 for i in range(1,20)])

    Our model gives us something that is 15 percent better in this case:
        sage: s2/s1
        1.154646361026...

    How does it compare overall?  To find out we do 100 simulations
    and for each we compute the percent that our model beats naively
    using the autocovariances of the sample:
        sage: y_out = finance.multifractal_cascade_random_walk_simulation(3700,0.02,0.01,0.01,1000,100)
        sage: s1 = []; s2 = []
        sage: for v in y_out:
        ...       s1.append(sum([(v[:-i].linear_forecast(F)-v[-i])^2 for i in range(1,20)]))
        ...       F2 = v[:-len(F)].linear_filter(len(F))
        ...       s2.append(sum([(v[:-i].linear_forecast(F2)-v[-i])^2 for i in range(1,20)]))
        ...

    We find that overall the model beats naive linear forecasting by 35 percent!
        sage: s = finance.TimeSeries([s2[i]/s1[i] for i in range(len(s1))])
        sage: s.mean()
        1.354073591877...
    """
    cdef TimeSeries c
    cdef Py_ssize_t i

    M = len(acvs)-1

    if M <= 0:
        raise ValueError, "M must be positive"

    if not isinstance(acvs, TimeSeries):
        c = TimeSeries(acvs)
    else:
        c = acvs

    # Also create the numpy vector v with entries c(1), ..., c(M).
    v = c[1:].numpy()

    # 2. Create the autocovariances numpy 2d array A whose i,j entry
    # is c(|i-j|).
    import numpy
    A = numpy.array([[c[abs(j-k)] for j in range(M)] for k in range(M)])

    # 3. Solve the equation A * x = v
    x = numpy.linalg.solve(A, v)

    # 4. The entries of x give the linear filter coefficients.
    return TimeSeries(x)
