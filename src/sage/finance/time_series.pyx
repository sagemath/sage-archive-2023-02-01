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

cdef extern from "math.h":
    double exp(double)
    double floor(double)
    double log(double)
    double sqrt(double)

cdef extern from "string.h":
    void* memcpy(void* dst, void* src, size_t len)

cdef extern from "stdlib.h":
    long random()
    long RAND_MAX

cdef inline float ranf():
    """
    Return random number uniformly distributed in 0..1.
    """
    return (<float> random())/RAND_MAX

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
        """
        cdef RealDoubleVectorSpaceElement z
        if isinstance(values, (int, long, Integer)):
            values = [0.0]*values
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
            right -- iterable that can be converted to a time series
        OUTPUT:
            a time series

        EXAMPLES:
            sage: v = finance.TimeSeries([1,2,3]); w = finance.TimeSeries([1,2])
            sage: v + w
            [1.0000, 2.0000, 3.0000, 1.0000, 2.0000]
            sage: v + xrange(4)
            [1.0000, 2.0000, 3.0000, 0.0000, 1.0000, 2.0000, 3.0000]
            sage: v = finance.TimeSeries([1,2,-5]); v
            [1.0000, 2.0000, -5.0000]
            sage: [1,5] + v
            [1.0000, 5.0000, 1.0000, 2.0000, -5.0000]
        """
        if not isinstance(right, TimeSeries):
            raise TypeError, "right must be a time series"
        if not isinstance(left, TimeSeries):
            raise TypeError, "right must be a time series"
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

    def diffs(self):
        """
        Return the new time series got by taking the differences of
        successive terms in the time series.  So if self is the time
        series $X_0, X_1, X_2, ...$, then this function outputs the
        series $X_1 - X_0, X_2 - X_1, ...$.  The output series has one
        less term than the input series.

        OUTPUT:
            a new time series.

        EXAMPLES:
            sage: v = finance.TimeSeries([5,4,1.3,2,8]); v
            [5.0000, 4.0000, 1.3000, 2.0000, 8.0000]
            sage: v.diffs()
            [-1.0000, -2.7000, 0.7000, 6.0000]
        """
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

        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length / k)  # in C / is floor division.
        for i from 0 <= i < self._length/k:
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
            v = self.scale_time(s).list()
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

    def exponential_moving_average(self, double alpha, double t0 = 0):
        """
        Return the exponential moving average time series over the
        last .  Assumes the input time series was constant
        with its starting value for negative time.  The t-th step of
        the output is the sum of the previous k-1 steps of self and
        the kth step divided by k.

        INPUT:
            alpha -- float; a smoothing factor with 0 <= alpha <= 1.

        OUTPUT:
            a time series with the same number of steps as self.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.exponential_moving_average(0,0)
            [0.0000, 1.0000, 1.0000, 1.0000, 1.0000]
            sage: v.exponential_moving_average(1,0)
            [0.0000, 1.0000, 1.0000, 1.0000, 2.0000]
            sage: v.exponential_moving_average(0.5,0)
            [0.0000, 1.0000, 1.0000, 1.0000, 1.5000]
        """
        if alpha < 0 or alpha > 1:
            raise ValueError, "alpha must be between 0 and 1"
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        if self._length == 0:
            return t
        t._values[0] = t0
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
            ([1, 5, 2],
             [(-5.0, 0.00033333333333285253),
              (0.00033333333333285253, 5.0006666666666657),
              (5.0006666666666657, 10.000999999999998)])
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

        for i from 0 <= i < self._length:
            cnts[<Py_ssize_t>floor((self._values[i] - mn)/step)] += 1

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
            normalize -- whether to normalize so the total area in the
                         bars of the histogram is 1.
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
        return s

    def numpy(self):
        """
        Return numpy 1-d array corresponding to this time series.

        OUTPUT:
            a numpy 1-d array

        EXAMPLES:
            sage: v = finance.TimeSeries([1,-3,4.5,-2])
            sage: v.numpy()
            array([ 1. , -3. ,  4.5, -2. ])
        """
        import numpy
        #TODO: make this faster by accessing raw memory?
        return numpy.array(self.list(), dtype=float)

    def randomize(self, distribution='uniform', loc=0, scale=1, **kwds):
        """
        INPUT:
            distribution -- 'uniform'
                            'normal'
                            'semicircle'
        """
        if distribution == 'uniform':
            self._randomize_uniform(loc, scale)
        elif distribution == 'normal':
            self._randomize_normal(loc, scale)
        elif distribution == 'semicircle':
            self._randomize_semicircle(loc)
        else:
            raise NotImplementedError

    def _randomize_uniform(self, double left, double right):
        if left >= right:
            raise ValueError, "left must be less than right"
        cdef Py_ssize_t k
        cdef double d = right - left
        for k from 0 <= k < self._length:
            self._values[k] = ranf() * d - left

    def _randomize_normal(self, double m, double s):
        """
        INPUT:
            m -- mean
            s -- standard deviation

        EXAMPLES:
        """
        # Ported from http://users.tkk.fi/~nbeijar/soft/terrain/source_o2/boxmuller.c
        # This the box muller algorithm.
        cdef double x1, x2, w, y1, y2
        cdef Py_ssize_t k
        for k from 0 <= k < self._length:
            while 1:
                x1 = 2.0 * ranf() - 1.0
                x2 = 2.0 * ranf() - 1.0
                w = x1 * x1 + x2 * x2
                if w < 1.0: break
            w = sqrt( (-2.0 * log( w ) ) / w )
            y1 = x1 * w
            y2 = x2 * w
            self._values[k] = m + y1*s
            k += 1
            if k < self._length:
                self._values[k] = m + y2*s

    def _randomize_semicircle(self, double center):
        cdef Py_ssize_t k
        cdef double x, y, s, d = 2, left = center - 1, z
        z = d*d
        s = 1.5707963267948966192  # pi/2
        for k from 0 <= k < self._length:
            while 1:
                x = ranf() * d - 1
                y = ranf() * s
                if y*y + x*x < 1:
                    break
            self._values[k] = x + center

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

