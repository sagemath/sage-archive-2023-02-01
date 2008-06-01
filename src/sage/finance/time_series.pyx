include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"

cdef extern from "math.h":
    double log(double)
    double exp(double)
    double sqrt(double)
    double floor(double)

cdef extern from "string.h":
    void* memcpy(void* dst, void* src, size_t len)

from sage.rings.integer import Integer

max_print = 10

cdef class TimeSeries:
    cdef double* _values
    cdef Py_ssize_t _length

    def __new__(self, values=None):
        """
        Create new empty uninitialized time series.

        EXAMPLES:
        """
        self._values = NULL

    def __init__(self, values):
        """
        Initialize new time series.

        INPUT:
            values -- integer (number of values) or an iterable of floats

        EXAMPLES:
        """
        if isinstance(values, (int, long, Integer)):
            values = [0.0]*values
        else:
            values = [float(x) for x in values]
        self._length = len(values)
        self._values = <double*> sage_malloc(sizeof(double) * self._length)
        if self._values == NULL:
            raise MemoryError
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            self._values[i] = values[i]

    def  __dealloc__(self):
        """
        Free up memory used by a time series.

        EXAMPLES:
        """
        if self._values:
            sage_free(self._values)

    def __repr__(self):
        return self._repr()

    def _repr(self, prec=4):
        """
        Print representation of a time series.

        EXAMPLES:
        """
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
        return self._length

    def __getitem__(self, i):
        cdef Py_ssize_t start, stop, step
        cdef TimeSeries t
        if PySlice_Check(i):
            start = 0 if (i.start is None) else i.start
            stop = self._length if (i.stop is None) else i.stop
            step = 1 if (i.step is None) else i.step
            if start < 0:
                start += self._length
                if start < 0: start = 0
            if stop < 0:
                stop += self._length
                if stop < 0: stop = 0
            t = new_time_series((stop-start)/step)
            if step > 1:
                for i from 0 <= i < (stop-start)/step:
                    t[i] = self._values[i*step+start]
            else:
                # do a memcopy
                memcpy(t._values, self._values + start, sizeof(double)*t._length)
            return t
        else:
            return self._values[i]

    def __setitem__(self, Py_ssize_t i, double x):
        self._values[i] = x

    def __copy__(self):
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        memcpy(t._values, self._values , sizeof(double)*self._length)
        return t

    def list(self):
        v = [0.0]*self._length
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            v[i] = self._values[i]
        return v

    def log(self):
        """
        Return new time series got by taking the logarithms of all the
        terms in the time series.
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
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = exp(self._values[i])
        return t

    def diffs(self):
        """
        Return the new time series got by taking the differences of
        successive terms in the time series.  So if self is the time
        series $X_0, X_1, X_2, ...$, then this function outputs
        the series $X_1 - X_0, X_2 - X_1, ...$.
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
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length / k)  # in C / is floor division.
        for i from 0 <= i < self._length/k:
            t._values[i] = self._values[i*k]
        return t


    def scale(self, double s):
        """
        Return new time series obtained by multiplying every value in the series by s.
        """
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = self._values[i] * s
        return t

    def sums(self, double s=0):
        """
        Return the new time series got by taking the running partial
        sums of the terms of this time series.

        INPUT:
            s -- starting value for partial sums
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            s += self._values[i]
            t._values[i] = s
        return t

    def plot(self, points=False, **kwds):
        from sage.plot.all import line, point
        v = self.list()
        w = list(enumerate(v))
        if points:
            L = point(w, **kwds)
        else:
            L = line(w, **kwds)
        L.ymin(min(v))
        L.ymax(max(v))
        L.xmin(0)
        L.xmax(len(v))
        return L

    def moving_average(self, Py_ssize_t k):
        """
        Return the moving average time series over the last k time units.
        Assumes the input time series was constant with its starting value
        for negative time.

        INPUT:
            k -- positive integer

        OUTPUT:
            a time series with the same number of steps as self.

        EXAMPLES:
            sage: v = finance.TimeSeries([1,1,1,2,3]); v
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.moving_average(0)
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.moving_average(1)
            [1.0000, 1.0000, 1.0000, 2.0000, 3.0000]
            sage: v.moving_average(2)
            [0.5000, 1.0000, 1.0000, 1.5000, 2.5000]
            sage: v.moving_average(3)
            [0.6667, 0.6667, 1.0000, 1.3333, 2.0000]
        """
        if k == 0:
            return self.__copy__()
        if k <= 0:
            raise ValueError, "k must be positive"
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        if self._length == 0:
            return t

        cdef double s = self._values[0] * (k-1)
        for i from 0 <= i < self._length:
            if i >= k-1:
                s -= self._values[i-k]
            else:
                s -= self._values[0]
            # I have a slight concern about accumulated rounding error given how
            # this algorithm adds and subtracts.
            s += self._values[i]
            t._values[i] = s/k
        return t

    def sum(self):
        cdef double s = 0
        cdef Py_ssize_t i
        for i from 0 <= i < self._length:
            s += self._values[i]
        return s

    def mean(self):
        return self.sum() / self._length

    def variance(self, bias=True):
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

    def standard_deviation(self, bias=True):
        return sqrt(self.variance(bias=bias))

    def min(self):
        """
        Return the smallest value in this time series. If this series
        has length 0 we raise a ValueError

        OUTPUT:
            float -- smallest value
            integer -- index of smallest value
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
        return s, j

    def max(self):
        """
        Return the largest value in this time series. If this series
        has length 0 we raise a ValueError

        OUTPUT:
            float -- largest value
            integer -- index of largest value
        """
        if self._length == 0:
            raise ValueError, "max() arg is an empty sequence"
        cdef Py_ssize_t i, j
        cdef double s = self._values[0]
        for i from 1 <= i < self._length:
            if self._values[i] > s:
                s = self._values[i]
                j = i
        return s, j

    def histogram(self, Py_ssize_t bins=50):
        if bins <= 0:
            raise ValueError, "bins must be positive"

        cdef double mn = self.min()[0], mx = self.max()[0]
        cdef double r = mx - mn + 0.001, step = r/bins

        v = [mn + j*step for j in range(bins)]
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

        counts = [cnts[i] for i in range(bins)]
        sage_free(cnts)

        return counts, v

    def plot_histogram(self, bins=50, **kwds):
        # NOT DONE -- needs to correctly make the x axis.
        from sage.plot.all import bar_chart
        counts, v = self.histogram(bins)
        return bar_chart(counts, **kwds)

    def numpy(self):
        import numpy
        #TODO: make this faster by accessing raw memory?
        return numpy.array(self.list(), dtype=float)


def new_time_series(length):
    """
    Return a new uninitialized time series of the given length.
    The entries of the time series are garbage.

    INPUT:
        length -- integer

    OUTPUT:
        TimeSeries
    """
    cdef TimeSeries t = PY_NEW(TimeSeries)
    t._length = length
    t._values = <double*> sage_malloc(sizeof(double)*length)
    return t

