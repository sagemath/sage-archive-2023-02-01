include "../ext/cdefs.pxi"
include "../ext/stdsage.pxi"

cdef extern from "math.h":
    double abs(double)
    double exp(double)
    double floor(double)
    double log(double)
    double sqrt(double)

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
            sage: v[:]
            [1.0000, -4.0000, 3.0000, -2.5000, -4.0000, 3.0000]
        """
        cdef Py_ssize_t start, stop, step
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
            t = new_time_series((stop-start)/step)
            if step > 1:
                for i from 0 <= i < (stop-start)/step:
                    t[i] = self._values[i*step+start]
            else:
                # do a memcopy
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
        """
        cdef Py_ssize_t i
        cdef TimeSeries t = new_time_series(self._length)
        for i from 0 <= i < self._length:
            t._values[i] = abs(self._values[i])
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

    def plot(self, points=False, **kwds):
        """
        Return a plot of this time series as a line or points through
        (i,T(i)), where i ranges over nonnegative integers up to the
        length of self.

        INPUT:
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

    def simple_moving_average(self, Py_ssize_t k):
        """
        Return the moving average time series over the last k time units.
        Assumes the input time series was constant with its starting value
        for negative time.  The t-th step of the output is the sum of
        the previous k-1 steps of self and the kth step divided by k.

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
            [0.5000, 1.0000, 1.0000, 1.5000, 2.5000]
            sage: v.simple_moving_average(3)
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

    def histogram(self, Py_ssize_t bins=50):
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
        cdef double r = mx - mn + 0.001, step = r/bins

        v = [(mn + j*step, mn + (j+1)*step) for j in range(bins)]
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
        """
        INPUT:
            bins -- positive integer (default: 50)
            **kwds -- passed to the bar_chart function
        OUTPUT:
            a histogram plot
        """
        from sage.plot.all import bar_chart, polygon
        counts, intervals = self.histogram(bins)
        s = 0
        for i, (x0,x1) in enumerate(intervals):
            s += polygon([(x0,0), (x0,counts[i]), (x1,counts[i]), (x1,0)], **kwds)
        return s
        #return bar_chart(counts, **kwds)

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

