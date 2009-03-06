#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "../ext/python_list.pxi"


cdef class TripleDict:
    """
    This is a hashtable specifically designed for (read) speed in
    the coercion model.

    It differs from a python dict in the following important ways:

       - All keys must be sequence of exactly three elements. All sequence
         types (tuple, list, etc.) map to the same item.
       - Comparison is done using the 'is' rather than '==' operator.

    There are special cdef set/get methods for faster access.
    It is bare-bones in the sense that not all dictionary methods are
    implemented.

    It is implemented as a list of lists (hereafter called buckets). The bucket
    is chosen according to a very simple hash based on the object pointer.
    and each bucket is of the form [k1, k2, k3, value, k1, k2, k3, value, ...]
    on which a linear search is performed.

    To spread objects evenly, the size should ideally be a prime, and certainly
    not divisible by 2.


    EXAMPLES:

        sage: from sage.structure.coerce_dict import TripleDict
        sage: L = TripleDict(31)
        sage: a = 'a'; b = 'b'; c = 'c'
        sage: L[a,b,c] = 1
        sage: L[a,b,c]
        1
        sage: L[c,b,a] = -1
        sage: list(L.iteritems())     # random order of output.
        [(('c', 'b', 'a'), -1), (('a', 'b', 'c'), 1)]
        sage: del L[a,b,c]
        sage: list(L.iteritems())
        [(('c', 'b', 'a'), -1)]
        sage: len(L)
        1
        sage: L.stats()             # min, avg, max (bucket length)
        (0, 0.032258064516129031, 1)
        sage: L.bucket_lens()       # random layout
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        sage: for i in range(1000):
        ...       L[i,i,i] = i
        sage: len(L)
        1001
        sage: L.stats()             # random
        (31, 32.29032258064516, 35)
        sage: L.bucket_lens()       # random layout
        [33, 34, 32, 32, 35, 32, 31, 33, 34, 31, 32, 34, 32, 31, 31, 32, 32, 31, 31, 33, 32, 32, 32, 33, 33, 33, 31, 33, 33, 32, 31]
        sage: L = TripleDict(101, L)
        sage: L.stats()             # random
        (8, 9.9108910891089117, 11)
        sage: L = TripleDict(3, L)
        sage: L.stats()             # random
        (291, 333.66666666666669, 410)
        sage: L[c,b,a]
        -1
        sage: L[a,b,c]
        Traceback (most recent call last):
        ...
        KeyError: ('a', 'b', 'c')
        sage: L[a]
        Traceback (most recent call last):
        ...
        KeyError: 'a'
        sage: L[a] = 1
        Traceback (most recent call last):
        ...
        KeyError: 'a'

    The following illistrates why even sizes are bad.
        sage: L = TripleDict(4, L)
        sage: L.stats()
        (0, 250.25, 1001)
        sage: L.bucket_lens()
        [1001, 0, 0, 0]


    AUTHOR:
       -- Robert Bradshaw, 2007-08
    """

    def __init__(self, size, data=None, threshold=0):
        """
        Create a special dict using triples for keys.

        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = 1
            sage: L[a,b,c]
            1
        """
        cdef int i
        self.threshold = threshold
        self.buckets = [[] for i from 0 <= i <  size]
        self._size = 0
        if data is not None:
            for k, v in data.iteritems():
                self[k] = v

    def __len__(self):
        """
        The number of items in self.

        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(37)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = 1
            sage: L[a,b,c] = -1 # re-assign
            sage: L[a,c,b] = 1
            sage: L[a,a,None] = None
            sage: len(L)
            3
        """
        return self._size

    def stats(self):
        """
        The distribution of items in buckets.

        OUTPUT:
            (min, avg, max)

        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(37)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.stats() # random
            (2, 2.7027027027027026, 4)

            sage: L = TripleDict(3007)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.stats() # random
            (0, 0.03325573661456601, 1)

            sage: L = TripleDict(1)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.stats()
            (100, 100.0, 100)
        """
        cdef Py_ssize_t size = self._size
        cdef Py_ssize_t cur, min = size, max = 0
        for bucket in self.buckets:
            if bucket:
                cur = len(bucket)/4
                if cur < min: min = cur
                if cur > max: max = cur
            else:
                min = 0
        return min, 1.0*size/len(self.buckets), max

    def bucket_lens(self):
        """
        The distribution of items in buckets.

        OUTPUT:
            A list of how many items are in each bucket.

        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(37)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.bucket_lens() # random
            [3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4]
            sage: sum(L.bucket_lens())
            100

            sage: L = TripleDict(1)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.bucket_lens()
            [100]
        """
        return [len(self.buckets[i])/4 for i from 0 <= i < len(self.buckets)]

    def _get_buckets(self):
        """
        The actual buckets of self, for debugging.

        EXAMPLE:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(3)
            sage: L[0,0,0] = None
            sage: L._get_buckets() # random
            [[0, 0, 0, None], [], []]
        """
        return self.buckets

    def __getitem__(self, k):
        """
        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = 1
            sage: L[a,b,c]
            1
        """
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError, k
        return self.get(k1, k2, k3)

    cdef get(self, k1, k2, k3):
        cdef Py_ssize_t h = (<Py_ssize_t><void *>k1 + 13*<Py_ssize_t><void *>k2 ^ 503*<Py_ssize_t><void *>k3)
        if h < 0: h = -h
        cdef Py_ssize_t i
        bucket = <object>PyList_GET_ITEM(self.buckets, h % PyList_GET_SIZE(self.buckets))
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
            if PyList_GET_ITEM(bucket, i) == <PyObject*>k1 and \
               PyList_GET_ITEM(bucket, i+1) == <PyObject*>k2 and \
               PyList_GET_ITEM(bucket, i+2) == <PyObject*>k3:
                return <object>PyList_GET_ITEM(bucket, i+3)
        raise KeyError, (k1, k2, k3)

    def __setitem__(self, k, value):
        """
        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = -1
            sage: L[a,b,c]
            -1
        """
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError, k
        self.set(k1, k2, k3, value)

    cdef set(self, k1, k2, k3, value):
        if self.threshold and self._size > len(self.buckets) * self.threshold:
            self.resize()
        cdef Py_ssize_t h = (<Py_ssize_t><void *>k1 + 13*<Py_ssize_t><void *>k2 ^ 503*<Py_ssize_t><void *>k3)
        if h < 0: h = -h
        cdef Py_ssize_t i
        bucket = <object>PyList_GET_ITEM(self.buckets, h % PyList_GET_SIZE(self.buckets))
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
            if PyList_GET_ITEM(bucket, i) == <PyObject*>k1 and \
               PyList_GET_ITEM(bucket, i+1) == <PyObject*>k2 and \
               PyList_GET_ITEM(bucket, i+2) == <PyObject*>k3:
                bucket[i+3] = value
                return
        bucket += [k1, k2, k3, value]
        self._size += 1

    def __delitem__(self, k):
        """
        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = -1
            sage: del L[a,b,c]
            sage: len(L)
            0
        """
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError, k
        cdef Py_ssize_t h = (<Py_ssize_t><void *>k1 + 13*<Py_ssize_t><void *>k2 ^ 503*<Py_ssize_t><void *>k3)
        if h < 0: h = -h
        cdef Py_ssize_t i
        bucket = <object>PyList_GET_ITEM(self.buckets, h % PyList_GET_SIZE(self.buckets))
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
            if PyList_GET_ITEM(bucket, i) == <PyObject*>k1 and \
               PyList_GET_ITEM(bucket, i+1) == <PyObject*>k2 and \
               PyList_GET_ITEM(bucket, i+2) == <PyObject*>k3:
                del bucket[i:i+4]
                self._size -= 1
                return
        raise KeyError, k

    def resize(self, int buckets=0):
        """
        Changes the number of buckets of self, while preserving the contents.

        If the number of buckes is 0 or not given, it resizes self to the
        smallest prime that is at least twice as large as self.

        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(8)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.bucket_lens() # random
            [50, 0, 0, 0, 50, 0, 0, 0]
            sage: L.resize(7) # random
            [15, 14, 14, 14, 14, 15, 14]
            sage: L.resize()
            sage: len(L.bucket_lens())
            17
        """
        if buckets == 0:
            buckets = next_odd_prime(2*len(self.buckets))
        cdef TripleDict new = TripleDict(buckets, self)
        self.buckets = new.buckets

    def iteritems(self):
        """
        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: L[1,2,3] = None
            sage: list(L.iteritems())
            [((1, 2, 3), None)]
        """
        return TripleDictIter(self)

    def __reduce__(self):
        """
        Note that we don't expect equality as this class concerns itself with
        object identy rather than object equality.

        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: L[1,2,3] = True
            sage: loads(dumps(L)) == L
            False
            sage: list(loads(dumps(L)).iteritems())
            [((1, 2, 3), True)]
        """
        return TripleDict, (len(self.buckets), dict(self.iteritems()), self.threshold)


cdef class TripleDictIter:
    def __init__(self, pairs):
        """
        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict, TripleDictIter
            sage: L = TripleDict(31)
            sage: L[1,2,3] = None
            sage: L.iteritems().next()
            ((1, 2, 3), None)
        """
        self.pairs = pairs
        self.buckets = iter(self.pairs.buckets)
    def __iter__(self):
        """
        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict, TripleDictIter
            sage: L = TripleDict(31)
            sage: L[1,2,3] = None
            sage: iter(L.iteritems()).next()
            ((1, 2, 3), None)
        """
        return self
    def __next__(self):
        """
        EXAMPLES:
            sage: from sage.structure.coerce_dict import TripleDict, TripleDictIter
            sage: L = TripleDict(31)
            sage: L[1,2,3] = None
            sage: L[3,2,1] = None
            sage: sorted(L.iteritems())
            [((1, 2, 3), None), ((3, 2, 1), None)]
        """
        while self.bucket_iter is None:
            self.bucket_iter = self.buckets.next()
        self.bucket_iter = iter(self.bucket_iter)
        try:
            k1 = self.bucket_iter.next()
            k2 = self.bucket_iter.next()
            k3 = self.bucket_iter.next()
            value = self.bucket_iter.next()
            return ((k1, k2, k3), value)
        except StopIteration:
            self.bucket_iter = None
            return self.next()



cdef long next_odd_prime(long n):
    if n % 2 == 0:
        n -= 1
    cdef long k
    while n > 0:
        n += 2
        k = 3
        while k*k <= n:
            if n % k == 0:
                break
            k += 2
        if k*k > n:
            return n
