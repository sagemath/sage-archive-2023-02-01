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
       - A value of None is treated as a non-existant value.

    In addition, there is the following difference which (unlike the above
    three) should probably be changed.

       - Its size is fixed at creation time, so no load-adjusting parameters
         are in place. One can re-size it by creating a new TripleDict
         passing in self as a parameter.

    There are special cdef set/get methods for faster access.
    It is bare-bones in the sense that not all dictionary methods are
    implemented.

    It is implemented as a list of lists (called buckets). The bucket
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
        sage: list(L.iteritems())
        [(('c', 'b', 'a'), -1), (('a', 'b', 'c'), 1)]
        sage: del L[a,b,c]
        sage: list(L.iteritems())
        [(('c', 'b', 'a'), -1)]
        sage: len(L)
        1
        sage: L.stats()             # random -- min, avg, max (bucket length)
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
        sage: L.stats()
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
        sage: L = TripleDict(4, L)
        sage: L.stats()
        (0, 250.25, 1001)
        sage: L.bucket_lens()
        [1001, 0, 0, 0]


    AUTHOR:
       -- Robert Bradshaw, 2007-08
    """

    def __init__(self, size, data=None):
        cdef int i
        self.buckets = [[] for i from 0 <= i <  size]
        if data is not None:
            for k, v in data.iteritems():
                self[k] = v

    def __len__(self):
        cdef Py_ssize_t size = 0
        for bucket in self.buckets:
            if bucket:
                size += len(bucket)/4
        return size

    def stats(self):
        cdef Py_ssize_t size = len(self)
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
        return [len(self.buckets[i])/4 for i from 0 <= i < len(self.buckets)]

    def _get_buckets(self):
        return self.buckets

    def __getitem__(self, k):
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError, k
        value = self.get(k1, k2, k3)
        if value is None:
            raise KeyError, k
        else:
            return value

    cdef get(self, k1, k2, k3):
        cdef Py_ssize_t h = (<Py_ssize_t>k1 + 13*<Py_ssize_t>k2 + 503*<Py_ssize_t>k3)
        if h < 0: h = -h
        cdef Py_ssize_t i
        bucket = <object>PyList_GET_ITEM(self.buckets, h % PyList_GET_SIZE(self.buckets))
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
            if PyList_GET_ITEM(bucket, i) == <PyObject*>k1 and \
               PyList_GET_ITEM(bucket, i+1) == <PyObject*>k2 and \
               PyList_GET_ITEM(bucket, i+2) == <PyObject*>k3:
                return <object>PyList_GET_ITEM(bucket, i+3)
        return None

    def __setitem__(self, k, value):
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError, k
        if value is None:
            # for consistency with get
            try:
                del self[k]
            except KeyError:
                pass
        else:
            self.set(k1, k2, k3, value)

    cdef set(self, k1, k2, k3, value):
        cdef Py_ssize_t h = (<Py_ssize_t>k1 + 13*<Py_ssize_t>k2 + 503*<Py_ssize_t>k3)
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

    def __delitem__(self, k):
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError, k
        cdef Py_ssize_t h = (<Py_ssize_t>k1 + 13*<Py_ssize_t>k2 + 503*<Py_ssize_t>k3)
        if h < 0: h = -h
        cdef Py_ssize_t i
        bucket = <object>PyList_GET_ITEM(self.buckets, h % PyList_GET_SIZE(self.buckets))
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
            if PyList_GET_ITEM(bucket, i) == <PyObject*>k1 and \
               PyList_GET_ITEM(bucket, i+1) == <PyObject*>k2 and \
               PyList_GET_ITEM(bucket, i+2) == <PyObject*>k3:
                del bucket[i:i+4]
                return
        raise KeyError, k

    def iteritems(self):
        return TripleDictIter(self)


cdef class TripleDictIter:
    def __init__(self, pairs):
        self.pairs = pairs
        self.buckets = iter(self.pairs.buckets)
    def __iter__(self):
        return self
    def __next__(self):
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



