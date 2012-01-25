#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#                     2012 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
"""
Containers for storing coercion data

This module provides :class:`TripleDict` and :class:`MonoDict`. These are
structures similar to ``WeakKeyDictionary`` in Python's weakref module,
and are optimized for lookup speed. The keys for :class:`TripleDict` consist
of triples (k1,k2,k3) and are looked up by identity rather than equality. The
keys are stored by weakrefs if possible. If any one of the components k1, k2,
k3 gets garbage collected, then the entry is removed from the :class:`TripleDict`.

Key components that do not allow for weakrefs are stored via a normal
refcounted reference. That means that any entry stored using a triple
(k1,k2,k3) so that none of the k1,k2,k3 allows a weak reference behaves
as an entry in a normal dictionary: Its existence in :class:`TripleDict`
prevents it from being garbage collected.

That container currently is used to store coercion and conversion maps
between two parents (:trac:`715`) and to store homsets of pairs of objects
of a category (:trac:`11521`). In both cases, it is essential that the parent
structures remain garbage collectable, it is essential that the data access
is faster than with a usual ``WeakKeyDictionary``, and we enforce the "unique
parent condition" in Sage (parent structures should be identical if they are
equal).

:class:`MonoDict` behaves similarly, but it takes a single item as a key. It
is used for caching the parents which allow a coercion map into a fixed other
parent (:trac:`12313`).

"""
include "../ext/python_list.pxi"

from weakref import KeyedRef, ref

############################################
# The following code is responsible for
# removing dead references from the cache
############################################

cdef class MonoDictEraser:
    """
    Erase items from a :class:`MonoDict` when a weak reference becomes
    invalid.

    This is of internal use only. Instances of this class will be passed as a
    callback function when creating a weak reference.

    EXAMPLES::

        sage: from sage.structure.coerce_dict import MonoDict
        sage: class A: pass
        sage: a = A()
        sage: M = MonoDict(11)
        sage: M[a] = 1
        sage: len(M)
        1
        sage: del a
        sage: import gc
        sage: n = gc.collect()
        sage: len(M)    # indirect doctest
        0

    AUTHOR:

    - Simon King (2012-01)
    """
    def __init__(self, D):
        """
        INPUT:

        A :class:`MonoDict`.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict, MonoDictEraser
            sage: D = MonoDict(11)
            sage: MonoDictEraser(D)
            <sage.structure.coerce_dict.MonoDictEraser object at ...>
        """
        self.D = ref(D)

    def __call__(self, r):
        """
        INPUT:

        A weak reference with key.

        When this is called with a weak reference ``r``, then each item
        containing ``r`` is removed from the associated :class:`MonoDict`.
        Normally, this only happens when a weak reference becomes invalid.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: class A: pass
            sage: a = A()
            sage: M = MonoDict(11)
            sage: M[a] = 1
            sage: len(M)
            1
            sage: del a
            sage: import gc
            sage: n = gc.collect()
            sage: len(M)    # indirect doctest
            0
        """
        # r is a (weak) reference (typically to a parent),
        # and it knows the stored key of the unique singleton r() had been part
        # of.
        # We remove that unique tuple from self.D -- if self.D is still there!
        cdef MonoDict D = self.D()
        if D is None:
            return
        cdef list buckets = D.buckets
        if buckets is None:
            return
        cdef size_t k = r.key
        cdef size_t h = k % PyList_GET_SIZE(buckets)
        cdef list bucket = <object>PyList_GET_ITEM(buckets,h)
        cdef Py_ssize_t i
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 2:
            if <size_t><object>PyList_GET_ITEM(bucket,i)==k:
                del bucket[i:i+2]
                D._size -= 1
                break
        try:
            D._refcache.__delitem__(k)
        except KeyError:
            pass

cdef class TripleDictEraser(MonoDictEraser):
    """
    Erases items from a :class:`TripleDict` when a weak reference becomes
    invalid.

    This is of internal use only. Instances of this class will be passed as a
    callback function when creating a weak reference.

    EXAMPLES::

        sage: from sage.structure.coerce_dict import TripleDict
        sage: class A: pass
        sage: a = A()
        sage: T = TripleDict(11)
        sage: T[a,ZZ,None] = 1
        sage: T[ZZ,a,1] = 2
        sage: T[a,a,ZZ] = 3
        sage: len(T)
        3
        sage: del a
        sage: import gc
        sage: n = gc.collect()
        sage: len(T)
        0

    AUTHOR:

    - Simon King (2012-01)
    """
    def __init__(self, D):
        """
        INPUT:

        A :class:`TripleDict`.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict, TripleDictEraser
            sage: D = TripleDict(11)
            sage: TripleDictEraser(D)
            <sage.structure.coerce_dict.TripleDictEraser object at ...>

        """
        self.D = ref(D)

    def __call__(self, r):
        """
        INPUT:

        A weak reference with key.

        When this is called with a weak reference ``r``, then each item
        containing ``r`` is removed from the associated :class:`TripleDict`.
        Normally, this only happens when a weak reference becomes invalid.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: class A: pass
            sage: a = A()
            sage: T = TripleDict(11)
            sage: T[a,ZZ,None] = 1
            sage: T[ZZ,a,1] = 2
            sage: T[a,a,ZZ] = 3
            sage: len(T)
            3
            sage: del a
            sage: import gc
            sage: n = gc.collect()
            sage: len(T)    # indirect doctest
            0
        """
        cdef TripleDict D = self.D()
        if D is None:
            return
        cdef list buckets = D.buckets
        if buckets is None:
            return
        # r is a (weak) reference (typically to a parent), and it knows the
        # stored key of the unique triple r() had been part of.
        # We remove that unique triple from self.D
        cdef size_t k1,k2,k3
        k1,k2,k3 = r.key
        cdef size_t h = (k1 + 13*k2 ^ 503*k3)
        cdef list bucket = <object>PyList_GET_ITEM(buckets, h % PyList_GET_SIZE(buckets))
        cdef Py_ssize_t i
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
            if <size_t><object>PyList_GET_ITEM(bucket, i)==k1 and \
               <size_t><object>PyList_GET_ITEM(bucket, i+1)==k2 and \
               <size_t><object>PyList_GET_ITEM(bucket, i+2)==k3:
                del bucket[i:i+4]
                D._size -= 1
                break
        try:
            D._refcache.__delitem__((k1,k2,k3))
        except KeyError:
            pass

cdef class MonoDict:
    """
    This is a hashtable specifically designed for (read) speed in
    the coercion model.

    It differs from a python dict in the following important way:

       - Comparison is done using the 'is' rather than '==' operator.
       - There are only weak references (if possible) to the keys.

    There are special cdef set/get methods for faster access.
    It is bare-bones in the sense that not all dictionary methods are
    implemented.

    It is implemented as a list of lists (hereafter called buckets). The bucket
    is chosen according to a very simple hash based on the object pointer,
    and each bucket is of the form [id(k1), value1, id(k2), value2, ...], on
    which a linear search is performed.

    To spread objects evenly, the size should ideally be a prime, and certainly
    not divisible by 2.

    EXAMPLES::

        sage: from sage.structure.coerce_dict import MonoDict
        sage: L = MonoDict(31)
        sage: a = 'a'; b = 'ab'; c = -15
        sage: L[a] = 1
        sage: L[b] = 2
        sage: L[c] = 3

    The key is expected to be a unique object. Hence, the item stored for ``c``
    can not be obtained by providing another equal number::

        sage: L[a]
        1
        sage: L[b]
        2
        sage: L[c]
        3
        sage: L[-15]
        Traceback (most recent call last):
        ...
        KeyError: -15

    Not all features of Python dictionaries are available, but iteration over
    the dictionary items is possible::

        sage: # for some reason the following fails in "make ptest"
        sage: # on some installations, see #12313 for details
        sage: sorted(L.iteritems()) # random layout
        [(-15, 3), ('a', 1), ('ab', 2)]
        sage: # the following seems to be more consistent
        sage: set(L.iteritems())
        set([('a', 1), ('ab', 2), (-15, 3)])
        sage: del L[c]
        sage: sorted(L.iteritems())
        [('a', 1), ('ab', 2)]
        sage: len(L)
        2
        sage: L.stats()             # random
        (0, 0.06451..., 1)
        sage: L.bucket_lens()       # random layout
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
        sage: for i in range(1000):
        ...       L[i] = i
        sage: len(L)
        1002
        sage: L.stats()             # random
        (26, 32.32258064516129, 37)
        sage: L.bucket_lens()       # random layout
        [32, 34, 31, 32, 33, 32, 32, 33, 34, 31, 32, 32, 31, 31, 32, 34, 31, 33, 34, 32, 32, 33, 33, 31, 33, 35, 32, 32, 32, 32, 31]
        sage: L = MonoDict(101, L)
        sage: L.stats()             # random
        (8, 9.92079207920792, 12)
        sage: L = MonoDict(3, L)
        sage: L.stats()             # random
        (0, 334.0, 985)
        sage: L['a']
        1
        sage: L['c']
        Traceback (most recent call last):
        ...
        KeyError: 'c'

    The following illustrates why even sizes are bad::

        sage: L = MonoDict(4, L)
        sage: L.stats()
        (0, 250.5, 1002)
        sage: L.bucket_lens()
        [1002, 0, 0, 0]

    Note that this kind of dictionary is also used for caching actions
    and coerce maps. In previous versions of Sage, the cache was by
    strong references and resulted in a memory leak in the following
    example. However, this leak was fixed by :trac:`715`, using
    weak references::

        sage: K = GF(1<<55,'t')
        sage: for i in range(50):
        ...     a = K.random_element()
        ...     E = EllipticCurve(j=a)
        ...     P = E.random_point()
        ...     Q = 2*P
        sage: import gc
        sage: n = gc.collect()
        sage: from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
        sage: LE = [x for x in gc.get_objects() if isinstance(x, EllipticCurve_finite_field)]
        sage: len(LE)    # indirect doctest
        1

    AUTHOR:

    - Simon King (2012-01)
    """
    def __init__(self, size, data=None, threshold=0):
        """
        Create a special dict using singletons for keys.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: a = 'a'
            sage: L[a] = 1
            sage: L[a]
            1
        """
        cdef Py_ssize_t i
        self.threshold = threshold
        self.buckets = [[] for i from 0 <= i < size]
        self._size = 0
        self.eraser = MonoDictEraser(self)
        self._refcache = {}
        if data is not None:
            for k, v in data.iteritems():
                self[k] = v

    def __len__(self):
        """
        The number of items in self.
        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(37)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a] = 1
            sage: L[a] = -1 # re-assign
            sage: L[b] = 1
            sage: L[c] = None
            sage: len(L)
            3
        """
        return self._size

    def stats(self):
        """
        The distribution of items in buckets.

        OUTPUT:

        - (min, avg, max)

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(37)
            sage: for i in range(100): L[i] = None
            sage: L.stats() # random
            (2, 2.7027027027027026, 4)

            sage: L = MonoDict(3007)
            sage: for i in range(100): L[i] = None
            sage: L.stats() # random
            (0, 0.03325573661456601, 1)

            sage: L = MonoDict(1)
            sage: for i in range(100): L[i] = None
            sage: L.stats()
            (100, 100.0, 100)
        """
        cdef Py_ssize_t size = self._size
        cdef Py_ssize_t cur, min = size, max = 0
        for bucket in self.buckets:
            if bucket:
                cur = len(bucket)/2
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

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(37)
            sage: for i in range(100): L[i] = None
            sage: L.bucket_lens() # random
            [3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4]
            sage: sum(L.bucket_lens())
            100

            sage: L = MonoDict(1)
            sage: for i in range(100): L[i] = None
            sage: L.bucket_lens()
            [100]
        """
        return [len(self.buckets[i])/2 for i from 0 <= i < len(self.buckets)]

    def _get_buckets(self):
        """
        The actual buckets of self, for debugging.

        EXAMPLE::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(3)
            sage: L[0,0,0] = None
            sage: L._get_buckets() # random
            [[0, 0, 0, None], [], []]
        """
        return self.buckets

    def __contains__(self, k):
        """
        Test if the dictionary contains a given key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: a = 'a'; b = 'ab'; c = 15
            sage: L[a] = 1
            sage: L[b] = 2
            sage: L[c] = 3
            sage: c in L         # indirect doctest
            True

        The keys are compared by identity, not by equality. Hence, we have::

            sage: c == 15
            True
            sage: 15 in L
            False
        """
        cdef size_t h = <size_t><void *>k
        cdef object r
        try:
            r = self._refcache[h]
        except KeyError:
            return False
        return not (isinstance(r,KeyedRef) and r() is None)

    def __getitem__(self, k):
        """
        Get the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: a = 'a'; b = 'b'; c = 15
            sage: L[a] = 1
            sage: L[b] = 2
            sage: L[c] = 3
            sage: L[c]                  # indirect doctest
            3

        Note that the keys are supposed to be unique::

            sage: c==15
            True
            sage: c is 15
            False
            sage: L[15]
            Traceback (most recent call last):
            ...
            KeyError: 15
        """
        cdef size_t h = <size_t><void *>k
        cdef object r
        try:
            r = self._refcache[h]
        except KeyError:
            raise KeyError, k
        if isinstance(r, KeyedRef) and r() is None:
            raise KeyError, k
        cdef Py_ssize_t i
        cdef list all_buckets = self.buckets
        cdef list bucket = <object>PyList_GET_ITEM(all_buckets, h % PyList_GET_SIZE(all_buckets))
        cdef object tmp
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 2:
            tmp = <object>PyList_GET_ITEM(bucket, i)
            if <size_t>tmp == h:
                return <object>PyList_GET_ITEM(bucket, i+1)
        raise KeyError, k

    def __setitem__(self, k, value):
        """
        Set the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: a = 'a'
            sage: L[a] = -1   # indirect doctest
            sage: L[a]
            -1
            sage: L[a] = 1
            sage: L[a]
            1
            sage: len(L)
            1
        """
        if self.threshold and self._size > len(self.buckets) * self.threshold:
            self.resize()
        cdef size_t h = <size_t><void *>k
        cdef Py_ssize_t i
        cdef list bucket = <object>PyList_GET_ITEM(self.buckets, h % PyList_GET_SIZE(self.buckets))
        cdef object tmp
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 2:
            tmp = <object>PyList_GET_ITEM(bucket, i)
            if <size_t>tmp == h:
                r = self._refcache[h]
                if isinstance(r, KeyedRef) and r() is None:
                    del bucket[i:i+2]
                    self._size -= 1
                    break
                bucket[i+1] = value
                return
        PyList_Append(bucket, h)
        PyList_Append(bucket, value)
        try:
            self._refcache[h] = KeyedRef(k,self.eraser,h)
        except TypeError:
            self._refcache[h] = k
        self._size += 1

    def __delitem__(self, k):
        """
        Delete the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: a = 15
            sage: L[a] = -1
            sage: len(L)
            1

        Note that the keys are unique, hence using a key that is equal but not
        identical to a results in an error::

            sage: del L[15]
            Traceback (most recent call last):
            ...
            KeyError: 15
            sage: a in L
            True
            sage: del L[a]
            sage: len(L)
            0
            sage: a in L
            False
        """
        cdef size_t h = <size_t><void *>k
        try:
            r = self._refcache[h]
        except KeyError:
            raise KeyError, k
        try:
            del self._refcache[h]
        except KeyError:
            raise KeyError, k
        if isinstance(r,KeyedRef) and r() is None:
            raise KeyError, k
        cdef Py_ssize_t i
        cdef object tmp
        cdef list bucket = <object>PyList_GET_ITEM(self.buckets, h % PyList_GET_SIZE(self.buckets))
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 2:
            tmp = <object>PyList_GET_ITEM(bucket, i)
            if <size_t>tmp == h:
                del bucket[i:i+2]
                self._size -= 1
                return

    def resize(self, int buckets=0):
        """
        Change the number of buckets of self, while preserving the contents.

        If the number of buckets is 0 or not given, it resizes self to the
        smallest prime that is at least twice as large as self.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(8)
            sage: for i in range(100): L[i] = None
            sage: L.bucket_lens() # random
            [50, 0, 0, 0, 50, 0, 0, 0]
            sage: L.resize(7)
            sage: L.bucket_lens() # random
            [15, 14, 14, 14, 14, 15, 14]
            sage: L.resize()
            sage: len(L.bucket_lens())
            17
        """
        if buckets == 0:
            buckets = next_odd_prime(2*len(self.buckets))
        cdef list old_buckets = self.buckets
        cdef list bucket
        cdef Py_ssize_t i
        cdef size_t k
        self.buckets = [[] for i from 0 <= i <  buckets]
        cdef object v
        for bucket in old_buckets:
            for i from 0 <= i < PyList_GET_SIZE(bucket) by 2:
                k = <size_t><object>PyList_GET_ITEM(bucket,i)
                v  = <object>PyList_GET_ITEM(bucket,i+1)
                self.buckets[k % buckets] += [k,v]

    def iteritems(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: L[1] = None
            sage: L[2] = True
            sage: list(sorted(L.iteritems()))
            [(1, None), (2, True)]
        """
        cdef list bucket
        cdef size_t i, h
        # We test whether the references are still valid.
        # However, we must not delete them, since we are
        # iterating.
        for bucket in self.buckets:
            for i from 0<=i<len(bucket) by 2:
                h = <object>PyList_GET_ITEM(bucket,i)
                try:
                    r = self._refcache[h]
                except KeyError:
                    # That should not happen, unless there is a nasty
                    # race condition. Anyway, the item isn't there,
                    # hence:
                    continue
                if isinstance(r, KeyedRef):
                    r = r()
                    if r is None:
                        continue
                yield r, <object>PyList_GET_ITEM(bucket,i+1)

    def __reduce__(self):
        """
        Note that we don't expect equality as this class concerns itself with
        object identity rather than object equality.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import MonoDict
            sage: L = MonoDict(31)
            sage: L[1] = True
            sage: loads(dumps(L)) == L
            False
            sage: list(loads(dumps(L)).iteritems())
            [(1, True)]
        """
        return MonoDict, (len(self.buckets), dict(self.iteritems()), self.threshold)

cdef class TripleDict(MonoDict):
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
    is chosen according to a very simple hash based on the object pointer,
    and each bucket is of the form [id(k1), id(k2), id(k3), value, id(k1),
    id(k2), id(k3), value, ...], on which a linear search is performed.

    To spread objects evenly, the size should ideally be a prime, and certainly
    not divisible by 2.

    EXAMPLES::

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
        (0, 0.03225806451612903, 1)
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

    The following illustrates why even sizes are bad (setting the threshold
    zero, so that no beneficial resizing happens)::

        sage: L = TripleDict(4, L, threshold=0)
        sage: L.stats()
        (0, 250.25, 1001)
        sage: L.bucket_lens()
        [1001, 0, 0, 0]

    Note that this kind of dictionary is also used for caching actions
    and coerce maps. In previous versions of Sage, the cache was by
    strong references and resulted in a memory leak in the following
    example. However, this leak was fixed by :trac:`715`, using weak
    references::

        sage: K = GF(1<<55,'t')
        sage: for i in range(50):
        ...     a = K.random_element()
        ...     E = EllipticCurve(j=a)
        ...     P = E.random_point()
        ...     Q = 2*P
        sage: import gc
        sage: n = gc.collect()
        sage: from sage.schemes.elliptic_curves.ell_finite_field import EllipticCurve_finite_field
        sage: LE = [x for x in gc.get_objects() if isinstance(x, EllipticCurve_finite_field)]
        sage: len(LE)    # indirect doctest
        1

    .. NOTE::

        The index `h` corresponding to the key [k1, k2, k3] is computed as a
        value of unsigned type size_t as follows:

        .. MATH::

            h = id(k1) + 13*id(k2) \oplus 503 id(k3)

        Indeed, although the PyList_GetItem function and corresponding
        PyList_GET_ITEM macro take a value of signed type Py_ssize_t as input
        for the index, they do not accept negative inputs as the higher level
        Python functions. Moreover, the above formula can overflow so that `h`
        might be considered as negative. Even though this value is taken
        modulo the size of the buckets' list before accessing the corresponding
        item, the Cython "%" operator behaves for values of type size_t and
        Py_ssize_t like the C "%" operator, rather than like the Python "%"
        operator as it does for values of type int. That is, it returns a
        result of the same sign as its input. Therefore, if `h` was defined as
        a signed value, we might access the list at a negative index and raise
        a segfault (and this has been observed on 32 bits systems, see
        :trac:`715` for details).

    AUTHORS:

    - Robert Bradshaw, 2007-08

    - Simon King, 2012-01
    """

    def __init__(self, size, data=None, threshold=0.7):
        """
        Create a special dict using triples for keys.

        EXAMPLES::

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
        self.eraser = TripleDictEraser(self)
        self._refcache = {}
        if data is not None:
            for (k1,k2,k3), v in data.iteritems():
                self.set(k1,k2,k3, v)

    def __len__(self):
        """
        The number of items in self.

        EXAMPLES::

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

        - (min, avg, max)

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(37)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.stats() # random
            (2, 2.7027027027027026, 4)

            sage: L = TripleDict(3007)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.stats() # random
            (0, 0.03325573661456601, 1)

        In order to have a test that isn't random, we use parameters
        that should not be used in real applications::

            sage: L = TripleDict(1, threshold=0)
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

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(37, threshold=0)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.bucket_lens() # random
            [3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4]
            sage: sum(L.bucket_lens())
            100

        In order to have a test that isn't random, we use parameters
        that should not be used in real applications::

            sage: L = TripleDict(1, threshold=0)
            sage: for i in range(100): L[i,i,i] = None
            sage: L.bucket_lens()
            [100]
        """
        return [len(self.buckets[i])/4 for i from 0 <= i < len(self.buckets)]

    def __contains__(self, k):
        """
        Test if the dictionary contains a given key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'ab'; c = 15
            sage: L[a,b,c] = 123
            sage: (a,b,c) in L         # indirect doctest
            True

        The keys are compared by identity, not by equality. Hence, we have::

            sage: c == 15
            True
            sage: (a,b,15) in L
            False
        """
        cdef object k1,k2,k3
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            return False
        cdef object r1,r2,r3
        try:
            r1,r2,r3 = self._refcache[<size_t><void *>k1,<size_t><void *>k2,<size_t><void *>k3]
        except KeyError:
            return False
        if (isinstance(r1,KeyedRef) and r1() is None) or \
           (isinstance(r2,KeyedRef) and r2() is None) or \
           (isinstance(r3,KeyedRef) and r3() is None):
            return False
        return True

    def __getitem__(self, k):
        """
        Get the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = 1
            sage: L[a,b,c]
            1
        """
        cdef object k1,k2,k3
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError, k
        return self.get(k1, k2, k3)

    cdef get(self, object k1, object k2, object k3):
        cdef size_t h1,h2,h3
        h1 = <size_t><void *>k1
        h2 = <size_t><void *>k2
        h3 = <size_t><void *>k3
        cdef object r1,r2,r3
        try:
            r1,r2,r3 = self._refcache[h1,h2,h3]
        except KeyError:
            raise KeyError, (k1,k2,k3)
        if (isinstance(r1,KeyedRef) and r1() is None) or \
           (isinstance(r2,KeyedRef) and r2() is None) or \
           (isinstance(r3,KeyedRef) and r3() is None):
            raise KeyError, (k1,k2,k3)
        cdef size_t h = (h1 + 13*h2 ^ 503*h3)
        cdef Py_ssize_t i
        cdef list all_buckets = self.buckets
        cdef list bucket = <object>PyList_GET_ITEM(all_buckets, h % PyList_GET_SIZE(all_buckets))
        cdef object tmp
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
            tmp = <object>PyList_GET_ITEM(bucket, i)
            if <size_t>tmp == <size_t><void *>k1:
                tmp = <object>PyList_GET_ITEM(bucket, i+1)
                if <size_t>tmp == <size_t><void *>k2:
                    tmp = <object>PyList_GET_ITEM(bucket, i+2)
                    if <size_t>tmp == <size_t><void *>k3:
                        return <object>PyList_GET_ITEM(bucket, i+3)
        raise KeyError, (k1, k2, k3)

    def __setitem__(self, k, value):
        """
        Set the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = -1
            sage: L[a,b,c]
            -1
        """
        cdef object k1,k2,k3
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError, k
        self.set(k1, k2, k3, value)

    cdef set(self, object k1, object k2, object k3, value):
        if self.threshold and self._size > len(self.buckets) * self.threshold:
            self.resize()
        cdef size_t h1 = <size_t><void *>k1
        cdef size_t h2 = <size_t><void *>k2
        cdef size_t h3 = <size_t><void *>k3
        cdef size_t h = (h1 + 13*h2 ^ 503*h3)
        cdef Py_ssize_t i
        cdef list bucket = <object>PyList_GET_ITEM(self.buckets, h % PyList_GET_SIZE(self.buckets))
        cdef object tmp
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
            tmp = <object>PyList_GET_ITEM(bucket, i)
            if <size_t>tmp == h1:
                tmp = <object>PyList_GET_ITEM(bucket, i+1)
                if <size_t>tmp == h2:
                    tmp = <object>PyList_GET_ITEM(bucket, i+2)
                    if <size_t>tmp == h3:
                        # Test whether the old references are still active
                        r1,r2,r3 = <tuple>(self._refcache[h1,h2,h3])
                        if (isinstance(r1,KeyedRef) and r1() is None) or \
                           (isinstance(r2,KeyedRef) and r2() is None) or \
                           (isinstance(r3,KeyedRef) and r3() is None):
                            del bucket [i:i+4]
                            self._size -= 1
                            break
                        bucket[i+3] = value
                        return
        PyList_Append(bucket, h1)
        PyList_Append(bucket, h2)
        PyList_Append(bucket, h3)
        PyList_Append(bucket, value)
        try:
            ref1 = KeyedRef(k1,self.eraser,(h1, h2, h3))
        except TypeError:
            ref1 = k1
        if k2 is not k1:
            try:
                ref2 = KeyedRef(k2,self.eraser,(h1, h2, h3))
            except TypeError:
                ref2 = k2
        else:
            ref2 = None
        if k3 is not k2 or k3 is not k1:
            try:
                ref3 = KeyedRef(k3,self.eraser,(h1, h2, h3))
            except TypeError:
                ref3 = k3
        else:
            ref3 = None
        self._refcache[h1,h2,h3] = (ref1,ref2,ref3)
        self._size += 1

    def __delitem__(self, k):
        """
        Delete the value corresponding to a key.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: a = 'a'; b = 'b'; c = 'c'
            sage: L[a,b,c] = -1
            sage: (a,b,c) in L
            True
            sage: del L[a,b,c]
            sage: len(L)
            0
            sage: (a,b,c) in L
            False
        """
        cdef object k1,k2,k3
        try:
            k1, k2, k3 = k
        except (TypeError,ValueError):
            raise KeyError, k
        try:
            r1,r2,r3 = self._refcache[<size_t><void *>k1,<size_t><void *>k2,<size_t><void *>k3]
        except KeyError:
            raise KeyError, k
        if (isinstance(r1,KeyedRef) and r1() is None) or \
           (isinstance(r2,KeyedRef) and r2() is None) or \
           (isinstance(r3,KeyedRef) and r3() is None):
            raise KeyError, k
        try:
            del self._refcache[<size_t><void *>k1,<size_t><void *>k2,<size_t><void *>k3]
        except KeyError:
            # This is to cope with a potential racing condition - if garbage
            # collection and weakref callback happens right between the
            # "if (isinstance(r1,..." and the "del", then the previously
            # existing entry might already be gone.
            raise KeyError, k
        cdef size_t h = (<size_t><void *>k1 + 13*<size_t><void *>k2 ^ 503*<size_t><void *>k3)
        cdef Py_ssize_t i
        cdef list bucket = <object>PyList_GET_ITEM(self.buckets, h % PyList_GET_SIZE(self.buckets))
        for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
            if <size_t><object>PyList_GET_ITEM(bucket, i) == <size_t><void *>k1 and \
               <size_t><object>PyList_GET_ITEM(bucket, i+1) == <size_t><void *>k2 and \
               <size_t><object>PyList_GET_ITEM(bucket, i+2) == <size_t><void *>k3:
                del bucket[i:i+4]
                self._size -= 1
                return
        raise KeyError, k

    def resize(self, int buckets=0):
        """
        Change the number of buckets of self, while preserving the contents.

        If the number of buckets is 0 or not given, it resizes self to the
        smallest prime that is at least twice as large as self.

        EXAMPLES::

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
        cdef list old_buckets = self.buckets
        cdef list bucket
        cdef Py_ssize_t i
        cdef size_t h
        self.buckets = [[] for i from 0 <= i <  buckets]
        cdef size_t k1,k2,k3
        cdef object v
        for bucket in old_buckets:
            for i from 0 <= i < PyList_GET_SIZE(bucket) by 4:
                k1 = <size_t><object>PyList_GET_ITEM(bucket, i)
                k2 = <size_t><object>PyList_GET_ITEM(bucket, i+1)
                k3 = <size_t><object>PyList_GET_ITEM(bucket, i+2)
                v  = <object>PyList_GET_ITEM(bucket, i+3)
                h = (k1 + 13*k2 ^ 503*k3)
                self.buckets[h % buckets] += [k1,k2,k3,v]

    def iteritems(self):
        """
        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: L[1,2,3] = None
            sage: list(L.iteritems())
            [((1, 2, 3), None)]
        """
        cdef list bucket
        cdef size_t i, h1,h2,h3
        # We test whether the references are still valid.
        # However, we must not delete them, since we are
        # iterating.
        for bucket in self.buckets:
            for i from 0<=i<len(bucket) by 4:
                h1,h2,h3 = bucket[i:i+3]
                try:
                    r1,r2,r3 = self._refcache[h1,h2,h3]
                except KeyError:
                    # That can only happen under a race condition.
                    # Anyway, it means the item is not there.
                    continue
                if isinstance(r1, KeyedRef):
                    r1 = r1()
                    if r1 is None:
                        continue
                if isinstance(r2, KeyedRef):
                    r2 = r2()
                    if r2 is None:
                        continue
                if isinstance(r3, KeyedRef):
                    r3 = r3()
                    if r3 is None:
                        continue
                yield (r1,r2,r3), <object>PyList_GET_ITEM(bucket,i+3)

    def __reduce__(self):
        """
        Note that we don't expect equality as this class concerns itself with
        object identity rather than object equality.

        EXAMPLES::

            sage: from sage.structure.coerce_dict import TripleDict
            sage: L = TripleDict(31)
            sage: L[1,2,3] = True
            sage: loads(dumps(L)) == L
            False
            sage: list(loads(dumps(L)).iteritems())
            [((1, 2, 3), True)]
        """
        return TripleDict, (len(self.buckets), dict(self.iteritems()), self.threshold)

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
