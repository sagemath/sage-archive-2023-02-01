"""
Fast and safe weak value dictionary

AUTHORS:

- Simon King (2013-10)
- Nils Bruin (2013-10)
- Julian Rueth (2014-03-16): improved handling of unhashable objects

Python's :mod:`weakref` module provides
:class:`~weakref.WeakValueDictionary`. This behaves similar to a dictionary,
but it does not prevent its values from garbage collection. Hence, it stores
the values by weak references with callback functions: The callback function
deletes a key-value pair from the dictionary, as soon as the value becomes
subject to garbage collection.

However, a problem arises if hash and comparison of the key depend on the
value that is being garbage collected::

    sage: import weakref
    sage: class Vals(object): pass
    sage: class Keys:
    ....:     def __init__(self, val):
    ....:         self.val = weakref.ref(val)
    ....:     def __hash__(self):
    ....:         return hash(self.val())
    ....:     def __eq__(self, other):
    ....:         return self.val() == other.val()
    sage: ValList = [Vals() for _ in range(10)]
    sage: D = weakref.WeakValueDictionary()
    sage: for v in ValList:
    ....:     D[Keys(v)] = v
    sage: len(D)
    10
    sage: del ValList, v
    Exception KeyError: (<__main__.Keys instance at ...>,) in <function remove at ...> ignored
    Exception KeyError: (<__main__.Keys instance at ...>,) in <function remove at ...> ignored
    Exception KeyError: (<__main__.Keys instance at ...>,) in <function remove at ...> ignored
    Exception KeyError: (<__main__.Keys instance at ...>,) in <function remove at ...> ignored
    ...
    sage: len(D) > 1
    True

Hence, there are scary error messages, and moreover the defunct items have not
been removed from the dictionary.

Therefore, Sage provides an alternative implementation
:class:`sage.misc.weak_dict.WeakValueDictionary`, using a callback that
removes the defunct item not based on hash and equality check of the key (this
is what fails in the example above), but based on comparison by identity. This
is possible, since references with callback function are distinct even if they
point to the same object. Hence, even if the same object ``O`` occurs as value
for several keys, each reference to ``O`` corresponds to a unique key. We see
no error messages, and the items get correctly removed::

    sage: ValList = [Vals() for _ in range(10)]
    sage: import sage.misc.weak_dict
    sage: D = sage.misc.weak_dict.WeakValueDictionary()
    sage: for v in ValList:
    ....:     D[Keys(v)] = v
    sage: len(D)
    10
    sage: del ValList
    sage: len(D)
    1
    sage: del v
    sage: len(D)
    0

Another problem arises when iterating over the items of a dictionary: If
garbage collection occurs during iteration, then the content of the dictionary
changes, and the iteration breaks for :class:`weakref.WeakValueDictionary`::

    sage: class Cycle:
    ....:     def __init__(self):
    ....:         self.selfref = self
    sage: C = [Cycle() for n in range(10)]
    sage: D = weakref.WeakValueDictionary(enumerate(C))
    sage: import gc
    sage: gc.disable()
    sage: del C[:5]
    sage: len(D)
    10

With :class:`~sage.misc.weak_dict.WeakValueDictionary`, the behaviour is
safer. Note that iteration over a WeakValueDictionary is non-deterministic,
since the lifetime of values (and hence the presence of keys) in the dictionary
may depend on when garbage collection occurs. The method implemented here
will at least postpone dictionary mutations due to garbage collection callbacks.
This means that as long as there is at least one iterator active on a dictionary,
none of its keys will be deallocated (which could have side-effects).
Which entries are returned is of course still dependent on when garbage
collection occurs. Note that when a key gets returned as "present" in the
dictionary, there is no guarantee one can actually retrieve its value: it may
have been garbage collected in the mean time.

Note that Sage's weak value dictionary is actually an instance of
:class:`dict`, in contrast to :mod:`weakref`'s weak value dictionary::

    sage: issubclass(weakref.WeakValueDictionary, dict)
    False
    sage: issubclass(sage.misc.weak_dict.WeakValueDictionary, dict)
    True

See :trac:`13394` for a discussion of some of the design considerations.
"""
########################################################################
#       Copyright (C) 2013 Simon King <simon.king@uni-jena.de>
#                          Nils Bruin <nbruin@sfu.ca>
#                          Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import weakref
from weakref import KeyedRef
from copy import deepcopy

from cpython.dict cimport *
from cpython.weakref cimport *
from cpython.list cimport *

from cpython cimport Py_XINCREF, Py_XDECREF

cdef extern from "Python.h":
    ctypedef struct PyDictEntry:
        Py_ssize_t me_hash
        PyObject* me_key
        PyObject* me_value
    ctypedef struct PyDictObject:
        Py_ssize_t ma_fill
        Py_ssize_t ma_used
        Py_ssize_t ma_mask
        PyDictEntry* ma_table
        PyDictEntry* (*ma_lookup)(PyDictObject *mp, PyObject *key, long hash) except NULL
        
    PyObject* Py_None
    #we need this redefinition because we want to be able to call
    #PyWeakref_GetObject with borrowed references. This is the recommended
    #strategy according to Cython/Includes/cpython/__init__.pxd
    PyObject* PyWeakref_GetObject(PyObject * wr)
    int PyList_SetItem(object list, Py_ssize_t index,PyObject * item) except -1
    #this one's just missing.
    long PyObject_Hash(object obj) except -1

cdef PyObject* PyDict_GetItemWithError(dict op, object key) except? NULL:
    cdef PyDictEntry* ep
    cdef PyDictObject* mp = <PyDictObject*><void *>op
    ep = mp.ma_lookup(mp, <PyObject*><void*>key, PyObject_Hash(key))
    if ep:
        return ep.me_value
    else:
        return NULL

#this routine extracts the "dummy" sentinel value that is used in dicts to mark
#"freed" slots. We need that to delete things ourselves.

cdef PyObject* init_dummy() except NULL:
    cdef dict D = dict()
    cdef PyDictObject* mp = <PyDictObject *><void *>D
    cdef size_t mask
    cdef PyDictEntry* ep0 = mp.ma_table
    cdef PyDictEntry* ep
    cdef size_t i
    global dummy

    D[0]=0; del D[0] #ensure that there is a "deleted" entry in the dict
    mask = mp.ma_mask
    #since our entry had hash 0, we should really succeed on the first iteration
    for i in range(mask+1):
        ep = &(ep0[i])
        if ep.me_key != NULL:
            return ep.me_key
    raise RuntimeError("Problem initializing dummy")

#note that dummy here is a borrowed reference. That's not a problem because
#we're never giving it up and dictobject.c is also holding a permanent reference
#to this object
cdef PyObject* dummy = init_dummy()

#this routine looks for the first entry in dict D with given hash of the
#key and given (identical!) value and deletes that entry.
cdef del_dictitem_by_exact_value(PyDictObject *mp, PyObject *value, long hash):
    """
    This is used in callbacks for the weak values of :class:`WeakValueDictionary`.

    INPUT:

    - ``PyDictObject *mp`` -- pointer to a dict
    - ``PyObject *value``  -- pointer to a value of the dictionary
    - ``long hash``        -- hash of the key by which the value is stored in the dict

    The hash bucket determined by the given hash is searched for the item
    containing the given value. If this item can't be found, the function is
    silently returning. Otherwise, the item is removed from the dict.

    TESTS:

    The following is an indirect doctest, as discussed on :trac:`13394`.
    ::

        sage: from sage.misc.weak_dict import WeakValueDictionary
        sage: V = [set(range(n)) for n in range(5)]
        sage: D = WeakValueDictionary(enumerate(V))

    The line ``V[k] = None`` triggers execution of the callback functions of
    the dict values. However, the actual deletion is postponed till after the
    iteration over the dictionary has finished. Hence, when the callbacks are
    executed, the values which the callback belongs to has already been
    overridded by a new value. Therefore, the callback does not delete the
    item::

        sage: for k in D.iterkeys():    # indirect doctest
        ....:     V[k] = None
        ....:     D[k] = ZZ
        sage: len(D)
        5
        sage: D[1]
        Integer Ring

    TESTS:

    The following shows that the deletion of deeply nested structures does not
    result in an error, by :trac:`15506`::

        sage: class A: pass
        sage: a = A(); prev = a
        sage: M = WeakValueDictionary()
        sage: for i in range(10^3+10): newA = A(); M[newA] = prev; prev = newA
        sage: del a

    """
    cdef size_t i
    cdef size_t perturb
    cdef size_t mask = <size_t> mp.ma_mask
    cdef PyDictEntry* ep0 = mp.ma_table
    cdef PyDictEntry* ep
    i = hash & mask
    ep = &(ep0[i])

    if ep.me_key == NULL:
        # key not found
        return

    perturb = hash
    while (<PyObject *>(ep.me_value) != value or ep.me_hash != hash):
        i = (i << 2) + i + perturb +1
        ep = &ep0[i & mask]
        if ep.me_key == NULL:
            # key not found
            return
        perturb = perturb >> 5 #this is the value of PERTURB_SHIFT

    T=PyList_New(2)
    PyList_SetItem(T,0,ep.me_key)
    if dummy == NULL:
        raise RuntimeError("dummy needs to be initialized")
    Py_XINCREF(dummy)
    ep.me_key = dummy
    PyList_SetItem(T,1,ep.me_value)
    ep.me_value = NULL
    mp.ma_used -= 1
    #We have transferred the to-be-deleted references to the list T
    #we now delete the list so that the actual decref happens through a
    #deallocation routine that uses the Python Trashcan macros to
    #avoid stack overflow in deleting deep structures.
    del T

def test_del_dictitem_by_exact_value(D, value, h):
    """
    This function helps testing some cdef function used to delete dictionary items.

    INPUT:

    - ``D`` -- a Python ``<dict>``.
    - ``value`` -- an object that is value ``D``.
    - ``h`` -- the hash of the key under which to find ``value`` in ``D``.

    The underlying cdef function deletes an item from ``D`` that is in the
    hash bucket determined by ``h`` and whose value is identic with
    ``value``. Of course, this only makes sense if the pairs ``(h, value)``
    corresponding to items in ``D`` are pair-wise distinct.

    If a matching item can not be found, the function does nothing and
    silently returns.

    TESTS:

    See :trac:`13394` for a discussion.
    ::

        sage: from sage.misc.weak_dict import test_del_dictitem_by_exact_value
        sage: B=1000
        sage: L=list(range(B))
        sage: D1=dict()
        sage: D2=dict()
        sage: for i in range(100000):        # long time
        ....:     ki=L[floor(random()*B)]
        ....:     vi=L[floor(random()*B)]
        ....:     D1[ki]=vi
        ....:     D2[ki]=vi
        ....:     ko=L[floor(random()*B)]
        ....:     if ko in D1:
        ....:         vo=D1[ko]
        ....:         del D1[ko]
        ....:         test_del_dictitem_by_exact_value(D2,vo,hash(ko))
        ....:     assert D1 == D2

    No action is taken if the item prescribed by key hash and value does not
    exist in the dictionary::

        sage: D = {1: ZZ}
        sage: test_del_dictitem_by_exact_value(D, ZZ, 2)
        sage: D
        {1: Integer Ring}
        sage: test_del_dictitem_by_exact_value(D, QQ, 1)
        sage: D
        {1: Integer Ring}

    """
    return del_dictitem_by_exact_value(<PyDictObject *>D, <PyObject *>value, h)

cdef class WeakValueDictEraser:
    """
    Erases items from a :class:`sage.misc.weak_dict.WeakValueDictionary` when
    a weak reference becomes invalid.

    This is of internal use only. Instances of this class will be passed as a
    callback function when creating a weak reference.

    EXAMPLES::

        sage: from sage.misc.weak_dict import WeakValueDictionary
        sage: v = frozenset([1])
        sage: D = WeakValueDictionary({1 : v})
        sage: len(D)
        1
        sage: del v
        sage: len(D)
        0

    AUTHOR:

     - Nils Bruin (2013-11)
    """
    cdef D
    def __init__(self, D):
        """
        INPUT:

        A :class:`sage.misc.weak_dict.WeakValueDictionary`.

        EXAMPLES::

            sage: v = frozenset([1])
            sage: D = sage.misc.weak_dict.WeakValueDictionary({ 1 : v })
            sage: len(D)
            1
            sage: del v
            sage: len(D) #indirect doctest
            0
        """
        self.D = PyWeakref_NewRef(D,None)
    def __call__(self, r):
        """
        INPUT:

        A weak reference with key.

        When this is called with a weak reference ``r``, then an entry from the
        dictionary pointed to by ``self.D`` is removed that has ``r`` as a value
        identically, stored under a key with hash ``r.key``. If no such key
        exists, or if the dictionary itself doesn't exist any more, then nothing
        happens.

        If the dictionary has an iterator active on it then the object is
        queued for removal when all iterators have concluded.

        EXAMPLES::

            sage: v = frozenset([1])
            sage: D = sage.misc.weak_dict.WeakValueDictionary({ 1 : v })
            sage: len(D)
            1
            sage: del v
            sage: len(D) #indirect doctest
            0
        """
        cdef WeakValueDictionary D = <object> PyWeakref_GetObject(<PyObject*> self.D)
        if D is None:
            return
        #The situation is the following:
        #in the underlying dictionary, we have stored a KeyedRef r
        #under a key k. The attribute r.key is the hash of k.
        if D._guard_level:
            D._pending_removals.append(r)
        else:
            del_dictitem_by_exact_value(<PyDictObject *>D, <PyObject *>r, r.key)

cdef class WeakValueDictionary(dict):
    """
    IMPLEMENTATION:

    The :class:`WeakValueDictionary` inherits from :class:`dict`. In its
    implementation, it stores weakrefs to the actual values under the keys.
    All access routines are wrapped to transparently place and remove these
    weakrefs.

    NOTE:

    In contrast to :class:`weakref.WeakValueDictionary` in Python's
    :mod:`weakref` module, the callback does not need to assume that the
    dictionary key is a valid Python object when it is called. There is no
    need to compute the hash or compare the dictionary keys. This is why
    the example below would not work with
    :class:`weakref.WeakValueDictionary`, but does work with
    :class:`sage.misc.weak_dict.WeakValueDictionary`.

    EXAMPLES::

        sage: import weakref
        sage: class Vals(object): pass
        sage: class Keys:
        ....:     def __init__(self, val):
        ....:         self.val = weakref.ref(val)
        ....:     def __hash__(self):
        ....:         return hash(self.val())
        ....:     def __eq__(self, other):
        ....:         return self.val() == other.val()
        sage: ValList = [Vals() for _ in range(10)]
        sage: import sage.misc.weak_dict
        sage: D = sage.misc.weak_dict.WeakValueDictionary()
        sage: for v in ValList:
        ....:     D[Keys(v)] = v
        sage: len(D)
        10
        sage: del ValList
        sage: len(D)
        1
        sage: del v
        sage: len(D)
        0

    TESTS:

    The following reflects the behaviour of the callback on weak dict values,
    as discussed on :trac:`13394`.  ::

        sage: from sage.misc.weak_dict import WeakValueDictionary
        sage: V = [set(range(n)) for n in range(5)]
        sage: D = WeakValueDictionary(enumerate(V))

    The line ``V[k] = None`` triggers execution of the callback functions of
    the dictionary values. However, the actual deletion is postponed till
    after the iteration over the dictionary has finished. Hence, when the
    callbacks are executed, the values which the callback belongs to has
    already been overridded by a new value. Therefore, the callback does not
    delete the item::

        sage: for k in D.iterkeys():    # indirect doctest
        ....:     V[k] = None
        ....:     D[k] = ZZ
        sage: len(D)
        5
        sage: D[1]
        Integer Ring

    The following is a stress test for weak value dictionaries::

        sage: class C(object):
        ....:     def __init__(self, n):
        ....:         self.n = n
        ....:     def __cmp__(self, other):
        ....:         return cmp(type(self),type(other)) or cmp(self.n, other.n)
        sage: B = 100
        sage: L = [None]*B
        sage: D1 = WeakValueDictionary()
        sage: D2 = WeakValueDictionary()
        sage: for i in range(10000):
        ....:     ki = floor(random()*B)
        ....:     vi = C(floor(random()*B))
        ....:     D1[ki] = vi
        ....:     D2[ki] = vi
        ....:     L[ki]  = vi
        ....:     del vi
        ....:     ko = floor(random()*B)
        ....:     if ko in D1:
        ....:         del D1[ko]
        ....:         L[ko] = None
        ....:     assert D1 == D2

    """
    cdef __weakref__
    cdef callback
    cdef int _guard_level
    cdef list _pending_removals

    def __init__(self, data=()):
        """
        Create a :class:`WeakValueDictionary` with given initial data.

        INPUT:

        Optional iterable of key-value pairs

        EXAMPLES::

            sage: L = [(p,GF(p)) for p in prime_range(10)]
            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: len(D)
            0
            sage: D = sage.misc.weak_dict.WeakValueDictionary(L)
            sage: len(D) == len(L)
            True

        """
        dict.__init__(self)
        self.callback = WeakValueDictEraser(self)
        self._guard_level = 0
        self._pending_removals = []
        try:
            data=data.iteritems()
        except AttributeError:
            pass
        for k,v in data:
            self[k] = v

    def __copy__(self):
        """
        Return a copy of this weak dictionary.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D[1] = QQ
            sage: D[2] = ZZ
            sage: D[None] = CC
            sage: E = copy(D)    # indirect doctest
            sage: set(E.items()) == set(D.items())
            True

        """
        return WeakValueDictionary(self.iteritems())

    def __deepcopy__(self, memo):
        """
        Return a copy of this dictionary using copies of the keys.

        NOTE:

        The values of the dictionary are not copied, since we can not copy the
        external strong references to the values, which are decisive for
        garbage collection.

        EXAMPLES::

            sage: class C(object): pass
            sage: V = [C(),C()]
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D[C()] = V[0]
            sage: D[C()] = V[1]
            sage: E = deepcopy(D)     # indirect doctest

        The keys are copied (in this silly example, the copies of the keys are
        actually not equal to the original keys)::

            sage: set(E.keys()) == set(D.keys())
            False

        However, the values are not copied::

            sage: set(E.values()) == set(D.values()) == set(V)
            True

        """
        out = WeakValueDictionary()
        for k,v in self.iteritems():
            out[deepcopy(k, memo)] = v
        return out

    def __repr__(self):
        """
        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: repr(sage.misc.weak_dict.WeakValueDictionary([(1,ZZ),(2,QQ)]))  # indirect doctest
            '<WeakValueDictionary at 0x...>'

        """
        return "<WeakValueDictionary at 0x%x>" % id(self)

    def __str__(self):
        """
        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: str(sage.misc.weak_dict.WeakValueDictionary([(1,ZZ),(2,QQ)]))  # indirect doctest
            '<WeakValueDictionary at 0x...>'

        """
        return "<WeakValueDictionary at 0x%x>" % id(self)

    def setdefault(self, k, default=None):
        """
        Return the stored value for a given key; return and store a default
        value if no previous value is stored.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: L = [(p,GF(p)) for p in prime_range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(L)
            sage: len(D)
            4

        The value for an existing key is returned and not overridden::

            sage: D.setdefault(5, ZZ)
            Finite Field of size 5
            sage: D[5]
            Finite Field of size 5

        For a non-existing key, the default value is stored and returned::

            sage: 4 in D
            False
            sage: D.setdefault(4, ZZ)
            Integer Ring
            sage: 4 in D
            True
            sage: D[4]
            Integer Ring
            sage: len(D)
            5

        TESTS:

        Check that :trac:`15956` has been fixed, i.e., a ``TypeError`` is
        raised for unhashable objects::

            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D.setdefault(matrix([]),ZZ)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        """
        cdef PyObject* wr = PyDict_GetItemWithError(self, k)
        if wr != NULL:
            out = PyWeakref_GetObject(wr)
            if out != Py_None:
                return <object>out
        PyDict_SetItem(self,k,KeyedRef(default,self.callback,PyObject_Hash(k)))
        return default

    def __setitem__(self, k, v):
        """
        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: ZZ in D
            False

        One can set new items::

            sage: D[ZZ] = QQ   # indirect doctest
            sage: D[ZZ]
            Rational Field
            sage: len(D)
            1
            sage: ZZ in D
            True

       One can also override existing items::

           sage: D[ZZ] = RLF
           sage: ZZ in D
           True
           sage: D[ZZ]
           Real Lazy Field
           sage: len(D)
           1

        TESTS:

        One may wonder whether it causes problems when garbage collection for
        a previously existing item happens *after* overriding the item. The
        example shows that it is not a problem::

            sage: class Cycle:
            ....:     def __init__(self):
            ....:         self.selfref = self
            sage: L = [Cycle() for _ in range(5)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: len(D)
            5
            sage: import gc
            sage: gc.disable()
            sage: del L
            sage: len(D)
            5
            sage: D[2] = ZZ
            sage: len(D)
            5
            sage: gc.enable()
            sage: _ = gc.collect()
            sage: len(D)
            1
            sage: D.items()
            [(2, Integer Ring)]

        Check that :trac:`15956` has been fixed, i.e., a ``TypeError`` is
        raised for unhashable objects::

            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D[matrix([])] = ZZ
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        """
        PyDict_SetItem(self,k,KeyedRef(v,self.callback,PyObject_Hash(k)))

    #def __delitem__(self, k):
    #we don't really have to override this method.

    def pop(self, k):
        """
        Return the value for a given key, and delete it from the dictionary.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: L = [GF(p) for p in prime_range(10^3)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: 20 in D
            True
            sage: D.pop(20)
            Finite Field of size 73
            sage: 20 in D
            False
            sage: D.pop(20)
            Traceback (most recent call last):
            ...
            KeyError: 20

        TESTS:

        Check that :trac:`15956` has been fixed, i.e., a ``TypeError`` is
        raised for unhashable objects::

            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D.pop(matrix([]))
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        """
        cdef PyObject* wr = PyDict_GetItemWithError(self, k)
        if wr == NULL:
            raise KeyError(k)
        cdef PyObject* outref = PyWeakref_GetObject(wr)
        if outref == Py_None:
            raise KeyError(k)
        #we turn the output into a new reference before deleting the item,
        #because the deletion can cause any kind of havoc.
        out = <object>outref
        del self[k]
        return out

    def popitem(self):
        """
        Return and delete some item from the dictionary.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D[1] = ZZ

        The dictionary only contains a single item, hence, it is clear which
        one will be returned::

            sage: D.popitem()
            (1, Integer Ring)

        Now, the dictionary is empty, and hence the next attempt to pop an
        item will fail with a ``KeyError``::

            sage: D.popitem()
            Traceback (most recent call last):
            ...
            KeyError: 'popitem(): weak value dictionary is empty'

        """
        for k,v in self.iteritems():
            del self[k]
            return k, v
        raise KeyError('popitem(): weak value dictionary is empty')

    def get(self, k, d=None):
        """
        Return the stored value for a key, or a default value for unkown keys.

        The default value defaults to None.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: L = [GF(p) for p in prime_range(10^3)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: 100 in D
            True
            sage: 200 in D
            False
            sage: D.get(100, "not found")
            Finite Field of size 547
            sage: D.get(200, "not found")
            'not found'
            sage: D.get(200) is None
            True

        TESTS:

        Check that :trac:`15956` has been fixed, i.e., a ``TypeError`` is
        raised for unhashable objects::

            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D.get(matrix([]))
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        """
        cdef PyObject * wr = PyDict_GetItemWithError(self, k)
        if wr == NULL:
            return d
        out = PyWeakref_GetObject(wr)
        if out == Py_None:
            return d
        else:
            return <object>out

    def __getitem__(self, k):
        """
        TESTS::

            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D[ZZ] = QQ
            sage: D[QQ]
            Traceback (most recent call last):
            ...
            KeyError: Rational Field
            sage: D[ZZ]     # indirect doctest
            Rational Field

        As usual, the dictionary keys are compared by ``==`` and not by
        identity::

            sage: D[10] = ZZ
            sage: D[int(10)]
            Integer Ring

        Check that :trac:`15956` has been fixed, i.e., a ``TypeError`` is
        raised for unhashable objects::

            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D[matrix([])]
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        """
        cdef PyObject* wr = PyDict_GetItemWithError(self, k)
        if wr == NULL:
            raise KeyError(k)
        out = PyWeakref_GetObject(wr)
        if out == Py_None:
            raise KeyError(k)
        return <object>out

    def __contains__(self, k):
        """
        Containment in the set of keys.

        TESTS::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: 3 in D     # indirect doctest
            True

        As usual, keys are compared by equality and not by identity::

            sage: int(3) in D
            True

        This is a weak value dictionary. Hence, the existence of the
        dictionary does not prevent the values from garbage collection,
        thereby removing the corresponding key-value pairs::

            sage: del L[3]
            sage: 3 in D
            False

        Check that :trac:`15956` has been fixed, i.e., a ``TypeError`` is
        raised for unhashable objects::

            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: matrix([]) in D
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable

        """
        cdef PyDictObject* mp=<PyDictObject*><void*>self
        cdef PyDictEntry* ep=mp.ma_lookup(mp,<PyObject*><void*>k, PyObject_Hash(k))
        return (ep.me_value != NULL) and (PyWeakref_GetObject(ep.me_value) != Py_None)

    #def __len__(self):
    #since GC is not deterministic, neither is the length of a WeakValueDictionary,
    #so we might as well just return the normal dictionary length.

    def iterkeys(self):
        """
        Iterate over the keys of this dictionary.

        .. WARNING::

            Iteration is unsafe, if the length of the dictionary changes
            during the iteration! This can also happen by garbage collection.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: del L[4]

        One item got deleted from the list ``L`` and hence the corresponding
        item in the dictionary got deleted as well. Therefore, the
        corresponding key 4 is missing in the list of keys::

            sage: list(sorted(D.iterkeys()))
            [0, 1, 2, 3, 5, 6, 7, 8, 9]

        """
        cdef PyObject *key, *wr
        cdef Py_ssize_t pos = 0
        try:
            self._enter_iter()
            while PyDict_Next(self, &pos, &key, &wr):
                #this check doesn't really say anything: by the time
                #the key makes it to the customer, it may have already turned
                #invalid. It's a cheap check, though.
                if PyWeakref_GetObject(wr)!=Py_None:
                    yield <object>key
        finally:
            self._exit_iter()

    def __iter__(self):
        """
        Iterate over the keys of this dictionary.

        .. WARNING::

            Iteration is unsafe, if the length of the dictionary changes
            during the iteration! This can also happen by garbage collection.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: del L[4]

        One item got deleted from the list ``L`` and hence the corresponding
        item in the dictionary got deleted as well. Therefore, the
        corresponding key 4 is missing in the list of keys::

            sage: sorted(list(D))    # indirect doctest
            [0, 1, 2, 3, 5, 6, 7, 8, 9]

        """
        return self.iterkeys()

    def keys(self):
        """
        The list of keys.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: del L[4]

        One item got deleted from the list ``L`` and hence the corresponding
        item in the dictionary got deleted as well. Therefore, the
        corresponding key 4 is missing in the list of keys::

            sage: sorted(D.keys())
            [0, 1, 2, 3, 5, 6, 7, 8, 9]

        """
        return list(self.iterkeys())

    def itervalues(self):
        """
        Iterate over the values of this dictionary.

        .. WARNING::

            Iteration is unsafe, if the length of the dictionary changes
            during the iteration! This can also happen by garbage collection.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals:
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __repr__(self):
            ....:         return "<%s>"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: L = [Vals(n) for n in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))

        We delete one item from ``D`` and we delete one item from the list
        ``L``. The latter implies that the corresponding item from ``D`` gets
        deleted as well. Hence, there remain eight values::

            sage: del D[2]
            sage: del L[5]
            sage: for v in sorted(D.itervalues()):
            ....:     print v
            <0>
            <1>
            <3>
            <4>
            <6>
            <7>
            <8>
            <9>

        """
        cdef PyObject *key, *wr
        cdef Py_ssize_t pos = 0
        try:
            self._enter_iter()
            while PyDict_Next(self, &pos, &key, &wr):
                out = PyWeakref_GetObject(wr)
                if out != Py_None:
                    yield <object>out
        finally:
            self._exit_iter()

    def values(self):
        """
        Return the list of values.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals:
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __repr__(self):
            ....:         return "<%s>"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: L = [Vals(n) for n in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))

        We delete one item from ``D`` and we delete one item from the list
        ``L``. The latter implies that the corresponding item from ``D`` gets
        deleted as well. Hence, there remain eight values::

            sage: del D[2]
            sage: del L[5]
            sage: sorted(D.values())
            [<0>, <1>, <3>, <4>, <6>, <7>, <8>, <9>]

        """
        return list(self.itervalues())

    def iteritems(self):
        """
        Iterate over the items of this dictionary.

        .. WARNING::

            Iteration is unsafe, if the length of the dictionary changes
            during the iteration! This can also happen by garbage collection.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals:
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __repr__(self):
            ....:         return "<%s>"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: class Keys(object):
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __hash__(self):
            ....:         if self.n%2:
            ....:             return 5
            ....:         return 3
            ....:     def __repr__(self):
            ....:         return "[%s]"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: L = [(Keys(n), Vals(n)) for n in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(L)

        We remove one dictionary item directly. Another item is removed by
        means of garbage collection. By consequence, there remain eight
        items in the dictionary::

            sage: del D[Keys(2)]
            sage: del L[5]
            sage: for k,v in sorted(D.iteritems()):
            ....:     print k, v
            [0] <0>
            [1] <1>
            [3] <3>
            [4] <4>
            [6] <6>
            [7] <7>
            [8] <8>
            [9] <9>

        """
        cdef PyObject *key, *wr
        cdef Py_ssize_t pos = 0
        try:
            self._enter_iter()
            while PyDict_Next(self, &pos, &key, &wr):
                out = PyWeakref_GetObject(wr)
                if out != Py_None:
                    yield <object>key, <object>out
        finally:
            self._exit_iter()

    def items(self):
        """
        The key-value pairs of this dictionary.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals:
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __repr__(self):
            ....:         return "<%s>"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: class Keys(object):
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __hash__(self):
            ....:         if self.n%2:
            ....:             return 5
            ....:         return 3
            ....:     def __repr__(self):
            ....:         return "[%s]"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: L = [(Keys(n), Vals(n)) for n in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(L)

        We remove one dictionary item directly. Another item is removed by
        means of garbage collection. By consequence, there remain eight
        items in the dictionary::

            sage: del D[Keys(2)]
            sage: del L[5]
            sage: sorted(D.items())
            [([0], <0>),
             ([1], <1>),
             ([3], <3>),
             ([4], <4>),
             ([6], <6>),
             ([7], <7>),
             ([8], <8>),
             ([9], <9>)]

        """
        return list(self.iteritems())

    cdef int _enter_iter(self) except -1:
        """
        Make sure that items of a weak value dictionary are not actually
        deleted, but only *marked* for deletion.

        TESTS::

            sage: from sage.misc.weak_dict import WeakValueDictionary
            sage: K = [frozenset([i]) for i in range(11)]
            sage: D = WeakValueDictionary((K[i],K[i+1]) for i in range(10))
            sage: k = K[10]
            sage: del K
            sage: i = D.iterkeys(); d = i.next(); del d
            sage: len(D.keys())
            10
            sage: del k
            sage: len(D.keys())
            9
            sage: del i
            sage: len(D.keys())
            0

        """
        self._guard_level += 1
        return 0

    cdef int _exit_iter(self) except -1:
        """
        Make sure that all items of a weak value dictionary that are marked
        for deletion are actually deleted, as soon as there is no iteration
        over the dictionary.

        TESTS::

            sage: from sage.misc.weak_dict import WeakValueDictionary
            sage: K = [frozenset([i]) for i in range(11)]
            sage: D = WeakValueDictionary((K[i],K[i+1]) for i in range(10))
            sage: k = K[10]
            sage: del K
            sage: i = D.iterkeys(); d = i.next(); del d
            sage: len(D.keys())
            10
            sage: del k
            sage: len(D.keys())
            9
            sage: del i
            sage: len(D.keys())
            0

        """
        self._guard_level -= 1
        #when the guard_level drops to zero, we try to remove all the
        #pending removals. Note that this could trigger another iterator
        #to become active, in which case we should back off.
        while (not self._guard_level) and self._pending_removals:
            self.callback(self._pending_removals.pop())
        return 0
