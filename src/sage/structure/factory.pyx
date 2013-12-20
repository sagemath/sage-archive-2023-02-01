r"""
Factory for unique objects

More general than :mod:`~sage.structure.unique_representation`, the
:class:`UniqueFactory` class can specify a subset of the arguments
that serve as the unique key. Typically, this is used to construct
objects that accept an optional ``check=[True|False]`` argument, but
whose result should be unique irregardless of said optional argument.
"""

#*****************************************************************************
#  Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#******************************************************************************

import types, copy_reg

from sage_object cimport SageObject

cdef sage_version
from sage.version import version as sage_version

sage_version = sage_version.split('.')
for i in range(len(sage_version)):
    try:
        sage_version[i] = int(sage_version[i])
    except ValueError:
        pass
sage_version = tuple(sage_version)

import sage.misc.weak_dict

cdef class UniqueFactory(SageObject):
    """
    This class is intended to make it easy to make unique, cached objects.

    It is based on the idea that the object is uniquely defined by a set of
    defining data (the key). There is also the possibility some non-defining
    data (extra args) which will be used in initial creation, but not affect
    the caching.

    The objects created are cached (using weakrefs) based on their key and
    returned directly rather than re-created if requested again. Pickling
    will return the same object for the same version of Sage, and distinct
    (but hopefully equal) objects for different versions of Sage.

    Typically one only needs to implement :meth:`create_key` and
    :meth:`create_object`.
    """

    cdef readonly _name
    cdef readonly _cache

    def __init__(self, name):
        """
        INPUT:

        - ``name`` -- string. A name in the global namespace referring
          to self or a fully qualified path name to self, which is
          used to locate the factory on unpickling.

        EXAMPLES::

            sage: from sage.structure.factory import UniqueFactory
            sage: fake_factory = UniqueFactory('ZZ')
            sage: loads(dumps(fake_factory))
            Integer Ring
            sage: fake_factory = UniqueFactory('sage.rings.all.QQ')
            sage: loads(dumps(fake_factory))
            Rational Field
        """
        self._name = name
        self._cache = sage.misc.weak_dict.WeakValueDictionary()

    def __reduce__(self):
        """
        EXAMPLES::

            sage: A = FiniteField(127)
            sage: A is loads(dumps(A)) # indirect doctest
            True
            sage: B = FiniteField(3^3,'b')
            sage: B is loads(dumps(B))
            True
            sage: C = FiniteField(2^16,'c')
            sage: C is loads(dumps(C))
            True
            sage: D = FiniteField(3^20,'d')
            sage: D is loads(dumps(D))
            True

        TESTS::

            sage: loads(dumps(FiniteField)) is FiniteField
            True
            sage: from sage.structure.test_factory import test_factory
            sage: loads(dumps(test_factory)) is test_factory
            True
        """
        return lookup_global, (self._name,)

    def __call__(self, *args, **kwds):
        """
        This is the method invoked to create objects. It first creates a key
        from the given parameters, then if an object with that key already
        exists returns it, and otherwise creates one and stores a weak reference
        to it in its dictionary.

        Do not override this method, override create_key and create_object and
        put the docstring in the body of the class.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: _ = test_factory(1,2,3); _
            Making object (1, 2, 3)
            <sage.structure.test_factory.A instance at ...>

        It already created one, so don't re-create::

            sage: test_factory(1,2,3)
            <sage.structure.test_factory.A instance at ...>
            sage: test_factory(1,2,3) is test_factory(1,2,3)
            True

        Of course, with a different key, a new object will be created::

            sage: test_factory(1,2,3) is test_factory(1,2,4)
            Making object (1, 2, 4)
            False
        """
        key, kwds = self.create_key_and_extra_args(*args, **kwds)
        version = self.get_version(sage_version)
        return self.get_object(version, key, kwds)

    cpdef get_object(self, version, key, extra_args):
        """
        Returns the object corresponding to key, creating it with extra_args
        if necessary (for example, it isn't in the cache or it is unpickling
        from an older version of Sage).

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: a = test_factory.get_object(3.0, 'a', {}); a
            Making object a
            <sage.structure.test_factory.A instance at ...>
            sage: test_factory.get_object(3.0, 'a', {}) is test_factory.get_object(3.0, 'a', {})
            True
            sage: test_factory.get_object(3.0, 'a', {}) is test_factory.get_object(3.1, 'a', {})
            Making object a
            False
            sage: test_factory.get_object(3.0, 'a', {}) is test_factory.get_object(3.0, 'b', {})
            Making object b
            False
        """
        try:
            return self._cache[version, key]
        except KeyError:
            pass
        obj = self.create_object(version, key, **extra_args)
        self._cache[version, key] = obj
        try:
            other_keys = self.other_keys(key, obj)
            for key in other_keys:
                try:
                    obj = self._cache[version, key]
                    break
                except KeyError:
                    pass
            for key in other_keys:
                self._cache[version, key] = obj
            obj._factory_data = self, version, key, extra_args
            if obj.__class__.__reduce__.__objclass__ is object:
                # replace the generic object __reduce__ to use this one
                obj.__reduce_ex__ = types.MethodType(generic_factory_reduce, obj)
        except AttributeError:
            pass
        return obj

    cpdef get_version(self, sage_version):
        """
        This is provided to allow more or less granular control over
        pickle versioning. Objects pickled in the same version of Sage
        will unpickle to the same rather than simply equal objects. This
        can provide significant gains as arithmetic must be performed on
        objects with identical parents. However, if there has been an
        incompatible change (e.g. in element representation) we want the
        version number to change so coercion is forced between the two
        parents.

        Defaults to the Sage version that is passed in, but courser
        granularity can be provided.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: test_factory.get_version((3,1,0))
            (3, 1, 0)
        """
        return sage_version

    def create_key_and_extra_args(self, *args, **kwds):
        r"""
        Return a tuple containing the key (uniquely defining data)
        and any extra arguments (empty by default).

        Defaults to :meth:`create_key`.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: test_factory.create_key_and_extra_args(1, 2, key=5)
            ((1, 2), {})
            sage: GF.create_key_and_extra_args(3, foo='value')
            ((3, None, None, None, "{'foo': 'value'}", 3, 1, True), {'foo': 'value'})
        """
        return self.create_key(*args, **kwds), {}

    def create_key(self, *args, **kwds):
        """
        Given the arguments and keywords, create a key that uniquely
        determines this object.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: test_factory.create_key(1, 2, key=5)
            (1, 2)
        """
        raise NotImplementedError

    def create_object(self, version, key, **extra_args):
        """
        Create the object from the key and extra arguments. This is only
        called if the object was not found in the cache.

        EXAMPLES::

            sage: from sage.structure.test_factory import test_factory
            sage: test_factory.create_object(0, (1,2,3))
            Making object (1, 2, 3)
            <sage.structure.test_factory.A instance at ...>
            sage: test_factory('a')
            Making object ('a',)
            <sage.structure.test_factory.A instance at ...>
            sage: test_factory('a') # NOT called again
            <sage.structure.test_factory.A instance at ...>
        """
        raise NotImplementedError

    cpdef other_keys(self, key, obj):
        """
        Sometimes during object creation, certain defaults are chosen which
        may result in a new (more specific) key. This allows the more specific
        key to be cached as well, and used for pickling.

        EXAMPLES::

            sage: key, _ = GF.create_key_and_extra_args(27, 'k'); key
            (27, ('k',), x^3 + 2*x + 1, None, '{}', 3, 3, True)
            sage: K = GF.create_object(0, key); K
            Finite Field in k of size 3^3
            sage: GF.other_keys(key, K)
            [(27, ('k',), x^3 + 2*x + 1, None, '{}', 3, 3, True),
             (27, ('k',), x^3 + 2*x + 1, 'givaro', '{}', 3, 3, True)]

            sage: K = GF(7^40, 'a')
            sage: loads(dumps(K)) is K
            True
        """
        return []

    cpdef reduce_data(self, obj):
        """
        The results of this function can be returned from
        :meth:`__reduce__`. This is here so the factory internals can
        change without having to re-write :meth:`__reduce__` methods
        that use it.

        EXAMPLE::

            sage: V = FreeModule(ZZ, 5)
            sage: factory, data = FreeModule.reduce_data(V)
            sage: factory(*data)
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
            sage: factory(*data) is V
            True

            sage: from sage.structure.test_factory import test_factory
            sage: a = test_factory(1, 2)
            Making object (1, 2)
            sage: test_factory.reduce_data(a)
            (<built-in function generic_factory_unpickle>,
             (<class 'sage.structure.test_factory.UniqueFactoryTester'>,
              (...),
              (1, 2),
              {}))

        Note that the ellipsis ``(...)`` here stands for the Sage
        version.
        """
        return generic_factory_unpickle, obj._factory_data


def generic_factory_unpickle(UniqueFactory factory, *args):
    """
    Method used for unpickling the object.

    The unpickling mechanism needs a plain Python function to call.
    It takes a factory as the first argument, passes the rest of the
    arguments onto the factory's :meth:`UniqueFactory.get_object` method.

    EXAMPLES::

        sage: V = FreeModule(ZZ, 5)
        sage: func, data = FreeModule.reduce_data(V)
        sage: func is sage.structure.factory.generic_factory_unpickle
        True
        sage: sage.structure.factory.generic_factory_unpickle(*data) is V
        True
    """
    return factory.get_object(*args)

def generic_factory_reduce(self, proto):
    """
    Used to provide a ``__reduce__`` method if one does not already exist.

    EXAMPLES::

        sage: V = QQ^6
        sage: sage.structure.factory.generic_factory_reduce(V, 1) == V.__reduce_ex__(1)
        True
    """
    if self._factory_data is None:
        raise NotImplementedError, "__reduce__ not implemented for %s" % type(self)
    else:
        return self._factory_data[0].reduce_data(self)

def lookup_global(name):
    """
    Used in unpickling the factory itself.

    EXAMPLES::

        sage: from sage.structure.factory import lookup_global
        sage: lookup_global('ZZ')
        Integer Ring
        sage: lookup_global('sage.rings.all.ZZ')
        Integer Ring
    """
    if '.' in name:
        module, name = name.rsplit('.', 1)
        all = __import__(module, fromlist=[name])
    else:
        import sage.all as all
    return getattr(all, name)


# To make the pickle jar happy:
from sage.structure.test_factory import test_factory
