r"""
Singleton categories
"""
#*****************************************************************************
#  Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#                     Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.constant_function import ConstantFunction
from sage.misc.lazy_attribute import lazy_attribute, lazy_class_attribute
from sage.categories.category import Category
from sage.structure.category_object cimport CategoryObject
from sage.structure.unique_representation import UniqueRepresentation

include "../ext/python_type.pxi"

# This helper class is used to implement Category_singleton.__contains__
# In particular, the docstring is what appears upon C.__contains__?
# for C a singleton category like Fields().

cdef class Category_contains_method_by_parent_class:
    """
    Returns whether ``x`` is an object in this category.

    More specifically, returns ``True`` if and only if ``x`` has a
    category which is a subcategory of this one.

    EXAMPLES::

        sage: ZZ in Sets()
        True
    """

    def __init__(self, category):
        """
        TESTS::

            sage: from sage.categories.category_singleton import Category_contains_method_by_parent_class
            sage: Category_contains_method_by_parent_class(Rings())
            <sage.categories.category_singleton.Category_contains_method_by_parent_class object at ...
        """
        self._parent_class_of_category = <type> category.parent_class

    def __call__(self, x):
        """
        EXAMPLES::

            sage: from sage.categories.category_singleton import Category_contains_method_by_parent_class
            sage: in_Fields = Category_contains_method_by_parent_class(Fields())
            sage: in_Fields(QQ)
            True
            sage: in_Fields(ZZ)
            False
            sage: in_Fields(1)              # Not a CategoryObject
            False
            sage: in_Fields(int(1))         # Not a SageObject
            False

        TESTS:

            The following used to segfault in a preliminary version of the
            code::

                sage: None in Rings()
                False
        """
        if x is None:
            return False
        cdef CategoryObject y
        try:
            y = x
            return PyType_IsSubtype(<type>((y._category or y.category()).parent_class), self._parent_class_of_category)
        except AttributeError:
            return False
        except TypeError: # this is for objects that aren't CategoryObjects
            try:
                return PyType_IsSubtype(<type>(x.category().parent_class), self._parent_class_of_category)
            except AttributeError:
                return False

cdef class FastHashable_class:
    """
    A class that has a fast hash method, returning a pre-assigned value.

    NOTE:

    This is for internal use only. The class has a cdef attribute ``_hash``,
    that needs to be assigned.

    TESTS::

        sage: issubclass(Rings, sage.categories.category_singleton.FastHashable_class)
        True

    """
    def __hash__(self):
        """
        TESTS::

            sage: hash(Rings()) == id(Rings)    # indirect doctest
            True

        """
        return self._hash


class Category_singleton(FastHashable_class,Category):
    """
    A base class for implementing singleton category

    A *singleton* category is a category whose class takes no
    parameters like ``Fields()`` or ``Rings()``. See also the
    `Singleton design pattern <http://en.wikipedia.org/wiki/Singleton_pattern>`_.

    This is a subclass of :class:`Category`, with a couple
    optimizations for singleton categories.

    The main purpose is to make the idioms:

        sage: QQ in Fields()
        True
        sage: ZZ in Fields()
        False

    as fast as possible, and in particular competitive to calling a
    constant Python method, in order to foster its systematic use
    throughout the Sage library. Such tests are time critical, in
    particular when creating a lot of polynomial rings over small
    fields like in the elliptic curve code.

    EXAMPLES::

        sage: from sage.categories.category_singleton import Category_singleton
        sage: class MyRings(Category):
        ...        super_categories = Rings.__dict__['super_categories']
        sage: class MyRingsSingleton(Category_singleton):
        ...        super_categories = Rings.__dict__['super_categories']

    We create three rings. One of them is contained in the usual category of
    rings, one in the category of "my rings" and the third in the category of
    "my rings singleton"::

        sage: R = QQ['x,y']
        sage: R1 = Parent(category = MyRings())
        sage: R2 = Parent(category = MyRingsSingleton())
        sage: R in MyRings()
        False
        sage: R1 in MyRings()
        True
        sage: R1 in MyRingsSingleton()
        False
        sage: R2 in MyRings()
        False
        sage: R2 in MyRingsSingleton()
        True

    One sees that containment tests for the singleton class is a lot faster
    than for a usual class::

        sage: timeit("R in MyRings()", number=10000)                  # not tested
        10000 loops, best of 3: 7.12 µs per loop
        sage: timeit("R1 in MyRings()", number=10000)                 # not tested
        10000 loops, best of 3: 6.98 µs per loop
        sage: timeit("R in MyRingsSingleton()", number=10000)         # not tested
        10000 loops, best of 3: 3.08 µs per loop
        sage: timeit("R2 in MyRingsSingleton()", number=10000)        # not tested
        10000 loops, best of 3: 2.99 µs per loop

    So this is an improvement, but not yet competitive with a pure
    Cython method:

        sage: timeit("R.is_ring()", number=10000)                     # not tested
        10000 loops, best of 3: 383 ns per loop

    However, it is competitive with a Python method. Actually it is faster,
    if one stores the category in a variable::

        sage: _Rings = Rings()
        sage: R3 = Parent(category = _Rings)
        sage: R3.is_ring.__module__
        'sage.categories.rings'
        sage: timeit("R3.is_ring()", number=10000)                    # not tested
        10000 loops, best of 3: 2.64 µs per loop
        sage: timeit("R3 in Rings()", number=10000)                   # not tested
        10000 loops, best of 3: 3.01 µs per loop
        sage: timeit("R3 in _Rings", number=10000)                    # not tested
        10000 loops, best of 3: 652 ns per loop

    This might not be easy to further optimize, since the time is
    consumed in many different spots::

        sage: timeit("MyRingsSingleton.__classcall__()", number=10000)# not tested
        10000 loops, best of 3: 306 ns per loop

        sage: X = MyRingsSingleton()
        sage: timeit("R in X  ", number=10000)                        # not tested
        10000 loops, best of 3: 699 ns per loop

        sage: c = MyRingsSingleton().__contains__
        sage: timeit("c(R)", number = 10000)                          # not tested
        10000 loops, best of 3: 661 ns per loop

    .. warning::

        A singleton concrete class `A` cannot have a subclass `B`
        (necessarily concrete). Otherwise, creating an instance `a` of
        `A` and an instance `b` of `B` would break the singleton
        principle: `A` would have two instances `a` and `b`.

        This implementation currently requires that:

        - a subclass of :class:`Category_singleton` should not have subclasses

        - :class:`Category_singleton` should not be instantiated.

        For example, subclassing our rings singleton produces wrong results::

            sage: class Disaster(MyRingsSingleton): pass
            sage: Disaster()   # Note that this is NOT "Category of disaster"!
            Category of my rings singleton
            sage: Disaster() is MyRingsSingleton()
            True

        For safety, we made sure that instanciating
        :class:`Category_singleton` does not break subclasses, by not
        caching the result::

            sage: Category_singleton() is Category_singleton()
            False
            sage: class NewSubclass(Category_singleton): pass
            sage: NewSubclass()
            Category of new subclass
            sage: NewSubclass() is NewSubclass()
            True

    TESTS::

        sage: import __main__
        sage: __main__.MyRings = MyRings
        sage: __main__.MyRingsSingleton = MyRingsSingleton
        sage: TestSuite(MyRingsSingleton()).run()
    """

    # That is just an optimized constant cached_method
    @lazy_class_attribute
    def __classcall__(object cls):
        """
        TESTS::

            sage: from sage.categories.category_singleton import Category_singleton
            sage: class MyRingsSingleton(Category_singleton):
            ...        super_categories = Rings.__dict__['super_categories']
            sage: MyRingsSingleton()
            Category of my rings singleton

        Note that when calling :class:`Category_singleton` directly
        then there is no cache used::

            sage: Category_singleton() is Category_singleton()
            False

        You must **never** subclass a subclass of :class:`Category_singleton`.
        Otherwise, the instance the sub-sub-class may become identical with
        the instance of the sub-class::

            sage: class MyStuff(Category_singleton): pass
            sage: MyStuff()
            Category of my stuff
            sage: class MySubStuff(MyStuff): pass
            sage: MySubStuff()
            Category of my stuff
        """
        cdef FastHashable_class obj
        if cls.__mro__[1] is not Category_singleton:
            # Actually this type error is invisible. But it makes sure that
            # the __classcall__ for Category_singleton is not overridden
            # when someone is calling it.
            raise TypeError, "%s is not a direct subclass of %s"%(cls,Category_singleton)
        obj = UniqueRepresentation.__classcall__(cls)
        obj._hash = <Py_ssize_t><void *>cls
        return ConstantFunction(obj)

#    @lazy_class_attribute
#    def __hash__(cls):
#        """
#        The hash of the unique instance of a category singleton is the address of its class.
#
#        EXAMPLES::
#
#            sage: hash(Rings()) == hash(Rings)     # indirect doctest
#            True
#
#        """
#        return ConstantFunction(id(cls))

    @lazy_class_attribute
    def __contains__(cls):
        """
        TESTS::

            sage: from sage.categories.category_singleton import Category_singleton
            sage: class MyRingsSingleton(Category_singleton):
            ...        super_categories = Rings.__dict__['super_categories']
            sage: ZZ in MyRingsSingleton()
            False
            sage: Parent(category=MyRingsSingleton()) in MyRingsSingleton()
            True
        """
        return Category_contains_method_by_parent_class(cls())
