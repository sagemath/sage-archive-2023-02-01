r"""
UniqueRepresentation
"""
#*****************************************************************************
#  Copyright (C) 2008 Nicolas M. Thiery <nthiery at users.sf.net>
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

from sage.misc.cachefunc import cached_function
from sage.misc.classcall_metaclass import ClasscallMetaclass

class UniqueRepresentation:
    """
    Classes derived from UniqueRepresentation inherit a unique
    representation behavior for their instances.

    EXAMPLES:

    The short story: to construct a class whose instances have a
    unique representation behavior one just have to do::

        sage: class MyClass(UniqueRepresentation):
        ...       # all the rest as usual
        ...       pass

    Everything below is for the curious or for advanced usage.

    .. rubric:: What is unique representation?

    Instances of a class have a *unique representation behavior* when
    several instances constructed with the same arguments share the
    same memory representation. For example, calling twice::

        sage: f = GF(7)
        sage: g = GF(7)

    to create the finite field of order 7 actually gives back the same
    object::

        sage: f == g
        True
        sage: f is g
        True

    This is a standard design pattern. Besides saving memory, it
    allows for sharing cached data (say representation theoretical
    information about a group) as well as for further optimizations
    (fast hashing, equality testing). This behaviour is typically
    desirable for parents and categories. It can also be useful for
    intensive computations where one wants to cache all the operations
    on a small set of elements (say the multiplication table of a
    small group), and access this cache as quickly as possible.

    The :class:`UniqueRepresentation` and :class:`UniqueFactory` classes
    provide two alternative implementations of this design pattern. Both
    implementations have their own merits. :class:`UniqueRepresentation` is
    very easy to use: a class just needs to derive from it, or make sure some
    of its super classes does. For basic usage . Also, it groups together the
    class and the factory in a single gadget; in the example above, one would
    want to do::

        sage: isinstance(f, GF)         # todo: not implemented
        True

    but this does not work, because GF is only the factory.

    On the other hand the :class:`UniqueRepresentation` class is more
    intrusive, as it imposes a behavior (and a metaclass) to all the
    subclasses. Its implementation is also more technical, which leads
    to some subtleties.

    EXAMPLES:

    We start with a simple class whose constructor takes a single
    value as argument (TODO: find a more meaningful example)::

        sage: class MyClass(UniqueRepresentation):
        ...       def __init__(self, value):
        ...           self.value = value
        ...

    Two coexisting instances of MyClass created with the same
    argument data are guaranteed to share the same identity::

        sage: x = MyClass(1)
        sage: y = MyClass(1)
        sage: x is y
        True
        sage: z = MyClass(2)
        sage: x is z
        False

    In particular, modifying any one of them modifies the other
    (reference effect)::

        sage: x.value = 3
        sage: x.value, y.value
        (3, 3)
        sage: y.value = 1
        sage: x.value, y.value
        (1, 1)

    Unless overridden by the derived class, equality testing is
    implemented by comparing identities, which is as fast as it can get::

        sage: x == y
        True
        sage: z = MyClass(2)
        sage: x == z, x is z
        (False, False)

    Similarly, the identity is used as hash function, which is also as
    fast as it can get. However this means that the hash function may
    change in between Sage sessions, or even within the same Sage
    session (see below). Subclasses should overload :meth:`__hash__`
    if this could be a problem.


    The arguments can consist of any combination of positional or
    keyword arguments, as taken by a usual :meth:`__init__`
    function. However, all values passed in should be hashable::

        sage: MyClass(value = [1,2,3])
        Traceback (most recent call last):
        ...
        TypeError: unhashable type: 'list'

    .. rubric:: Argument preprocessing

    Sometimes, one wants to do some preprocessing on the arguments, to
    put them in some canonical form. The following example illustrates
    how to achieve this; it takes as argument any iterable, and
    canonicalizes it into a tuple (which is hashable!)::

        sage: class MyClass2(UniqueRepresentation):
        ...       @staticmethod
        ...       def __classcall__(cls, iterable):
        ...           t = tuple(iterable)
        ...           return super(MyClass2, cls).__classcall__(cls, t)
        ...
        ...       def __init__(self, value):
        ...           self.value = value
        ...
        sage: x = MyClass2([1,2,3])
        sage: y = MyClass2(tuple([1,2,3]))
        sage: z = MyClass2(i for i in [1,2,3])
        sage: x.value
        (1, 2, 3)
        sage: x is y, y is z
        (True, True)

    A similar situation arises when the constructor accepts default
    values for some of its parameters. Alas, the obvious
    implementation does not work::

        sage: class MyClass3(UniqueRepresentation):
        ...       def __init__(self, value = 3):
        ...           self.value = value
        ...
        sage: MyClass3(3) is MyClass3()
        False

    Instead, one should do::

        sage: class MyClass3(UniqueRepresentation):
        ...       @staticmethod
        ...       def __classcall__(cls, value = 3):
        ...           return super(MyClass3, cls).__classcall__(cls, value)
        ...
        ...       def __init__(self, value):
        ...           self.value = value
        ...
        sage: MyClass3(3) is MyClass3()
        True

    A bit of explanation is in order. First, the call
    ``MyClass2([1,2,3])`` triggers a call to
    ``MyClass2.__classcall__(MyClass2, [1,2,3])``. This is an extension of
    the standard Python behavior, needed by :class:`UniqueRepresentation`,
    and implemented by the ``ClasscallMetaclass``. Then,
    ``MyClass2.__classcall__`` does the desired transformations on the
    arguments. Finally, it uses ``super`` to call the default
    implementation of ``__classcall__`` provided by
    :class:`UniqueRepresentation`. This one in turn handles the caching
    and, if needed, constructs and initializes a new object in the
    class using :meth:`__new__` and :meth:`__init__` as usual.

    Constraints:

     - :meth:`__classcall__` is a staticmethod (like, implicitly, :meth:`__new__`)

     - the preprocessing on the arguments should be
       idempotent. Namely, If :meth:`MyClass2.__classcall__` calls
       ``UniqueRepresentation.__classcall__(<some_arguments>)``, then
       it should accept <some_arguments> as its own input, and pass it
       down unmodified to :meth:`UniqueRepresentation.__classcall__`.
     - :meth:`MyClass2.__classcall__` should return the result of
       :meth:`UniqueRepresentation.__classcall__` without modifying it.

    Other than that :meth:`MyClass2.__classcall__` may play any tricks,
    like acting as a Factory and returning object from other classes.

    .. rubric:: Unique representation and mutability

    :class:`UniqueRepresentation` is primarily intended for implementing
    objects which are (at least semantically) immutable. This is in
    particular assumed by the default implementations of ``copy`` and
    ``deepcopy``::

        sage: copy(x) is x
        True
        sage: from copy import deepcopy
        sage: deepcopy(x) is x
        True

    Using :class:`UniqueRepresentation` on mutable objects may lead to
    subtle behavior::

        sage: t = MyClass(3)
        sage: z = MyClass(2)
        sage: t.value = 2

    Now x and z have the same data structure, but are not considered
    as equal::

        sage: t.value == z.value
        True
        sage: t == z
        False

    .. rubric:: More on unique representation and identity

    :class:`UniqueRepresentation` is implemented by mean of a cache. This
    cache uses weak references so that, when all other references to,
    say, ``MyClass(1)`` have been deleted, the instance is actually
    deleted from memory. A later call to ``MyClass(1)`` reconstructs the
    instance from scratch, *most likely with a different id*.

    TODO: add an example illustrating this behavior


    .. rubric:: Unique representation and pickling

    The default Python pickling implementation (by reconstructing an
    object from its class and dictionary, see "The pickle protocol" in
    the Python Library Reference) does not preserves unique
    representation, as Python has no chance to know whether and where
    the same object already exists.

    :class:`UniqueRepresentation` tries to ensure appropriate pickling by
    implementing a :meth:`__reduce__` method returning the arguments
    passed to the constructor::

        sage: import __main__             # Fake MyClass being defined in a python module
        sage: __main__.MyClass = MyClass
        sage: x = MyClass(1)
        sage: loads(dumps(x)) is x
        True

    :class:`UniqueRepresentation` uses the :meth:`__reduce__` pickle
    protocol rather than :meth:`__getnewargs__` because the later does
    not handle keyword arguments::

        sage: x = MyClass(value = 1)
        sage: x.__reduce__()
        (<function unreduce at ...>, (<class '__main__.MyClass'>, (), {'value': 1}))
        sage: x is loads(dumps(x))
        True

    Caveat: the default implementation of :meth:`__reduce__` in
    :class:`UniqueRepresentation` requires to store the constructor's
    arguments in the instance dictionary upon construction:

        sage: x.__dict__
        {'_reduction': (<class '__main__.MyClass'>, (), {'value': 1}), 'value': 1}

    It is often easy in a derived subclass to reconstruct the
    constructors arguments from the instance data structure. When this
    is the case, :meth:`__reduce__` should be overridden; automagically
    the arguments won't be stored anymore:

        sage: class MyClass3(UniqueRepresentation):
        ...       def __init__(self, value):
        ...           self.value = value
        ...
        ...       def __reduce__(self):
        ...           return (MyClass3, (self.value,))
        ...
        sage: import __main__; __main__.MyClass3 = MyClass3  # Fake MyClass3 being defined in a python module
        sage: x = MyClass3(1)
        sage: loads(dumps(x)) is x
        True
        sage: x.__dict__
        {'value': 1}

    .. rubric:: Migrating classes to ``UniqueRepresentation`` and unpickling

    We check that, when migrating a class to ``UniqueRepresentation``,
    older pickle can still be reasonably unpickled. Let us create a
    (new style) class, and pickle one of its instances::

        sage: class MyClass4(object):
        ...       def __init__(self, value):
        ...           self.value = value
        ...
        sage: import __main__; __main__.MyClass4 = MyClass4  # Fake MyClass4 being defined in a python module
        sage: pickle = dumps(MyClass4(1))

    It can be unpickled::

        sage: y = loads(pickle)
        sage: y.value
        1

    Now, we upgrade the class to derive from ``UniqueRepresentation``::

        sage: class MyClass4(UniqueRepresentation, object):
        ...       def __init__(self, value):
        ...           self.value = value
        sage: import __main__; __main__.MyClass4 = MyClass4  # Fake MyClass4 being defined in a python module
        sage: __main__.MyClass4 = MyClass4

    The pickle can still be unpickled::

        sage: y = loads(pickle)
        sage: y.value
        1

    Note however that, for the reasons explained above, unique
    representation is not guaranteed in this case::

        sage: y is MyClass4(1)
        False

    Todo: illustrate how this can be fixed on a case by case basis.


    Now, we redo the same test for a class deriving from SageObject::

        sage: class MyClass4(SageObject):
        ...       def __init__(self, value):
        ...           self.value = value
        sage: import __main__; __main__.MyClass4 = MyClass4  # Fake MyClass4 being defined in a python module
        sage: pickle = dumps(MyClass4(1))

        sage: class MyClass4(UniqueRepresentation, SageObject):
        ...       def __init__(self, value):
        ...           self.value = value
        sage: __main__.MyClass4 = MyClass4
        sage: y = loads(pickle)
        sage: y.value
        1

    Caveat: unpickling instances of a formerly old-style class is not supported yet by default::

        sage: class MyClass4:
        ...       def __init__(self, value):
        ...           self.value = value
        sage: import __main__; __main__.MyClass4 = MyClass4  # Fake MyClass4 being defined in a python module
        sage: pickle = dumps(MyClass4(1))

        sage: class MyClass4(UniqueRepresentation, SageObject):
        ...       def __init__(self, value):
        ...           self.value = value
        sage: __main__.MyClass4 = MyClass4
        sage: y = loads(pickle)  # todo: not implemented
        sage: y.value            # todo: not implemented
        1



    .. rubric:: Rationale for the current implementation

    :class:`UniqueRepresentation` and derived classes use the
    ``ClasscallMetaclass`` of the standard Python type. The following
    example explains why.

    We define a variant of ``MyClass`` where the calls to :meth:`__init__`
    are traced::

        sage: class MyClass(UniqueRepresentation):
        ...       def __init__(self, value):
        ...           print "initializing object"
        ...           self.value = value
        ...

    Let us create an object twice::

        sage: x = MyClass(1)
        initializing object
        sage: z = MyClass(1)

    As desired the __init__ method was only called the first time,
    which is an important feature.

    As far as we can tell, this is not achievable while just using
    :meth:`__new__` and :meth:`__init__` (as defined by type; see
    Section "Basic Customization" in the Python Reference
    Manual). Indeed, :meth:`__init__` is called systematically on the
    result of :meth:`__new__` whenever the result is an instance of
    the class.

    Another difficulty is that argument preprocessing (as in the
    example above) cannot be handled by :meth:`__new__`, since the
    unprocessed arguments will be passed down to :meth:`__init__`.

    .. rubric: Mixing super types and super classes

    TESTS:

    For the record, this test did fail with previous implementation
    attempts::

        sage: class bla(UniqueRepresentation, SageObject):
        ...       pass
        ...
        sage: b = bla()

    """

    __metaclass__ = ClasscallMetaclass

    @cached_function # automatically a staticmethod
    def __classcall__(cls, *args, **options):
        """
        Constructs a new object of this class or reuse an existing one.

        See also :class:`UniqueRepresentation` for a discussion.

        EXAMPLES::

            sage: x = UniqueRepresentation()
            sage: y = UniqueRepresentation()
            sage: x is y
            True
        """
        instance = type.__call__(cls, *args, **options)
        assert isinstance( instance, cls )
        if instance.__class__.__reduce__ == UniqueRepresentation.__reduce__:
            instance._reduction = (cls, args, options)
        return instance

    # Should be cythoned
    def __eq__(self, other):
        """
        Test if ``self`` and ``other` are equal by comparing their
        identity.

        See also :class:`UniqueRepresentation` for a discussion.

        EXAMPLES::

            sage: x = UniqueRepresentation()
            sage: y = UniqueRepresentation()
            sage: x == y
            True
            sage: x is y
            True
            sage: x == 3
            False
        """
        return self is other

    # Should be cythoned
    def __hash__(self):
        """
        Returns the hash value of ``self``.

        See also :class:`UniqueRepresentation` for a discussion.

        EXAMPLES::

            sage: x = UniqueRepresentation()
            sage: y = UniqueRepresentation()
            sage: hash(x) # random
            74153040
            sage: type(hash(x))
            <type 'int'>
            sage: hash(x) == hash(y)
            True
            sage: class bla(UniqueRepresentation, SageObject): pass
            sage: x = bla(); hx = hash(x)
            sage: x.rename("toto")
            sage: hx == hash(x)
            True
        """
        return id(self)

    def __reduce__(self):
        """
        Returns the argument that have been passed to :meth:`__new__`
        to construct this object, as per the pickle protocol.

        See also :class:`UniqueRepresentation` for a discussion.

        EXAMPLES::

            sage: x = UniqueRepresentation()
            sage: x.__reduce__()
            (<function unreduce at ...>, (<class 'sage.structure.unique_representation.UniqueRepresentation'>, (), {}))
        """
        return (unreduce, self._reduction)

    def __copy__(self):
        """
        Returns self, as a semantic copy of self

        This assume that the object is semantically immutable.

        EXAMPLES::

            sage: x = UniqueRepresentation()
            sage: x is copy(x)
            True
        """
        return self

    def __deepcopy__(self, memo):
        """
        Returns self, as a semantic deep copy of self

        This assume that the object is semantically immutable.

        EXAMPLES::

            sage: from copy import deepcopy
            sage: x = UniqueRepresentation()
            sage: x is deepcopy(x)
            True
        """
        return self

def unreduce(cls, args, keywords):
    """
    Calls a class on the given arguments::

        sage: sage.structure.unique_representation.unreduce(Integer, (1,), {})
        1

    Todo: should reuse something preexisting ...
    """
    return cls(*args, **keywords)
