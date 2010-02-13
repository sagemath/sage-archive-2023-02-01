r"""
Special Methods for Classes.
"""
#*****************************************************************************
#  Copyright (C) 2009    Nicolas M. Thiery <nthiery at users.sf.net>
#  Copyright (C) 2010    Florent Hivert <Florent.Hivert at univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from nested_class import NestedClassMetaclass

class ClasscallMetaclass(NestedClassMetaclass):
    """
    A metaclass providing support for special methods for classes.

    From the Section ``Special method names`` of the Python Reference
    Manual: 'a class ``cls`` can implement certain operations on its
    instances that are invoked by special syntax (such as arithmetic
    operations or subscripting and slicing) by defining methods with
    special names'. The purpose of this metaclass is to allow the
    class ``cls`` to implement analogues of those special methods for
    the operations on the class itself.

    Currently, the following special methods are supported:

     - ``.__classcall__`` (and ``.__classcall_private__``) for
       customizing ``cls(...)`` (analogue of ``.__call__``).

     - ``.__classget__`` for customizing the binding behavior in
       ``foo.cls`` (analogue of ``.__get__``).

    See the documentation of :meth:`.__call__`` and of
    :meth:`.__get__`` for the description of the respective protocols.

    TODO: find a good name for this metaclass.

    AUTHORS:

     - Nicolas M. Thiery (2009-10) first implementation of __classcall__ and __classget__
     - Florent Hivert (2010-01): implementation of __classcall_private__, doc
    """

    def __get__(cls, instance, owner):
        """
        This method implements instance binding behavior for nested classes.

        Suppose that a class ``Outer`` contains a nested class ``cls`` which
        is an instance of this metaclass. For any object ``obj`` of ``cls``,
        this method implements a instance binding behavior for ``obj.cls`` by
        delegating it to ``cls.__classget__(Outer, obj, owner)`` if available.
        Otherwise, ``obj.cls`` results in ``cls``, as usual.

        Similarily, a class binding as in ``Outer.cls`` is delegated
        to ``cls.__classget__(Outer, None, owner)`` if available and
        to ``cls`` if not.

        For technical details, and in particular the description of
        the ``owner`` argument, see the Section ``Implementing
        Descriptors`` in the Python reference manual.

        EXAMPLES:

        We show how to implement a nested class ``Outer.Inner`` with a
        binding behavior, as if it was a method of ``Outer``: namely,
        for ``obj`` an instance of ``Outer``, calling
        ``obj.Inner(...)`` is equivalent to ``Outer.Inner(obj, ...)``::

            sage: import functools
            sage: from sage.misc.nested_class import NestedClassMetaclass
            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class Outer:
            ...       __metaclass__ = NestedClassMetaclass # workaround for python pickling bug
            ...
            ...       class Inner(object):
            ...           __metaclass__ = ClasscallMetaclass
            ...           @staticmethod
            ...           def __classget__(cls, instance, owner):
            ...               print "calling __classget__(%s, %s, %s)"%(
            ...                          cls, instance, owner)
            ...               if instance is None:
            ...                   return cls
            ...               return functools.partial(cls, instance)
            ...           def __init__(self, instance):
            ...               self.instance = instance
            sage: obj = Outer()
            sage: bar = obj.Inner()
            calling __classget__(<class '__main__.Outer.Inner'>, <__main__.Outer object at 0x...>, <class '__main__.Outer'>)
            sage: bar.instance == obj
            True

        Calling ``Outer.Inner`` returns the (unbinded) class as usual::

            sage: Inner = Outer.Inner
            calling __classget__(<class '__main__.Outer.Inner'>, None, <class '__main__.Outer'>)
            sage: Inner
            <class '__main__.Outer.Inner'>
            sage: type(bar) is Inner
            True

        ..warning:: Inner has to be a new style class (i.e. a subclass of object).

        ..warning:: calling ``obj.Inner`` does no longer return a class::

            sage: bind = obj.Inner
            calling __classget__(<class '__main__.Outer.Inner'>, <__main__.Outer object at 0x...>, <class '__main__.Outer'>)
            sage: bind
            <functools.partial object at 0x...>
        """
        if hasattr(cls, "__classget__"):
            return cls.__classget__(cls, instance, owner)
        else:
            return cls

    def __call__(cls, *args, **options):
        """
        This method implements ``cls(<some arguments>)``.

        Let ``cls`` be a class in ``ClasscallMetaclass``, and consider
        a call of the form:

            ``cls(<some arguments>)``

        If ``cls`` defines a method ``__classcall_private__``, then
        this results in a call to::

         - ``cls.__classcall_private__(cls, <some arguments>)``

        Otherwise, if ``cls`` has a method ``__classcall__``, then instead
        the following is called:

         - ``cls.__classcall__(cls, <some arguments>)``

        If neither of these two methods are implemented, then the standard
        ``type.__call__(cls, <some arguments>)`` is called, which in turn
        uses :meth:`__new__` and :meth:`__init__` as usual (see Section
        "Basic Customization" in the Python Reference Manual).

        EXAMPLES::

            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class Foo(object):
            ...       __metaclass__ = ClasscallMetaclass
            ...       @staticmethod
            ...       def __classcall__(cls):
            ...           print "calling classcall"
            ...           return type.__call__(cls)
            ...       def __new__(cls):
            ...           print "calling new"
            ...           return super(Foo, cls).__new__(cls)
            ...       def __init__(self):
            ...           print "calling init"
            ...
            sage: Foo()
            calling classcall
            calling new
            calling init
            <__main__.Foo object at ...>

       This behavior is inherited::

            sage: class Bar(Foo): pass
            sage: Bar()
            calling classcall
            calling new
            calling init
            <__main__.Bar object at ...>

       We now show the usage of :meth:`__classcall_private__`::

            sage: class FooNoInherits(object):
            ...       __metaclass__ = ClasscallMetaclass
            ...       @staticmethod
            ...       def __classcall_private__(cls):
            ...           print "calling private classcall"
            ...           return type.__call__(cls)
            ...
            sage: FooNoInherits()
            calling private classcall
            <__main__.FooNoInherits object at ...>

            sage: class BarNoInherits(FooNoInherits): pass
            sage: BarNoInherits()
            <__main__.BarNoInherits object at ...>

       We now show the usage of both::

            sage: class Foo2(object):
            ...       __metaclass__ = ClasscallMetaclass
            ...       @staticmethod
            ...       def __classcall_private__(cls):
            ...           print "calling private classcall"
            ...           return type.__call__(cls)
            ...       @staticmethod
            ...       def __classcall__(cls):
            ...           print "calling classcall with %s"%cls
            ...           return type.__call__(cls)
            ...
            sage: Foo2()
            calling private classcall
            <__main__.Foo2 object at ...>

            sage: class Bar2(Foo2): pass
            sage: Bar2()
            calling classcall with <class '__main__.Bar2'>
            <__main__.Bar2 object at ...>


        ..rubric: Discussion

        Typical applications include the implementation of factories or of
        unique representation (see :class:`UniqueRepresentation`). Such
        features are traditionaly implemented by either using a wrapper
        function, or fiddling with :meth:`__new__`.

        The benefit, compared with fiddling directly with :meth:`__new__`
        is a clear separation of the three distinct roles:

         - :meth:`cls.__classcall__`: what cls(<...>) does
         - :meth:`cls.__new__`: memory allocation for a *new* instance
         - :meth:`cls.__init__`: initialization of a newly created instance

        The benefit, compared with using a wrapper function, is that the
        user interface has a single handle for the class::

            sage: x = Partition([3,2,2])
            sage: isinstance(x, Partition)		# todo: not implemented

        instead of::

            sage: isinstance(x, sage.combinat.partition.Partition_class)
            True

        Another difference is that :meth:`__classcall__` is inherited by
        subclasses, which may be desirable, or not. If not, one should
        instead define the method :meth:`__classcall_private__` which will
        not be called for subclasses. Specifically, if a class ``cls``
        defines both methods ``__classcall__`` and
        ``__classcall_private__`` then, for any subclass ``sub`` of ``cls``:

         - ``cls(<args>)`` will call ``cls.__classcall_private__(cls, <args>)``
         - ``sub(<args>)`` will call ``cls.__classcall__(sub, <args>)``



        TESTS:

        We check that the idiom ``method_name in cls.__dict__`` works
        for extension types::

           sage: "_sage_" in SageObject.__dict__, "_sage_" in Parent.__dict__
           (True, False)
        """
        if '__classcall_private__' in cls.__dict__:
            return cls.__classcall_private__(cls, *args, **options)
        elif hasattr(cls, "__classcall__"):
            return cls.__classcall__(cls, *args, **options)
        else:
            return type.__call__(cls, *args, **options)
