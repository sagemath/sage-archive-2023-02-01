r"""
ClasscallMetaclass
"""
#*****************************************************************************
#  Copyright (C) 2008    Nicolas M. Thiery <nthiery at users.sf.net>
#  Copyright (C) 2009    Florent Hivert <Florent.Hivert at univ-rouen.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from nested_class import NestedClassMetaclass

class ClasscallMetaclass(NestedClassMetaclass):
    """
    A (trivial) metaclass for customizing class calls via a static method

    For a class ``cls`` in this metaclass, the call ``cls(<some arguments>)``
    triggers in this priority order the following calls if the respective
    methods exists:

     - ``cls.__classcall_private__(cls, <some arguments>)``
     - ``cls.__classcall__(cls, <some arguments>)``

    If neither of these two methods are implemented, then
    ``type.__call__(cls, <some arguments>)`` is called as usual, which
    in turns uses :meth:`__new__` and :meth:`__init__` as usual (see
    Section "Basic Customization" in the Python Reference Manual).

    See ``sage.misc.classcall_metaclass.ClasscallMetaclass.__call__?``
    for an example.

    Typical applications include the implementation of factories or of
    unique representation (see :class:`UniqueRepresentation`). This is
    usually done by either using a wrapper function, or fiddling with
    :meth:`__new__`.

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
    subclasses. If this is not desirable, one should define the method
    :meth:`__classcall_private__` which will not be called by
    subclass. Specifically, If a class ``cls`` define both methods
    ``__classcall__`` and ``__classcall_private__`` then, for any subclass
    ``sub`` of ``cls``

     - ``cls(<args>)`` will call ``cls.__classcall_private__(cls, <args>)``
     - ``sub(<args>)`` will call ``cls.__classcall__(sub, <args>)``

    AUTHORS:

     - Nicolas M. Thiery (2009-04) first release
     - Florent Hivert (2010-01) added __classcall_private__
    """

    def __call__(cls, *args, **options):
        """
        This method implements ``cls(<some arguments>)``.

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
        """
        if '__classcall_private__' in cls.__dict__:
            return cls.__classcall_private__(cls, *args, **options)
        elif hasattr(cls, "__classcall__"):
            return cls.__classcall__(cls, *args, **options)
        else:
            return type.__call__(cls, *args, **options)
