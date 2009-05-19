r"""
ClasscallMetaclass
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

class ClasscallMetaclass(type):
    """
    A (trivial) metaclass for customizing class calls via a static method

    For a class ``cls`` in this metaclass, ``cls(<some arguments>)``
    triggers a call to ``cls.__classcall__(cls, <some arguments>)``.
    If ``cls.__classcall__`` is not implemented, then
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
    subclasses. Depending on the context this may, or may not, be
    desirable.

    Todo: any better name for ``__classcall__``?
    """

    def __call__(cls, *args, **options):
        """
        This method implements ``cls(<some arguments>)``.

        EXAMPLES::

            sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
            sage: class Foo(object):
            ...       __metaclass__ = sage.structure.unique_representation.ClasscallMetaclass
            ...       @staticmethod
            ...       def __classcall__(cls):
            ...           print "calling call_cls"
            ...           return type.__call__(cls)
            ...       def __new__(cls):
            ...           print "calling new"
            ...           return super(Foo, cls).__new__(cls)
            ...       def __init__(self):
            ...           print "calling init"
            ...
            sage: Foo()
            calling call_cls
            calling new
            calling init
            <__main__.Foo object at ...>
        """
        try:
            return cls.__classcall__(cls, *args, **options)
        except AttributeError:
            return type.__call__(cls, *args, **options)
