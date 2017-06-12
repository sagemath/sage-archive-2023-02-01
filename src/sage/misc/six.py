# -*- encoding: utf-8 -*-
"""
Python 2 and 3 Compatibility
"""

from __future__ import absolute_import

from six import *


def with_metaclass(meta, *bases):
    """
    Create a base class with a metaclass.

    TESTS::

        sage: from sage.misc.six import with_metaclass
        sage: class Meta(type): pass
        sage: class X(with_metaclass(Meta)): pass
        sage: type(X) is Meta
        True
        sage: issubclass(X, object)
        True
        sage: class Base(object): pass
        sage: class X(with_metaclass(Meta, Base)): pass
        sage: type(X) is Meta
        True
        sage: issubclass(X, Base)
        True
        sage: class Base2(object): pass
        sage: class X(with_metaclass(Meta, Base, Base2)): pass
        sage: type(X) is Meta
        True
        sage: issubclass(X, Base)
        True
        sage: issubclass(X, Base2)
        True
        sage: X.__mro__ == (X, Base, Base2, object) or X.__mro__
        True

    Check that :trac:`18503` is fixed, i.e. that with_metaclass
    works with cdef'ed metaclasses::

        sage: from sage.misc.classcall_metaclass import ClasscallMetaclass
        sage: class X(with_metaclass(ClasscallMetaclass)): pass
        sage: type(X) is ClasscallMetaclass
        True
        sage: X.__mro__ == (X, object) or X.__mro__
        True

    Check a fix for :trac:`16074`::

        sage: from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
        sage: from sage.modules.with_basis.morphism import ModuleMorphismByLinearity
        sage: from sage.structure.unique_representation import UniqueRepresentation
        sage: class ExteriorAlgebraDifferential(with_metaclass(
        ....:        InheritComparisonClasscallMetaclass,
        ....:        ModuleMorphismByLinearity, UniqueRepresentation
        ....:     )):
        ....:     pass
    """
    # This requires a bit of explanation: with_metaclass() returns
    # a dummy class which will setup the correct metaclass when this
    # dummy class is used as base class.
    #
    # In detail: let T = with_metaclass(meta, *bases). When the user does
    # "class X(with_metaclass(meta, *bases))", Python will first determine
    # the metaclass of X from its bases. So, the metaclass will be
    # type(T).
    # Next, Python will call type(T)("X", (T,), methods) and it is this
    # call that we want to override. So we need to define a __call__
    # method in the metaclass of type(T), which needs a metametaclass,
    # called "_WithMetaclass".
    # The metaclass type(T) is returned by _WithMetaclass.get(...) and
    # the dummy() method creates an instance of this metaclass, which
    # is a regular class.
    return _WithMetaclass.get(meta, bases).dummy()


class _WithMetaclass(type):
    """
    Metaclasses to be used in ``with_metaclass``.
    """
    @classmethod
    def get(cls, meta, bases):
        """
        Return a metaclass which, if an instance of it is used as base
        class, will create a class with metaclass ``meta`` and bases
        ``bases``.
        """
        return cls("metaclass", (type,), dict(meta=meta, bases=bases))

    def dummy(self):
        """
        Return a dummy instance of this metaclass.
        """
        return self.__new__(self, "dummy", (object,), {})

    def __call__(self, name, _unused_bases, d):
        """
        Create the eventual class with metaclass ``self.meta`` and bases
        ``self.bases``.
        """
        return self.meta(name, self.bases, d)


def u(x):
    r"""
    Convert `x` to unicode, assuming UTF-8 encoding.

    Python2 behaviour:

    If input is unicode, returns the input.

    If input is str (assumed to be utf-8 encoded), convert to unicode.

    Python3 behaviour:

    If input is str, returns the input.

    If input is bytes (assumed to be utf-8 encoded), convert to unicode.

    EXAMPLES::

        sage: from sage.misc.six import u
        sage: u("500 €")
        u'500 \u20ac'
        sage: u(u"500 \u20ac")
        u'500 \u20ac'
    """
    if isinstance(x, text_type):  # py2 unicode and py3 str
        return x
    if isinstance(x, bytes):
        return x.decode("utf-8")
    raise TypeError('input has no conversion to unicode')
