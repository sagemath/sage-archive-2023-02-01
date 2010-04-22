"""
Lazy format strings
"""
from copy import copy

class LazyFormat(str):
    """
    Lazy format string which only calls ``__repr__()`` and ``__str__()`` of
    their parameters if printed.

    The class :class:`LazyFormat` allows to create format strings which calls
    their argument's ``__repr__`` only if needed.  Otherwise it behaves as
    usual format string.

    EXAMPLES:

    Under normal circumstances, :class:`Lazyformat` strings behaves as usual::

        sage: from sage.misc.lazy_format import LazyFormat
        sage: LazyFormat("<%s>.+.")%"toto"
        <toto>.+.
        sage: LazyFormat("<%s>.+.<%s>")%("toto", 3)
        <toto>.+.<3>

    Lets build an object with a broken ``__repr__`` method::

        sage: class IDontLikeBeingPrinted(object):
        ...    def __repr__(self):
        ...        raise ValueError, "Don't ever try to print me !"

    There is no error when binding a lazy format with the broken object::

        sage: lf = LazyFormat("<%s>")%IDontLikeBeingPrinted()

    The error occurs only when printed::

        sage: lf
        Traceback (most recent call last):
        ...
        ValueError: Don't ever try to print me !


    .. rubric:: Common Usage:

    A priori, since they are mostly called during user interaction, there is
    no particular need to write fast ``__repr__`` methods, and indeed there
    are some objects ``x`` whose call ``x.__repr__()`` is quite time
    expensive. However, there are several use case where it is actually called
    and when the results is simply discarded. In particular writing ``"%s"%x``
    actually calls ``x.__repr__()`` even the results is not printed. A typical
    use case is during tests when there is a tester which tests an assertion
    and prints an error message on failure::

        sage: QQ._tester().assertTrue(0 in QQ,
        ...                "%s doesn't contain 0"%QQ)

    I the previous case ``QQ.__repr__()`` has been called. To show this we
    replace QQ in the format string argument with out broken object::

        sage: QQ._tester().assertTrue(True,
        ...                "%s doesn't contain 0"%IDontLikeBeingPrinted())
        Traceback (most recent call last):
        ...
        ValueError: Don't ever try to print me !

    There is no need to call ``IDontLikeBeingPrinted().__repr__()``, but itx
    can't be avoided with usual strings. Note that with the usual assert, the
    call is not performed::

        sage: assert True, "%s is wrong"%IDontLikeBeingPrinted()

    We now check that :class:`LazyFormat` indeed solve the assertion problem::

        sage: QQ._tester().assertTrue(True,
        ...               LazyFormat("%s is wrong")%IDontLikeBeingPrinted())
        sage: QQ._tester().assertTrue(False,
        ...               LazyFormat("%s is wrong")%IDontLikeBeingPrinted())
        Traceback (most recent call last):
        ...
        AssertionError: <unprintable AssertionError object>
    """

    def __mod__(self, args):
        """
        Binds the lazy format with its parameters

        EXAMPLES::

            sage: from sage.misc.lazy_format import LazyFormat
            sage: form = LazyFormat("<%s>");
            sage: form
            unbound LazyFormat("<%s>")
            sage: form%"params"
            <params>
        """
        if hasattr(self, "_args"): # self is already bound...
            self = copy(self)
        self._args = args
        return self

    def __repr__(self):
        """
        TESTS::

            sage: from sage.misc.lazy_format import LazyFormat
            sage: form = LazyFormat("<%s>");
            sage: form
            unbound LazyFormat("<%s>")
            sage: print form
            unbound LazyFormat("<%s>")
            sage: form%"toto"
            <toto>
            sage: print form%"toto"
            <toto>
        """
        try:
            args = self._args
        except AttributeError:
            return "unbound LazyFormat(\""+self+"\")"
        else:
            return str.__mod__(self, args)

    __str__ = __repr__
