def python(z, locals=None):
    """
    Return the closest Python/Sage equivalent of the given PARI object.

    This is deprecated, use the ``python`` method of :class:`gen`
    instead.

    TESTS::

        sage: from sage.libs.pari.gen_py import python
        sage: python(pari(3))
        doctest:...: DeprecationWarning: gen_py.python is deprecated, use sage.libs.pari.gen.gentoobj or the .python() method instead
        See http://trac.sagemath.org/19888 for details.
        3
    """
    from sage.misc.superseded import deprecation
    deprecation(19888, 'gen_py.python is deprecated, use sage.libs.pari.gen.gentoobj or the .python() method instead')
    from sage.libs.pari.gen import gentoobj
    return gentoobj(z, locals)
