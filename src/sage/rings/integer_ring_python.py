def iterator(self):
    r"""
    Iterate over all integers: 0 1 -1 2 -2 3 -3 ...

    This is the iterator for the integer ring `\ZZ`.

    EXAMPLES::

        sage: for n in ZZ: # indirect doctest
        ...    if n < 3: print n
        ...    else: break
        0
        1
        -1
        2
        -2
    """
    yield self(0)
    n = self(1)
    while True:
        yield n
        yield -n
        n += 1

