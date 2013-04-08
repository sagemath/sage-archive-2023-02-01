We explicitly raise an interrupt and expect it::

    sage: raise KeyboardInterrupt
    Traceback (most recent call last):
    ...
    KeyboardInterrupt

We now raise an interrupt without expecting it, which should lead to
an ordinary doctest failure::

    sage: raise KeyboardInterrupt
