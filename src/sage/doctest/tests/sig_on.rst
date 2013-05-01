Test a bad value of ``sig_on_count``::

    sage: cython('sig_on()')

The following test should succeed as usual::

    sage: 3**12 - 2**19
    7153
