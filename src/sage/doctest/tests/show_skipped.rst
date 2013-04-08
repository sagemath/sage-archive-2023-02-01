We test that optional tests are run correctly and the
``--show-skipped`` option displays the correctly::

    sage: 2 + 3
    5
    sage: 1 + 1 # optional - gap
    2
    sage: 1 - 1 # known bug
    17
    sage: 4 / 2 # long time
    2
    sage: 1 - 2 # optional - bug
    16
    sage: 8 + 1 # optional
    9
    sage: 1 / 0 # not tested
    I'm on fire!
