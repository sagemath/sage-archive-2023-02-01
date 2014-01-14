We test that optional tests are run correctly and the
``--show-skipped`` option displays the correctly.
We also test case insensitivity::

    sage: 2 + 3
    5
    sage: 1 + 1 # Optional --- Gap
    2
    sage: 1 - 1 # Known Bug
    17
    sage: 4 / 2 # Long Time
    2
    sage: 1 - 2 # Optional - Bug
    16
    sage: 8 + 1 # Optional
    9
    sage: 1 / 0 # Not Tested
    I'm on fire!
