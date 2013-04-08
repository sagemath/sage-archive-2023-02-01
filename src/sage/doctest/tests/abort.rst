This causes an "Unhandled SIGABRT..." which should be handled by
the doctester::

    sage: from sage.tests.interrupt import *
    sage: unguarded_abort()
