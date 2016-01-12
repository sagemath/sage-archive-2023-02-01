"""
Moved to sage.libs.eclib.all

TESTS::

    sage: from sage.libs.cremona.all import get_precision
    doctest:...: DeprecationWarning: the module sage.libs.cremona.all has moved to sage.libs.eclib.all
    See http://trac.sagemath.org/19818 for details.
    sage: get_precision()
    50
"""

from sage.misc.superseded import deprecation
deprecation(19818, "the module sage.libs.cremona.all has moved to sage.libs.eclib.all")

from sage.libs.eclib.all import *
