"""
Moved to sage.libs.cremona.mwrank

TESTS::

    sage: from sage.libs.mwrank.mwrank import get_precision
    doctest:...: DeprecationWarning: the module sage.libs.mwrank.mwrank has moved to sage.libs.cremona.mwrank
    See http://trac.sagemath.org/19818 for details.
    sage: get_precision()
    50
"""

from sage.misc.superseded import deprecation
deprecation(19818, "the module sage.libs.mwrank.mwrank has moved to sage.libs.cremona.mwrank")

from sage.libs.cremona.mwrank import *
