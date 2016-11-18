"""
Moved to sage.libs.eclib.interface

TESTS::

    sage: from sage.libs.mwrank.interface import mwrank_EllipticCurve
    doctest:...: DeprecationWarning: the module sage.libs.mwrank.interface has moved to sage.libs.eclib.interface
    See http://trac.sagemath.org/19818 for details.
"""

from sage.misc.superseded import deprecation
deprecation(19818, "the module sage.libs.mwrank.interface has moved to sage.libs.eclib.interface")

from sage.libs.eclib.interface import *
