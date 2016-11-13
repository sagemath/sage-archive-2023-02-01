"""
Moved to sage.libs.eclib.newforms

TESTS::

    sage: from sage.libs.cremona.newforms import ECModularSymbol
    doctest:...: DeprecationWarning: the module sage.libs.cremona.newforms has moved to sage.libs.eclib.newforms
    See http://trac.sagemath.org/19818 for details.
    sage: E = EllipticCurve('11a')
    sage: ECModularSymbol(E,1)
    Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
"""

from sage.misc.superseded import deprecation
deprecation(19818, "the module sage.libs.cremona.newforms has moved to sage.libs.eclib.newforms")

from sage.libs.eclib.newforms import *
