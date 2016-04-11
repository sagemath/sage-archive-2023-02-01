"""
Moved to sage.libs.eclib.homspace

TESTS::

    sage: from sage.libs.cremona.homspace import ModularSymbols
    doctest:...: DeprecationWarning: the module sage.libs.cremona.homspace has moved to sage.libs.eclib.homspace
    See http://trac.sagemath.org/19818 for details.
    sage: ModularSymbols(225)
    Cremona Modular Symbols space of dimension 61 for Gamma_0(225) of weight 2 with sign 0
"""

from sage.misc.superseded import deprecation
deprecation(19818, "the module sage.libs.cremona.homspace has moved to sage.libs.eclib.homspace")

from sage.libs.eclib.homspace import *
