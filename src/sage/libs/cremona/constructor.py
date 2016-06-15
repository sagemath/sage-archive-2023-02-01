"""
Moved to sage.libs.eclib.constructor

TESTS::

    sage: from sage.libs.cremona.constructor import CremonaModularSymbols
    doctest:...: DeprecationWarning: the module sage.libs.cremona.constructor has moved to sage.libs.eclib.constructor
    See http://trac.sagemath.org/19818 for details.
    sage: CremonaModularSymbols(225)
    Cremona Modular Symbols space of dimension 61 for Gamma_0(225) of weight 2 with sign 0
"""

from sage.misc.superseded import deprecation
deprecation(19818, "the module sage.libs.cremona.constructor has moved to sage.libs.eclib.constructor")

from sage.libs.eclib.constructor import *
