"""
Big O for various types (power series, p-adics, etc.)
"""

import sage.rings.arith as arith
import laurent_series_ring_element
import sage.rings.padics.factory as padics_factory
import sage.rings.padics.padic_generic_element as padic_generic_element
import power_series_ring_element
import integer
import rational
from sage.rings.integer_mod import Mod

def O(x):
    if isinstance(x, power_series_ring_element.PowerSeries):
        return x.parent()(0, x.degree())

    elif isinstance(x, laurent_series_ring_element.LaurentSeries):
        return laurent_series_ring_element.LaurentSeries(x.parent(), x.valuation_zero_part(),
                             x.valuation()).add_bigoh(x.degree())

    elif isinstance(x, (int,long,integer.Integer,rational.Rational)):  # p-adic number
        if x <= 0:
            raise ArithmeticError, "x must be a prime power >= 2"
        F = arith.factor(x)
        if len(F) != 1:
            raise ArithmeticError, "x must be prime power"
        p, r = F[0]
        return padics_factory.Zp(p, prec = r, type = 'capped-rel')(0, absprec = r)

    elif isinstance(x, padic_generic_element.pAdicGenericElement):
         return x.parent()(0, absprec = x.valuation())
    raise ArithmeticError, "O(x) not defined"

