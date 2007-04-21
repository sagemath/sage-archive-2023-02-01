"""
Big O for various types (power series, p-adics, etc.)
"""

import arith
import laurent_series_ring_element_pyx as laurent_series_ring_element
import padics.qp
import padics.padic_generic
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
        return padics.qp.Qp(p, prec = r, type = 'capped-rel')(Mod(0, x))

    elif isinstance(x, padics.padic_generic.pAdicGeneric):
         return padics.qp.Qp(p, prec = x.parent().precision_cap(), type = 'capped-rel')(x.parent().prime_pow(x.valuation()))
    raise ArithmeticError, "O(x) not defined"

