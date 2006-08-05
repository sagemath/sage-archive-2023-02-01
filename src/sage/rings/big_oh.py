"""
Big O for various types (power series, p-adics, etc.)
"""

import arith
import laurent_series_ring_element
import padic_field
import power_series_ring_element
import integer

def O(x):
    if isinstance(x, power_series_ring_element.PowerSeries):
        return x.parent()(0, x.degree())

    elif isinstance(x, laurent_series_ring_element.LaurentSeries):
        return laurent_series_ring_element.LaurentSeries(x.parent(), x.unit_part(),
                             x.valuation()).add_bigoh(x.degree())

    elif isinstance(x, (int,long,integer.Integer)):  # p-adic number
        if x <= 1:
            raise ArithmeticError, "x must be a prime power >= 2"
        F = arith.factor(x)
        if len(F) != 1:
            raise ArithmeticError, "x must be prime power"
        p, r = F[0]
        return padic_field.Qp(p)(0,r)

    raise ArithmeticError, "O(x) not defined"

