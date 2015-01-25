"""
Symbolic Series

Symbolic series are special kinds of symbolic expressions that are
constructed via the
:meth:`Expression.series <sage.symbolic.expression.Expression.series>`
method.
They usually have an ``Order()`` term unless the series representation
is exact, see
:meth:`~sage.symbolic.series.SymbolicSeries.is_terminating_series`.

For series over general rings see
:class:`power series <sage.rings.power_series_poly.PowerSeries_poly>`
and
:class:`Laurent series<sage.rings.laurent_series_ring_element.LaurentSeries>`.

EXAMPLES:

We expand a polynomial in `x` about 0, about `1`, and also truncate
it back to a polynomial::

    sage: var('x,y')
    (x, y)
    sage: f = (x^3 - sin(y)*x^2 - 5*x + 3); f
    x^3 - x^2*sin(y) - 5*x + 3
    sage: g = f.series(x, 4); g
    3 + (-5)*x + (-sin(y))*x^2 + 1*x^3
    sage: g.truncate()
    x^3 - x^2*sin(y) - 5*x + 3
    sage: g = f.series(x==1, 4); g
    (-sin(y) - 1) + (-2*sin(y) - 2)*(x - 1) + (-sin(y) + 3)*(x - 1)^2 + 1*(x - 1)^3
    sage: h = g.truncate(); h
    (x - 1)^3 - (x - 1)^2*(sin(y) - 3) - 2*(x - 1)*(sin(y) + 1) - sin(y) - 1
    sage: h.expand()
    x^3 - x^2*sin(y) - 5*x + 3

We compute another series expansion of an analytic function::

    sage: f = sin(x)/x^2
    sage: f.series(x,7)
    1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
    sage: f.series(x==1,3)
    (sin(1)) + (cos(1) - 2*sin(1))*(x - 1) + (-2*cos(1) + 5/2*sin(1))*(x - 1)^2 + Order((x - 1)^3)
    sage: f.series(x==1,3).truncate().expand()
    -2*x^2*cos(1) + 5/2*x^2*sin(1) + 5*x*cos(1) - 7*x*sin(1) - 3*cos(1) + 11/2*sin(1)

Following the GiNaC tutorial, we use John Machin's amazing
formula `\pi = 16 \mathrm{tan}^{-1}(1/5) - 4 \mathrm{tan}^{-1}(1/239)`
to compute digits of `\pi`. We expand the arc tangent around 0 and insert
the fractions 1/5 and 1/239.

::

    sage: x = var('x')
    sage: f = atan(x).series(x, 10); f
    1*x + (-1/3)*x^3 + 1/5*x^5 + (-1/7)*x^7 + 1/9*x^9 + Order(x^10)
    sage: float(16*f.subs(x==1/5) - 4*f.subs(x==1/239))
    3.1415926824043994
"""
########################################################################
#       Copyright (C) 2015 Ralf Stephan <ralf@ark.in-berlin.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
########################################################################

include "sage/ext/interrupt.pxi"
include "sage/ext/stdsage.pxi"
include "sage/ext/cdefs.pxi"
include "sage/ext/python.pxi"

from ginac cimport *
from sage.symbolic.expression cimport Expression, new_Expression_from_GEx

cdef class SymbolicSeries(Expression):
    def __init__(self, SR):
        Expression.__init__(self, SR, 0)
        self._parent = SR

    def is_series(self):
        """
        Return True.

        EXAMPLES::

            sage: exp(x).series(x,10).is_series()
            True
        """
        return True

    def is_terminating_series(self):
        """
        Return True if ``self`` is without order term.

        A series is terminating if it can be represented exactly,
        without requiring an order term.

        OUTPUT:

        Boolean. Whether ``self`` has no order term.

        EXAMPLES::

            sage: (x^5+x^2+1).series(x,10)
            1 + 1*x^2 + 1*x^5
            sage: (x^5+x^2+1).series(x,10).is_series()
            True
            sage: (x^5+x^2+1).series(x,10).is_terminating_series()
            True
            sage: SR(5).is_terminating_series()
            False
            sage: exp(x).series(x,10).is_terminating_series()
            False
        """
        return g_is_a_terminating_series((<Expression>self)._gobj)

    def truncate(self):
        """
        Given a power series or expression, return the corresponding
        expression without the big oh.

        INPUT:

        - ``self`` -- a series as output by the :meth:`series` command.

        OUTPUT:

        A symbolic expression.

        EXAMPLES::

            sage: f = sin(x)/x^2
            sage: f.truncate()
            sin(x)/x^2
            sage: f.series(x,7)
            1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
            sage: f.series(x,7).truncate()
            -1/5040*x^5 + 1/120*x^3 - 1/6*x + 1/x
            sage: f.series(x==1,3).truncate().expand()
            -2*x^2*cos(1) + 5/2*x^2*sin(1) + 5*x*cos(1) - 7*x*sin(1) - 3*cos(1) + 11/2*sin(1)
        """
        return new_Expression_from_GEx(self._parent, series_to_poly(self._gobj))

    def default_variable(self):
        """
        Return the expansion variable of this symbolic series.

        EXAMPLES::

            sage: s=(1/(1-x)).series(x,3); s
            1 + 1*x + 1*x^2 + Order(x^3)
            sage: s.default_variable()
            x
        """
        cdef GEx x = g_series_var(self._gobj)
        cdef Expression ex = new_Expression_from_GEx(self._parent, x)
        return ex

    def coefficients(self, x=None, sparse=True):
        r"""
        Return the coefficients of this symbolic series as a polynomial in x.

        INPUT:

        -  ``x`` -- optional variable.

        OUTPUT:

        Depending on the value of ``sparse``,

        - A list of pairs ``(expr, n)``, where ``expr`` is a symbolic
          expression and ``n`` is a power (``sparse=True``, default)

        - A list of expressions where the ``n``-th element is the coefficient of
          ``x^n`` when self is seen as polynomial in ``x`` (``sparse=False``).

        EXAMPLES::

            sage: s=(1/(1-x)).series(x,6); s
            1 + 1*x + 1*x^2 + 1*x^3 + 1*x^4 + 1*x^5 + Order(x^6)
            sage: s.coefficients()
            [[1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [1, 5]]
            sage: s.coefficients(x, sparse=False)
            [1, 1, 1, 1, 1, 1]
            sage: x,y = var("x,y")
            sage: s=(1/(1-y*x-x)).series(x,3); s
            1 + (y + 1)*x + ((y + 1)^2)*x^2 + Order(x^3)
            sage: s.coefficients(x, sparse=False)
            [1, y + 1, (y + 1)^2]

        """
        if x is None:
            x = self.default_variable()
        l = [[self.coefficient(x, d), d] for d in xrange(self.degree(x))]
        if sparse is True:
            return l
        else:
            from sage.rings.integer_ring import ZZ
            if any(not c[1] in ZZ for c in l):
                raise ValueError("Cannot return dense coefficient list with noninteger exponents.")
            val = l[0][1]
            if val < 0:
                raise ValueError("Cannot return dense coefficient list with negative valuation.")
            deg = l[-1][1]
            ret = [ZZ(0)] * int(deg+1)
            for c in l:
                ret[c[1]] = c[0]
            return ret

