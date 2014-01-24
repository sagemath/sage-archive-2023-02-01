# -*- coding: utf-8 -*-
#*****************************************************************************
#       Copyright (C) 2014 Ralf Stephan <gtrwst9@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL) v2.0
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/gpl-2.0.html
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.power_series_ring_element import PowerSeries
"""
Linear Recurrences


This class provides linear recurrence objects, in particular homogenous linear
recurrences with constant coefficients (abbreviated "linrec" in the following).
Linrecs contain the recurrence degree, the coefficients, and the starting values.
They can be created from these values, from their generating function (which
is always a polynomial fraction), or from a power series by guessing.


AUTHORS:

-Ralf Stephan (2014): initial version

REFERENCES
    .. [GK82] Greene, Daniel H.; Knuth, Donald E. (1982), "2.1.1 Constant coefficients – A) Homogeneous equations", Mathematics for the Analysis of Algorithms (2nd ed.), Birkhäuser, p. 17.
    .. [SZ94] Bruno Salvy and Paul Zimmermann. — Gfun: a Maple package for the manipulation of generating and holonomic functions in one variable. — Acm transactions on mathematical software, 20.2:163-177, 1994.
"""

class LinearRecurrence(SageObject):

    def __init__(self, deg, c, a):

        """
    Create an (integer) homogenous linear recurrence with constant coefficients,
    given its degree, its coefficients, and its starting values.

    INPUT:

    - ``deg`` -- degree, a nonnegative integer

    - ``c`` -- coefficients, a list of integers with size deg

    - ``a`` -- start values, a list of integers with size deg


    OUTPUT:

    - A LinearRecurrence object

    EXAMPLES::

        sage: r = LinearRecurrence(2,[1,1],[2,1])     # the Lucas sequence
        sage: r
        Homogenous linear recurrence with constant coefficients of degree 2: a(n+2) = a(n) + a(n+1), starting a_0=2, a_1=1.

        """
        if not isinstance(deg, Integer):
            raise ValueError("Wrong type for recurrence degree.")
        if deg<0:
            raise ValueError("Degree must be nonnegative.")
        self._deg = deg
        if not isinstance(c, List):
            raise ValueError("Wrong type for recurrence coefficient list.")
        self._c = c
        if not isinstance(a, List):
            raise ValueError("Wrong type for recurrence start value list.")
        self._a = a

    def __init__(self, obj):
        """
    Convert another object to a linrec. If obj is a polynomial or polynomial fraction,
    take it as an ordinary generating function.

    - ``obj`` -- object to convert

    OUTPUT:

    - A LinearRecurrence object

    EXAMPLES::

        sage: r = LinearRecurrence(x/(1-x)^2)         # the natural numbers
        sage: s = LinearRecurrence(2,[2,-1],[0,1])
        sage: r == s
        True
    """

        self._ogf = obj
        if isinstance(obj, Polynomial):
            self._deg = 0
            self._c = []
            self._a = obj.list()
        elif not isinstance(obj, FractionFieldElement):
            raise ValueError("Cannot convert to linrec.")
        
        self._deg = obj.denominator().degree()
        s = PowerSeries(obj, prec=1+deg).list()
        self._a = s
        self._c = s

    def __repr__(self):
        """
        Give string representation of the class.

        EXAMPLES::

        sage: r = LinearRecurrence(2,[1,1],[2,1])     # the Lucas sequence
        sage: r
        Homogenous linear recurrence with constant coefficients of degree 2: a(n+2) = a(n) + a(n+1), starting a_0=2, a_1=1.

        """
        return 'Homogenous linear recurrence with constant coefficients of degree ' + str(self.deg)

    def _latex_(self):
       r"""
       Return the LaTeX representation of the linrec.

       EXAMPLES::

           sage: powers_of_two = LinearRecurrence(1,[2],[1])
           sage: latex(powers_of_two)
           '\\frac{1}{2}'
       """
       return '' #'\\big\\{a_{n\ge0}|a_{n+%s}=\sum_{i=0}^{%s}c_ia_{n+i}, c=\\{%s\\}, a_{n<$s}=\\{%s\\}\\big\\}'% (latex(self.deg), latex(self.deg-1), latex(self.c}, latex(self.deg), latex(self.a)

    def __eq__(self, other):
        """
    Compare two linrecs. We do not attempt to compare values, which can get
    complicated. Instead we compare generating functions lazily.

    EXAMPLES::

        sage: r = LinearRecurrence(2,[1,1],[2,1])
        sage: s = LinearRecurrence(1,[-1],[1])
        sage: r == s
        False
        sage: r = LinearRecurrence(1,[-1],[1])
        sage: s = LinearRecurrence(1,[-1],[1])
        sage: r == s
        True
        sage: r = LinearRecurrence(1,[-1],[1])
        sage: s = LinearRecurrence(1,[-1],[1,-1])
        sage: r == s
        True

        """

        return self.ogf() == other.ogf()

    def __call__(self, n):
        """
        Give the nth term of the sequence generated by the linrec. 

        INPUT:

        - ``n`` -- a nonnegative integer

        EXAMPLES::

            sage: r = LinearRecurrence(2,[3,3],[2,1])
            sage: r(2)
            9
            sage: r(101)
            16158686318788579168659644539538474790082623100896663971001

        """
        F = matrix(R, [[0,1],[self.c,self.b]])            # F*[u_{n}, u_{n+1}]^T = [u_{n+1}, u_{n+2}]^T (T indicates transpose).
        v = matrix(R, [[self.u0],[self.u1]])
        return list(F**n*v)[0][0]

    def ogf(self):
        """
        Return the ordinary generating function associated with the linrec.
        This is always a polynomial fraction.

        EXAMPLES::
        
            sage: R = LinearRecurrence(1,[2],[1])   # powers of 2
            sage: R.ogf()
            1/(1-2*x)
        """

        return self._ogf

    def series(self, n):
        """
        Return a power series generated by the linrec with precision n.
        
        INPUT:

        - ``n`` -- a nonnegative integer

        EXAMPLES::

            sage: r = LinearRecurrence(2,[2,-1],[0,1])
            sage: r.series(5)
            x + 2*x^2 + 3*x^3 + 4*x^4 + O(x^5)
        """

        R = PowerSeriesRing(ZZ, 'x')
        return R(self._ogf, prec=n)

    def seq(self, n):
        """
        Return a list containing the sequence generated by the linrec up to
        the term with index n-1. The first term has index 0.
        
        INPUT:

        - ``n`` -- a nonnegative integer

        EXAMPLES::

            sage: r = LinearRecurrence(2,[2,-1],[0,1])
            sage: r.seq(5)
            [0, 1, 2, 3, 4]
        """

        R = PowerSeriesRing(ZZ, 'x')
        return R(self._ogf, prec=n).padded_list(n)

    def egf(self):
        """
        Return the exponential generating function associated with this linrec as
        a symbolic expression (using the exponential function).

        EXAMPLES::
        
            sage: r = LinearRecurrence(1/(1-2*x))        # not implemented
            sage: R.egf()
            exp(2*x)
        """

        return ''

    @staticmethod
    def guess(list):
        """
        Return the minimal linrec that generates the sequence. Assume the first
        value has index 0.

        INPUT:

        - ``list`` -- list of integers

        EXAMPLES::

            sage: LinearRecurrence.guess([1,2,4,8,16,32])   # not implemented
            [1, [2], [1]]
        """
        
        return LinearRecurrence.random_object()

    @staticmethod
    def random_object():
        """
        Return random linrec.
        
        EXAMPLES::

            sage: LinearRecurrence.guess([1,2,4,8,16,32])   # not implemented
            [1, [2], [1]]
        """
        
        return LinearRecurrence(0)

