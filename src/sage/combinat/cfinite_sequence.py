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

from sage.rings.all import Integers, Integer, PolynomialRing, PowerSeriesRing, Rationals, Rational

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.matrix.constructor import Matrix
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.rings.fraction_field_element import is_FractionFieldElement
#from sage.rings.polynomial.polynomial_element import Polynomial
#from sage.rings.fraction_field_element import FractionFieldElement
from sage.rings.power_series_ring_element import PowerSeries
"""
C-Finite Sequences

This class provides C-finite sequence objects, also called homogenous linear
recurrences with constant coefficients.
CFiniteSequences contain their ordinary generating function (which
is always a polynomial fraction) as main defining object,
but also equivalently the recurrence degree, the coefficients, and the starting values.
They can be created from both of these representations, or from a power series by guessing.


AUTHORS:

-Ralf Stephan (2014): initial version

REFERENCES
    .. [GK82] Greene, Daniel H.; Knuth, Donald E. (1982), "2.1.1 Constant coefficients – A) Homogeneous equations", Mathematics for the Analysis of Algorithms (2nd ed.), Birkhäuser, p. 17.
    .. [SZ94] Bruno Salvy and Paul Zimmermann. — Gfun: a Maple package for the manipulation of generating and holonomic functions in one variable. — Acm transactions on mathematical software, 20.2:163-177, 1994.
"""

class CFiniteSequence(SageObject):

    def __init__(self, ogf):

        """
    Create an (integer) homogenous linear recurrence with constant coefficients,
    given its ordinary generating function.

    INPUT:

    - ``ogf`` -- the ordinary generating function, a fraction of polynomials


    OUTPUT:

    - A CFiniteSequence object

    EXAMPLES::

        sage: R.<x> = ZZ[]
        sage: CFiniteSequence((2-x)/(1-x-x^2))     # the Lucas sequence
        Homogenous linear recurrence with constant coefficients of degree 2: a(n+2) = 1*a(n) + 1*a(n+1), starting a(0) = 2, a(1) = 1
        sage: CFiniteSequence(x/(1-x)^3)           # triangular numbers
        Homogenous linear recurrence with constant coefficients of degree 3: a(n+3) = 3*a(n) - 3*a(n+1) + 1*a(n+2), starting a(0) = 0, a(1) = 1, a(2) = 3
        """

        if isinstance(ogf, Integer):
            self._ogf = ogf
            self._deg = 0
            self._c = []
            self._a = [ogf]
        if is_Polynomial(ogf):
            self._ogf = ogf
            self._deg = 0
            self._c = []
            self._a = ogf.list()
        elif is_FractionFieldElement(ogf):
            num = ogf.numerator()
            den = ogf.denominator()
            f = den.constant_coefficient()
            num = num/f
            den = den/f
            if den == 1:
                self._deg = 0
            else:
                self._deg = den.degree()
            R = PowerSeriesRing(Rationals(), 'x')
            self._a = (R(num, prec=self._deg) / R(den, prec=self._deg)).padded_list()
            self._c = [-den.list()[i] for i in range(1,self._deg+1)]
            self._ogf = num/den
        else: raise ValueError("Cannot convert to CFiniteSequence.")

    @classmethod
    def from_attributes(cls, deg, coefficients, values):
        """
    Create an (integer) homogenous linear recurrence with constant coefficients,
    given its degree, its coefficients, and its starting values.

    INPUT:

    - ``deg`` -- degree, a nonnegative integer; or an object to convert

    - ``coefficients`` -- a list of integers with size deg

    - ``values`` -- start values, a list of integers with size deg

    OUTPUT:

    - A CFiniteSequence object

    EXAMPLES::

        sage: r=CFiniteSequence.from_attributes(2,[1,1],[0,1])   # Fibonacci numbers
        sage: r
        Homogenous linear recurrence with constant coefficients of degree 2: a(n+2) = 1*a(n) + 1*a(n+1), starting a(0) = 0, a(1) = 1
        sage: CFiniteSequence.from_attributes(2,[-1,2],[0,1])    # natural numbers
        Homogenous linear recurrence with constant coefficients of degree 2: a(n+2) = -1*a(n) + 2*a(n+1), starting a(0) = 0, a(1) = 1
        """

        if deg < 0:
            raise ValueError("Degree must be nonnegative.")
        if not isinstance(coefficients, list):
            raise ValueError("Wrong type for recurrence coefficient list.")
        if not isinstance(values, list):
            raise ValueError("Wrong type for recurrence start value list.")
        if len(values) > deg:
            raise ValueError("More than deg start values not implemented.")

        R = PolynomialRing(Rationals(), 'x')
        x = R.gen()
        den = -1 + sum([x**(n+1)*coefficients[n] for n in range(deg)])
        num = -values[0] + sum([x**n * (-values[n] + sum([values[k]*coefficients[n-1-k] for k in range(n)])) for n in range(1,deg)])
        return cls(num/den)

    def __repr__(self):
        """
        Give string representation of the class.
        """
        
        if self._deg == 0:
            cstr = 'a(n) = 0'
        else:
            cstr = 'a(n+' + str(self._deg) + ') = ' + str(self._c[0]) + '*a(n)'
            for i in range(1, self._deg):
                if self._c[i] < 0:
                    cstr = cstr + ' - ' + str(-(self._c[i])) + '*a(n+' + str(i) + ')'
                elif self._c[i] > 0:
                    cstr = cstr + ' + ' + str(self._c[i]) + '*a(n+' + str(i) + ')'
                else:
                    cstr = cstr + ' + ' + str(self._c[i]) + '*a(n+' + str(i) + ')'
        astr = ', starting a(0) = ' + str(self._a[0])
        for i in range(1, self._deg):
            astr = astr + ', a(' + str(i) + ')' + ' = ' + str(self._a[i])
        return 'Homogenous linear recurrence with constant coefficients of degree ' + str(self._deg) + ': ' + cstr + astr

    def __eq__(self, other):
        """
    Compare two CFiniteSequences. We do not attempt to compare values, which can get
    complicated. Instead we compare generating functions lazily.

    EXAMPLES::

        sage: r = CFiniteSequence.from_attributes(2,[1,1],[2,1])
        sage: s = CFiniteSequence.from_attributes(1,[-1],[1])
        sage: r == s
        False
        sage: R.<x> = ZZ[]
        sage: r = CFiniteSequence.from_attributes(1,[-1],[1])
        sage: s = CFiniteSequence(1/(1+x))
        sage: r == s
        True

        """

        return self.ogf() == other.ogf()

    def __call__(self, n):
        """
        Give the nth term of the sequence generated by the CFiniteSequence. 

        INPUT:

        - ``n`` -- a nonnegative integer

        EXAMPLES::

            sage: r = CFiniteSequence.from_attributes(2,[3,3],[2,1])
            sage: r(2)
            9
            sage: r(101)
            16158686318788579168659644539538474790082623100896663971001
            sage: R.<x> = ZZ[]
            sage: r = CFiniteSequence(1/(1-x))
            sage: r(5)
            1
            sage: r = CFiniteSequence(x)
            sage: r(0)
            0
            sage: r(1)
            1
            sage: r = CFiniteSequence(0)
            sage: r(66)
            0

        """
        if self._deg == 0:
            if len(self._a) > n: return self._a[n]
            else:                return 0
        A = Matrix(Rationals(), 1, self._deg, self._c)
        B = Matrix.identity(Rationals(), self._deg-1)
        C = Matrix(Rationals(), self._deg-1, 1, 0)
        V = Matrix(Rationals(), self._deg, 1, self._a[::-1])
        M = Matrix.block([[A],[B,C]],subdivide=False)
        return list(M**(n-1)*V)[0][0]

    def ogf(self):
        """
        Return the ordinary generating function associated with the CFiniteSequence.
        This is always a polynomial fraction.

        EXAMPLES::
        
            sage: R = CFiniteSequence.from_attributes(1,[2],[1])   # powers of 2
            sage: R.ogf()
            1/(-2*x + 1)
        """

        return self._ogf

    def series(self, n):
        """
        Return a power series generated by the CFiniteSequence with precision n.
        
        INPUT:

        - ``n`` -- a nonnegative integer

        EXAMPLES::

            sage: r = CFiniteSequence.from_attributes(2,[2,-1],[0,1])
            sage: r.series(5)
            x + 2*x^2 + 3*x^3 + 4*x^4 + O(x^5)
        """

        R = PowerSeriesRing(Rationals(), 'x')
        return R(self._ogf.numerator(), prec=n) / R(self._ogf.denominator(), prec=n)

    def seq(self, n):
        """
        Return a list containing the sequence generated by the CFiniteSequence up to
        the term with index n-1. The first term has index 0.
        
        INPUT:

        - ``n`` -- a nonnegative integer

        EXAMPLES::

            sage: r = CFiniteSequence.from_attributes(2,[2,-1],[0,1])
            sage: r.seq(5)
            [0, 1, 2, 3, 4]
        """

        R = PowerSeriesRing(Rationals(), 'x')
        return (R(self._ogf.numerator(), prec=n) / R(self._ogf.denominator(), prec=n)).padded_list(n)

    @staticmethod
    def guess(list):
        """
        Return the minimal CFiniteSequence that generates the sequence. Assume the first
        value has index 0.

        INPUT:

        - ``list`` -- list of integers

        EXAMPLES::

            sage: CFiniteSequence.guess([1,2,4,8,16,32])   # not implemented
            [1, [2], [1]]
        """
        
        return CFiniteSequence.random_object()

    @staticmethod
    def random_object():
        """
        Return random CFiniteSequence.
        
        EXAMPLES::

            sage: CFiniteSequence.random_object()   # not implemented
            [1, [2], [1]]
        """
        
        return CFiniteSequence(0)

"""
.. TODO::

        sage: CFiniteSequence(x/(1-x))         # specific output not implemented
        Constant sequence a(n) = 1, starting a(0) = 0
        sage: CFiniteSequence(x^3/(1-x-x^2))       # gf>1 output not implemented
        Homogenous linear recurrence with constant coefficients of degree 2: a(n+2) = 1*a(n) + 1*a(n+1), starting a(0) = 0, a(1) = 0, a(2) = 0, a(3) = 1
        sage: CFiniteSequence.from_attributes(0,[],[1])          # specific output not implemented
        Constant sequence a(n) = 0, starting a(0) = 1
        sage: r=CFiniteSequence.from_attributes(2,[1,1],[0,0,0,1])   # gf>1 from attr not implemented
        sage: r.ogf()      # not implemented
        x^3/(-x^2 - x + 1)
        sage: r.egf()      # not implemented
        exp(2*x)
        sage: latex(r)        # not implemented
        \big\{a_{n\ge0}\big|a_{n+2}=\sum_{i=0}^{1}c_ia_{n+i}, c=\{1,1\}, a_{n<2}=\{0,0,0,1\}\big\}
        sage: r = CFiniteSequence.from_attributes(1,[-1],[1])       # not implemented
        sage: s = CFiniteSequence.from_attributes(1,[-1],[1,-1])    # not implemented
        sage: r == s                                                # not implemented
        True
        """
