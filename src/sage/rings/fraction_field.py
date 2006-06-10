"""
Fraction Field of Integral Domains

AUTHOR: William Stein (with input from David Joyner, David Kohel, and
Joe Wetherell)
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import integral_domain
import field
import fraction_field_element
import sage.misc.latex as latex

def FractionField(R):
    """
    Create the fraction field of the integral domain R.

    INPUT:
        R -- an integral domain

    EXAMPLES:
    We create some example fraction fields.
        sage: FractionField(IntegerRing())
        Rational Field
        sage: FractionField(PolynomialRing(RationalField()))
        Fraction Field of Univariate Polynomial Ring in x over Rational Field
        sage: FractionField(PolynomialRing(IntegerRing()))
        Fraction Field of Univariate Polynomial Ring in x over Integer Ring
        sage: FractionField(MPolynomialRing(RationalField(),2))
        Fraction Field of Polynomial Ring in x0, x1 over Rational Field

    Dividing elements often implicitly creates elements of the fraction field.
        sage: x = PolynomialRing(RationalField()).gen()
        sage: f = x/(x+1)
        sage: g = x**3/(x+1)
        sage: f/g
        1/x^2
        sage: g/f
        x^2

    """
    if not isinstance(R, integral_domain.IntegralDomain):
        raise TypeError, "R must be an integral domain."
    return R.fraction_field()

def is_FractionField(x):
    return isinstance(x, FractionField_generic)

class FractionField_generic(field.Field):
    """
    The fraction field of an integral domain.
    """
    def __init__(self, R):
        """
        Create the fraction field of the integral domain R.
        INPUT:
            R -- an integral domain
        EXAMPLES:
            sage: K, x = FractionField(PolynomialRing(QQ)).objgen()
            sage: K
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
        """
        self.__R = R

    def is_field(self):
        return True

    def base_ring(self):
        return self.__R.base_ring().fraction_field()

    def characteristic(self):
        return self.ring().characteristic()

    def __repr__(self):
        return "Fraction Field of %s"%self.ring()

    def _latex_(self):
        return "\\mbox{\\rm Frac}(%s)"%latex.latex(self.ring())

    def ring(self):
        return self.__R

    def __call__(self, x, y=1, coerce=True):
        if isinstance(x, fraction_field_element.FractionFieldElement):
            if x.parent() == self:
                return x
            else:
                R = self.ring()
                return fraction_field_element.FractionFieldElement(self, R(x.numerator()), R(x.denominator()))
        if coerce:
            R = self.ring()
            x = R(x); y = R(y)
        return fraction_field_element.FractionFieldElement(self, x, y, coerce=False)

    def _coerce_(self, x):
        if isinstance(x, fraction_field_element.FractionFieldElement) \
           and x.parent() == self:
            return x
        if x.parent() == self.ring():
            return self(x)
        raise TypeError, "no canonical coercion of x(=%s) into %s"%(x,self)

    def __cmp__(self, other):
        if not isinstance(other, FractionField_generic):
            return -1
        if self.ring() == other.ring():
            return 0
        return 1

    def ngens(self):
        return self.ring().ngens()

    def gen(self, i=0):
        return fraction_field_element.FractionFieldElement(self, self.ring().gen(i))

