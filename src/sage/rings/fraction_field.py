"""
Fraction Field of Integral Domains

AUTHOR: William Stein (with input from David Joyner, David Kohel, and
        Joe Wetherell)

TESTS:
    sage: F = FractionField(IntegerRing())
    sage: F == loads(dumps(F))
    True

    sage: F = FractionField(PolynomialRing(RationalField(),'x'))
    sage: F == loads(dumps(F))
    True

    sage: F == FractionField(PolynomialRing(IntegerRing(),'x'))
    sage: F == loads(dumps(F))
    True

    sage: FractionField(MPolynomialRing(RationalField(),2,'x'))
    sage: F == loads(dumps(F))
    True


"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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

import ring
import integral_domain
import field
import fraction_field_element
import sage.misc.latex as latex

from sage.structure.parent_base import ParentWithBase

def FractionField(R, names=None):
    """
    Create the fraction field of the integral domain R.

    INPUT:
        R -- an integral domain
        names -- ignored

    EXAMPLES:
    We create some example fraction fields.
        sage: FractionField(IntegerRing())
        Rational Field
        sage: FractionField(PolynomialRing(RationalField(),'x'))
        Fraction Field of Univariate Polynomial Ring in x over Rational Field
        sage: FractionField(PolynomialRing(IntegerRing(),'x'))
        Fraction Field of Univariate Polynomial Ring in x over Integer Ring
        sage: FractionField(MPolynomialRing(RationalField(),2,'x'))
        Fraction Field of Polynomial Ring in x0, x1 over Rational Field

    Dividing elements often implicitly creates elements of the fraction field.
        sage: x = PolynomialRing(RationalField(), 'x').gen()
        sage: f = x/(x+1)
        sage: g = x**3/(x+1)
        sage: f/g
        1/x^2
        sage: g/f
        x^2

    The input must be an integral domain.
        sage: Frac(Integers(4))
        Traceback (most recent call last):
        ...
        TypeError: R must be an integral domain.
    """
    if not ring.is_Ring(R):
        raise TypeError, "R must be a ring"
    if not R.is_integral_domain():
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
            sage: Frac(QQ['x'])
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
        """
        self.__R = R
        ParentWithBase.__init__(self, R)

    def is_field(self):
        """
        Returns True, since the fraction field is a field.

        EXAMPLES:
            sage: Frac(ZZ).is_field()
            True
        """
        return True

    def base_ring(self):
        """
        Return the base ring of self; this is the base ring of the ring which
        this fraction field is the fraction field of.

        EXAMPLES:
            sage: R = Frac(ZZ['t'])
            sage: R.base_ring()
            Integer Ring
        """
        return self.__R.base_ring()

    def characteristic(self):
        """
        Return the characteristic of this fraction field.

        EXAMPLES:
            sage: R = Frac(ZZ['t'])
            sage: R.base_ring()
            Integer Ring
            sage: R = Frac(ZZ['t']); R.characteristic()
            0
            sage: R = Frac(GF(5)['w']); R.characteristic()
            5
        """
        return self.ring().characteristic()

    def __repr__(self):
        return "Fraction Field of %s"%self.ring()

    def _latex_(self):
        return "\\mbox{\\rm Frac}(%s)"%latex.latex(self.ring())

    def ring(self):
        """
        Return the ring that this is the fraction field of.

        EXAMPLES:
            sage: R = Frac(QQ['x,y'])
            sage: R
            Fraction Field of Polynomial Ring in x, y over Rational Field
            sage: R.ring()
            Polynomial Ring in x, y over Rational Field
        """
        return self.__R

    def __call__(self, x, y=1, coerce=True):
        if isinstance(x, fraction_field_element.FractionFieldElement):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                return fraction_field_element.FractionFieldElement(self, x.numerator(), x.denominator())
            else:
                R = self.ring()
                return fraction_field_element.FractionFieldElement(self, R(x.numerator()), R(x.denominator()))
        if coerce:
            R = self.ring()
            x = R(x); y = R(y)
        return fraction_field_element.FractionFieldElement(self, x, y, coerce=False)

    def _coerce_impl(self, x):
        """
        Return the canonical coercion of x into this fraction field, or raise a TypeError.

        The rings that canonically coerce to the fraction field are
           * the fraction field itself
           * any fraction field that of the form Frac(S) where S canonically
             coerces to this ring.
           * any ring that canonically coerces to the ring R such that this
             fraction field is Frac(R)

        """
        try:
            P = x.parent()
            if is_FractionField(P):
                R = P.ring()
                S = self.ring()
                if S.has_coerce_map_from(R):
                    return self(x)
            else:
                S = self.ring()
                if S.has_coerce_map_from(P):
                    return self(x)
                raise TypeError
        except AttributeError:
            pass
        return self._coerce_try(x, [self.ring()])

    def __cmp__(self, other):
        if not isinstance(other, FractionField_generic):
            return cmp(type(self), type(other))
        return cmp(self.ring(), other.ring())

    def ngens(self):
        """
        This is the same as for the parent object.

        EXAMPLES:
            sage: R = Frac(PolynomialRing(QQ,'z',10)); R
            Fraction Field of Polynomial Ring in z0, z1, z2, z3, z4, z5, z6, z7, z8, z9 over Rational Field
            sage: R.ngens()
            10
        """
        return self.ring().ngens()

    def gen(self, i=0):
        """
        Return the ith generator of self.

        EXAMPLES:
            sage: R = Frac(PolynomialRing(QQ,'z',10)); R
            Fraction Field of Polynomial Ring in z0, z1, z2, z3, z4, z5, z6, z7, z8, z9 over Rational Field
            sage: R.0
            z0
            sage: R.gen(3)
            z3
            sage: R.3
            z3
        """
        x = self.ring().gen(i)
        one = self.ring()(1)
        r = fraction_field_element.FractionFieldElement(self, x, one,
                                                           coerce=False, reduce=False)
        return r
