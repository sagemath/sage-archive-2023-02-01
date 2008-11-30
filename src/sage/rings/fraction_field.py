"""
Fraction Field of Integral Domains

AUTHOR: William Stein (with input from David Joyner, David Kohel, and
        Joe Wetherell)

EXAMPLES:
Quotienting is a constructor for an element of the fraction field:
    sage: R.<x> = QQ[]
    sage: (x^2-1)/(x+1)
    x - 1
    sage: parent((x^2-1)/(x+1))
    Fraction Field of Univariate Polynomial Ring in x over Rational Field


The GCD is not taken (since it doesn't converge sometimes) in the inexact case.
    sage: Z.<z> = CC[]
    sage: I = CC.gen()
    sage: (1+I+z)/(z+0.1*I)
    (1.00000000000000*z + 1.00000000000000 + 1.00000000000000*I)/(1.00000000000000*z + 0.100000000000000*I)
    sage: (1+I*z)/(z+1.1)
    (1.00000000000000*I*z + 1.00000000000000)/(1.00000000000000*z + 1.10000000000000)


TESTS:
    sage: F = FractionField(IntegerRing())
    sage: F == loads(dumps(F))
    True

    sage: F = FractionField(PolynomialRing(RationalField(),'x'))
    sage: F == loads(dumps(F))
    True

    sage: F = FractionField(PolynomialRing(IntegerRing(),'x'))
    sage: F == loads(dumps(F))
    True

    sage: F = FractionField(PolynomialRing(RationalField(),2,'x'))
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
        sage: FractionField(PolynomialRing(RationalField(),2,'x'))
        Fraction Field of Multivariate Polynomial Ring in x0, x1 over Rational Field

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
    """
    Tests whether or not x inherits from FractionField_generic.

    EXAMPLES:
        sage: from sage.rings.fraction_field import is_FractionField
        sage: is_FractionField(Frac(ZZ['x']))
        True
        sage: is_FractionField(QQ)
        False
    """
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
            sage: Frac(QQ['x,y']).variable_names()
            ('x', 'y')
        """
        self.__R = R
        ParentWithBase.__init__(self, R)
        self._assign_names(R._names)

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
        """
        EXAMPLES:
            sage: str(Frac(ZZ['x']))
            'Fraction Field of Univariate Polynomial Ring in x over Integer Ring'
        """
        return "Fraction Field of %s"%self.ring()

    def _latex_(self):
        """
        EXAMPLES:
            sage: latex(Frac(GF(7)['x,y,z']))
            \mbox{\rm Frac}(\mathbf{F}_{7}[x, y, z])
        """
        return "\\mbox{\\rm Frac}(%s)"%latex.latex(self.ring())

    def _magma_init_(self, magma):
        """
        Return a string representation of self in the given magma instance.

        EXAMPLES:
            sage: QQ['x'].fraction_field()._magma_init_(magma)            # optional - magma
            'SageCreateWithNames(FieldOfFractions(SageCreateWithNames(PolynomialRing(RationalField()),["x"])),["x"])'
            sage: GF(9,'a')['x,y,z'].fraction_field()._magma_init_(magma) # optional - magma
            'SageCreateWithNames(FieldOfFractions(SageCreateWithNames(PolynomialRing(_sage_ref...,3,"grevlex"),["x","y","z"])),["x","y","z"])'

        _magma_init_ gets called implicitly below.
            sage: magma(QQ['x,y'].fraction_field())                  # optional - magma
            Multivariate rational function field of rank 2 over Rational Field
            Variables: x, y
            sage: magma(ZZ['x'].fraction_field())                    # optional - magma
            Univariate rational function field over Integer Ring
            Variables: x

        Verify that conversion is being properly cached:
            sage: k = Frac(QQ['x,z'])                                # optional - magma
            sage: magma(k) is magma(k)                               # optional - magma
            True
        """
        s = 'FieldOfFractions(%s)'%self.ring()._magma_init_(magma)
        return magma._with_names(s, self.variable_names())

    def ring(self):
        """
        Return the ring that this is the fraction field of.

        EXAMPLES:
            sage: R = Frac(QQ['x,y'])
            sage: R
            Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field
            sage: R.ring()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self.__R

    def __call__(self, x, y=None, check=True):
        """
        Try to create an element out of x. Optionally, a denominator can be provided as well.

        EXAMPLES:
            sage: R.<x> = ZZ[]
            sage: K = Frac(R)
            sage: K(3)
            3
            sage: K(3,2)
            3/2
            sage: K(x,x^2-1)
            x/(x^2 - 1)
            sage: K(2/3)
            2/3
            sage: K(3/x)
            3/x
            sage: K(x + 1/2)
            (2*x + 1)/2
        """
        if isinstance(x, fraction_field_element.FractionFieldElement):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                return fraction_field_element.FractionFieldElement(self, x.numerator(), x.denominator())
            else:
                R = self.ring()
                return fraction_field_element.FractionFieldElement(self, R(x.numerator()), R(x.denominator()))
        if check:
            R = self.ring()
            if y is None:
                try:
                    x, y = R(x), R(1)
                except:
                    if hasattr(x, 'numerator') and hasattr(x, 'denominator'):
                        x, y = R(x.numerator()), R(x.denominator())
                    else:
                        raise
            else:
                x = R(x); y = R(y)
        else:
            if y is None:
                y = R(1)
        return fraction_field_element.FractionFieldElement(self, x, y,
                                            coerce=False, reduce = self.is_exact())

    def is_exact(self):
        """
        EXAMPLES:
            sage: Frac(ZZ['x']).is_exact()
            True
            sage: Frac(CDF['x']).is_exact()
            False
        """
        try:
            return self.__is_exact
        except AttributeError:
            r = self.ring().is_exact()
            self.__is_exact = r
        return r

    def construction(self):
        """
        EXAMPLES:
            sage: Frac(ZZ['x']).construction()
            (FractionField, Univariate Polynomial Ring in x over Integer Ring)
            sage: K = Frac(GF(3)['t'])
            sage: f, R = K.construction()
            sage: f(R)
            Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 3
            sage: f(R) == K
            True
        """
        from sage.categories.pushout import FractionField
        return FractionField(), self.ring()

    def _coerce_impl(self, x):
        """
        Return the canonical coercion of x into this fraction field, or raise a TypeError.

        The rings that canonically coerce to the fraction field are
           * the fraction field itself
           * any fraction field that of the form Frac(S) where S canonically
             coerces to this ring.
           * any ring that canonically coerces to the ring R such that this
             fraction field is Frac(R)

        EXAMPLES:
            sage: K = Frac(QQ['x'])
            sage: K.coerce(4)
            4
            sage: K.coerce(ZZ['x'].0)
            x
            sage: K.coerce(Frac(ZZ['x']).0)
            x
            sage: K.coerce(RR['x'].0)
            Traceback (most recent call last):
            ...
            TypeError: no cannonical coercion from Univariate Polynomial Ring in x over Real Field with 53 bits of precision to Fraction Field of Univariate Polynomial Ring in x over Rational Field
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
        except AttributeError:
            pass
        return self._coerce_try(x, [self.ring()])

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: Frac(ZZ['x']) == Frac(ZZ['x'])
            True
            sage: Frac(ZZ['x']) == Frac(QQ['x'])
            False
            sage: Frac(ZZ['x']) == Frac(ZZ['y'])
            False
            sage: Frac(ZZ['x']) == QQ['x']
            False
        """
        if not isinstance(other, FractionField_generic):
            return cmp(type(self), type(other))
        return cmp(self.ring(), other.ring())

    def ngens(self):
        """
        This is the same as for the parent object.

        EXAMPLES:
            sage: R = Frac(PolynomialRing(QQ,'z',10)); R
            Fraction Field of Multivariate Polynomial Ring in z0, z1, z2, z3, z4, z5, z6, z7, z8, z9 over Rational Field
            sage: R.ngens()
            10
        """
        return self.ring().ngens()

    def gen(self, i=0):
        """
        Return the ith generator of self.

        EXAMPLES:
            sage: R = Frac(PolynomialRing(QQ,'z',10)); R
            Fraction Field of Multivariate Polynomial Ring in z0, z1, z2, z3, z4, z5, z6, z7, z8, z9 over Rational Field
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
