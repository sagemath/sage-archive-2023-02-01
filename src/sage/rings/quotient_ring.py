"""
Quotient Rings

AUTHOR:
    -- William Stein
"""

#*****************************************************************************
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import quotient_ring_element
import sage.misc.latex as latex
import commutative_ring
import ideal
from sage.structure.all import Generators

def QuotientRing(R, I):
    if not isinstance(R, commutative_ring.CommutativeRing):
        raise TypeError, "R must be a commutative ring."
    if not isinstance(I, ideal.Ideal_generic) or I.ring() != R:
        I = R.ideal(I)
    try:
        if I.is_principal():
            return R.quotient_by_principal_ideal(I.gen())
    except (AttributeError, NotImplementedError):
        pass
    if isinstance(R, QuotientRing_generic):
        pi = R.cover()
        S = pi.domain()
        G = [pi.lift(x) for x in I.gens()]
        I_lift = S.ideal(G)
        J = R.defining_ideal()
        return QuotientRing_generic(S, I_lift + J)

    return QuotientRing_generic(R, I)

def is_QuotientRing(x):
    return isinstance(x, QuotientRing_generic)

class QuotientRing_generic(commutative_ring.CommutativeRing, Generators):
    """
    The quotient ring of $R$ by the ideal $I$.

    EXAMPLES:
        sage: R.<x> = PolynomialRing(Z)
        sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
        sage: S = R.quotient_ring(I); S
        Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 1, x^2 + 3*x + 4)

        sage: R, (x,y) = PolynomialRing(Q, 2, 'xy').objgens()
        sage: S, (a,b) = (R/(x^2 + y^2)).objgens('ab')
        sage: a^2 + b^2 == 0
        True
        sage: S(0) == a^2 + b^2
        True

    EXAMPLE: Quotient of quotient

    A quotient of a quotient is just the quotient of the
    original top ring by the sum of two ideals.
        sage: R, (x,y) = PolynomialRing(Q, 2, 'xy').objgens()
        sage: S, (a,b) = (R/(1 + y^2)).objgens('ab')
        sage: T, (c,d) = (S/(a, )).objgens('cd')
        sage: T
        Quotient of Polynomial Ring in x, y over Rational Field by the ideal (x, 1 + y^2)
        sage: T.gens()
        (0, d)
    """
    def __init__(self, R, I, names=None):
        """
        Create the quotient ring of R by the ideal I.

        INPUT:
            R -- a commutative ring
            I -- an ideal
        """
        self.__R = R
        self.__I = I
        self.assign_names(names)

    def _repr_(self):
        return "Quotient of %s by the ideal %s"%(self.cover_ring(), self.defining_ideal()._repr_short())

    def _latex_(self):
        return "%s/%s"%(latex.latex(self.cover_ring()), latex.latex(self.defining_ideal()))

    def cover(self):
        r"""
        The covering ring homomorphism $R \to R/I$, equipped with a section.

        EXAMPLES:
            sage: R = Z/(3*Z)
            sage: pi = R.cover()
            sage: pi
            Ring morphism:
              From: Integer Ring
              To:   Ring of integers modulo 3
              Defn: Natural quotient map
            sage: pi(5)
            2
            sage:
            sage: l = pi.lift()

        EXAMPLES:
            sage: R, (x,y)  = QQ['x,y'].objgens()
            sage: Q = R/(x^2,y^2)
            sage: pi = Q.cover()
            sage: pi(x^3+y)
            y
            sage: l = pi.lift(x+y^3)
            sage: l
            x
            sage: l = pi.lift(); l
            Set-theoretic ring morphism:
              From: Quotient of Polynomial Ring in x, y over Rational Field by the ideal (y^2, x^2)
              To:   Polynomial Ring in x, y over Rational Field
              Defn: Choice of lifting map
            sage: l(x+y^3)
            x
        """
        try:
            return self.__cover
        except AttributeError:
            import morphism
            pi = morphism.RingHomomorphism_cover(self.__R, self)
            #pi = self.__R.homomorphism(self.__R.gens(), self)       # needlessly slow (!)
            lift = self.lift()
            pi._set_lift(lift)
            self.__cover = pi
            return self.__cover

    def lift(self):
        """
        Return the lifting map to the cover.

        EXAMPLES:
            sage: R, (x,y) = PolynomialRing(Q, 2, 'xy').objgens()
            sage: S = R.quotient(x^2 + y^2, names=['xbar', 'ybar'])
            sage: pi = S.cover(); pi
            Ring morphism:
              From: Polynomial Ring in x, y over Rational Field
              To:   Quotient of Polynomial Ring in x, y over Rational Field by the ideal (y^2 + x^2)
              Defn: Natural quotient map
            sage: L = S.lift(); L
            Set-theoretic ring morphism:
              From: Quotient of Polynomial Ring in x, y over Rational Field by the ideal (y^2 + x^2)
              To:   Polynomial Ring in x, y over Rational Field
              Defn: Choice of lifting map
            sage: L(S.0)
            x
            sage: L(S.1)
            y

        Note that some reduction may be applied so that the lift of a reduction
        need not equal the original element.
            sage: z = pi(x^3 + 2*y^2); z
            2*ybar^2 - xbar*ybar^2
            sage: L(z)
            2*y^2 - x*y^2
            sage: L(z) == x^3 + 2*y^2
            False
        """
        try:
            return self.__lift
        except AttributeError:
            from morphism import RingMap_lift
            self.__lift = RingMap_lift(self, self.__R)
            return self.__lift

    def characteristic(self):
        raise NotImplementedError

    def defining_ideal(self):
        return self.__I

    def is_field(self):
        return self.defining_ideal().is_maximal()

    def is_integral_domain(self):
        return self.defining_ideal().is_prime()

    def cover_ring(self):
        return self.__R

    def __call__(self, x, coerce=True):
        if isinstance(x, quotient_ring_element.QuotientRingElement):
            if x.parent() is self:
                return x
            x = x.lift()
        if coerce:
            R = self.cover_ring()
            x = R(x)
        return quotient_ring_element.QuotientRingElement(self, x)

    def __cmp__(self, other):
        r"""
        Only quotients by the \emph{same} (in "is") ring
        and same ideal are considered equal.
        """
        if not isinstance(other, QuotientRing_generic):
            return -1
        if self.cover_ring() is other.cover_ring() and self.defining_ideal() is other.defining_ideal():
            return 0
        return 1

    def ngens(self):
        return self.cover_ring().ngens()

    def gen(self, i=0):
        return quotient_ring_element.QuotientRingElement(self, self.__R.gen(i))


