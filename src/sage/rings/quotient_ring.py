"""
Quotient Rings

AUTHOR:
    -- William Stein

TESTS:
    sage: R = PolynomialRing(ZZ,'x')
    sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
    sage: S = R.quotient_ring(I);
    sage: S == loads(dumps(S))
    True


"""

################################################################################
#
#   SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
################################################################################

import quotient_ring_element
import sage.misc.latex as latex
import commutative_ring
import ideal
import multi_polynomial_ideal
import sage.structure.parent_gens
from sage.interfaces.all import singular as singular_default, is_SingularElement

def QuotientRing(R, I, names=None):
    if not isinstance(R, commutative_ring.CommutativeRing):
        raise TypeError, "R must be a commutative ring."
    if names is None:
        names = tuple([x + 'bar' for x in R.variable_names()])
    else:
        names = sage.structure.parent_gens.normalize_names(R.ngens(), names)
    if not isinstance(I, ideal.Ideal_generic) or I.ring() != R:
        I = R.ideal(I)
    try:
        if I.is_principal():
            return R.quotient_by_principal_ideal(I.gen(), names)
    except (AttributeError, NotImplementedError):
        pass
    if isinstance(R, QuotientRing_generic):
        pi = R.cover()
        S = pi.domain()
        G = [pi.lift(x) for x in I.gens()]
        I_lift = S.ideal(G)
        J = R.defining_ideal()
        return QuotientRing_generic(S, I_lift + J, names)

    return QuotientRing_generic(R, I, names)

def is_QuotientRing(x):
    return isinstance(x, QuotientRing_generic)

class QuotientRing_generic(commutative_ring.CommutativeRing, sage.structure.parent_gens.ParentWithGens):
    """
    The quotient ring of $R$ by the ideal $I$.

    EXAMPLES:
        sage: R.<x> = PolynomialRing(ZZ,'x')
        sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
        sage: S = R.quotient_ring(I); S
        Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 1, x^2 + 3*x + 4)

        sage: R.<x,y> = PolynomialRing(QQ)
        sage: S.<a,b> = R.quo(x^2 + y^2)
        sage: a^2 + b^2 == 0
        True
        sage: S(0) == a^2 + b^2
        True

    EXAMPLE: Quotient of quotient

    A quotient of a quotient is just the quotient of the
    original top ring by the sum of two ideals.
        sage: R.<x,y> = PolynomialRing(QQ,2)
        sage: S.<a,b> = R.quo(1 + y^2)
        sage: T.<c,d> = S.quo(a)
        sage: T
        Quotient of Polynomial Ring in x, y over Rational Field by the ideal (x, 1 + y^2)
        sage: T.gens()
        (0, d)
    """
    def __init__(self, R, I, names):
        """
        Create the quotient ring of R by the ideal I.

        INPUT:
            R -- a commutative ring
            I -- an ideal
        """
        self.__R = R
        self.__I = I
        sage.structure.parent_gens.ParentWithGens.__init__(self, R.base_ring(), names)

    def _repr_(self):
        return "Quotient of %s by the ideal %s"%(self.cover_ring(), self.defining_ideal()._repr_short())

    def _latex_(self):
        return "%s/%s"%(latex.latex(self.cover_ring()), latex.latex(self.defining_ideal()))

    def cover(self):
        r"""
        The covering ring homomorphism $R \to R/I$, equipped with a section.

        EXAMPLES:
            sage: R = ZZ.quo(3*ZZ)
            sage: pi = R.cover()
            sage: pi
            Ring morphism:
              From: Integer Ring
              To:   Ring of integers modulo 3
              Defn: Natural quotient map
            sage: pi(5)
            2
            sage: l = pi.lift()

        EXAMPLES:
            sage: R.<x,y>  = PolynomialRing(QQ)
            sage: Q = R.quo( (x^2,y^2) )
            sage: pi = Q.cover()
            sage: pi(x^3+y)
            ybar
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
            pi = morphism.RingHomomorphism_cover(self.__R.Hom(self))
            lift = self.lift()
            pi._set_lift(lift)
            self.__cover = pi
            return self.__cover

    def lift(self):
        """
        Return the lifting map to the cover.

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S = R.quotient(x^2 + y^2)
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
        r"""
        If this function returns True then self is definitely an
        integral domain.  If it returns False, then either self is
        definitely not an integral domain or this function was unable
        to determine whether or not self is an integral domain.

        Use \code{self.defining_ideal().is_prime()} to find out for
        sure whether this quotient ring is really not an integral
        domain, of if SAGE is unable to determine the answer.

        EXAMPLES:
            sage: R = Integers(8)
            sage: R.is_integral_domain()
            False
        """
        try:
            return self.defining_ideal().is_prime()
        except NotImplementedError:
            return False

    def cover_ring(self):
        return self.__R

    def ideal(self, gens, coerce=False):
        if not self.__R._has_singular:
            # pass through
            MPolynomialRing_generic.ideal(self,gens,coerce)
        if is_SingularElement(gens):
            gens = list(gens)
            coerce = True
        elif not isinstance(gens, (list, tuple)):
            gens = [gens]
        if coerce:
            gens = [self(x) for x in gens]  # this will even coerce from singular ideals correctly!
        return multi_polynomial_ideal.MPolynomialIdeal(self, gens, coerce=False)

    def _can_convert_to_singular(self):
        return self.__R._can_convert_to_singular()

    def __call__(self, x, coerce=True):
        if isinstance(x, quotient_ring_element.QuotientRingElement):
            if x.parent() is self:
                return x
            x = x.lift()
        if is_SingularElement(x):
            #self._singular_().set_ring()
            x = x.sage_poly(self)
            return x
        if coerce:
            R = self.cover_ring()
            x = R(x)
        return quotient_ring_element.QuotientRingElement(self, x)

    def _coerce_impl(self, x):
        """
        Return the coercion of x into this quotient ring.

        The rings that coerce into the quotient ring canonically, are:

           * this ring.
           * anything that coerces into the ring of which this is the quotient

        EXAMPLES:
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = R.quotient(x^2 + y^2)
            sage: S._coerce_(0)
            0
            sage: S._coerce_(2/3)
            2/3
            sage: S._coerce_(a^2 - b)
            -1*b - b^2
            sage: S._coerce_(GF(7)(3))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of element into self
        """
        return self._coerce_try(x, [self.cover_ring()])

    def __cmp__(self, other):
        r"""
        Only quotients by the \emph{same} ring and same ideal (with
        the same generators!!) are considered equal.
        """
        if not isinstance(other, QuotientRing_generic):
            return cmp(type(self), type(other))
        return cmp((self.cover_ring(), self.defining_ideal().gens()),
                   (other.cover_ring(), other.defining_ideal().gens()))

    def ngens(self):
        return self.cover_ring().ngens()

    def gen(self, i=0):
        return quotient_ring_element.QuotientRingElement(self, self.__R.gen(i))


    def _singular_(self, singular=singular_default):
        """
        Returns the Singular quotient ring of self if the base ring is
        coercable to Singular.

        If a valid singular representation is found it is used
        otherwise a new 'qring' is created.

        INPUT:
            singular -- Singular instance (default: the default Singular instance)

        \note{This method also sets the current ring in Singular to self}
        """

        try:
            Q = self.__singular
            if not (Q.parent() is singular):
                raise ValueError
            Q._check_valid()
            return Q
        except (AttributeError, ValueError):
            return self._singular_init_(singular)

    def _singular_init_(self,singular=singular_default):
        """
        Returns a newly created Singular quotient ring matching self
        if the base ring is coecable to Singular.

        See self._singular_

        """
        self.__R._singular_().set_ring()
        self.__singular = singular("%s"%self.__I._singular_().name(),"qring")
        return self.__singular

