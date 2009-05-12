"""
Quotient Rings

AUTHORS:

- William Stein

TESTS::

    sage: R.<x> = PolynomialRing(ZZ)
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
import sage.rings.polynomial.multi_polynomial_ideal
import sage.structure.parent_gens
from sage.interfaces.all import singular as singular_default, is_SingularElement

def QuotientRing(R, I, names=None):
    r"""
    Creates a quotient ring of the ring R by the ideal I. Variables are
    labeled by names. (If the quotient ring is a quotient of a
    polynomial ring.). If names isn't given, 'bar' will be appended to
    the variable names in R.

    INPUTS:

    - ``R`` - a commutative ring

    - ``I`` - an ideal of R

    - ``names`` - a list of
      strings to be used as names for the variables in the quotient ring
      R/I

    OUTPUTS: R/I - the quotient ring R mod the ideal I

    EXAMPLES:

    Some simple quotient rings with the integers::

        sage: R = QuotientRing(ZZ,7*ZZ); R
        Quotient of Integer Ring by the ideal (7)
        sage: R.gens()
        (1,)
        sage: 1*R(3); 6*R(3); 7*R(3)
        3
        4
        0

    ::

        sage: S = QuotientRing(ZZ,ZZ.ideal(8)); S
        Quotient of Integer Ring by the ideal (8)
        sage: 2*S(4)
        0

    With polynomial rings: (note that the variable name of the quotient
    ring can be specified as shown below)

    ::

        sage: R.<xx> = QuotientRing(QQ[x], QQ[x].ideal(x^2 + 1)); R
        Univariate Quotient Polynomial Ring in xx over Rational Field with modulus x^2 + 1
        sage: R.gens(); R.gen()
        (xx,)
        xx
        sage: for n in range(4): xx^n
        1
        xx
        -1
        -xx

    ::

        sage: S = QuotientRing(QQ[x], QQ[x].ideal(x^2 - 2)); S
        Univariate Quotient Polynomial Ring in xbar over Rational Field with
        modulus x^2 - 2
        sage: xbar = S.gen(); S.gen()
        xbar
        sage: for n in range(3): xbar^n
        1
        xbar
        2

    Sage coerces objects into ideals when possible::

        sage: R = QuotientRing(QQ[x], x^2 + 1); R
        Univariate Quotient Polynomial Ring in xbar over Rational Field with
        modulus x^2 + 1

    By Noether's homomorphism theorems, the quotient of a quotient ring
    in R is just the quotient of R by the sum of the ideals. In this
    example, we end up modding out the ideal (x) from the ring
    QQ[x,y]::

        sage: R.<x,y> = PolynomialRing(QQ,2)
        sage: S.<a,b> = QuotientRing(R,R.ideal(1 + y^2))
        sage: T.<c,d> = QuotientRing(S,S.ideal(a))
        sage: T
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x, y^2 + 1)
        sage: R.gens(); S.gens(); T.gens()
        (x, y)
        (a, b)
        (0, d)
        sage: for n in range(4): d^n
        1
        d
        -1
        -d
    """
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
    if R.is_field():
        if I.is_zero():
            return R
        else:
            return 0
    return QuotientRing_generic(R, I, names)

def is_QuotientRing(x):
    """
    Tests whether or not x inherits from QuotientRing_generic.

    EXAMPLES::

        sage: from sage.rings.quotient_ring import is_QuotientRing
        sage: R.<x> = PolynomialRing(ZZ,'x')
        sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
        sage: S = R.quotient_ring(I)
        sage: is_QuotientRing(S)
        True
        sage: is_QuotientRing(R)
        False
    """
    return isinstance(x, QuotientRing_generic)

class QuotientRing_generic(commutative_ring.CommutativeRing, sage.structure.parent_gens.ParentWithGens):
    """
    The quotient ring of `R` by the ideal `I`.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(ZZ,'x')
        sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
        sage: S = R.quotient_ring(I); S
        Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 3*x + 4, x^2 + 1)

    ::

        sage: R.<x,y> = PolynomialRing(QQ)
        sage: S.<a,b> = R.quo(x^2 + y^2)
        sage: a^2 + b^2 == 0
        True
        sage: S(0) == a^2 + b^2
        True

    EXAMPLE: Quotient of quotient

    A quotient of a quotient is just the quotient of the original top
    ring by the sum of two ideals.

    ::

        sage: R.<x,y> = PolynomialRing(QQ,2)
        sage: S.<a,b> = R.quo(1 + y^2)
        sage: T.<c,d> = S.quo(a)
        sage: T
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x, y^2 + 1)
        sage: T.gens()
        (0, d)
    """
    def __init__(self, R, I, names):
        """
        Create the quotient ring of R by the ideal I.

        INPUT:


        -  ``R`` - a commutative ring

        -  ``I`` - an ideal


        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: R.quotient_ring(x^2 + y^2)
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
        """
        self.__R = R
        self.__I = I
        sage.structure.parent_gens.ParentWithGens.__init__(self, R.base_ring(), names)

    def construction(self):
        """
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ,'x')
            sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
            sage: R.quotient_ring(I).construction()
            (QuotientFunctor, Univariate Polynomial Ring in x over Integer Ring)

        TESTS::

            sage: F, R = Integers(5).construction()
            sage: F(R)
            Ring of integers modulo 5
            sage: F, R = GF(5).construction()
            sage: F(R)
            Finite Field of size 5
        """
        from sage.categories.pushout import QuotientFunctor
        # Is there a better generic way to distinguish between things like Z/pZ as a field and Z/pZ as a ring?
        from sage.rings.field import Field
        return QuotientFunctor(self.__I, as_field=isinstance(self, Field)), self.__R

    def _repr_(self):
        """
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ,'x')
            sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
            sage: R.quotient_ring(I)._repr_()
            'Quotient of Univariate Polynomial Ring in x over Integer Ring by the ideal (x^2 + 3*x + 4, x^2 + 1)'
        """
        return "Quotient of %s by the ideal %s"%(self.cover_ring(), self.defining_ideal()._repr_short())

    def _latex_(self):
        """
        EXAMPLES::

            sage: R.<x> = PolynomialRing(ZZ,'x')
            sage: I = R.ideal([4 + 3*x + x^2, 1 + x^2])
            sage: R.quotient_ring(I)._latex_()
            '\\Bold{Z}[x]/\\left(x^{2} + 3x + 4, x^{2} + 1\\right)\\Bold{Z}[x]'
        """
        return "%s/%s"%(latex.latex(self.cover_ring()), latex.latex(self.defining_ideal()))

    def cover(self):
        r"""
        The covering ring homomorphism `R \to R/I`, equipped with a
        section.

        EXAMPLES::

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

        EXAMPLES::

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
              From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2, y^2)
              To:   Multivariate Polynomial Ring in x, y over Rational Field
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

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S = R.quotient(x^2 + y^2)
            sage: pi = S.cover(); pi
            Ring morphism:
              From: Multivariate Polynomial Ring in x, y over Rational Field
              To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
              Defn: Natural quotient map
            sage: L = S.lift(); L
            Set-theoretic ring morphism:
              From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
              To:   Multivariate Polynomial Ring in x, y over Rational Field
              Defn: Choice of lifting map
            sage: L(S.0)
            x
            sage: L(S.1)
            y

        Note that some reduction may be applied so that the lift of a
        reduction need not equal the original element.

        ::

            sage: z = pi(x^3 + 2*y^2); z
            -xbar*ybar^2 + 2*ybar^2
            sage: L(z)
            -x*y^2 + 2*y^2
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
        r"""
        Return the characteristic of the quotient ring.

        TODO: Not yet implemented!

        EXAMPLES::

            sage: Q = QuotientRing(ZZ,7*ZZ)
            sage: Q.characteristic()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def defining_ideal(self):
        r"""
        Returns the ideal generating this quotient ring.

        EXAMPLES:

        In the integers::

            sage: Q = QuotientRing(ZZ,7*ZZ)
            sage: Q.defining_ideal()
            Principal ideal (7) of Integer Ring

        An example involving a quotient of a quotient. By Noether's
        homomorphism theorems, this is actually a quotient by a sum of two
        ideals::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = QuotientRing(R,R.ideal(1 + y^2))
            sage: T.<c,d> = QuotientRing(S,S.ideal(a))
            sage: S.defining_ideal()
            Ideal (y^2 + 1) of Multivariate Polynomial Ring in x, y over Rational Field
            sage: T.defining_ideal()
            Ideal (x, y^2 + 1) of Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self.__I

    def is_field(self):
        r"""
        Returns True if the quotient ring is a field. Checks to see if the
        defining ideal is maximal.

        TESTS:

        Requires the ``is_maximal`` function to be
        implemented::

            sage: Q = QuotientRing(ZZ,7*ZZ)
            sage: Q.is_field()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return self.defining_ideal().is_maximal()

    def is_integral_domain(self):
        r"""
        If this function returns True then self is definitely an integral
        domain. If it returns False, then either self is definitely not an
        integral domain or this function was unable to determine whether or
        not self is an integral domain.

        Use ``self.defining_ideal().is_prime()`` to find out
        for sure whether this quotient ring is really not an integral
        domain, of if Sage is unable to determine the answer.

        EXAMPLES::

            sage: R = Integers(8)
            sage: R.is_integral_domain()
            False
            sage: R.<a,b,c> = ZZ['a','b','c']
            sage: I = R.ideal(a,b)
            sage: Q = R.quotient_ring(I)
            sage: Q.is_integral_domain()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return self.defining_ideal().is_prime()

    def cover_ring(self):
        r"""
        Returns the cover ring of the quotient ring: that is, the original
        ring R from which we modded out an ideal, I.

        TODO: PolynomialQuotientRings_field objects don't have a
        ``cover_ring`` function.

        EXAMPLES::

            sage: Q = QuotientRing(ZZ,7*ZZ)
            sage: Q.cover_ring()
            Integer Ring

        ::

            sage: Q = QuotientRing(QQ[x], x^2 + 1)
            sage: Q.cover_ring()
            Traceback (most recent call last):
            ...
            AttributeError: 'PolynomialQuotientRing_field' object has no attribute 'cover_ring'
        """
        return self.__R

    def ideal(self, *gens, **kwds):
        """
        Return the ideal of self with the given generators.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2+y^2)
            sage: S.ideal()
            Ideal (0) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
            sage: S.ideal(x+y+1)
            Ideal (xbar + ybar + 1) of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
        """
        if len(gens) == 1:
            gens = gens[0]
        from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomialRing_libsingular
        if not isinstance(self.__R,MPolynomialRing_libsingular) and not self.__R._has_singular:
            # pass through
            MPolynomialRing_generic.ideal(self,gens,**kwds)
        if is_SingularElement(gens):
            gens = list(gens)
            coerce = True
        elif not isinstance(gens, (list, tuple)):
            gens = [gens]
        if kwds.has_key('coerce') and kwds['coerce']:
            gens = [self(x) for x in gens]  # this will even coerce from singular ideals correctly!
        return sage.rings.polynomial.multi_polynomial_ideal.MPolynomialIdeal(self, gens, **kwds)

    def __call__(self, x, coerce=True):
        """
        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2+y^2)
            sage: S(x)
            xbar
            sage: S(x^2 + y^2)
            0
        """
        if isinstance(x, quotient_ring_element.QuotientRingElement):
            if x.parent() is self:
                return x
            x = x.lift()
        if is_SingularElement(x):
            #self._singular_().set_ring()
            x = quotient_ring_element.QuotientRingElement(self, x.sage_poly(self.cover_ring()))
            return x
        if coerce:
            R = self.cover_ring()
            x = R(x)
        return quotient_ring_element.QuotientRingElement(self, x)

    def _coerce_impl(self, x):
        """
        Return the coercion of x into this quotient ring.

        The rings that coerce into the quotient ring canonically, are:

        - this ring

        - anything that coerces into the ring of which this is the
          quotient

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = R.quotient(x^2 + y^2)
            sage: S._coerce_(0)
            0
            sage: S._coerce_(2/3)
            2/3
            sage: S._coerce_(a^2 - b)
            -b^2 - b
            sage: S._coerce_(GF(7)(3))
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion of element into self
        """
        return self._coerce_try(x, [self.cover_ring()])

    def __cmp__(self, other):
        r"""
        Only quotients by the *same* ring and same ideal (with the same
        generators!!) are considered equal.

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2 + y^2)
            sage: S == R.quotient_ring(x^2 + y^2)
            True

        The ideals `(x^2 + y^2)` and `(-x^2-y^2)` are
        equal, but since the generators are different, the corresponding
        quotient rings are not equal::

            sage: R.ideal(x^2+y^2) == R.ideal(-x^2 - y^2)
            True
            sage: R.quotient_ring(x^2 + y^2) == R.quotient_ring(-x^2 - y^2)
            False
        """
        if not isinstance(other, QuotientRing_generic):
            return cmp(type(self), type(other))
        return cmp((self.cover_ring(), self.defining_ideal().gens()),
                   (other.cover_ring(), other.defining_ideal().gens()))

    def ngens(self):
        r"""
        Returns the number of generators for this quotient ring.

        TODO: Note that ``ngens`` counts 0 as a generator. Does
        this make sense? That is, since 0 only generates itself and the
        fact that this is true for all rings, is there a way to "knock it
        off" of the generators list if a generator of some original ring is
        modded out?

        EXAMPLES::

            sage: R = QuotientRing(ZZ,7*ZZ)
            sage: R.gens(); R.ngens()
            (1,)
            1

        ::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = QuotientRing(R,R.ideal(1 + y^2))
            sage: T.<c,d> = QuotientRing(S,S.ideal(a))
            sage: T
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x, y^2 + 1)
            sage: R.gens(); S.gens(); T.gens()
            (x, y)
            (a, b)
            (0, d)
            sage: R.ngens(); S.ngens(); T.ngens()
            2
            2
            2
        """
        return self.cover_ring().ngens()

    def gen(self, i=0):
        r"""
        Returns the ith generator for this quotient ring.

        EXAMPLES::

            sage: R = QuotientRing(ZZ,7*ZZ)
            sage: R.gen(0)
            1

        ::

            sage: R.<x,y> = PolynomialRing(QQ,2)
            sage: S.<a,b> = QuotientRing(R,R.ideal(1 + y^2))
            sage: T.<c,d> = QuotientRing(S,S.ideal(a))
            sage: T
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x, y^2 + 1)
            sage: R.gen(0); R.gen(1)
            x
            y
            sage: S.gen(0); S.gen(1)
            a
            b
            sage: T.gen(0); T.gen(1)
            0
            d
        """
        return self(self.__R.gen(i))


    def _singular_(self, singular=singular_default):
        """
        Returns the Singular quotient ring of self if the base ring is
        coercable to Singular.

        If a valid singular representation is found it is used otherwise a
        new 'qring' is created.

        INPUT:


        -  ``singular`` - Singular instance (default: the
           default Singular instance)


        .. note::

           This method also sets the current ring in Singular to self

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2+y^2)
            sage: S._singular_()
            //   characteristic : 0
            //   number of vars : 2
            //        block   1 : ordering dp
            //                  : names    x y
            //        block   2 : ordering C
            // quotient ring from ideal
            _[1]=x2+y2
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
        Returns a newly created Singular quotient ring matching self if the
        base ring is coecable to Singular.

        See self._singular_

        EXAMPLES::

            sage: R.<x,y> = PolynomialRing(QQ)
            sage: S = R.quotient_ring(x^2+y^2)
            sage: T = S._singular_init_()
            sage: parent(S)
            <class 'sage.rings.quotient_ring.QuotientRing_generic'>
            sage: parent(T)
            Singular
        """
        self.__R._singular_().set_ring()
        self.__singular = singular("%s"%self.__I._singular_().name(),"qring")
        return self.__singular

    def _magma_init_(self, magma):
        r"""
        Return string that evaluates to Magma version of this quotient
        ring. This is called implicitly when doing conversions to Magma.

        INPUT:


        -  ``magma`` - a magma instance


        EXAMPLE::

            sage: P.<x,y> = PolynomialRing(GF(2))
            sage: Q = P.quotient(sage.rings.ideal.FieldIdeal(P))
            sage: magma(Q)                                         # optional - magma
            Affine Algebra of rank 2 over GF(2)
            Graded Reverse Lexicographical Order
            Variables: x, y
            Quotient relations:
            [
            x^2 + x,
            y^2 + y
            ]
        """
        R = magma(self.__R)
        I = magma(self.__I.gens())
        return "quo<%s|%s>"%(R.name(), I._ref())
