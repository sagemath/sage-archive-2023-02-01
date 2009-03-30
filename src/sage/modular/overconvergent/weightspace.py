r"""
The space of `p`-adic weights

A `p`-adic weight is a continuous character `\mathbb{Z}_p^\times \to
\mathbb{C}_p^\times`. These are the `\mathbb{C}_p`-points of a rigid space over
`\mathbb{Q}_p`, which is isomorphic to a disjoint union of copies (indexed by
`(\mathbb{Z}/p\mathbb{Z})^\times`) of the open unit `p`-adic disc.

Sage supports both "classical points", which are determined by the data of a
Dirichlet character modulo `p^m` for some `m` and an integer `k` (corresponding
to the character `z \mapsto z^k \chi(z)`) and "non-classical points" which are
determined by the data of an element of `(\mathbb{Z}/p\mathbb{Z})^\times` and
an element `w \in \mathbb{C}_p` with `|w - 1| < 1`.
"""

from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.modular.dirichlet import DirichletGroup, DirichletCharacter, trivial_character
from sage.rings.all import ZZ, QQ, divisors, IntegerModRing, Qp
from sage.misc.misc import sxrange
from sage.rings.padics.padic_generic_element import pAdicGenericElement
from sage.misc.misc import verbose
import weakref

_wscache = {}
def WeightSpace_constructor(p):
    r"""
    Construct the p-adic weight space for the given prime p.

    EXAMPLES::

        sage: WeightSpace(3) # indirect doctest
        Space of 3-adic weight-characters
        sage: WeightSpace(10)
        Traceback (most recent call last):
        ...
        ValueError: p must be prime
    """
    if _wscache.has_key(p):
        m = _wscache[p]()
        if m is not None:
            return m
    m = WeightSpace_class(p)
    _wscache[p] = weakref.ref(m)
    return m

class WeightSpace_class(Parent):
    r"""
    The space of `p`-adic weight-characters `\mathcal{W} = {\rm
    Hom}(\mathbb{Z}_p^\times, \mathbb{C}_p^\times)`. This isomorphic to a
    disjoint union of `(p-1)` open discs of radius 1 (or 2 such discs if `p =
    2`), with the parameter on the open disc corresponding to the image of `1 +
    p` (or 5 if `p = 2`)

    TESTS::

        sage: W = WeightSpace(3)
        sage: W is loads(dumps(W))
        True
    """

    def __init__(self, p):
        r"""
        Initialisation function.

        EXAMPLE::

            sage: WeightSpace(17)
            Space of 17-adic weight-characters
        """
        p = ZZ(p)
        if not p.is_prime():
            raise ValueError, "p must be prime"
        self.p = p
        self.param = Qp(p)((p == 2 and 5) or (p + 1))

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: WeightSpace(17)._repr_()
            'Space of 17-adic weight-characters'
        """
        return "Space of %s-adic weight-characters" % self.p

    def __reduce__(self):
        r"""
        Used for pickling.

        EXAMPLE::

            sage: WeightSpace(3).__reduce__()
            (<function WeightSpace_constructor at ...>, (3,))
        """
        return (WeightSpace_constructor, (self.p,))

    def __call__(self, arg1, arg2 = None, algebraic=True):
        r"""
        Create an element of this space.

        If ``algebraic = True`` (the default), create a locally algebraic
        character. The arguments should be `(k, \chi)` with `k \in \mathbb{Z}`
        and `\chi` a Dirichlet character of `p`-power conductor defined over a
        `p`-adic field; this corresponds to the weight-character `x \mapsto x^k
        \chi(x)`. If `\chi` is omitted, it defaults to the trivial character.

        If ``algebraic = False``, create a general character. The arguments are
        now (t, w) where `t \in \mathbb{Z}/(p-1)\mathbb{Z}` and `w \in
        \mathbb{C}_p` with `|w - 1| < 1`. This corresponds to the character
        `\kappa` satisfying `\kappa(\mu) = \mu^t` where `\mu` is a `(p-1)`-st
        root of unity, and `\kappa(1 + p) = w`.

        EXAMPLES::

            sage: W = WeightSpace(17)
            sage: W(4)
            4
            sage: W(4, DirichletGroup(17, Qp(17)).0)
            (4, 17, [3 + 13*17 + ... + O(17^20)])
            sage: W(1 + O(17^5), 4, algebraic = False)
            [1 + O(17^5), 4]
        """

        if isinstance(arg1, WeightCharacter):
            if arg1.parent().p != self.p:
                raise TypeError
            else:
                arg1._parent = self
                return arg1

        if algebraic:
            return AlgebraicWeight(self, arg1, arg2)
        else:
            return ArbitraryWeight(self, arg1, arg2)

class WeightCharacter(Element):
    r"""
    Abstract base class representing an element of the p-adic weight space
    `Hom(\mathbb{Z}_p^\times, \mathbb{C}_p^\times)`.
    """

    # This should probably derive from Morphism or even from
	# AbelianGroupMorphism; but SAGE doesn't know about the abelian group
	# Z_p^*, so Hom(Z_p^*, C_p^*) is a bit beyond it!

    def __init__(self, parent):
        r"""
        Initialisation function.

        EXAMPLE::

            sage: WeightSpace(17)(0)
            0
        """

        Element.__init__(self, parent)
        self.prime = self.parent().p

    def is_even(self):
        r"""
        Return True if this weight-character sends -1 to +1.

        EXAMPLE::

            sage: WeightSpace(17)(0).is_even()
            True
            sage: WeightSpace(17)(11).is_even()
            False
            sage: WeightSpace(17)(1 + 17 + O(17^20), 3, False).is_even()
            False
            sage: WeightSpace(17)(1 + 17 + O(17^20), 4, False).is_even()
            True
        """
        if self(-1) == -1:
            return False
        else:
            return True

    def pAdicEisensteinSeries(self, ring, prec=20):
        r"""
        Calculate the q-expansion of the p-adic Eisenstein series of given
        weight-character, normalised so the constant term is 1.

        EXAMPLE::

            sage: kappa = WeightSpace(3)(3, DirichletGroup(3,QQ).0)
            sage: kappa.pAdicEisensteinSeries(QQ[['q']], 20)
            1 - 9*q + 27*q^2 - 9*q^3 - 117*q^4 + 216*q^5 + 27*q^6 - 450*q^7 + 459*q^8 - 9*q^9 - 648*q^10 + 1080*q^11 - 117*q^12 - 1530*q^13 + 1350*q^14 + 216*q^15 - 1845*q^16 + 2592*q^17 + 27*q^18 - 3258*q^19 + O(q^20)
        """
        if not self.is_even():
            raise ValueError, "Eisenstein series not defined for odd weight-characters"
        q = ring.gen()
        s = ring(1) + 2*self.one_over_Lvalue() * sum([sum([self(d)/d for d in divisors(n)]) * q**n for n in xrange(1, prec)])
        return s.add_bigoh(prec)

    def values_on_gens(self):
        r"""
        If `\kappa` is this character, calculate the values `(\kappa(r), t)`
        where `r` is `1 + p` (or 5 if `p = 2`) and `t` is the unique element of
        `\mathbb{Z}/(p-1)\mathbb{Z}` such that `\kappa(\mu) = \mu^t` for `\mu`
        a (p-1)st root of unity. (If `p = 2`, we take `t` to be 0 or 1
        according to whether `\kappa` is odd or even.) These two values
        uniquely determine the character `\kappa`.

        EXAMPLES::

            sage: W=WeightSpace(11); W(2).values_on_gens()
            (1 + 2*11 + 11^2 + O(11^20), 2)
            sage: W(2, DirichletGroup(11, QQ).0).values_on_gens()
            (1 + 2*11 + 11^2 + O(11^20), 7)
            sage: W(1 + 2*11 + O(11^5), 4, algebraic = False).values_on_gens()
            (1 + 2*11 + O(11^5), 4)
        """

        return ( self(self.parent().param), self.teichmuller_type())

    def is_trivial(self):
        r"""
        Return True if and only if this is the trivial character.

        EXAMPLES::

            sage: WeightSpace(11)(2).is_trivial()
            False
            sage: WeightSpace(11)(2, DirichletGroup(11, QQ).0).is_trivial()
            False
            sage: WeightSpace(11)(0).is_trivial()
            True
        """
        if self.values_on_gens() == (1, 0):
            return True
        else:
            return False

    def __cmp__(self, other):
        r"""
        Compare self to other.

        EXAMPLES::

            sage: W=WeightSpace(11)
            sage: W(2) == W(3)
            False
            sage: W(2, DirichletGroup(11, QQ).0) == W(2)
            False
            sage: W(2, DirichletGroup(11, QQ).0) == W(144 + O(11^20), 7, False)
            True
        """
        if not isinstance(other, WeightCharacter):
            return cmp(type(self), type(other))
        else:
            return cmp(self.values_on_gens(), other.values_on_gens())

    def Lvalue(self):
        r"""
        Return the value of the p-adic L-function of `\mathbb{Q}`, which can be
        regarded as a rigid-analytic function on weight space, evaluated at
        this character.

        EXAMPLES::

            sage: W = WeightSpace(11)
            sage: sage.modular.overconvergent.weightspace.WeightCharacter(W).Lvalue()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def one_over_Lvalue(self):
        r"""
        Return the reciprocal of the p-adic L-function evaluated at this
        weight-character. If the weight-character is odd, then the L-function
        is zero, so an error will be raised.

        EXAMPLES::

            sage: WeightSpace(11)(4).one_over_Lvalue()
            -12/133
            sage: WeightSpace(11)(3, DirichletGroup(11, QQ).0).one_over_Lvalue()
            -1/6
            sage: WeightSpace(11)(3).one_over_Lvalue()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Rational division by zero
        """
        if self.is_trivial():
            return 0
        else:
            return 1/self.Lvalue()

class AlgebraicWeight(WeightCharacter):
    r"""
    A point in weight space corresponding to a locally algebraic character, of
    the form `x \mapsto \chi(x) x^k` where `k` is an integer and `\chi` is a
    Dirichlet character modulo `p^n` for some `n`.

    TESTS::

        sage: w = WeightSpace(23)(12, DirichletGroup(23, QQ).0) # exact
        sage: w == loads(dumps(w))
        True
        sage: w = WeightSpace(23)(12, DirichletGroup(23, Qp(23)).0) # inexact
        sage: w == loads(dumps(w))
        True
        sage: w is loads(dumps(w)) # elements are not globally unique
        False
    """

    def __init__(self, parent, k, chi=None):
        r"""
        Create a locally algebraic weight-character.

        EXAMPLES::

            sage: WeightSpace(29)(13, DirichletGroup(29, Qp(29)).0)
            (13, 29, [2 + 2*29 + ... + O(29^20)])
        """
        WeightCharacter.__init__(self, parent)
        k = ZZ(k)
        self.k = k
        if chi is None:
            chi = trivial_character(self.prime, QQ)
        n = ZZ(chi.conductor())
        if n == 1:
            n = self.prime
        if not n.is_power_of(self.prime):
            raise ValueError, "Character must have %s-power conductor" % p
        self.chi = DirichletGroup(n, chi.base_ring())(chi)

    def __call__(self, x):
        r"""
        Evaluate this character at an element of `\mathbb{Z}_p^\times`.

        EXAMPLES:

        Exact answers are returned when this is possible::

            sage: kappa = WeightSpace(29)(13, DirichletGroup(29, QQ).0)
            sage: kappa(1)
            1
            sage: kappa(0)
            0
            sage: kappa(12)
            -106993205379072
            sage: kappa(-1)
            -1
            sage: kappa(13 + 4*29 + 11*29^2 + O(29^3))
            9 + 21*29 + 27*29^2 + O(29^3)

        When the character chi is defined over a p-adic field, the results returned are inexact::

            sage: kappa = WeightSpace(29)(13, DirichletGroup(29, Qp(29)).0^14)
            sage: kappa(1)
            1 + O(29^20)
            sage: kappa(0)
            0
            sage: kappa(12)
            17 + 11*29 + 7*29^2 + 4*29^3 + 5*29^4 + 2*29^5 + 13*29^6 + 3*29^7 + 18*29^8 + 21*29^9 + 28*29^10 + 28*29^11 + 28*29^12 + 28*29^13 + 28*29^14 + 28*29^15 + 28*29^16 + 28*29^17 + 28*29^18 + 28*29^19 + O(29^20)
            sage: kappa(12) == -106993205379072
            True
            sage: kappa(-1) == -1
            True
            sage: kappa(13 + 4*29 + 11*29^2 + O(29^3))
            9 + 21*29 + 27*29^2 + O(29^3)
        """
        if isinstance(x, pAdicGenericElement):
            if x.parent().prime() != self.prime:
                raise TypeError, "x must be an integer or a %s-adic integer" % self.prime
            if self.prime**(x.precision_absolute()) < self.chi.conductor():
                raise Exception, "Precision too low"
            xint = x.lift()
        else:
            xint = x
        if (xint % self.prime == 0): return 0
        return self.chi(xint) * x**self.k

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: WeightSpace(17)(2)._repr_()
            '2'
            sage: WeightSpace(17)(2, DirichletGroup(17, QQ).0)._repr_()
            '(2, 17, [-1])'
            sage: WeightSpace(17)(2, DirichletGroup(17, QQ).0^2)._repr_()
            '2'
        """
        if self.chi.is_trivial():
            return "%s" % self.k
        else:
            return "(%s, %s, %s)" % (self.k, self.chi.modulus(), self.chi)

    def teichmuller_type(self):
        r"""
        Return the Teichmuller type of this weight-character `\kappa`, which is
        the unique `t \in \mathbb{Z}/(p-1)\mathbb{Z}` such that `\kappa(\mu) =
        \mu^t` for \mu a `(p-1)`-st root of 1.

        For `p = 2` this doesn't make sense, but we still want the Teichmuller
        type to correspond to the index of the component of weight space in
        which `\kappa` lies, so we return 1 if `\kappa` is odd and 0 otherwise.

        EXAMPLE::

            sage: WeightSpace(11)(2, DirichletGroup(11,QQ).0).teichmuller_type()
            7
            sage: WeightSpace(29)(13, DirichletGroup(29, Qp(29)).0).teichmuller_type()
            14
            sage: WeightSpace(2)(3, DirichletGroup(4,QQ).0).teichmuller_type()
            0
        """
        # Special case p == 2
        if self.prime == 2:
            if self.is_even():
                return IntegerModRing(2)(0)
            else:
                return IntegerModRing(2)(1)
        m = IntegerModRing(self.prime).multiplicative_generator()
        x = [y for y in IntegerModRing(self.chi.modulus()) if y == m and y**(self.prime - 1) == 1]
        if len(x) != 1: raise ArithmeticError
        x = x[0]
        f = IntegerModRing(self.prime)(self.chi(x)).log(m)
        return IntegerModRing(self.prime - 1)(self.k + f)

    def Lvalue(self):
        r"""
        Return the value of the p-adic L-function of `\mathbb{Q}` evaluated at
        this weight-character. If the character is `x \mapsto x^k \chi(x)`
        where `k > 0` and `\chi` has conductor a power of `p`, this is an
        element of the number field generated by the values of `\chi`, equal to
        the value of the complex L-function `L(1-k, \chi)`. If `\chi` is
        trivial, it is equal to `(1 - p^{k-1})\zeta(1-k)`.

        At present this is not implemented in any other cases, except the
        trivial character (for which the value is `\infty`).

        TODO: Implement this more generally using the Amice transform machinery
        in sage/schemes/elliptic_curves/padic_lseries.py, which should clearly
        be factored out into a separate class.

        EXAMPLES::

            sage: WeightSpace(7)(4).Lvalue() == (1 - 7^3)*zeta__exact(-3)
            True
            sage: WeightSpace(7)(5, DirichletGroup(7, Qp(7)).0^4).Lvalue()
            0
            sage: WeightSpace(7)(6, DirichletGroup(7, Qp(7)).0^4).Lvalue()
            1 + 2*7 + 7^2 + 3*7^3 + 3*7^5 + 4*7^6 + 2*7^7 + 5*7^8 + 2*7^9 + 3*7^10 + 6*7^11 + 2*7^12 + 3*7^13 + 5*7^14 + 6*7^15 + 5*7^16 + 3*7^17 + 6*7^18 + O(7^19)
        """
        if self.k > 0:
            return -self.chi.bernoulli(self.k)/self.k
        if self.is_trivial():
            return Infinity
        else:
            raise NotImplementedError, "Don't know how to compute value of this L-function"

class ArbitraryWeight(WeightCharacter):

    def __init__(self, parent, w, t):
        r"""
        Create the element of p-adic weight space in the given component
        mapping 1 + p to w. Here w must be an element of a p-adic field, with
        finite precision.

        EXAMPLE::

            sage: WeightSpace(17)(1 + 17^2 + O(17^3), 11, False)
            [1 + 17^2 + O(17^3), 11]
        """
        WeightCharacter.__init__(self, parent)

        self.t = ZZ(t) % (self.prime > 2 and (self.prime - 1) or 2)
		# do we store w precisely?
        if (w - 1).valuation() <= 0:
            raise ValueError, "Must send generator to something nearer 1"
        self.w = w

    def _repr_(self):
        r"""String representation of this character.

        EXAMPLES::

            sage: WeightSpace(97)(1 + 2*97 + O(97^20), 12, False)._repr_()
            '[1 + 2*97 + O(97^20), 12]'
        """
        return "[%s, %s]" % (self.w, self.t)

    def __call__(self, x):
        r"""
        Evaluate this character at an element of `\mathbb{Z}_p^\times`.

        EXAMPLES::

            sage: kappa = WeightSpace(23)(1 + 23^2 + O(23^20), 4, False)
            sage: kappa(2)
            16 + 7*23 + 7*23^2 + 16*23^3 + 23^4 + 20*23^5 + 15*23^7 + 11*23^8 + 12*23^9 + 8*23^10 + 22*23^11 + 16*23^12 + 13*23^13 + 4*23^14 + 19*23^15 + 6*23^16 + 7*23^17 + 11*23^19 + O(23^20)
            sage: kappa(-1)
            1 + O(23^20)
            sage: kappa(23)
            0
            sage: kappa(2 + 2*23 + 11*23^2 + O(23^3))
            16 + 7*23 + O(23^3)
        """

        if not isinstance(x, pAdicGenericElement):
            x = Qp(self.prime)(x)
        if x.valuation() != 0:
            return 0

        teich = x.parent().teichmuller(x)
        xx = x / teich
        if (xx - 1).valuation() <= 0:
            raise ArithmeticError
        verbose("Normalised element is %s" % xx)

        e = xx.log() / self.parent().param.log()
        verbose("Exponent is %s" % e)

        return teich**(self.t) * (self.w.log() * e).exp()

    def teichmuller_type(self):
        r"""
        Return the Teichmuller type of this weight-character `\kappa`, which is
        the unique `t \in \mathbb{Z}/(p-1)\mathbb{Z}` such that `\kappa(\mu) =
        \mu^t` for \mu a `(p-1)`-st root of 1.

        For `p = 2` this doesn't make sense, but we still want the Teichmuller
        type to correspond to the index of the component of weight space in
        which `\kappa` lies, so we return 1 if `\kappa` is odd and 0 otherwise.

        EXAMPLES::

            sage: WeightSpace(17)(1 + 3*17 + 2*17^2 + O(17^3), 8, False).teichmuller_type()
            8
            sage: WeightSpace(2)(1 + 2 + O(2^2), 1, False).teichmuller_type()
            1
        """
        return self.t

