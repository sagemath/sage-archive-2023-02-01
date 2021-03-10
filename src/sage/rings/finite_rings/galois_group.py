r"""
Galois groups of Finite Fields
"""

from sage.groups.abelian_gps.abelian_group import AbelianGroup_class
from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
from sage.groups.galois_group import _GaloisMixin
from sage.rings.integer_ring import ZZ
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic, FrobeniusEndomorphism_finite_field

class GaloisGroup_GFElement(AbelianGroupElement):
    def as_hom(self):
        r"""
        Return the automorphism of the finite field corresponding to this element.

        EXAMPLES::

            sage: GF(3^6).galois_group()([4]).as_hom()
            Frobenius endomorphism z6 |--> z6^(3^4) on Finite Field in z6 of size 3^6
        """
        n = self.exponents()[0]
        return self.parent()._field.frobenius_endomorphism(n)

    def __call__(self, x):
        r"""
        Return the action of this automorphism on an element `x` of the finite field.

        EXAMPLES::

            sage: k.<a> = GF(3^6)
            sage: g = k.galois_group()([4])
            sage: g(a) == a^(3^4)
            True
        """
        return self.as_hom()(x)

    def fixed_field(self):
        r"""
        The fixed field of this automorphism.

        EXAMPLES::

            sage: k.<a> = GF(3^12)
            sage: g = k.galois_group()([8])
            sage: k0, embed = g.fixed_field()
            sage: k0.cardinality()
            81
            sage: embed.domain() is k0
            True
            sage: embed.codomain() is k
            True
        """
        return self.as_hom().fixed_field()

# We don't inherit from sage.groups.galois_group.GaloisGroup since we want to use Sage's AbelianGroup instead
class GaloisGroup_GF(AbelianGroup_class, _GaloisMixin):
    r"""
    The Galois group of a finite field.
    """
    Element = GaloisGroup_GFElement

    def __init__(self, finite_field):
        r"""
        Create a Galois group.

        TESTS::

            sage: TestSuite(GF(9).galois_group()).run()
        """
        self._field = finite_field
        AbelianGroup_class.__init__(self, (finite_field.degree(),), names="Frob")

    def order(self, algorithm=None, recompute=False):
        r"""
        Return the order of this Galois group, which is just the degree of the extension since finite fields are Galois.

        EXAMPLES::

            sage: GF(9).galois_group().order()
            2
        """
        return self._field.degree()

    def _repr_(self):
        r"""
        String representation of this Galois group

        EXAMPLES::

            sage: GF(9).galois_group()
            Galois group C2 of GF(3^2)
        """
        return "Galois group C{0} of GF({1}^{0})".format(self._field.degree(), self._field.characteristic())

    def is_galois(self):
        r"""
        All finite fields are Galois.

        For compatibility with Galois groups of number fields.

        EXAMPLES::

            sage: GF(9).galois_group().is_galois()
            True
        """
        return True

    def transitive_number(self, algorithm=None, recompute=False):
        r"""
        Returns the transitive number for the action of the roots of the defining polynomial.

        EXAMPLES::

            sage: GF(2^8).galois_group().transitive_number()
            1
            sage: GF(3^32).galois_group().transitive_number()
            33
            sage: GF(2^60).galois_group().transitive_number()
            Traceback (most recent call last):
            ...
            NotImplementedError: transitive database only computed up to degree 47
        """
        d = self.order()
        if d > 47:
            raise NotImplementedError("transitive database only computed up to degree 47")
        elif d == 32:
            # I don't know why this case is special, but you can check this in Magma (GAP only goes up to 22)
            return ZZ(33)
        else:
            return ZZ(1)

    def signature(self):
        r"""
        Return 1 if contained in the alternating group, -1 otherwise.

        EXAMPLES::

            sage: GF(3^2).galois_group().signature()
            -1
            sage: GF(3^3).galois_group().signature()
            1
        """
        return ZZ(1) if (self._field.degree() % 2) else ZZ(-1)

    @lazy_attribute
    def _gcdata(self):
        r"""
        Return the Galois closure (ie, the finite field itself) together with the identity

        EXAMPLES::

            sage: GF(3^2).galois_group()._gcdata
            (Finite Field in z2 of size 3^2,
             Identity endomorphism of Finite Field in z2 of size 3^2)
        """
        k = self._field
        return k, k.Hom(k).identity()

    def _element_constructor_(self, x, check=True):
        r"""
        Create an element of this Galois group from ``x``.

        INPUT:

        - ``x`` -- one of the following (`G` is this Galois group):

          - the integer 1, denoting the identity of `G`;

          - an element of `G`;

          - a list of length 1, giving the exponent of Frobenius

          - a permutation of the right length that defines an element of `G`,
            or anything that coerces into such a permutation;

          - an automorphism of the finite field.

        EXAMPLES::

            sage: k = GF(3^3)
            sage: G = k.galois_group()
            sage: G(1)
            1
            sage: G([2])
            Frob^2
            sage: G(G.gens()[0])
            Frob
            sage: G([(1,3,2)])
            Frob^2
            sage: G(k.hom(k.gen()^3, k))
            Frob
            sage: G(k.frobenius_endomorphism())
            Frob
        """
        if x == 1:
            return self.element_class(self, [0])

        k = self._field
        d = k.degree()
        n = None
        if isinstance(x, GaloisGroup_GFElement) and x.parent() is self:
            n = x.exponents()[0]
        elif isinstance(x, FiniteFieldHomomorphism_generic):
            if check and not (x.domain() is k and x.codomain() is k):
                raise ValueError("Not an automorphism of the correct finite field")
            a = k.gen()
            b = x(a)
            q = k.base_ring().cardinality()
            n = 0
            while n < d:
                if a == b:
                    break
                n += 1
                a = a**q
            else:
                raise RuntimeError("Automorphism was not a power of Frobenius")
        elif isinstance(x, FrobeniusEndomorphism_finite_field):
            if check and not x.domain() is k:
                raise ValueError("Not an automorphism of the correct finite field")
            n = x.power()
        elif isinstance(x, list) and len(x) == 1 and x[0] in ZZ:
            n = x[0]
        else:
            g = self.permutation_group()(x)
            n = g(1) - 1
        return self.element_class(self, [n])


