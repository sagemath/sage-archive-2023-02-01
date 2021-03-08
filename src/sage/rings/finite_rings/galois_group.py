r"""
Galois groups of Finite Fields
"""

from sage.groups.galois_group import GaloisGroup as GaloisGroup_base
from sage.rings.integer_ring import ZZ
from sage.misc.lazy_attribute import lazy_attribute

class GaloisGroup_GF(GaloisGroup_base):
    r"""
    The Galois group of a finite field.
    """
    def __init__(self, finite_field):
        r"""
        Create a Galois group.

        TESTS::

            sage: TestSuite(GF(9).galois_group()).run()
        """
        super().__init__(finite_field, algorithm=None, names=None, gc_numbering=False)

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

    @lazy_attribute
    def _gens(self):
        r"""
        Computes the generators
        """