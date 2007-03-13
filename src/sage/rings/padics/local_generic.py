"""
Superclass for p-adic and power series rings.
"""

#*****************************************************************************
#       Copyright (C) 2007 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.rings.ring
import sage.structure.parent_gens


class LocalGeneric(sage.rings.ring.CommutativeRing):
    def __init__(self, prec, names):
        self._prec = prec
        sage.structure.parent_gens.ParentWithGens.__init__(self, self, (names,), normalize=False)

    def __call__(self, x):
        raise NotImplementedError

    def __cmp__(self, other):
        raise NotImplementedError

    def __contains__(self, x):
        raise NotImplementedError

    def _coerce_impl(self, x):
        raise NotImplementedError

    def _repr_(self, do_latex = False, mode = None):
        return "Generic Local Ring"

    def _latex_(self):
        return self._repr_(do_latex = True)

    def characteristic(self):
        raise NotImplementedError

    def precision_cap(self):
        r"""
        Returns the precision cap for self.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            integer -- self's precision cap

        EXAMPLES:
            sage: R = Zp(3, 10,'fixed-mod'); R.precision_cap()
                10

        NOTES:
            This will have different meanings depending on the type of local ring.  For fixed modulus rings, all elements are considered modulo self.prime()^self.precision_cap().  For rings with an absolute cap (i.e. the class pAdicRingCappedAbsolute), each element has a precision that is tracked and is bounded above by self.precision_cap().  That element self.prime()^precision.  Rings with relative caps (i.e. the class pAdicRingCappedRelative) are the same except that the precision is the precision of the unit part of each element.  For lazy rings, this gives the initial precision to which elements are computed.
        """
        return self._prec

    def print_mode(self):
        raise NotImplementedError

    def is_atomic_repr(self):
        r"""
        Return False, since we want p-adics to be printed with parentheses around them
        when they are coefficients, e.g., in a polynomial.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            boolean -- whether self's representation is atomic, i.e., False

        """
        return False

    def is_exact(self):
        r"""
        Returns whether this p-adic ring is exact, i.e. False.

        INPUT:
            self -- a p-adic rinng

        OUTPUT:
            boolean -- whether self is exact, i.e. False.
        """
        return False

    def residue_characteristic(self):
        return self.residue_class_field().characteristic()

    def residue_class_field(self):
        raise NotImplementedError

    def defining_polynomial(self, var = 'x'):
        r"""
        Returns the defining polynomial of this local ring

        INPUT:
            self -- a local ring
            var -- string (default: 'x') the name of the variable

        OUTPUT:
            polynomial -- the defining polynomial of this ring as an extension over its ground ring
        """
        from sage.rings.polynomial_ring import PolynomialRing
        x = PolynomialRing(self, var).gen()
        return x - 1

    def ground_ring(self):
        r"""
        Returns self.

        Will be overridden by extensions.

        INPUT:
            self -- a local ring

        OUTPUT:
            the ground ring of self, i.e., itself
        """
        return self

    def ground_ring_of_tower(self):
        r"""
        Returns self.

        Well be overridden by extensions.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            the ground ring of the tower for self, i.e., itself
        """
        return self

    def degree(self):
        r"""
        Returns the degree of self over the ground ring, i.e. 1.

        INPUT:
            self -- a local ring

        OUTPUT:
            integer -- the degree of this ring, i.e., 1
        """
        return 1

    def ramification_index(self):
        r"""
        Returns the ramification index over the ground ring, i.e. 1.

        INPUT:
            self -- a local ring

        OUTPUT:
            integer -- the ramification index of this ring, i.e., 1
        """
        return 1

    def e(self):
        r"""
        Returns the ramification index over the ground ring, i.e. 1.

        INPUT:
            self -- a local ring

        OUTPUT:
            integer -- the ramification index of this ring, i.e., 1
        """
        return ramification_index(self)

    def inertia_degree(self):
        r"""
        Returns the inertia degree over the ground ring, i.e. 1.

        INPUT:
            self -- a local ring

        OUTPUT:
            integer -- the inertia degree of this ring, i.e., 1
        """
        return 1

    def f(self):
        r"""
        Returns the inertia degree over the ground ring, i.e. 1.

        INPUT:
            self -- a local ring

        OUTPUT:
            integer -- the inertia degree of this ring, i.e., 1
        """
        return inertia_degree(self)

    def inertia_subring(self):
        r"""
        Returns the inertia subring, i.e. self.

        INPUT:
            self -- a local ring

        OUTPUT:
            the inertia subring of self, i.e., itself
        """
        return self

    def get_extension(self):
        r"""
        Returns the trivial extension of self.
        """
        raise NotImplementedError

    def uniformizer(self):
        r"""
        Returns a uniformizer.
        """
        raise NotImplementedError

    def uniformiser(self):
        r"""
        Returns a uniformiser.
        """
        return self.uniformizer()

    def has_root_of_unity(self, n):
        raise NotImplementedError

    def hom(self, ring):
        r"""
        Returns the set of homomorphisms from $\Z_p$ into another ring.
        """
        raise NotImplementedError

    def is_isomorphic(self, ring):
        raise NotImplementedError

    def random_element(self):
        raise NotImplementedError

    def unit_group(self):
        raise NotImplementedError

    def unit_group_gens(self):
        raise NotImplementedError

    def principal_unit_group(self):
        raise NotImplementedError

    def zeta(self, n = None):
        raise NotImplementedError

    def zeta_order(self):
        raise NotImplementedError

    def krull_dimension(self):
        raise NotImplementedError

    def is_finite(self):
        r"""
        Returns whether this ring is finite, i.e. False.

        INPUT:
            self -- a p-adic ring

        OUTPUT:
            boolean -- whether self is finite, i.e., False

        EXAMPLES:
            sage: R = Zp(3, 10,'fixed-mod'); R.is_finite()
            False
        """
        return False
