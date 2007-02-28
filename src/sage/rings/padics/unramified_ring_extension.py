from sage.rings.padics.padic_ring_generic import pAdicRingGeneric
from sage.rings.padics.padic_ring_generic import pAdicRingBaseGeneric
from sage.rings.polynomial_quotient_ring import *
from sage.rings.integer import Integer
from sage.rings.finite_field import GF
from sage.rings.polynomial_ring import PolynomialRing
from sage.rings.padics.unramified_ring_extension_element import UnramifiedRingExtensionElement
import sage.rings.polynomial_quotient_ring_element
import sage.rings.padics.padic_ring_generic
import sage.rings.polynomial_element

Polynomial = sage.rings.polynomial_element.Polynomial
pAdicRingBaseGeneric = sage.rings.padics.padic_ring_generic.pAdicRingBaseGeneric
PQRElement = sage.rings.polynomial_quotient_ring_element.PolynomialQuotientRingElement

class UnramifiedRingExtension(pAdicRingGeneric, PolynomialQuotientRing_domain):
    r"""
    Base class for extensions of Z_p.

    As of right now, implemented in the simplest way possible:
    we just take a polynomial and a ring, and return the polynomial
    quotient ring, after checking that the defining polynomial is
    irreducible.

    TODO:
    We should use the Round4 algorithm to break the given
    polynomial into a totally ramified and unramified part, and
    create two subextensions, so that the whole extension is
    given as a tower.
    """

    def __init__(self, poly, prec = None, print_mode = None):
        #need to check stuff about poly: irreducible, monic, unramified, valid coefficient ring
        R = poly.base_ring()
        if prec is None:
            prec = R.precision_cap()
        if print_mode is None:
            print_mode = R.get_print_mode()
        pAdicRingGeneric.__init__(self, R.prime(), prec, print_mode)
        PolynomialQuotientRing_domain.__init__(self, poly.parent(), poly)

    def __call__(self, x, prec = None):
        return UnramifiedRingExtensionElement(self, x, prec)
    #        r"""
    #        This is currently unimplemented, as understanding how to
    #        embed one ring in another when appropriate is nontrivial.
    #        """
    #        raise NotImplementedError, "not smart enough to embed rings"

    def __cmp__(self, other):
        raise NotImplementedError

    def __contains__(self, x):
        if isinstance(x, PQRElement) and x.parent() is self:
            return True
        if isinstance(x, (int, long, Integer)):
            return True
        if isinstance(x, pAdicRingBaseGeneric) and x.parent().prime() == self.prime() and (not isinstance(x, pAdicRingFixedModElement) or isinstance(self.ground_ring_of_tower(), pAdicRingFixedMod)):
            return True
        if isinstance(x, Polynomial) and x in self.polynomial_ring():
            return True
        #have not yet added support for coercion from other unramified extension rings
        return False

    def _coerce_impl(self, x):
        if self.__contains__(x):
            return self.__call__(x)
        else:
            raise TypeError, "cannot coerce %s into %s"%(x, self)

    def _repr_(self, do_latex = False):
        return "Unramified Extension of %s in %s defined by %s"%(
            self.base_ring(), self.variable_name(), self.modulus())

    def ngens(self):
        return 1

    def gen(self, n = 0):
        if n == 0:
            try:
                return self._gen
            except AttributeError:
                self._gen = UnramifiedRingExtensionElement(self, self.polynomial_ring().gen(), check = False, construct = True)
                return self._gen


    def defining_polynomial(self):
        return self.modulus()

    def ground_ring(self):
        return self.base_ring()

    def ground_ring_of_tower(self):
        if isinstance(self.base_ring(), pAdicRingBaseGeneric):
            return self.base_ring()
        else:
            return self.base_ring().ground_ring_of_tower()

    def ramification_index(self):
        return 1

    def inertia_degree(self):
        return self.modulus().degree()

    def inertia_subring(self):
        raise NotImplementedError

    def get_extension(self):
        raise NotImplementedError

    def is_isomorphic(self, ring):
        raise NotImplementedError

    def residue_class_field(self):
        return GF(self.prime_pow(self.modulus().degree()), name = self.polynomial_ring().variable_name(), modulus = PolynomialRing(self.base_ring().residue_class_field(), self.polynomial_ring().variable_name())(self.modulus()))

    def residue_system(self):
        return [self(i) for i in self.residue_class_field()]

    def teichmuller(self, x, prec = None):
        if prec is None:
            prec = self.precision_cap()
        q = self.prime_pow(self.modulus().degree())
        x = self(x, prec)
        xnew = x**q
        while x != xnew:
            x = xnew
            xnew = x**q
        return x

    def teichmuller_system(self):
        return [self.teichmuller(i) for i in self.residue_class_field() if i != 0]

    def absolute_discriminant(self):
        raise NotImplementedError

    def discriminant(self, K=None):
        raise NotImplementedError

    def automorphisms(self):
        raise NotImplementedError

    def galois_group(self):
        raise NotImplementedError

    def is_abelian(self):
        raise NotImplementedError

    def is_normal(self):
        return True

    def uniformizer(self):
        return self(self.prime())

    def has_pth_root(self):
        raise NotImplementedError

    def has_root_of_unity(self, n):
        raise NotImplementedError

    def hasGNB(self):
        raise NotImplementedError

    def hom(self, ring):
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

    def fraction_field(self):
        r"""
        Returns the fraction field of this extension, which is just
        the extension of base.fraction_field() determined by the
        same polynomial.
        """

        return self.base_ring().fraction_field().extension(self._poly)
