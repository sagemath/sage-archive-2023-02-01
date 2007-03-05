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
pAdicGeneric = sage.rings.padics.padic_ring_generic.pAdicRingGeneric
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

    def __init__(self, poly, prec = None, print_mode = None, check = True):
        if check:
            if not (isinstance(poly, Polynomial) and isinstance(poly.base_ring(), pAdicRingGeneric)):
                raise TypeError, "Defining must be a polynomial with p-adic base ring"
            if poly.leading_coefficient() != 1:
                raise ValueError, "Defining polynomial must be monic"
            if poly[0].valuation() > 0:
                raise ValueError, "Defining polynomial must be irreducible and unramified"
            F = PolynomialRing(f.base_ring().residue_class_field(), f.parent().variable_name())(f).factor()
            if len(F) != 1 or F[0][1] != 1:
                raise ValueError, "Defining polynomial must be irreducible over the residue field"
        R = poly.base_ring()
        if prec is None:
            prec = R.precision_cap()
        if print_mode is None:
            print_mode = R.get_print_mode()
        PolynomialQuotientRing_domain.__init__(self, poly.parent(), poly)
        pAdicRingGeneric.__init__(self, R.prime(), prec, print_mode)

    def __call__(self, x, prec = None):
        return UnramifiedRingExtensionElement(self, x, prec)
    #        r"""
    #        This is currently unimplemented, as understanding how to
    #        embed one ring in another when appropriate is nontrivial.
    #        """
    #        raise NotImplementedError, "not smart enough to embed rings"

    def __cmp__(self, other):
        if isinstance(other, UnramifiedRingExtension):
            if self.prime() < other.prime():
                return -1
            elif self.prime() > other.prime():
                return 1
            elif cmp(type(self.base_ring()), type(other.base_ring())) != 0:
                return cmp(type(self.base_ring()), type(other.base_ring()))
            elif self.base_ring().precision_cap() < other.base_ring().precision_cap():
                return -1
            elif self.base_ring().precision_cap() > other.base_ring().precision_cap():
                return 1
            else:
                return cmp(self.defining_polynomial(), other.defining_polynomial())
        else:
            return cmp(type(self), type(other))

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

    def ramification_index(self, K = None):
        if K is None:
            return 1
        else:
            raise NotImplementedError

    def inertia_degree(self, K = None):
        if K is None:
            return self.modulus().degree()
        else:
            raise NotImplementedError

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
        r"""
        Return the absolute discriminant of self over Zp.
        """
        return 1

    def discriminant(self, K=None):
        if K is None:
            return 1
        else:
            raise NotImplementedError

    def automorphisms(self):
        raise NotImplementedError

    def galois_group(self):
        r"""
        Returns the galois group of self's fraction field over Qp.
        """
        ##
        ## If K is a number field, then K.galois_group() can return
        ## other variants, i.e. via Pari or KASH. We could consider
        ## doing this.
        ##
        from sage.groups.perm_gps.permgroup import CyclicPermutationGroup
        return CyclicPermutationGroup(self.modulus().degree())

    def gen(self):
        return self.polynomial_ring().gen()

    def is_abelian(self):
        return True

    def is_normal(self):
        return True

    def uniformizer(self):
        return self(self.ground_ring().uniformizer())

    def _uniformizer_sym(self, do_latex = False):
        return self.ground_ring()._uniformizer_sym(do_latex)

    def has_pth_root(self):
        return (self.prime() == 2)

    def has_root_of_unity(self, n):
        ##
        ## I wouldn't mind if someone wanted to check to make sure
        ## I've got the right formula for the 2-adic case.
        ##
        if (self.prime() == 2):
            return n.divides(2*(self.residue_class_field().order()-1))
        else:
            return n.divides(self.residue_class_field().order() - 1)

    def hasGNB(self):
        raise NotImplementedError

    def hom(self, ring):
        raise NotImplementedError

    def random_element(self):
        return reduce(lambda x,y: x+y,map(lambda a,b:a*b,[self.base_ring().random_element() for _ in range(self.modulus().degree())],[self.gen()**i for i in range(self.modulus().degree())]),0)

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

    def set_print_mode(self, mode):
        return self.base_ring().set_print_mode(mode)

    def get_print_mode(self):
        return self.base_ring().get_print_mode()

    def fraction_field(self):
        r"""
        Returns the fraction field of this extension, which is just
        the extension of base.fraction_field() determined by the
        same polynomial.
        """
        raise NotImplementedError
        #return self.base_ring().fraction_field().extension(self._poly)
