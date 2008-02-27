import padic_ring_base_generic
import sage.rings.finite_field
import sage.rings.polynomial.polynomial_ring
import sage.rings.polynomial.polynomial_quotient_ring as pqr
#import sage.rings.padics.unramified_extension_generic_element
#import sage.rings.polynomial.polynomial_quotient_ring_element
#import sage.rings.padics.padic_ring_generic
#import sage.rings.polynomial.polynomial_element
#import sage.rings.padics.unramified_extension_generic
#import sage.rings.padics.padic_ring_extension_generic
#import sage.rings.padics.padic_ring_fixed_mod
#import sage.rings.padics.unramified_ring_extension_element
import sage.rings.padics.padic_generic
import padic_extension_generic

#pAdicRingBaseGeneric = padic_ring_base_generic.pAdicRingBaseGeneric
pAdicExtensionGeneric = padic_extension_generic.pAdicExtensionGeneric
#UnramifiedExtensionGenericElement = sage.rings.padics.unramified_extension_generic_element.UnramifiedExtensionGenericElement
PolynomialRing = sage.rings.polynomial.polynomial_ring.PolynomialRing
GF = sage.rings.ring.FiniteField
#pAdicGeneric = sage.rings.padics.padic_generic.pAdicGeneric

class UnramifiedExtensionGeneric(pAdicExtensionGeneric):
    """
    An unramified extension of Qp or Zp.
    """
    def __init__(self, poly, prec, print_mode, names, element_class):
        base = poly.base_ring()
        #if base.is_field():
        #    self._PQR = pqr.PolynomialQuotientRing_field(poly.parent(), poly, name = names)
        #else:
        #    self._PQR = pqr.PolynomialQuotientRing_domain(poly.parent(), poly, name = names)
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)

    def _repr_(self, do_latex = False):
        return "Unramified Extension of %s in %s defined by %s"%(
            self.ground_ring(), self.variable_name(), self.modulus())

    def ramification_index(self, K = None):
        if K is None:
            return 1
        elif K is self:
            return 1
        else:
            raise NotImplementedError

    e = ramification_index

    def inertia_degree(self, K = None):
        if K is None:
            return self.modulus().degree()
        elif K is self:
            return 1
        else:
            raise NotImplementedError

    residue_class_degree = inertia_degree

    f = inertia_degree

    def inertia_subring(self):
        return self

    maximal_unramified_subextension = inertia_subring

    def extension(self, *args, **kwds):
        raise NotImplementedError

    ext = extension

    def get_extension(self):
        raise NotImplementedError

    def residue_class_field(self):
        #should eventually take advantage of finite field \code{extension} or finite field \code{unramified_extension_of_degree} over the automatic coercion base.
        return GF(self.prime_pow(self.modulus().degree()), name = self.variable_name(), modulus = PolynomialRing(self.ground_ring().residue_class_field(), self.polynomial_ring().variable_name())(self.modulus()))

    residue_field = residue_class_field

    def discriminant(self, K=None):
        if K is None or K is self:
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

    def is_abelian(self):
        return True

    def is_normal(self):
        return True

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "only one generator"
        return self([0,1])

    def uniformizer_pow(self, n):
        if n is infinity:
            return self(0)
        return self(self.prime_pow.pow_Integer(n))

    uniformiser_pow = uniformizer_pow

    def uniformizer(self):
        return self(self.ground_ring().uniformizer())

    uniformiser = uniformizer

    def _uniformizer_print(self):
        return self.ground_ring()._uniformizer_print()

    def _unram_print(self):
        return self.variable_name()

    def has_pth_root(self):
        return self.ground_ring().has_pth_root()

    def has_root_of_unity(self, n):
        if (self.prime() == 2):
            return n.divides(2*(self.residue_class_field().order()-1))
        else:
            return n.divides(self.residue_class_field().order() - 1)
