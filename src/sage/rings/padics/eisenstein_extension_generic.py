import sage.rings.polynomial.polynomial_quotient_ring as pqr
import sage.rings.padics.padic_generic
#import sage.rings.padics.eisenstein_extension_generic_element
import sage.rings.padics.padic_ring_generic
import sage.rings.padics.padic_field_generic
import padic_extension_generic

pAdicGeneric = sage.rings.padics.padic_generic.pAdicGeneric
#PolynomialQuotientRing_domain = sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_domain
#EisensteinExtensionGenericElement = sage.rings.padics.eisenstein_extension_generic_element.EisensteinExtensionGenericElement
#pAdicRingBaseGeneric = sage.rings.padics.padic_ring_generic.pAdicRingBaseGeneric
#pAdicFieldBaseGeneric = sage.rings.padics.padic_field_generic.pAdicFieldBaseGeneric
pAdicExtensionGeneric = padic_extension_generic.pAdicExtensionGeneric

class EisensteinExtensionGeneric(pAdicExtensionGeneric):
    def __init__(self, poly, prec, print_mode, names, element_class):
        #base = poly.base_ring()
        #if base.is_field():
        #    self._PQR = pqr.PolynomialQuotientRing_field(poly.parent(), poly, name = names)
        #else:
        #    self._PQR = pqr.PolynomialQuotientRing_domain(poly.parent(), poly, name = names)
        #if isinstance(names, tuple):
        #    names = names[0]
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
        #self._precompute()

    def _repr_(self, do_latex = False):
        return "Eisenstein Extension of %s in %s defined by %s"%(
            self.ground_ring(), self.variable_name(), self.modulus())

    def ramification_index(self, K = None):
        if K is None or K is self.ground_ring():
            return self.modulus().degree()
        else:
            raise NotImplementedError

    def inertia_degree(self, K = None):
        if K is None or K is self.ground_ring():
            return 1
        else:
            raise NotImplementedError

    def inertia_subring(self):
        raise NotImplementedError

    def extension(self):
        raise NotImplementedError

    ext = extension

    def residue_class_field(self):
        return self.ground_ring().residue_class_field()

    def discriminant(self, K=None):
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
        raise NotImplementedError

    def is_abelian(self):
        raise NotImplementedError

    def is_normal(self):
        raise NotImplementedError

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "only one generator"
        return self([0,1])

    def uniformizer_pow(self, n):
        ans = self(0)
        if n is not infinity:
            ans._set_uniformizer_pow(n)
        return ans

    def uniformizer(self):
        return self.gen()

    def _uniformizer_print(self):
        return self.variable_name()

    def has_pth_root(self):
        raise NotImplementedError

    def has_root_of_unity(self, n):
        raise NotImplementedError
