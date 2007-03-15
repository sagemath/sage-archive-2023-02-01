from sage.rings.padics.padic_ring_generic import pAdicRingBaseGeneric
import sage.rings.finite_field
import sage.rings.polynomial_ring
import sage.rings.padics.unramified_extension_element
#import sage.rings.polynomial_quotient_ring_element
#import sage.rings.padics.padic_ring_generic
#import sage.rings.polynomial_element
#import sage.rings.padics.unramified_extension_generic
#import sage.rings.padics.padic_ring_extension_generic
#import sage.rings.padics.padic_ring_fixed_mod
#import sage.rings.padics.unramified_ring_extension_element
import sage.rings.polynomial_quotient_ring
import sage.rings.padics.padic_generic

UnramifiedExtensionElement = sage.rings.padics.unramified_extension_element.UnramifiedExtensionElement
PolynomialRing = sage.rings.polynomial_ring.PolynomialRing
GF = sage.rings.finite_field.GF
pAdicGeneric = sage.rings.padics.padic_generic.pAdicGeneric
PolynomialQuotientRing_domain = sage.rings.polynomial_quotient_ring.PolynomialQuotientRing_domain

class UnramifiedExtensionGeneric(pAdicGeneric, PolynomialQuotientRing_domain):
    def __init__(self, poly, names, prec, print_mode):
        R = poly.base_ring()
        PolynomialQuotientRing_domain.__init__(self, poly.parent(), poly, name = names)
        pAdicGeneric.__init__(self, R.prime(), prec, print_mode, names)

    def __contains__(self, x):
        if isinstance(x, UnramifiedExtensionElement) and (x.parent() is self or x.parent().fraction_field() is self):
            return True
        if self.base_ring().__contains__(x):
            return True
        #have not yet added support for coercion from appropriate unramified subextensions of the base ring
        return False

    def _coerce_impl(self, x):
        if self.__contains__(x):
            return self.__call__(x)
        else:
            raise TypeError, "cannot coerce %s into %s"%(x, self)

    def __cmp__(self, other):
        if isinstance(other, type(self)):
            basecmp = self.base_ring().__cmp__(other.base_ring())
            if basecmp == 0:
                return cmp(self.defining_polynomial(), other.defining_polynomial())
            else:
                return basecmp
        else:
            return cmp(type(self), type(other))

    def _repr_(self, do_latex = False):
        return "Unramified Extension of %s in %s defined by %s"%(
            self.base_ring(), self.variable_name(), self.modulus())

    def ngens(self):
        return 1

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
        #should eventually take advantage of finite field \code{extension} or finite field \code{unramified_extension_of_degree} over the automatic coercion base.
        return GF(self.prime_pow(self.modulus().degree()), name = self.polynomial_ring().variable_name(), modulus = PolynomialRing(self.base_ring().residue_class_field(), self.polynomial_ring().variable_name())(self.modulus()))

    def residue_system(self):
        return [self(i) for i in self.residue_class_field()]

    def teichmuller(self, x, prec = None):
        if prec is None:
            prec = self.precision_cap()
        q = self.residue_class_field().order()
        x = self(x, prec)
        xnew = x**q
        while x != xnew:
            x = xnew
            xnew = x**q
        return x

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


    def fraction_field(self):
        r"""
        Returns the fraction field of this extension, which is just
        the extension of base.fraction_field() determined by the
        same polynomial.
        """
        if self.is_field():
            return self
        K = self.base_ring().fraction_field()
        return K.extension(self.polynomial_ring().base_extend(K)(self.defining_polynomial()), names = self.variable_name())


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

    def uniformizer(self):
        return self(self.ground_ring().uniformizer())

    def _uniformizer_print(self):
        return self.base_ring()._uniformizer_print()

    def has_pth_root(self):
        return self.base_ring().has_pth_root()

    def has_root_of_unity(self, n):
        ##
        ## I wouldn't mind if someone wanted to check to make sure
        ## I've got the right formula for the 2-adic case.
        ##
        ## Also, this is wrong when the base is not Zp
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

