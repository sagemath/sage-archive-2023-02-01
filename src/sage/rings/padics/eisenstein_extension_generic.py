import sage.rings.polynomial_quotient_ring as pqr
import sage.rings.padics.padic_generic
#import sage.rings.padics.eisenstein_extension_generic_element
import sage.rings.padics.padic_ring_generic
import sage.rings.padics.padic_field_generic
import padic_extension_generic

pAdicGeneric = sage.rings.padics.padic_generic.pAdicGeneric
#PolynomialQuotientRing_domain = sage.rings.polynomial_quotient_ring.PolynomialQuotientRing_domain
#EisensteinExtensionGenericElement = sage.rings.padics.eisenstein_extension_generic_element.EisensteinExtensionGenericElement
#pAdicRingBaseGeneric = sage.rings.padics.padic_ring_generic.pAdicRingBaseGeneric
#pAdicFieldBaseGeneric = sage.rings.padics.padic_field_generic.pAdicFieldBaseGeneric
pAdicExtensionGeneric = padic_extension_generic.pAdicExtensionGeneric

class EisensteinExtensionGeneric(pAdicExtensionGeneric):
    def __init__(self, poly, prec, print_mode, names, element_class):
        base = poly.base_ring()
        if base.is_field():
            self._PQR = pqr.PolynomialQuotientRing_field(poly.parent(), poly, name = names)
        else:
            self._PQR = pqr.PolynomialQuotientRing_domain(poly.parent(), poly, name = names)
        #if isinstance(names, tuple):
        #    names = names[0]
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
        self._precompute()

    def _precompute(self):
        """
        Precomputes some commonly used information about the ring.

        self._translate_p_pi[m, n] is the coefficient of pi^n in the expansion of p^m if pi is the uniformizer of self, p is the uniformizer of self.ground_ring()
        """
        #self._translate_p_pi = {}
        #cap = self.precision_cap()
        #e = self.e()
        #p = self.ground_ring().uniformizer()
        ##first we initialize all entries to 0
        #zero = self.ground_ring()(0)

        #v = [((m,n),zero) for m in range(cap//e) for n in range(cap)]
        #self._translate_p_pi = dict(v)

        #since p^0 = pi^0 = 1, the only entry for m = 0 is the following
        #self._translate_p_pi[0, 0] = self.ground_ring()(1)
        #we now compute the first row.  ppolylist is the list associated to the polynomial obtained by solving the defining polynomial of self for p, using the fact that that polynomial is Eisenstein.
        ppolylist = (self.ground_ring()(~(self.defining_polynomial()[0] / self.ground_ring().uniformizer()))*(self.defining_polynomial()[0] - self.defining_polynomial())).list()
        #We make a duplicate.  shift1list satisfies p = shift1list[1]*pi + ... + shift1list[e - 1] * pi^(e-1) + pi^e
        self._shift1list = ppolylist

        #For our use of the shift1list in shifting, we want it over by 1, so we pop.
        self._shift1list.pop(0)
        #REMEMBER TO MOVE THIS BACK TO THE END
        return

        #We add zeros to ppolylist because we're going to try to make all of the entries lifts of finite field elements.
        ppolylist.extend([self.ground_ring()(0)]*(cap - len(ppolylist)))
        #We now modify ppolylist, propogating changes so that p is always equal to the analogous sum that shift1list satisfies above.
        for k in range(1, cap):
            rsh = ppolylist[k] << 1
            ppolylist[k] = ppolylist[k] % p
            for i in range(1, min(e, cap - k)):
                ppolylist[k + i] += rsh * self._shift1list[i]
        #ppolylist should now satisfy p = ppolylist[1]*pi + ... + ppolylist[cap - 1] * pi^(cap - 1)
        #Now that we have it, we create a polynomial (mod pi^cap).  We will be raising this poly to powers to fill in _translate_p_pi.
        ppoly = self.polynomial_ring().quotient(self.polynomial_ring().gen()**cap, self.polynomial_ring().variable_name())(self.polynomial_ring()(ppolylist))
        # make a copy to raise to powers
        ppowpoly = ppoly
        # fill in _translate_p_pi
        for m in range(1, cap // e):
            for n in range(e * m, cap):
                self._translate_p_pi[m, n] = ppowpoly[n]
            ppowpoly *= ppoly

    def _repr_(self, do_latex = False):
        return "Eisenstein Extension of %s in %s defined by %s"%(
            self.ground_ring(), self.variable_name(), self.modulus())

    def ramification_index(self, K = None):
        if K is None:
            return self.modulus().degree()
        else:
            raise NotImplementedError

    def inertia_degree(self, K = None):
        if K is None:
            return 1
        else:
            raise NotImplementedError

    def inertia_subring(self):
        raise NotImplementedError

    def get_extension(self):
        raise NotImplementedError

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

    def uniformizer(self):
        return self.gen()

    def _uniformizer_print(self):
        return self.variable_name()

    def has_pth_root(self):
        raise NotImplementedError

    def has_root_of_unity(self, n):
        raise NotImplementedError
