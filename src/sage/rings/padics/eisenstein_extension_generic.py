import sage.rings.polynomial_quotient_ring
import sage.rings.padics.padic_generic

pAdicGeneric = sage.rings.padics.padic_generic.pAdicGeneric
PolynomialQuotientRing_domain = sage.rings.polynomial_quotient_ring.PolynomialQuotientRing_domain

class EisensteinExtensionGeneric(pAdicGeneric, PolynomialQuotientRing_domain):
    def __init__(self, poly, prec, print_mode, names):
        R = poly.base_ring()
        PolynomialQuotientRing_domain.__init__(self, poly.parent(), poly, name = names)
        pAdicGeneric.__init__(self, R.prime(), prec, print_mode, names)
        self._precompute()

    def _precompute(self):
        """
        Precomputes some commonly used information about the ring.

        self._translate_p_pi[m, n] is the coefficient of pi^n in the expansion of p^m if pi is the uniformizer of self, p is the uniformizer of self.base_ring()
        """
        self._translate_p_pi = {}
        cap = self.precision_cap()
        e = self.e()
        p = self.base_ring().uniformizer()
        #first we initialize all entries to 0
        zero = self.base_ring()(0)

        v = [((m,n),zero) for m in range(cap//e) for n in range(cap)]
        self._translate_p_pi = dict(v)

        #for m in range(cap // e):
        #    for n in range(cap):
        #        self._translate_p_pi[m, n] = zero

        #since p^0 = pi^0 = 1, the only entry for m = 0 is the following
        self._translate_p_pi[0, 0] = self.base_ring()(1)
        #we now compute the first row.  ppolylist is the list associated to the polynomial obtained by solving the defining polynomial of self for p, using the fact that that polynomial is Eisenstein.
        ppolylist = (self.base_ring()(~(self.defining_polynomial()[0] / self.base_ring().uniformizer()))*(self.defining_polynomial()[0] - self.defining_polynomial())).list()
        #We make a duplicate.  shift1list satisfies p = shift1list[1]*pi + ... + shift1list[e - 1] * pi^(e-1) + pi^e
        self._shift1list = ppolylist
        #We add zeros to ppolylist because we're going to try to make all of the entries lifts of finite field elements.
        ppolylist.extend([self.base_ring()(0)]*(cap - len(ppolylist)))
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

        #For our use of the shift1list in shifting, we want it over by 1, so we pop.
        self._shift1list.pop(0)


