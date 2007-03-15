import sage.rings.polynomial_quotient_ring
import sage.rings.padics.padic_generic
import sage.rings.padics.eisenstein_extension_element
import sage.rings.padics.padic_ring_generic
import sage.rings.padics.padic_field_generic

pAdicGeneric = sage.rings.padics.padic_generic.pAdicGeneric
PolynomialQuotientRing_domain = sage.rings.polynomial_quotient_ring.PolynomialQuotientRing_domain
EisensteinExtensionElement = sage.rings.padics.eisenstein_extension_element.EisensteinExtensionElement
pAdicRingBaseGeneric = sage.rings.padics.padic_ring_generic.pAdicRingBaseGeneric
pAdicFieldBaseGeneric = sage.rings.padics.padic_field_generic.pAdicFieldBaseGeneric

class EisensteinExtensionGeneric(pAdicGeneric, PolynomialQuotientRing_domain):
    def __init__(self, poly, prec, print_mode, names):
        R = poly.base_ring()
        PolynomialQuotientRing_domain.__init__(self, poly.parent(), poly, name = names)
        if isinstance(names, tuple):
            names = names[0]
        print names
        print self.variable_name()
        pAdicGeneric.__init__(self, R.prime(), prec, print_mode, names)
        self._precompute()

    def _precompute(self):
        """
        Precomputes some commonly used information about the ring.

        self._translate_p_pi[m, n] is the coefficient of pi^n in the expansion of p^m if pi is the uniformizer of self, p is the uniformizer of self.base_ring()
        """
        #self._translate_p_pi = {}
        #cap = self.precision_cap()
        #e = self.e()
        #p = self.base_ring().uniformizer()
        ##first we initialize all entries to 0
        #zero = self.base_ring()(0)

        #v = [((m,n),zero) for m in range(cap//e) for n in range(cap)]
        #self._translate_p_pi = dict(v)

        #since p^0 = pi^0 = 1, the only entry for m = 0 is the following
        #self._translate_p_pi[0, 0] = self.base_ring()(1)
        #we now compute the first row.  ppolylist is the list associated to the polynomial obtained by solving the defining polynomial of self for p, using the fact that that polynomial is Eisenstein.
        ppolylist = (self.base_ring()(~(self.defining_polynomial()[0] / self.base_ring().uniformizer()))*(self.defining_polynomial()[0] - self.defining_polynomial())).list()
        #We make a duplicate.  shift1list satisfies p = shift1list[1]*pi + ... + shift1list[e - 1] * pi^(e-1) + pi^e
        self._shift1list = ppolylist

        #For our use of the shift1list in shifting, we want it over by 1, so we pop.
        self._shift1list.pop(0)
        #REMEMBER TO MOVE THIS BACK TO THE END
        return

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



    def __contains__(self, x):
        if isinstance(x, EisensteinExtensionElement) and (x.parent() is self or x.parent().fraction_field() is self):
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
        return "Eisenstein Extension of %s in %s defined by %s"%(
            self.base_ring(), self.variable_name(), self.modulus())

    def ngens(self):
        return 1

    def defining_polynomial(self):
        return self.modulus()

    def ground_ring(self):
        return self.base_ring()

    def ground_ring_of_tower(self):
        if isinstance(self.base_ring(), (pAdicRingBaseGeneric, pAdicFieldBaseGeneric)):
            return self.base_ring()
        else:
            return self.base_ring().ground_ring_of_tower()

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

    def is_isomorphic(self, ring):
        raise NotImplementedError

    def residue_class_field(self):
        return self.base_ring().residue_class_field()

    def residue_system(self):
        return [self(i) for i in self.residue_class_field()]

    def teichmuller(self, x, prec = None):
        if prec is None:
            prec = self.precision_cap()
        q = self.residue_class_field.order()
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
        raise NotImplementedError

    def discriminant(self, K=None):
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
        raise NotImplementedError

    def is_abelian(self):
        return True

    def is_normal(self):
        return True

    def uniformizer(self):
        return self.gen()

    def has_pth_root(self):
        raise NotImplementedError

    def has_root_of_unity(self, n):
        raise NotImplementedError

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

    def _uniformizer_print(self):
        return self.variable_name()

