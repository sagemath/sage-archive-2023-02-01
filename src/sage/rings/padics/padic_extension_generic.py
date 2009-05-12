import padic_generic
#import padic_extension_generic_element
import padic_base_generic
from sage.structure.element import Element

pAdicGeneric = padic_generic.pAdicGeneric
pAdicBaseGeneric = padic_base_generic.pAdicBaseGeneric
#pAdicExtensionGenericElement = padic_extension_generic_element.pAdicExtensionGenericElement

class pAdicExtensionGeneric(pAdicGeneric):
    def __init__(self, poly, prec, print_mode, names, element_class):
        #type checking done in factory
        self._given_poly = poly
        R = poly.base_ring()
        # We'll deal with the different names better later.
        # Using a tuple here is mostly needed for more general extensions
        # (ie not eisenstein or unramified)
        print_mode['unram_name'] = names[2]
        print_mode['ram_name'] = names[3]
        print_mode['var_name'] = names[0]
        names = names[0]
        pAdicGeneric.__init__(self, R, R.prime(), prec, print_mode, names, element_class)

    def __reduce__(self):
        """
        For pickling.

        This function is provided because prime_pow needs to be set before _printer, so the standard unpickling fails.
        """
        from sage.rings.padics.factory import ExtensionFactory
        return ExtensionFactory, (self.base_ring(), self._pre_poly, self.precision_cap(), self.print_mode(), None, self.variable_name())

    def __contains__(self, x):
        if isinstance(x, Element) and (x.parent() is self or x.parent().fraction_field() is self):
            return True
        if self.ground_ring().__contains__(x):
            return True
        #have not yet added support for more coercion.  See pAdicTodo on sage/home/padicgroup
        return False

    def __cmp__(self, other):
        """
        Returns 0 if self == other, and 1 or -1 otherwise.

        We consider two p-adic rings or fields to be equal if they are equal mathematically, and also have the same precision cap and printing parameters.

        EXAMPLES:
        sage: R.<a> = Qq(27)
        sage: S.<a> = Qq(27,print_mode='val-unit')
        sage: R == S
        False
        sage: S.<a> = Qq(27,type='capped-rel')
        sage: R == S
        True
        sage: R is S
        True
        """
        c = cmp(type(self), type(other))
        if c != 0:
            return c
        groundcmp = self.ground_ring().__cmp__(other.ground_ring())
        if groundcmp != 0:
            return groundcmp
        c = cmp(self.defining_polynomial(), other.defining_polynomial())
        if c != 0:
            return c
        c = cmp(self.precision_cap(), other.precision_cap())
        if c != 0:
            return c
        return self._printer.cmp_modes(other._printer)

    def gen(self, n=0):
        raise NotImplementedError

    def defining_polynomial(self):
        return self._given_poly

    def modulus(self):
        return self._given_poly

    def ground_ring(self):
        return self._given_poly.base_ring()

    def ground_ring_of_tower(self):
        if isinstance(self.ground_ring(), pAdicBaseGeneric):
            return self.ground_ring()
        else:
            return self.ground_ring().ground_ring_of_tower()

    def is_isomorphic(self, ring):
        raise NotImplementedError

    def polynomial_ring(self):
        return self._given_poly.parent()

    #def teichmuller(self, x, prec = None):
    #    if prec is None:
    #        prec = self.precision_cap()
    #    x = self(x, prec)
    #    if x.valuation() > 0:
    #        return self(0)
    #    q = self.residue_class_field().order()
    #    u = 1 / self(1 - q, prec)
    #    delta = u * (1 - x ** (q - 1))
    #    xnew = x - x*delta*(1 - q * delta)
    #    while x != xnew:
    #        x = xnew
    #        delta = u*(1-x**(q-1))
    #        xnew = x - x*delta*(1-q*delta)
    #    return x

    def absolute_discriminant(self):
        r"""
        Return the absolute discriminant of self over Zp.
        """
        raise NotImplementedError

    def fraction_field(self, print_mode=None):
        r"""
        Returns the fraction field of this extension, which is just
        the extension of base.fraction_field() determined by the
        same polynomial.

        INPUT:
        print_mode - a dictionary containing print options.  Defaults to the same options as this ring.
        OUTPUT:
        the fraction field of self.

        EXAMPLES:
        sage: U.<a> = Zq(17^4, 6, print_mode='val-unit', print_max_terse_terms=3)
        sage: U.fraction_field()
        Unramified Extension of 17-adic Field with capped relative precision 6 in a defined by (1 + O(17^6))*x^4 + (O(17^6))*x^3 + (7 + O(17^6))*x^2 + (10 + O(17^6))*x + (3 + O(17^6))
        sage: U.fraction_field({"pos":False}) == U.fraction_field()
        False
        """
        if self.is_field() and print_mode is None:
            return self
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'mode': print_mode}
        for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
            if not print_mode.has_key(option):
                print_mode[option] = self._printer.dict()[option]
        K = self.ground_ring().fraction_field(print_mode)
        if self.is_lazy():
            return K.extension(self.polynomial_ring().base_extend(K)(self.defining_polynomial()), prec = self.precision_cap(), print_mode = print_mode, halt = self.halting_parameter(), names = self.variable_name())
        else:
            return K.extension(self.polynomial_ring().base_extend(K)(self.defining_polynomial()), prec = self.precision_cap(), print_mode = print_mode, names = self.variable_name())

    def integer_ring(self, print_mode=None):
        r"""
        Returns the ring of integers of self, which is just the
        extension of base.integer_ring() determined by the same
        polynomial.
        """
        #Currently does not support fields with non integral defining polynomials.  This should change when the padic_general_extension framework gets worked out.
        if not self.is_field() and print_mode is None:
            return self
        if print_mode is None:
            print_mode = {}
        elif isinstance(print_mode, str):
            print_mode = {'mode': print_mode}
        for option in ['mode', 'pos', 'ram_name', 'unram_name', 'var_name', 'max_ram_terms', 'max_unram_terms', 'max_terse_terms', 'sep', 'alphabet']:
            if not print_mode.has_key(option):
                print_mode[option] = self._printer.dict()[option]
        K = self.ground_ring().integer_ring(print_mode)
        if self.is_lazy():
            return K.extension(self.polynomial_ring().change_ring(K)(self.defining_polynomial()), prec = self.precision_cap(), print_mode = self.print_mode(), halt = self.halting_parameter(), names = self.variable_name())
        else:
            return K.extension(self.polynomial_ring().change_ring(K)(self.defining_polynomial()), prec = self.precision_cap(), print_mode = self.print_mode(), names = self.variable_name())

    #def hasGNB(self):
    #    raise NotImplementedError

    def random_element(self):
        return reduce(lambda x,y: x+y,map(lambda a,b:a*b,[self.ground_ring().random_element() for _ in range(self.modulus().degree())],[self.gen()**i for i in range(self.modulus().degree())]),0)

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

