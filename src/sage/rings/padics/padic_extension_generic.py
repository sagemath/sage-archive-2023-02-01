import padic_generic
import padic_extension_generic_element
import padic_base_generic

pAdicGeneric = padic_generic.pAdicGeneric
pAdicBaseGeneric = padic_base_generic.pAdicBaseGeneric
pAdicExtensionGenericElement = padic_extension_generic_element.pAdicExtensionGenericElement

class pAdicExtensionGeneric(pAdicGeneric):
    def __init__(self, poly, prec, print_mode, names, element_class):
        #type checking done in factory
        self._given_poly = poly
        R = poly.base_ring()
        pAdicGeneric.__init__(self, R.prime(), prec, print_mode, names, element_class)

    def __contains__(self, x):
        if isinstance(x, pAdicExtensionGenericElement) and (x.parent() is self or x.parent().fraction_field() is self):
            return True
        if self.ground_ring().__contains__(x):
            return True
        #have not yet added support for more coercion.  See pAdicTodo on sage/home/padicgroup
        return False

    def __cmp__(self, other):
        if isinstance(other, type(self)):
            groundcmp = self.ground_ring().__cmp__(other.ground_ring())
            if groundcmp == 0:
                return cmp(self.defining_polynomial(), other.defining_polynomial())
            else:
                return groundcmp
        else:
            return cmp(type(self), type(other))

    def gen(self, n=0):
        if n != 0:
            raise IndexError, "extension has only one generator"
        #this won't work for padic_general extensions
        return self(self._PQR.polynomial_ring().gen())

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

    def teichmuller(self, x, prec = None):
        if prec is None:
            prec = self.precision_cap()
        x = self(x, prec)
        if x.residue(1) == 0:
            return self(0)
        q = self.residue_class_field().order()
        u = 1 / self(1 - q, prec)
        delta = u * (1 - x ** (q - 1))
        xnew = x - x*delta*(1 - q * delta)
        while x != xnew:
            x = xnew
            delta = u*(1-x**(q-1))
            xnew = x - x*delta*(1-q*delta)
        return x

    def absolute_discriminant(self):
        r"""
        Return the absolute discriminant of self over Zp.
        """
        raise NotImplementedError

    def fraction_field(self):
        r"""
        Returns the fraction field of this extension, which is just
        the extension of base.fraction_field() determined by the
        same polynomial.
        """
        if self.is_field():
            return self
        K = self.ground_ring().fraction_field()
	if self.is_lazy():
            return K.extension(self.polynomial_ring().base_extend(K)(self.defining_polynomial()), prec = self.precision_cap(), print_mode = self.print_mode(), halt = self.halting_parameter(), names = self.variable_name())
        else:
            return K.extension(self.polynomial_ring().base_extend(K)(self.defining_polynomial()), prec = self.precision_cap(), print_mode = self.print_mode(), names = self.variable_name())

    def integer_ring(self):
        r"""
        Returns the ring of integers of self, which is just the
        extension of base.integer_ring() determined by the same
        polynomial.
        """
        #Currently does not support fields with non integral defining polynomials.  This should change when the padic_general_extension framework gets worked out.
        if not self.is_field():
            return self
        K = self.ground_ring().integer_ring()
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

