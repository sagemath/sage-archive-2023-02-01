import padic_generic
import sage.rings.infinity
from sage.rings.padics.pow_computer import PowComputer

infinity = sage.rings.infinity.infinity

class pAdicBaseGeneric(padic_generic.pAdicGeneric):
    def __init__(self, p, prec, print_mode, names, element_class):
        self.prime_pow = PowComputer(p, max(min(prec - 1, 30), 1), prec, self.is_field())
        padic_generic.pAdicGeneric.__init__(self, self, p, prec, print_mode, names, element_class)

    def __reduce__(self):
        """
        For pickling.

        This function is provided because prime_pow needs to be set before _printer, so the standard unpickling fails.
        """
        from sage.rings.padics.factory import Zp, Qp
        if self.is_field():
            if self.is_capped_relative():
                return Qp, (self.prime(), self.precision_cap(), 'capped-rel', self.print_mode(), 40, self.variable_name())
            else:
                return Qp, (self.prime(), self.precision_cap(), 'lazy', self.print_mode(), self.halting_parameter(), self.variable_name())
        else:
            if self.is_capped_relative():
                return Zp, (self.prime(), self.precision_cap(), 'capped-rel', self.print_mode(), 40, self.variable_name())
            elif self.is_capped_absolute():
                return Zp, (self.prime(), self.precision_cap(), 'capped-abs', self.print_mode(), 40, self.variable_name())
            elif self.is_fixed_mod():
                return Zp, (self.prime(), self.precision_cap(), 'fixed-mod', self.print_mode(), 40, self.variable_name())
            else:
                return Zp, (self.prime(), self.precision_cap(), 'lazy', self.print_mode(), self.halting_parameter(), self.variable_name())

    def is_isomorphic(self, ring):
        r"""
        Returns whether self and ring are isomorphic, i.e. whether ring is an implementation of $\Z_p$ for the same prime as self.

        INPUT:
            self -- a p-adic ring
            ring -- a ring

        OUTPUT:
            boolean -- whether ring is an implementation of $\Z_p$ for the same prime as self.
        """
        return is_instance(ring, pAdicBaseGeneric) and self.prime() == ring.prime() and self.is_field() == ring.is_field()

    def uniformizer_pow(self, n):
        return self(self.prime_pow(n))

    def _uniformizer_print(self):
        return self.variable_name()
