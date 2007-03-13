import sage.rings.polynomial_element_generic
import sage.rings.padics.padic_ring
import sage.rings.padics.padic_ring_element
import sage.rings.padics.padic_field
import sage.rings.padics.padic_field_element

class pAdicPolynomial(sage.rings.polynomial_element_generic.Polynomial_generic_dense):
    def _mul_(self, right):
        """
        Override the default Karatsuba because it leaks precision
        """
        if right == 0 or self == 0:
            return self.polynomial(0)
        return self._mul_generic(right)

    def factor(self):
        """
        Factors self.

        Algorithm -- Round4
        """

def _round4_(f):
    pass
#Can assume f is square-free, monic


def _hensel_lift(f, g, h):
    pass

def _hensel_internal(f, g, h, s, t, p, n): #f = gh (mod p^n), sg + ht = 1 (
    pass
