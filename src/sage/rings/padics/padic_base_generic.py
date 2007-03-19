import padic_generic

class pAdicBaseGeneric(padic_generic.pAdicGeneric):
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

