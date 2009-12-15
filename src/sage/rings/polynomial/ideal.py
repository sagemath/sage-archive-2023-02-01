

from sage.rings.ideal import Ideal_pid

class Ideal_1poly_field(Ideal_pid):
    """
    An ideal in a univariate polynomial ring over a field.
    """
    def residue_class_degree(self):
        """
        Returns the degree of the generator of this ideal.

        This function is included for compatibility with ideals in rings of integers of number fields.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: P = R.ideal(t^4 + t + 1)
            sage: P.residue_class_degree()
            4
        """
        return self.gen().degree()
