

def is_MPolynomial(x):
    return isinstance(x, MPolynomial)

cdef class MPolynomial(CommutativeRingElement):

    ####################
    # Some standard conversions
    ####################
    def __int__(self):
        if self.degree() == 0:
            return int(self.constant_coefficient())
        else:
            raise TypeError

    def __long__(self):
        if self.degree() == 0:
            return long(self.constant_coefficient())
        else:
            raise TypeError

    def __float__(self):
        if self.degree() == 0:
            return float(self.constant_coefficient())
        else:
            raise TypeError

    def _mpfr_(self, R):
        if self.degree() == 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _complex_mpfr_field_(self, R):
        if self.degree() == 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _complex_double_(self, R):
        if self.degree() == 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _real_double_(self, R):
        if self.degree() == 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _rational_(self):
        if self.degree() == 0:
            from sage.rings.rational import Rational
            return Rational(repr(self))
        else:
            raise TypeError

    def _integer_(self):
        if self.degree() == 0:
            from sage.rings.integer import Integer
            return Integer(repr(self))
        else:
            raise TypeError
