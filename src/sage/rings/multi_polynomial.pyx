

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
            from rational import Rational
            return Rational(repr(self))
        else:
            raise TypeError

    def _integer_(self):
        if self.degree() == 0:
            from integer import Integer
            return Integer(repr(self))
        else:
            raise TypeError

    def subs(self,in_dict=None,**kwds):
        """
        This is an alias for self.fix().

        INPUT:
            in_dict -- (optional) dictionary of inputs
            **kwds  -- named parameters

        OUTPUT:
            new MPolynomial

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f((5,y))
            25*y^2 + y + 30
            sage: f.subs({x:5})
            25*y^2 + y + 30
            sage: f.subs(x=5)
            25*y^2 + y + 30
        """
        return self.fix(in_dict,**kwds)

    def substitute(self,in_dict=None,**kwds):
        """
        This is an alias for self.fix().

        INPUT:
            in_dict -- (optional) dictionary of inputs
            **kwds  -- named parameters

        OUTPUT:
            new MPolynomial

        EXAMPLES:
            sage: x, y = MPolynomialRing(ZZ,2,'xy').gens()
            sage: f = x^2 + y + x^2*y^2 + 5
            sage: f((5,y))
            25*y^2 + y + 30
            sage: f.substitute({x:5})
            25*y^2 + y + 30
            sage: f.substitute(x=5)
            25*y^2 + y + 30
         """
        return self.fix(in_dict,**kwds)
