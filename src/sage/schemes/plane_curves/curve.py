from sage.misc.all import latex

from sage.rings.all import is_Field

from sage.schemes.generic.algebraic_scheme import (
    AlgebraicScheme_subscheme, AlgebraicScheme_subscheme_projective)

from sage.schemes.generic.divisor_group import DivisorGroup

from sage.schemes.generic.divisor import Divisor_curve

class Curve_generic(AlgebraicScheme_subscheme):
    def _repr_(self):
        return "%s Curve over %s defined by %s"%(
            self._repr_type(), self.base_ring(), self.defining_polynomial())

    def _repr_type(self):
        return "Generic"

    def _latex_(self):
        """
        EXAMPLES:
            sage: x,y,z = PolynomialRing(QQ, 3, names='x,y,z').gens()
            sage: C = Curve(y^2*z - x^3 - 17*x*z^2 + y*z^2)
            sage: latex(C)
            - x^{3} + y^{2} z - 17 x z^{2} + y z^{2}

            sage: A2 = AffineSpace(2, QQ, names=['x','y'])
            sage: x, y = A2.coordinate_ring().gens()
            sage: C = Curve(y^2 - x^3 - 17*x + y)
            sage: latex(C)
            - x^{3} + y^{2} - 17 x + y
        """
        return latex(self.defining_polynomial())

    def defining_polynomial(self):
        return self.defining_polynomials()[0]

    def divisor_group(self, K=None):
        """
        EXAMPLES:
            sage: x,y,z = PolynomialRing(QQ, 3, names='x,y,z').gens()
            sage: C  = Curve(y^2*z - x^3 - 17*x*z^2 + y*z^2)
            sage: Cp = Curve(y^2*z - x^3 - 17*x*z^2 + y*z^2)
            sage: C.divisor_group() is Cp.divisor_group()
            True
        """
        if K is None:
            K = self.base_ring()
        elif not is_Field(K):
            raise TypeError, "Argument K (=%s) must be a field"%K
        # TODO: check that there exists a canonical map self.base_ring -> K
        # TODO: allow a morphism of rings
        return DivisorGroup(self, K)

    def divisor(self, v, check=True, reduce=True):
        return Divisor_curve(v, check=check, reduce=reduce, parent=self.divisor_group())

    def genus(self):
        """
        The geometric genus of the curve.
        """
        return self.geometric_genus()

    def geometric_genus(self):
        r"""
        The geometric genus of the curve, which is by definition the
        genus of the normalization of the projective closure of the
        curve over the algebraic closure of the base field; the base
        field must be a prime field.

        \note{Calls Singular's genus command.}

        EXAMPLE:
        Examples of projective curves.
            sage: P2 = ProjectiveSpace(2, GF(5), names=['x','y','z'])
            sage: x, y, z = P2.coordinate_ring().gens()
            sage: C = Curve(y^2*z - x^3 - 17*x*z^2 + y*z^2)
            sage: C.geometric_genus()
                  1
            sage: C = Curve(y^2*z - x^3)
            sage: C.geometric_genus()
                  0
            sage: C = Curve(x^10 + y^7*z^3 + z^10)
            sage: C.geometric_genus()
                  3

        Examples of affine curves.
            sage: x, y = PolynomialRing(GF(5), 2, 'xy').gens()
            sage: C = Curve(y^2 - x^3 - 17*x + y)
            sage: C.geometric_genus()
            1
            sage: C = Curve(y^2 - x^3)
            sage: C.geometric_genus()
            0
            sage: C = Curve(x^10 + y^7 + 1)
            sage: C.geometric_genus()
            3

        """
        try:
            return self.__genus
        except AttributeError:
            self.__genus = self.defining_ideal().genus()
            return self.__genus

    def union(self, other):
        from constructor import Curve
        return Curve(AlgebraicScheme_subscheme.union(self, other))

class Curve_generic_projective(Curve_generic, AlgebraicScheme_subscheme_projective):
    pass
