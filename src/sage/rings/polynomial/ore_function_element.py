from sage.structure.richcmp import richcmp, op_EQ, op_NE
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex

from sage.categories.map import Map
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from sage.structure.element import AlgebraElement


class OreFunction(AlgebraElement):
    def __init__(self, parent, numerator, denominator=None, simplify=True):
        AlgebraElement.__init__(self, parent)
        ring = parent.ring()
        numerator = ring(numerator)
        if denominator is None:
            denominator = ring.one()
        else:
            denominator = ring(denominator)
            if not denominator:
                raise ZeroDivisionError("denominator must be nonzero")
        # We normalize the fraction
        if numerator:
            if simplify and parent._simplification:
                D = numerator.left_gcd(denominator, monic=False)
                numerator, _ = numerator.left_quo_rem(D)
                denominator, _ = denominator.left_quo_rem(D)
            s = denominator.leading_coefficient()
            if s != 1:
                s = ~s
                numerator = s*numerator
                denominator = s*denominator
            self._numerator = numerator
            self._denominator = denominator
        else:
            self._numerator = ring.zero()
            self._denominator = ring.one()

    def _repr_(self):
        if not self._numerator:
            return "0"
        if self._denominator == 1:
            return str(self._numerator)
        if self._denominator.is_monomial():
            s = "%s^(-%s)" % (self.parent().variable_name(), self._denominator.degree())
        else:
            s = "(%s)^(-1)" % self._denominator
        if self._numerator == 1:
            return s
        if self._numerator._is_atomic():
            return "%s * %s" % (s, self._numerator)
        else:
            return "%s * (%s)" % (s, self._numerator)

    def _latex_(self):
        if not self._numerator:
            return "0"
        if self._denominator == 1:
            return latex(self._numerator)
        if self._denominator.is_monomial():
            s = "%s^{-%s}" % (self.parent().latex_variable_names()[0], self._denominator.degree())
        else:
            s = "\\left(%s\\right)^{-1}" % self._denominator
        if self._numerator == 1:
            return s
        if self._numerator._is_atomic():
            return "%s \\cdot %s" % (s, latex(self._numerator))
        else:
            return "%s \\cdot \\left(%s\\right)" % (s, latex(self._numerator))

    def __hash__(self):
        return hash((self._numerator, self._denominator))

    def _richcmp_(self, other, op):
        if self.parent()._simplification:
            return richcmp((self._numerator, self._denominator), (other._numerator, other._denominator), op) 
        if op == op_EQ or op == op_NE:
            _, U, V = self._denominator.left_xlcm(other._denominator)
            return richcmp(U * self._numerator, V * other._numerator, op)
        return NotImplemented

    def left_denominator(self):
        return self._denominator

    def right_numerator(self):
        return self._numerator

    @cached_method
    def _reverse_fraction(self):
        _, denominator, numerator = self._numerator.right_xlcm(self._denominator, monic=False)
        d = denominator.degree()
        s = ~(denominator.leading_coefficient())
        morphism = self.parent().twisting_morphism(-d)
        if morphism is not None:
            s = morphism(s)
        numerator = numerator * s
        denominator = denominator * s
        return numerator, denominator

    def right_denominator(self):
        return self._reverse_fraction()[1]

    def left_numerator(self):
        return self._reverse_fraction()[0]

    def is_zero(self):
        return self._numerator.is_zero()

    def _add_(self, other):
        denominator, U, V = self._denominator.left_xlcm(other._denominator)
        numerator = U * self._numerator + V * other._numerator
        return self.parent()(numerator, denominator)

    def _sub_(self, other):
        denominator, U, V = self._denominator.left_xlcm(other._denominator)
        numerator = U * self._numerator - V * other._numerator
        return self.parent()(numerator, denominator)

    def _neg_(self):
        return self.parent()(-self._numerator, self._denominator, simplify=False)

    def _mul_(self, other):
        if self.is_zero():
            return self
        if other.is_zero():
            return other
        L, U, V = self._numerator.left_xlcm(other._denominator, monic=False)
        denominator = U * self._denominator
        numerator = V * other._numerator
        return self.parent()(numerator, denominator)

    def _div_(self, other):
        if not other._numerator:
            raise ZeroDivisionError("cannot divide by zero")
        L, U, V = self._numerator.left_xlcm(other._numerator, monic=False)
        denominator = U * self._denominator
        numerator = V * other._denominator
        return self.parent()(numerator, denominator)

    def __invert__(self):
        if not self._numerator:
            raise ZeroDivisionError("cannot divide by zero")
        return self.parent()(self._denominator, self._numerator)

    def hilbert_shift(self, s):
        numerator = self._numerator.hilbert_shift(s)
        denominator = self._denominator.hilbert_shift(s)
        parent = numerator.parent().fraction_field()
        return parent(numerator, denominator, simplify=False)


class ConstantOreFunctionSection(Map):
    def _call_(self, x):
        if x.degree() <= 0:
            return x[0]
        else:
            raise TypeError("not a constant polynomial")

class OreFunctionBaseringInjection(Morphism):
    def __init__(self, domain, codomain):
        assert codomain.base_ring() is domain, \
            "the domain of the injection must be the base ring of the Ore function field"
        Morphism.__init__(self, Hom(domain,codomain))
        self._an_element = codomain.gen()
        self._repr_type_str = "Ore Function base injection"

    def an_element(self):
        return self._an_element

    def _call_(self, e):
        codomain = self.codomain()
        try:
            return codomain._element_constructor_(e)
        except AttributeError:
            return codomain(e)

    def section(self):
        return ConstantOrePolynomialSection(self._codomain, self.domain())
