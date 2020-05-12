from sage.misc.cachefunc import cached_method
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
                raise ZeroDivisionError("division by zero")

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
        if self._denominator._is_atomic():
            s = "%s^(-1)" % self._denominator
        else:
            s = "(%s)^(-1)" % self._denominator
        if self._numerator == 1:
            return s
        if self._numerator._is_atomic():
            return "%s * %s" % (s, self._numerator)
        else:
            return "%s * (%s)" % (s, self._numerator)

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
