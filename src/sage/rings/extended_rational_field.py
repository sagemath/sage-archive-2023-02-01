import sage.rings.rational_field
import sage.rings.rational
import sage.rings.integer
import sage.rings.infinity
import sage.structure.element

Rational = sage.rings.rational.Rational
RationalField = sage.rings.rational_field.RationalField

_obj = {}
class _uniq(object):
    def __new__(cls):
        if _obj.has_key(0):
            return _obj[0]
        O = field.Field.__new__(cls)
        _obj[0] = O
        return O


class ExtendedRationalField_class(_uniq, RationalField):
    def __init__(self):
        ParentWithGens.__init__(self, self)
        self._assign_names(('x'),normalize=False)

    def _repr_(self):
        return "Extended Rational Field"

    def _latex_(self):
        return "\\mathbf{Q}\\cup\\{\\pm\\infty\\}"

    def __call__(self, x, base = 0):
        if isinstance(x, sage.rings.infinity.MinusInfinity):
            return self.gen(2)
        if isinstance(x, sage.structure.element.InfinityElement):
            return self.gen(1)
        if isinstance(x, sage.rings.infinity.FiniteNumber):
            if x == 0:
                y = Rational(0)
                y._parent = self
                return y
            raise TypeError, "cannot coerce unknown finite number into the extended rationals"
        y = ExtendedRational(x, base)
        y._parent = self
        return y

    def _coerce_impl(self, x):
        if isinstance(x, (int, long, sage.rings.integer.Integer, Rational)):
            return self(x)
        if isinstance(x, (RationalPlusInfinity, RationalMinusInfinity, IntegerPlusInfinity, IntegerMinusInfinity)):
            return self(x)
        raise TypeError, "no implicit coercion of element to the rational numbers"

    def _is_valid_homomorphism(self, codomain, im_gens):
        raise NotImplementedError

    def __iter__(self):
        yield self(0)
        yield self(1)
        yield self(-1)
        yield self.gen(1)
        yield self.gen(2)
        from integer_ring import IntegerRing
        for n in IntegerRing():
            m = abs(n)
            for d in abs(n).coprime_integers(m):
                yield n/d
                yield d/n

    def complex_embedding(self, prec=53):
        raise NotImplementedError

    def gens(self):
        return (self(1), self.gen(1), self.gen(2), )

    def gen(self, n=0):
        if n == 0:
            return self(1)
        elif n == 1:
            try:
                return self.gen1
            except AttributeError:
                self.gen1 = sage.rings.infinity.PlusInfinity(self)
                return self.gen1
        elif n == 2:
            try:
                return self.gen2
            except AttributeError:
                self.gen2 = sage.rings.infinity.MinusInfinity(self)
                return self.gen1
        else:
            raise IndexError, "n must be 0, 1 or 2"

    def is_prime_field(self):
        return False

    def ngens(self):
        return 3

    def numberfield(self, poly_var, nf_var):
        raise NotImplementedError

ExtendRationalField = ExtendedRationalField_class()

class ExtendedRational(Rational):
    def __init__(self, x = None, base = 0):
        self._parent = ExtendedRationalField
        Rational.__init__(self, x, base)


