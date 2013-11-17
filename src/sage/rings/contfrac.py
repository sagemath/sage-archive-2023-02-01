r"""
Continued Fractions

Sage implements the field ``ContinuedFractionField`` (or ``CFF``
for short) of finite simple continued fractions.  This is really
isomorphic to the field `\QQ` of rational numbers, but with different
printing and semantics.  It should be possible to use this field in
most cases where one could use `\QQ`, except arithmetic is slower.

The ``continued_fraction(x)`` command returns an element of
``CFF`` that defines a continued fraction expansion to `x`.  The
command ``continued_fraction(x,bits)`` computes the continued
fraction expansion of an approximation to `x` with given bits of
precision.  Use ``show(c)`` to see a continued fraction nicely
typeset, and ``latex(c)`` to obtain the typeset version, e.g., for
inclusion in a paper.

EXAMPLES:
We create some example elements of the continued fraction field::

    sage: c = continued_fraction([1,2]); c
    [1, 2]
    sage: c = continued_fraction([3,7,15,1,292]); c
    [3, 7, 15, 1, 292]
    sage: c = continued_fraction(pi); c
    [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
    sage: c.value()
    80143857/25510582
    sage: QQ(c)
    80143857/25510582
    sage: RealField(200)(QQ(c) - pi)
    -5.7908701643756732744264903067012149647564522968979302505514e-16

We can also create matrices, polynomials, vectors, etc., over the continued
fraction field.

::

    sage: a = random_matrix(CFF, 4)
    sage: a
    [    [-1, 2] [-1, 1, 94]      [0, 2]       [-12]]
    [       [-1]      [0, 2]  [-1, 1, 3]   [0, 1, 2]]
    [    [-3, 2]         [0]   [0, 1, 2]        [-1]]
    [        [1]        [-1]      [0, 3]         [1]]
    sage: f = a.charpoly()
    sage: f
    [1]*x^4 + ([-2, 3])*x^3 + [14, 1, 1, 1, 9, 1, 8]*x^2 + ([-13, 4, 1, 2, 1, 1, 1, 1, 1, 2, 2])*x + [-6, 1, 5, 9, 1, 5]
    sage: f(a)
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    sage: vector(CFF, [1/2, 2/3, 3/4, 4/5])
    ([0, 2], [0, 1, 2], [0, 1, 3], [0, 1, 4])

AUTHORS:

- Niles Johnson (2010-08): Trac #3893: ``random_element()`` should pass on ``*args`` and ``**kwds``.

"""
from sage.structure.element         import FieldElement
from sage.structure.parent_gens     import ParentWithGens
from sage.libs.pari.pari_instance   import pari

from field                          import Field
from rational_field                 import QQ
from integer_ring                   import ZZ
from infinity                       import infinity
from real_mpfr                      import is_RealNumber, RealField
from real_double                    import RDF
from arith                          import (continued_fraction_list,
                                            convergent, convergents)


class ContinuedFractionField_class(Field):
    """
    The field of all finite continued fraction of real numbers.

    EXAMPLES::

        sage: CFF
        Field of all continued fractions

    The continued fraction field inherits from the base class
    :class:`sage.rings.ring.Field`. However it was initialised
    as such only since trac ticket #11900::

        sage: CFF.category()
        Category of fields

    """
    def __init__(self):
        """
        EXAMPLES::

            sage: ContinuedFractionField()
            Field of all continued fractions

        TESTS::

            sage: CFF._repr_option('element_is_atomic')
            False
        """
        Field.__init__(self, self)
        self._assign_names(('x'),normalize=False)

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: CFF == ContinuedFractionField()
            True
            sage: CFF == CDF
            False
            sage: loads(dumps(CFF)) == CFF
            True
        """
        return cmp(type(self), type(right))

    def __iter__(self):
        """
        EXAMPLES::

            sage: i = 0
            sage: for a in CFF:
            ...    print a
            ...    i += 1
            ...    if i > 5: break
            ...
            [0]
            [1]
            [-1]
            [0, 2]
            [-1, 2]
            [2]
        """
        for n in QQ:
            yield self(n)


    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(CFF)
            \Bold{CFF}
        """
        return "\\Bold{CFF}"

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        Return whether or not the map to codomain by sending the
        continued fraction [1] of self to im_gens[0] is a
        homomorphism.

        EXAMPLES::

            sage: CFF._is_valid_homomorphism_(ZZ,[ZZ(1)])
            False
            sage: CFF._is_valid_homomorphism_(CFF,[CFF(1)])
            True
        """
        try:
            return im_gens[0] == codomain._coerce_(self(1))
        except TypeError:
            return False

    def _repr_(self):
        """
        EXAMPLES::

            sage: CFF
            Field of all continued fractions
        """
        return "Field of all continued fractions"

    def _coerce_impl(self, x):
        """
        Anything that implicitly coerces to the rationals or a real
        field, implicitly coerces to the continued fraction field.

        EXAMPLES:

        The additions below call _coerce_impl implicitly::

            sage: a = CFF(3/5); a
            [0, 1, 1, 2]
            sage: a + 2/5
            [1]
            sage: 2/5 + a
            [1]
            sage: 1.5 + a
            [2, 10]
            sage: a  + 1.5
            [2, 10]
        """
        if is_RealNumber(x):
            return self(x)
        return self._coerce_try(x, [QQ, RDF])

    def __call__(self, x, bits=None, nterms=None):
        """
        INPUT:

            - `x` -- a number

            - ``bits`` -- integer (optional) the number of bits of the
              input number to use when computing the continued fraction.

        EXAMPLES::

            sage: CFF(1.5)
            [1, 2]
            sage: CFF(e)
            [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1]
            sage: CFF(pi)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
            sage: CFF([1,2,3])
            [1, 2, 3]
            sage: CFF(15/17)
            [0, 1, 7, 2]
            sage: c2 = loads(dumps(CFF))
            sage: c2(CFF(15/17)).parent() is c2
            True

        We illustrate varying the bits parameter::

            sage: CFF(pi)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
            sage: CFF(pi, bits=20)
            [3, 7]
            sage: CFF(pi, bits=80)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1]
            sage: CFF(pi, bits=100)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1, 84, 2, 1, 1, 15, 3]

        And varying the nterms parameter::

            sage: CFF(pi, nterms=3)
            [3, 7, 15]
            sage: CFF(pi, nterms=10)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1]
            sage: CFF(pi, bits=10, nterms=10)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1]
        """
        return ContinuedFraction(self, x, bits, nterms)

    def __len__(self):
        """
        EXAMPLES::

            sage: len(CFF)
            Traceback (most recent call last):
            ...
            TypeError: len() of unsized object
        """
        raise TypeError, 'len() of unsized object'

    def gens(self):
        """
        EXAMPLES::

            sage: CFF.gens()
            ([1],)
        """
        return (self(1), )

    def gen(self, n=0):
        """
        EXAMPLES::

            sage: CFF.gen()
            [1]
            sage: CFF.0
            [1]
        """

        if n == 0:
            return self(1)
        else:
            raise IndexError, "n must be 0"

    def degree(self):
        """
        EXAMPLES::

            sage: CFF.degree()
            1
        """
        return 1

    def ngens(self):
        """
        EXAMPLES::

            sage: CFF.ngens()
            1
        """
        return 1

    def is_field(self, proof = True):
        """
        Return True, since the continued fraction field is a field.

        EXAMPLES::

            sage: CFF.is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Return False, since the continued fraction field is not finite.

        EXAMPLES::

            sage: CFF.is_finite()
            False
        """
        return False

    def characteristic(self):
        """
        Return 0, since the continued fraction field has characteristic 0.

        EXAMPLES::

            sage: c = CFF.characteristic(); c
            0
            sage: parent(c)
            Integer Ring
        """
        return ZZ(0)

    def order(self):
        """
        EXAMPLES::

            sage: CFF.order()
            +Infinity
        """
        return infinity

    def random_element(self, *args, **kwds):
        """
        EXAMPLES::

           sage: CFF.random_element(10,10)
           [0, 4]

        Passes extra positional or keyword arguments through::

           sage: [CFF.random_element(den_bound=10, num_bound=2) for x in range(4)]
           [[-1, 1, 3], [0, 7], [0, 3], [0, 4]]



        """
        return self(QQ.random_element(*args, **kwds))


class ContinuedFraction(FieldElement):
    """
    A continued fraction object.

    EXAMPLES::

        sage: continued_fraction(pi)
        [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
        sage: CFF(pi)
        [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
    """
    def __init__(self, parent, x, bits=None, nterms=None):
        """
        EXAMPLES::

            sage: sage.rings.contfrac.ContinuedFraction(CFF,[1,2,3,4,1,2])
            [1, 2, 3, 4, 1, 2]
            sage: sage.rings.contfrac.ContinuedFraction(CFF,[1,2,3,4,-1,2])
            Traceback (most recent call last):
            ...
            ValueError: each entry except the first must be positive
        """
        FieldElement.__init__(self, parent)
        if isinstance(x, ContinuedFraction):
            self._x = list(x._x)
        elif isinstance(x, (list, tuple)):
            x = [ZZ(a) for a in x]
            for i in range(1,len(x)):
                if x[i] <= 0:
                    raise ValueError, "each entry except the first must be positive"
            self._x = list(x)
        else:
            self._x = [ZZ(a) for a in continued_fraction_list(x, bits=bits, nterms=nterms)]

    def __getitem__(self, n):
        """
        Returns `n`-th term of the continued fraction.

        OUTPUT:
            - an integer or a a continued fraction

        EXAMPLES::

            sage: a = continued_fraction(pi); a
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
            sage: a[4]
            292
            sage: a[-1]
            14
            sage: a[2:5]
            [15, 1, 292]
            sage: a[:3]
            [3, 7, 15]
            sage: a[4:]
            [292, 1, 1, 1, 2, 1, 3, 1, 14]
            sage: a[4::2]
            [292, 1, 2, 3, 14]
        """
        if isinstance(n, slice):
            start, stop, step = n.indices(len(self))
            return ContinuedFraction(self.parent(), self._x[start:stop:step])
        else:
            return self._x[n]

    def _repr_(self):
        """
        EXAMPLES::

            sage: a = continued_fraction(pi); a
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
            sage: a.rename('continued fraction of pi')
            sage: a
            continued fraction of pi
        """
        return str(self._x)

    def convergents(self):
        """
        Return a list of rational numbers, which are the partial
        convergents of this continued fraction.

        OUTPUT:

            - list of rational numbers

        EXAMPLES::

            sage: a = CFF(pi, bits=34); a
            [3, 7, 15, 1, 292]
            sage: a.convergents()
            [3, 22/7, 333/106, 355/113, 103993/33102]
            sage: a.value()
            103993/33102
            sage: a[:-1].value()
            355/113
        """
        return convergents(self._x)

    def convergent(self, n):
        """
        Return the `n`-th partial convergent to self.

        OUTPUT:

            rational number

        EXAMPLES::

            sage: a = CFF(pi, bits=34); a
            [3, 7, 15, 1, 292]
            sage: a.convergents()
            [3, 22/7, 333/106, 355/113, 103993/33102]
            sage: a.convergent(0)
            3
            sage: a.convergent(1)
            22/7
            sage: a.convergent(4)
            103993/33102
        """
        return convergent(self._x, n)

    def __len__(self):
        """
        Return the number of terms in this continued fraction.

        EXAMPLES::

            sage: len(continued_fraction([1,2,3,4,5]) )
            5
        """
        return len(self._x)

    def pn(self, n):
        """
        Return the numerator of the `n`-th partial convergent, computed
        using the recurrence.

        EXAMPLES::

            sage: c = continued_fraction(pi); c
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
            sage: c.pn(0), c.qn(0)
            (3, 1)
            sage: len(c)
            13
            sage: c.pn(12), c.qn(12)
            (80143857, 25510582)
        """
        if n < -2:
            raise ValueError, "n must be at least -2"
        if n > len(self._x):
            raise ValueError, "n must be at most %s"%len(self._x)
        try:
            return self.__pn[n+2]
        except AttributeError:
            self.__pn = [0, 1, self._x[0]]
            self.__qn = [1, 0, 1]
        except IndexError:
            pass
        for k in range(len(self.__pn), n+3):
            self.__pn.append(self._x[k-2]*self.__pn[k-1] + self.__pn[k-2])
            self.__qn.append(self._x[k-2]*self.__qn[k-1] + self.__qn[k-2])
        return self.__pn[n+2]

    def qn(self, n):
        """
        Return the denominator of the `n`-th partial convergent, computed
        using the recurrence.

        EXAMPLES::

            sage: c = continued_fraction(pi); c
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14]
            sage: c.qn(0), c.pn(0)
            (1, 3)
            sage: len(c)
            13
            sage: c.pn(12), c.qn(12)
            (80143857, 25510582)
        """
        if n < -2:
            raise ValueError, "n must be at least -2"
        if n > len(self._x):
            raise ValueError, "n must be at most %s"%len(self._x)
        try:
            return self.__qn[n+2]
        except (AttributeError, IndexError):
            pass
        self.pn(n)
        return self.__qn[n+2]

    def _rational_(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: a._rational_()
            -17/389
            sage: QQ(a)
            -17/389
        """
        try:
            return self.__rational
        except AttributeError:
            r = convergents(self._x)[-1]
            self.__rational =r
            return r

    def value(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: a.value()
            -17/389
            sage: QQ(a)
            -17/389
        """
        return self._rational_()

    def numerator(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: a.numerator()
            -17
        """
        return self._rational_().numerator()

    def denominator(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: a.denominator()
            389
        """
        return self._rational_().denominator()

    def __int__(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: int(a)
            -1
        """
        return int(self._rational_())

    def __long__(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: long(a)
            -1L
        """
        return long(self._rational_())

    def __float__(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: float(a)
            -0.043701799485861184
        """
        return float(self._rational_())

    def _add_(self, right):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: c = a+b; c
            [-1, 1, 23, 3, 5]
            sage: c.value()
            -16/389
        """
        return ContinuedFraction(self.parent(),
                     self._rational_() + right._rational_())

    def _sub_(self, right):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: c = a - b; c
            [-1, 1, 20, 1, 1, 1, 1, 3]
            sage: c.value()
            -18/389
        """
        return ContinuedFraction(self.parent(),
                     self._rational_() - right._rational_())

    def _mul_(self, right):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: c = a * b; c
            [-1, 1, 8900, 4, 4]
            sage: c.value(), (-1/389)*(17/389)
            (-17/151321, -17/151321)
        """
        return ContinuedFraction(self.parent(),
                     self._rational_() * right._rational_())

    def _div_(self, right):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: c = a / b; c
            [-17]
            sage: c.value(), (17/389) / (-1/389)
            (-17, -17)
        """
        return ContinuedFraction(self.parent(),
                     self._rational_() / right._rational_())

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: b = CFF(1/389)
            sage: a < b
            True
            sage: QQ(a) < QQ(b)
            True
            sage: QQ(a)
            -17/389
            sage: QQ(b)
            1/389
        """
        return cmp(self._rational_(), right._rational_())

    def _latex_(self):
        """
        EXAMPLES::

            sage: a = CFF(-17/389)
            sage: latex(a)
            -1+ \frac{\displaystyle 1}{\displaystyle 1+ \frac{\displaystyle 1}{\displaystyle 21+ \frac{\displaystyle 1}{\displaystyle 1+ \frac{\displaystyle 1}{\displaystyle 7+ \frac{\displaystyle 1}{\displaystyle 2}}}}}
        """
        v = self._x
        if len(v) == 0:
            return '0'
        s = str(v[0])
        for i in range(1,len(v)):
            s += '+ \\frac{\\displaystyle 1}{\\displaystyle %s'%v[i]
        s += '}'*(len(v)-1)
        return s

    def sqrt(self, prec=53, all=False):
        """
        Return continued fraction approximation to square root of the
        value of this continued fraction.

        INPUT:

            - `prec` -- integer (default: 53) precision of square root
              that is approximated

            - `all` -- bool (default: False); if True, return all square
              roots of self, instead of just one.

        EXAMPLES::

            sage: a = CFF(4/19); a
            [0, 4, 1, 3]
            sage: b = a.sqrt(); b
            [0, 2, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1]
            sage: b.value()
            4508361/9825745
            sage: float(b.value()^2 - a)
            -5.451492525672688e-16
            sage: b = a.sqrt(prec=100); b
            [0, 2, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5]
            sage: b^2
            [0, 4, 1, 3, 7849253184229368265220252099, 1, 3]
            sage: a.sqrt(all=True)
            [[0, 2, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1],
             [-1, 1, 1, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1]]
            sage: a = CFF(4/25).sqrt(); a
            [0, 2, 2]
            sage: a.value()
            2/5
        """
        r = self._rational_()
        if r < 0:
            raise ValueError, "self must be positive"
        X = r.sqrt(all=all, prec=prec)
        if not all:
            return ContinuedFraction(self.parent(), X)
        else:
            return [ContinuedFraction(self.parent(), x) for x in X]

    def list(self):
        """
        Return copy of the underlying list of this continued fraction.

        EXAMPLES::

            sage: a = CFF(e); v = a.list(); v
            [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1]
            sage: v[0] = 5
            sage: a
            [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1]
        """
        return list(self._x)

    def __hash__(self):
        """
        Return hash of self, which is the same as the hash of the value
        of self, as a rational number.

        EXAMPLES::

            sage: a = CFF(e)
            sage: hash(a)
            19952398
            sage: hash(QQ(a))
            19952398
        """
        return hash(self._rational_())

    def __invert__(self):
        """
        Return the multiplicative inverse of self.

        EXAMPLES::

            sage: a = CFF(e)
            sage: b = ~a; b
            [0, 2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 2]
            sage: b*a
            [1]
        """
        return ContinuedFraction(self.parent(),
                     self._rational_().__invert__())

    def __pow__(self, n):
        """
        Return self to the power of `n`.

        EXAMPLES::

            sage: a = CFF([1,2,3]); a
            [1, 2, 3]
            sage: a^3
            [2, 1, 10, 1, 4, 1, 4]
            sage: QQ(a)^3 == QQ(a^3)
            True
            sage: a^(-3)
            [0, 2, 1, 10, 1, 4, 1, 4]
            sage: QQ(a)^(-3) == QQ(a^(-3))
            True
        """
        return ContinuedFraction(self.parent(),
                     self._rational_()**n)

    def __neg__(self):
        """
        Return additive inverse of self.

        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: -a
            [0, 22, 1, 7, 2]
            sage: QQ(-a)
            17/389
        """
        return ContinuedFraction(self.parent(),
                     self._rational_().__neg__())

    def __abs__(self):
        """
        Return absolute value of self.

        EXAMPLES::

            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: abs(a)
            [0, 22, 1, 7, 2]
            sage: QQ(abs(a))
            17/389
        """
        return ContinuedFraction(self.parent(),
                     self._rational_().__abs__())

    def is_one(self):
        """
        Return True if self is one.

        EXAMPLES::

            sage: continued_fraction(1).is_one()
            True
            sage: continued_fraction(2).is_one()
            False
        """
        return self._rational_().is_one()

    def __nonzero__(self):
        """
        Return False if self is zero.

        EXAMPLES::

            sage: continued_fraction(0).is_zero()
            True
            sage: continued_fraction(1).is_zero()
            False
        """
        return not self._rational_().is_zero()

    def _pari_(self):
        """
        Return PARI list corresponding to this continued fraction.

        EXAMPLES::

            sage: c = continued_fraction(0.12345); c
            [0, 8, 9, 1, 21, 1, 1]
            sage: pari(c)
            [0, 8, 9, 1, 21, 1, 1]
        """
        return pari(self._x)

    def _interface_init_(self, I=None):
        """
        Return list representation for other systems corresponding to
        this continued fraction.

        EXAMPLES::

            sage: c = continued_fraction(0.12345); c
            [0, 8, 9, 1, 21, 1, 1]
            sage: gp(c)
            [0, 8, 9, 1, 21, 1, 1]
            sage: gap(c)
            [ 0, 8, 9, 1, 21, 1, 1 ]
            sage: maxima(c)
            [0,8,9,1,21,1,1]
        """
        return str(self._x)

    def additive_order(self):
        """
        Return the additive order of this continued fraction,
        which we defined to be the additive order of its value.

        EXAMPLES::

            sage: CFF(-1).additive_order()
            +Infinity
            sage: CFF(0).additive_order()
            1
        """
        return self.value().additive_order()

    def multiplicative_order(self):
        """
        Return the multiplicative order of this continued fraction,
        which we defined to be the multiplicative order of its value.

        EXAMPLES::

            sage: CFF(-1).multiplicative_order()
            2
            sage: CFF(1).multiplicative_order()
            1
            sage: CFF(pi).multiplicative_order()
            +Infinity
        """
        return self.value().multiplicative_order()


CFF = ContinuedFractionField_class()

def ContinuedFractionField():
    """
    Return the (unique) field of all continued fractions.

    EXAMPLES::

        sage: ContinuedFractionField()
        Field of all continued fractions
    """
    return CFF

def continued_fraction(x, bits=None, nterms=None):
    """
    Return the truncated continued fraction expansion of the real number
    `x`, computed with an interval floating point approximation of `x`
    to the given number of bits of precision. The returned continued
    fraction is a list-like object, with a value method and partial
    convergents method.

    If bits is not given, then use the number of valid bits of
    precision of `x`, if `x` is a floating point number, or 53 bits
    otherwise. If nterms is given, the precision is increased until
    the specified number of terms can be computed, if possible.

    INPUT:

        - `x` -- number

        - ``bits`` -- None (default) or a positive integer

        - ``nterms`` -- None (default) or a positive integer

    OUTPUT:

        - a continued fraction

    EXAMPLES::

        sage: v = continued_fraction(sqrt(2)); v
        [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: v = continued_fraction(sqrt(2), nterms=22); v
        [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        sage: type(v)
        <class 'sage.rings.contfrac.ContinuedFraction'>
        sage: parent(v)
        Field of all continued fractions
        sage: v.value()
        131836323/93222358
        sage: RR(v.value()) == RR(sqrt(2))
        True
        sage: v.convergents()
        [1, 3/2, 7/5, 17/12, 41/29, 99/70, 239/169,...131836323/93222358]
        sage: [RR(x) for x in v.convergents()]
        [1.00000000000000, 1.50000000000000, 1.40000000000000, 1.41666666666667, ...1.41421356237310]
        sage: continued_fraction(sqrt(2), 10)
        [1, 2, 2]
        sage: v.numerator()
        131836323
        sage: v.denominator()
        93222358
        sage: [v.pn(i) for i in range(10)]
        [1, 3, 7, 17, 41, 99, 239, 577, 1393, 3363]
        sage: [v.qn(i) for i in range(10)]
        [1, 2, 5, 12, 29, 70, 169, 408, 985, 2378]

    Here are some more examples::

        sage: continued_fraction(e, bits=20)
        [2, 1, 2, 1, 1, 4, 1, 1]
        sage: continued_fraction(e, bits=30)
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8]
        sage: continued_fraction(RealField(200)(e))
        [2, 1, 2, 1, 1, 4, 1, 1, 6, ...36, 1, 1, 38, 1, 1]

    Initial rounding can result in incorrect trailing digits::

        sage: continued_fraction(RealField(39)(e))
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 2]
        sage: continued_fraction(RealIntervalField(39)(e))
        [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10]
    """
    return CFF(x, bits=bits, nterms=nterms)




