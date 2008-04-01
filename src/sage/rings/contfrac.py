r"""
Continued Fractions

\sage implements the field code{ContinuedFractionField} (or \code{CFF}
for short) of finite simple continued fractions.  This is really
isomorphic to the field $\QQ$ of rational numbers, but with different
printing and semantics.  It should be possible to use this field in
most cases where one could use $\QQ$, except arithmetic is slower.

The \code{continued\_fraction(x)} command returns an element of
\code{CFF} that defines a continued fraction expansion to $x$.  The
command \code{continued\_fraction(x,bits)} computes the continued
fraction expansion of an approximation to $x$ with given bits of
precision.  Use \code{show(c)} to see a continued fraction nicely
typeset, and \code{latex(c)} to obtain the typeset version, e.g., for
inclusion in a paper.

EXAMPLES:
We create some example elements of the continued fraction field.

    sage: c = continued_fraction([1,2]); c
    [1, 2]
    sage: c = continued_fraction([3,7,15,1,292]); c
    [3, 7, 15, 1, 292]
    sage: c = continued_fraction(pi); c
    [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
    sage: c.value()
    245850922/78256779
    sage: QQ(c)
    245850922/78256779
    sage: RealField(200)(QQ(c) - pi)
    -7.8179366199075435400152113059910891481153981448107195930950e-17

We can also create matrices, polynomials, vectors, etc., over the continued
fraction field.
    sage: a = random_matrix(CFF, 4)
    sage: a
    [[0, 2]    [1]    [1] [0, 2]]
    [  [-1]    [1]   [-1]   [-2]]
    [   [2]    [1]   [-1]    [1]]
    [  [-1]    [1]   [-2] [0, 2]]
    sage: f = a.charpoly()
    sage: f
    [1]*x^4 + ([-1])*x^3 + [3, 1, 3]*x^2 + [3]*x + [-18]
    sage: f(a)
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    [[0] [0] [0] [0]]
    sage: vector(CFF, [1/2, 2/3, 3/4, 4/5])
    ([0, 2], [0, 1, 2], [0, 1, 3], [0, 1, 4])

"""
from sage.structure.element         import FieldElement
from sage.structure.parent_gens     import ParentWithGens
from sage.libs.pari.all             import pari

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

    EXAMPLES:
        sage: CFF
        Field of all continued fractions
        sage: ContinuedFractionField()
        Field of all continued fractions
    """
    def __init__(self):
        ParentWithGens.__init__(self, self)
        self._assign_names(('x'),normalize=False)

    def __cmp__(self, right):
        """
        EXAMPLES:
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
        EXAMPLES:
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
        EXAMPLES:
            sage: latex(CFF)
            \mathbf{CFF}
        """
        return "\\mathbf{CFF}"

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            return im_gens[0] == codomain._coerce_(self(1))
        except TypeError:
            return False

    def _repr_(self):
        """
        EXAMPLES:
            sage: CFF
            Field of all continued fractions
        """
        return "Field of all continued fractions"

    def _coerce_impl(self, x):
        """
        Anything that implicitly coerces to the rationals or a real
        field, implicitly coerces to the continued fraction field.

        EXAMPLES:
        The additions below call _coerce_impl implicitly:
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

    def __call__(self, x, bits=None):
        """
        INPUT:
            x -- a number
            bits -- integer (optional) the number of bits of the
                 input number to use when computing the continued fraction.

        EXAMPLES:
            sage: CFF(1.5)
            [1, 2]
            sage: CFF(e)
            [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 11]
            sage: CFF(pi)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
            sage: CFF([1,2,3])
            [1, 2, 3]
            sage: CFF(15/17)
            [0, 1, 7, 2]
            sage: c2 = loads(dumps(CFF))
            sage: c2(CFF(15/17)).parent() is c2
            True

        We illustrate varying the bits parameter:
            sage: CFF(pi)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
            sage: CFF(pi, bits=5)
            [3, 8]
            sage: CFF(pi, bits=80)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1, 84, 1]
            sage: CFF(pi, bits=100)
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1, 84, 2, 1, 1, 15, 3, 14]
        """
        if not bits is None:
            x = RealField(bits)(x)
        return ContinuedFraction(self, x)

    def __len__(self):
        """
        EXAMPLES:
            sage: len(CFF)
            Traceback (most recent call last):
            ...
            TypeError: len() of unsized object
        """
        raise TypeError, 'len() of unsized object'

    def gens(self):
        """
        EXAMPLES:
            sage: CFF.gens()
            ([1],)
        """
        return (self(1), )

    def gen(self, n=0):
        """
        EXAMPLES:
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
        EXAMPLES:
            sage: CFF.degree()
            1
        """
        return 1

    def ngens(self):
        """
        EXAMPLES:
            sage: CFF.ngens()
            1
        """
        return 1

    def is_field(self):
        """
        Return True, since the continued fraction field is a field.

        EXAMPLES:
            sage: CFF.is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Return False, since the continued fraction field is not finite.

        EXAMPLES:
            sage: CFF.is_finite()
            False
        """
        return False

    def is_atomic_repr(self):
        """
        Return False, since continued fractions are not atomic.

        EXAMPLES:
            sage: CFF.is_atomic_repr()
            False
        """
        return False

    def characteristic(self):
        """
        Return 0, since the continued fraction field has characteristic 0.

        EXAMPLES:
            sage: c = CFF.characteristic(); c
            0
            sage: parent(c)
            Integer Ring
        """
        return ZZ(0)

    def order(self):
        """
        EXAMPLES:
            sage: CFF.order()
            +Infinity
        """
        return infinity

    def random_element(self, num_bound=2, den_bound=2):
        """
        EXAMPLES:
           sage: CFF.random_element(10,10)
           [0, 4]
        """
        return self(QQ.random_element(num_bound = num_bound, den_bound = den_bound))


class ContinuedFraction(FieldElement):
    def __init__(self, parent, x):
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
            self._x = [ZZ(a) for a in continued_fraction_list(x)]

    def __getitem__(self, n):
        """
        Returns n-th term of the continued fraction.

        OUTPUT:
            an integer

        EXAMPLES:
            sage: a = continued_fraction(pi); a
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
            sage: a[4]
            292
            sage: a[-2]
            14
        """
        return self._x[n]

    def __getslice__(self, i, j):
        """
        OUTPUT:
            a continued fraction

        EXAMPLES:
            sage: a = continued_fraction(pi); a
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
            sage: a[2:5]
            [15, 1, 292]
            sage: a[:3]
            [3, 7, 15]
            sage: a[4:]
            [292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
        """
        return ContinuedFraction(self.parent(), self._x[i:j])

    def _repr_(self):
        """
        EXAMPLES:
            sage: a = continued_fraction(pi); a
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
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
            list of rational numbers

        EXAMPLES:
            sage: a = CFF(pi, bits=33); a
            [3, 7, 15, 1, 292, 2]
            sage: a.convergents()
            [3, 22/7, 333/106, 355/113, 103993/33102, 208341/66317]
            sage: a.value()
            208341/66317
            sage: a[:-1].value()
            103993/33102
        """
        return convergents(self._x)

    def convergent(self, n):
        """
        Return the n-th partial convergent to self.

        OUTPUT:
            rational number

        EXAMPLES:
            sage: a = CFF(pi, bits=33); a
            [3, 7, 15, 1, 292, 2]
            sage: a.convergents()
            [3, 22/7, 333/106, 355/113, 103993/33102, 208341/66317]
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

        EXAMPLES:
            sage: len(continued_fraction([1,2,3,4,5]) )
            5
        """
        return len(self._x)

    def pn(self, n):
        """
        Return the number of the n-th partial convergent, computed
        using the recurrence.

        EXAMPLES:
            sage: c = continued_fraction(pi); c
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
            sage: c.pn(0), c.qn(0)
            (3, 1)
            sage: len(c)
            14
            sage: c.pn(13), c.qn(13)
            (245850922, 78256779)
        """
        n = ZZ(n)
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
        Return the denominator of the n-th partial convergent, computed
        using the recurrence.

        EXAMPLES:
            sage: c = continued_fraction(pi); c
            [3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 3]
            sage: c.pn(0), c.qn(0)
            (3, 1)
            sage: len(c)
            14
            sage: c.pn(13), c.qn(13)
            (245850922, 78256779)
        """
        n = ZZ(n)
        if n < -2:
            raise ValueError, "n must be >= -2"
        if len(self.__qn) < n+3:
            self.pn(n)
        return self.__qn[n+2]

    def _rational_(self):
        """
        EXAMPLES:
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
        EXAMPLES:
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
        EXAMPLES:
            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: a.numerator()
            -17
        """
        return self._rational_().numerator()

    def denominator(self):
        """
        EXAMPLES:
            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: a.denominator()
            389
        """
        return self._rational_().denominator()

    def __int__(self):
        """
        EXAMPLES:
            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: int(a)
            -1
        """
        return int(self._rational_())

    def __long__(self):
        """
        EXAMPLES:
            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: long(a)
            -1L
        """
        return long(self._rational_())

    def __float__(self):
        """
        EXAMPLES:
            sage: a = CFF(-17/389); a
            [-1, 1, 21, 1, 7, 2]
            sage: float(a)
            -0.043701799485861177
        """
        return float(self._rational_())

    def _add_(self, right):
        """
        EXAMPLES:
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
        EXAMPLES:
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
        EXAMPLES:
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
        EXAMPLES:
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
        EXAMPLES:
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
        EXAMPLES:
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
            prec -- integer (default: 53) precision of square root
                    that is approximated
            all -- bool (default: False); if True, return all square
                   roots of self, instead of just one.

        EXAMPLES:
            sage: a = CFF(4/19); a
            [0, 4, 1, 3]
            sage: b = a.sqrt(); b
            [0, 2, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1, 15, 2]
            sage: b.value()
            146231375/318703893
            sage: float(b.value()^2 - a)
            4.0935373134057017e-17
            sage: b = a.sqrt(prec=100); b
            [0, 2, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 5]
            sage: b^2
            [0, 4, 1, 3, 49545773063556658177372134479, 1, 3, 4]
            sage: a.sqrt(all=True)
            [[0, 2, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1, 15, 2],
             [-1, 1, 1, 5, 1, 1, 2, 1, 16, 1, 2, 1, 1, 5, 4, 5, 1, 1, 2, 1, 15, 2]]
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
        Return copy of the underlying list of this continued
        fraction.

        EXAMPLES:
            sage: a = CFF(e); v = a.list(); v
            [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 11]
            sage: v[0] = 5
            sage: a
            [2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 11]
        """
        return list(self._x)

    def __hash__(self):
        """
        Return hash of self, which is the same as the hash of the value
        of self, as a rational number.

        EXAMPLES:
            sage: a = CFF(e)
            sage: hash(a)
            340185673
            sage: hash(QQ(a))
            340185673
        """
        return hash(self._rational_())

    def __invert__(self):
        """
        Return the multiplicative inverse of self.

        EXAMPLES:
            sage: a = CFF(e)
            sage: b = ~a; b
            [0, 2, 1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 11]
            sage: b*a
            [1]
        """
        return ContinuedFraction(self.parent(),
                     self._rational_().__invert__())

    def __pow__(self, n):
        """
        Return self to the power of n.

        EXAMPLES:
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

        EXAMPLES:
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

        EXAMPLES:
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

        EXAMPLES:
            sage: continued_fraction(1).is_one()
            True
            sage: continued_fraction(2).is_one()
            False
        """
        return self._rational_().is_one()

    def __nonzero__(self):
        """
        Return False if self is zero.

        EXAMPLES:
            sage: continued_fraction(0).is_zero()
            True
            sage: continued_fraction(1).is_zero()
            False
        """
        return not self._rational_().is_zero()

    def _pari_(self):
        """
        Return PARI list corresponding to this continued fraction.

        EXAMPLES:
            sage: c = continued_fraction(0.12345); c
            [0, 8, 9, 1, 21, 1, 1, 4, 1]
            sage: pari(c)
            [0, 8, 9, 1, 21, 1, 1, 4, 1]
        """
        return pari(self._x)

    def _interface_init_(self):
        """
        Return list representation for other systems corresponding to
        this continued fraction.

        EXAMPLES:
            sage: c = continued_fraction(0.12345); c
            [0, 8, 9, 1, 21, 1, 1, 4, 1]
            sage: gp(c)
            [0, 8, 9, 1, 21, 1, 1, 4, 1]
            sage: gap(c)
            [ 0, 8, 9, 1, 21, 1, 1, 4, 1 ]
            sage: maxima(c)
            [0,8,9,1,21,1,1,4,1]
        """
        return str(self._x)


CFF = ContinuedFractionField_class()

def ContinuedFractionField():
    return CFF

def continued_fraction(x, bits=None):
    return CFF(x, bits=bits)




