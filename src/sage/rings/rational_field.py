r"""
Field $\Q$ of Rational Numbers.

The class \class{RationalField} represents the field $\Q$ of
(arbitrary precision) rational numbers.  Each rational number is an
instance of the class \class{Rational}.
"""

import random
import field
import ring
import sage.rings.rational
import sage.structure.factorization
import complex_field
import infinity

_obj = {}
class _uniq(object):
    def __new__(cls):
        if _obj.has_key(0):
            return _obj[0]
        O = field.Field.__new__(cls)
        _obj[0] = O
        return O

class RationalField(_uniq, field.Field):
    r"""
    The class \class{RationalField} represents the field $\Q$ of
    rational numbers.
    """
    def __init__(self):
        """
        We create the rational numbers $\\Q$, and call a few functions:

            sage: Q = RationalField(); Q
            Rational Field
            sage: Q.characteristic()
            0
            sage: Q.is_field()
            True
            sage: Q.zeta()
            -1

        We next illustrate arithmetic in $\\Q$.

            sage: Q('49/7')
            7
            sage: type(Q('49/7'))
            <type 'sage.rings.rational.Rational'>
            sage: a = Q('19/374'); b = Q('17/371'); print a, b
            19/374 17/371
            sage: a + b
            13407/138754
            sage: b + a
            13407/138754
            sage: a * b
            19/8162
            sage: b * a
            19/8162
            sage: a - b
            691/138754
            sage: b - a
            -691/138754
            sage: a / b
            7049/6358
            sage: b / a
            6358/7049
            sage: b < a
            True
            sage: a < b
            False

        Next finally illustrate arithmetic with automatic coercion.
        The types that coerce into the rational field include
        \\code{str, int, long, Integer}.

            sage: a + Q('17/371')
            13407/138754
            sage: a * 374
            19
            sage: 374 * a
            19
            sage: a/19
            1/374
            sage: a + 1
            393/374
        """

    def __repr__(self):
        return "Rational Field"

    def _latex_(self):
        return "\\mathbf{Q}"

    def __len__(self):
        raise TypeError, 'len() of unsized object'

    def __call__(self, x, base=0):
        """
        Coerce x into the field of rational numbers.

        EXAMPLES:
            sage: a = long(901824309821093821093812093810928309183091832091)
            sage: b = QQ(a); b
            901824309821093821093812093810928309183091832091
            sage: QQ(b)
            901824309821093821093812093810928309183091832091
            sage: QQ(int(93820984323))
            93820984323
            sage: QQ(ZZ(901824309821093821093812093810928309183091832091))
            901824309821093821093812093810928309183091832091
            sage: QQ('-930482/9320842317')
            -930482/9320842317
            sage: QQ((-930482, 9320842317))
            -930482/9320842317
            sage: QQ([9320842317])
            9320842317
            sage: QQ(pari(39029384023840928309482842098430284398243982394))
            39029384023840928309482842098430284398243982394
            sage: QQ('sage')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sage to a rational

        Coercion from the reals to the rational is done by default
        using continued fractions.

            sage: QQ(RR(3929329/32))
            3929329/32
            sage: QQ(-RR(3929329/32))
            -3929329/32
            sage: QQ(RR(1/7)) - 1/7
            0

        If you specify an optional second base argument, then
        the string representation of the float is used.
            sage: QQ(23.2, 2)
            6530219459687219/281474976710656
            sage: 6530219459687219.0/281474976710656
            23.19999999999999929
            sage: QQ(23.2, 10)
            116/5

        Here's a nice example involving elliptic curves:
            sage: E = EllipticCurve('11a')
            sage: L = E.Lseries_at1(300)[0]; L
            0.25384186085600002
            sage: O = E.omega(); O
            1.269209304279553421688794613    # 32-bit
            1.269209304279553421688794616754547304  # 64-bit
            sage: t = L/O; t
            0.20000000000007040
            sage: QQ(t)
            1/5
        """
        if isinstance(x, sage.rings.rational.Rational):
            return x
        return sage.rings.rational.Rational(x, base)

    def _coerce_(self, x):
        if isinstance(x, sage.rings.rational.Rational):
            return x
        elif isinstance(x, (int, long, sage.rings.integer.Integer)):
            return self(x)
        raise TypeError

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def __iter__(self):
        yield self(0)
        yield self(1)
        yield self(-1)
        from integer_ring import IntegerRing
        for n in IntegerRing():
            m = abs(n)
            for d in abs(n).coprime_integers(m):
                yield n/d
                yield d/n

    def complex_embedding(self, prec=53):
        CC = complex_field.ComplexField(prec)
        return self.hom([CC(1)])

    def gens(self):
        return (self(1), )

    def gen(self, n=0):
        if n == 0:
            return self(1)
        else:
            raise IndexError, "n must be 0"

    def degree(self):
        return 1

    def ngens(self):
        return 1

    def is_field(self):
        """
        Return True, since the rational field is a field.
        """
        return True

    def is_finite(self):
        """
        Return False, since the rational field is not finite.
        """
        return False

    def is_prime_field(self):
        return True

    def is_atomic_repr(self):
        return True

    def characteristic(self):
        """
        Return 0, since the rational field has characteristic 0.

        EXAMPLES:
            sage: c = QQ.characteristic(); c
            0
            sage: parent(c)
            Integer Ring
        """
        return sage.rings.integer.Integer(0)

    def number_field(self):
        from sage.rings.number_field.all import NumberField
        x = sage.rings.polynomial_ring.PolynomialRing(self).gen()
        return NumberField(x-1)

    def order(self):
        """
        EXAMPLES:
            sage: QQ.order()
            Infinity
        """
        return infinity.infinity

    def random_element(self, num_bound=1, den_bound=1):
        """
        EXAMPLES:
            sage: QQ.random_element(10,10)
            -5/3
        """
        return self("%s/%s"%(random.randrange(-num_bound-1, num_bound+1), \
                             random.randrange(1,den_bound+1)))

    def __cmp__(self, other):
        if isinstance(other, RationalField):
            return 0
        if ring.is_Ring(other):
            if other.characteristic() == 0 and field.is_Field(other):
                return -1
            else:
                return 1
        return -1

    def zeta(self, n=2):
        if n == 1:
            return sage.rings.rational.Rational(1)
        elif n == 2:
            return sage.rings.rational.Rational(-1)
        else:
            raise ValueError, "no n-th root of unity in rational field"

    #################################
    ## Coercions to interfaces
    #################################
    def _gap_init_(self):
        """
        EXAMPLES:
            sage: gap(QQ)
            Rationals
        """
        return 'Rationals'

    def _magma_init_(self):
        """
        EXAMPLES:
            sage: magma(QQ)                       # optional
            Rational Field
        """
        return 'RationalField()'


QQ = RationalField()
Q = QQ

def is_RationalField(x):
    return isinstance(x, RationalField)

def frac(n,d):
    return sage.rings.rational.Rational(n)/sage.rings.rational.Rational(d)

def factor(x):
    assert isinstance(x, sage.rings.rational.Rational)
    return x.numerator().factor() * \
           sage.structure.factorization.Factorization([(p,-e) for p, e in x.denominator().factor()])

