r"""
Field $\Q$ of Rational Numbers.

The class \class{RationalField} represents the field $\Q$ of
(arbitrary precision) rational numbers.  Each rational number is an
instance of the class \class{Rational}.

TEST:
   sage: Q = RationalField()
   sage: Q == loads(dumps(Q))
   True

"""

import random
import field
import ring
import sage.rings.rational
import sage.structure.factorization
import infinity

ZZ = None

from sage.structure.parent_gens import ParentWithGens

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
        ParentWithGens.__init__(self, self)
        self._assign_names(('x'),normalize=False)

    def __hash__(self):
        return -11115808

    def _repr_(self):
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
            23.199999999999999
            sage: a = 23.2; a
            23.2000000000000
            sage: QQ(a, 10)
            116/5

        Here's a nice example involving elliptic curves:
            sage: E = EllipticCurve('11a')
            sage: L = E.Lseries_at1(300)[0]; L
            0.253841860855911
            sage: O = E.omega(); O
            1.269209304279553421688794616754547305219492241830608667967136921230408338613     # 32-bit
            1.26920930427955342168879461675454730521949224183060866796713692123040833861277772269036230592151260731164529627832128743728170032847684397649271401057075        # 64-bit
            sage: t = L/O; t
            0.200000000000000
            sage: QQ(RealField(45)(t))
            1/5
        """
        if isinstance(x, sage.rings.rational.Rational):
            return x
        return sage.rings.rational.Rational(x, base)

    def construction(self):
        from sage.categories.pushout import FractionField
        import integer_ring
        return FractionField(), integer_ring.ZZ

    def completion(self, p, prec, extras = {}):
        if p == infinity.Infinity:
            from sage.rings.real_mpfr import RealField
            return RealField(prec, **extras) # TODO: handle RDF/RQDF/rounding via extras, probably via an extra utility function
        else:
            from sage.rings.padics.factory import Qp
            return Qp(p, prec, **extras)

    def _coerce_impl(self, x):
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            return self(x)
        raise TypeError, 'no implicit coercion of element to the rational numbers'

    def coerce_map_from_impl(self, S):
        global ZZ
        if ZZ is None:
            import integer_ring
            ZZ = integer_ring.ZZ
        if S is ZZ:
            return sage.rings.rational.Z_to_Q()
        else:
            return field.Field.coerce_map_from_impl(self, S)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def __iter__(self):
        r"""
        Creates an iterator that generates the rational numbers without
        repetition. It uses the sequence defined by $a_0=0$ and
        $a_{n+1}=\frac{1}{2\lfloor a_n\rfloor+1-a_n}$ and generates the
        sequence $$a_0,a_1,-a_1,a_2,-a_2,\ldots$$

        EXAMPLES:

        This example creates a list consisting of the first 10 terms
        generated by this function.

            sage: import itertools
            sage: [a for a in itertools.islice(Rationals(),10)]
            [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3/2]

        NOTES:
            A proof of the correctness of this formula is attributed to
            Sam Vandervelde and Don Zagier [A002487], but a better
            reference for the origin of this formula would be welcome.

            REFERENCES:
                 [A002487] Sloane's OLEIS,
                 http://www.research.att.com/~njas/sequences/A002487

        AUTHORS:
            - Nils Bruin (2007-02-20)
        """

        from sage.rings.arith import integer_floor as floor

        n=self(0)
        yield n
        while True:
          n=1/(2*floor(n)+1-n)
          yield n
          yield -n

    def complex_embedding(self, prec=53):
        import complex_field
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

    def is_subring(self, K):
        if K.is_field():
            return K.characteristic() == 0
        if K.characteristic() != 0:
            return False
        raise NotImplementedError

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

    def number_field(self, poly_var='x', nf_var='a'):
        from sage.rings.number_field.all import NumberField
        x = sage.rings.polynomial.polynomial_ring.PolynomialRing(self, poly_var).gen()
        return NumberField(x-1, nf_var)

    def order(self):
        """
        EXAMPLES:
            sage: QQ.order()
            +Infinity
        """
        return infinity.infinity

    def _an_element_impl(self):
        return sage.rings.rational.Rational((1,2))

    def random_element(self, num_bound=None, den_bound=None, distribution=None):
        """
        EXAMPLES:
            sage: QQ.random_element(10,10)
            -5/3
        """
        global ZZ
        if ZZ is None:
            import integer_ring
            ZZ = integer_ring.ZZ
        if num_bound == None:
            return self((ZZ.random_element(distribution=distribution),
                         ZZ.random_element(distribution=distribution)))
        else:
            if num_bound == 0:
                num_bound = 2
            if den_bound is None:
                den_bound = num_bound
                if den_bound < 1:
                    den_bound = 2
            return self((ZZ.random_element(-num_bound, num_bound+1, distribution=distribution),
                         ZZ.random_element(1, den_bound+1, distribution=distribution)))
    def zeta(self, n=2):
        """
        Return a root of unity in self.

        INPUT:
            n -- integer (default: 2) order of the root of unity

        EXAMPLES:
            sage: QQ.zeta()
            -1
            sage: QQ.zeta(2)
            -1
            sage: QQ.zeta(1)
            1
            sage: QQ.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: no n-th root of unity in rational field
        """
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

