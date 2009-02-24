r"""
Field `\mathbb{Q}` of Rational Numbers.

The class ``RationalField`` represents the field
`\mathbb{Q}` of (arbitrary precision) rational numbers.
Each rational number is an instance of the class
``Rational``.

Interactively, an instance of ``RationalField`` is
available as ``QQ``.

::

    sage: QQ
    Rational Field

Values of various types can be converted to rational numbers by
using the ``__call__`` method of
``RationalField`` (that is, by treating
``QQ`` as a function).

::

    sage: RealField(9).pi()
    3.1
    sage: QQ(RealField(9).pi())
    22/7
    sage: QQ(RealField().pi())
    245850922/78256779
    sage: QQ(35)
    35
    sage: QQ('12/347')
    12/347
    sage: QQ(exp(pi*I))
    -1
    sage: x = polygen(ZZ)
    sage: QQ((3*x)/(4*x))
    3/4

TEST::

    sage: Q = RationalField()
    sage: Q == loads(dumps(Q))
    True
    sage: RationalField() is RationalField()
    True
"""

import random
import field
import ring
import rational
import integer
import infinity
ZZ = None

from sage.structure.parent_gens import ParentWithGens
import sage.rings.number_field.number_field_base as number_field_base


_obj = {}
class _uniq(object):
    def __new__(cls):
        if _obj.has_key(0):
            return _obj[0]
        O = number_field_base.NumberField.__new__(cls)
        _obj[0] = O
        return O

class RationalField(_uniq, number_field_base.NumberField):
    r"""
    The class ``RationalField`` represents the field
    `\mathbb{Q}` of rational numbers.

    EXAMPLES::

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

    Coercion from the reals to the rational is done by default using
    continued fractions.

    ::

        sage: QQ(RR(3929329/32))
        3929329/32
        sage: QQ(-RR(3929329/32))
        -3929329/32
        sage: QQ(RR(1/7)) - 1/7
        0

    If you specify an optional second base argument, then the string
    representation of the float is used.

    ::

        sage: QQ(23.2, 2)
        6530219459687219/281474976710656
        sage: 6530219459687219.0/281474976710656
        23.199999999999999
        sage: a = 23.2; a
        23.2000000000000
        sage: QQ(a, 10)
        116/5

    Here's a nice example involving elliptic curves::

        sage: E = EllipticCurve('11a')
        sage: L = E.lseries().at1(300)[0]; L
        0.253841860855911
        sage: O = E.period_lattice().omega(); O
        1.26920930427955
        sage: t = L/O; t
        0.200000000000000
        sage: QQ(RealField(45)(t))
        1/5

    Elements from the extended rational field can be forced back into
    the rational field.

    ::

        sage: E = ExtendedRationalField
        sage: QQ(E(2))
        2
        sage: type(_)
        <type 'sage.rings.rational.Rational'>
    """

    def __init__(self):
        r"""
        We create the rational numbers `\mathbb{Q}`, and call a few
        functions::

            sage: Q = RationalField(); Q
            Rational Field
            sage: Q.characteristic()
            0
            sage: Q.is_field()
            True
            sage: Q.zeta()
            -1

        We next illustrate arithmetic in `\mathbb{Q}`.

        ::

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

        Next finally illustrate arithmetic with automatic coercion. The
        types that coerce into the rational field include ``str, int,
        long, Integer``.

        ::

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

        TESTS::

            sage: QQ.variable_name()
            'x'
            sage: QQ.variable_names()
            ('x',)
        """
        ParentWithGens.__init__(self, self)
        self._assign_names(('x',),normalize=False) # ???
        self._populate_coercion_lists_(element_constructor=rational.Rational, init_no_parent=True)

    def __hash__(self):
        return -11115808

    def _repr_(self):
        return "Rational Field"

    def _latex_(self):
        return "\\mathbf{Q}"

    def __len__(self):
        raise TypeError, 'len() of unsized object'

    def construction(self):
        from sage.categories.pushout import FractionField
        import integer_ring
        return FractionField(), integer_ring.ZZ

    def completion(self, p, prec, extras = {}):
        if p == infinity.Infinity:
            from sage.rings.real_mpfr import create_RealField
            return create_RealField(prec, **extras)
        else:
            from sage.rings.padics.factory import Qp
            return Qp(p, prec, **extras)

    def _coerce_map_from_(self, S):
        """
        EXAMPLES::

            sage: f = QQ.coerce_map_from(ZZ); f
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
            sage: f(3)
            3
            sage: f(3^99) - 3^99
            0
            sage: f = QQ.coerce_map_from(int); f
            Native morphism:
              From: Set of Python objects of type 'int'
              To:   Rational Field
            sage: f(44)
            44

        ::

            sage: QQ.coerce_map_from(long)
            Composite map:
              From: Set of Python objects of type 'long'
              To:   Rational Field
              Defn:   Native morphism:
                      From: Set of Python objects of type 'long'
                      To:   Integer Ring
                    then
                      Natural morphism:
                      From: Integer Ring
                      To:   Rational Field
        """
        global ZZ
        if ZZ is None:
            import integer_ring
            ZZ = integer_ring.ZZ
        if S is ZZ:
            return rational.Z_to_Q()
        elif S is int:
            return rational.int_to_Q()
        elif ZZ.has_coerce_map_from(S):
            return rational.Z_to_Q() * ZZ.coerce_map_from(S)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def __iter__(self):
        r"""
        Creates an iterator that generates the rational numbers without
        repetition, in order of the height.

        See also QQ.range_by_height?

        EXAMPLES:

        The first 17 rational numbers, ordered by height::

            sage: import itertools
            sage: lst = [a for a in itertools.islice(Rationals(),17)]
            sage: lst
            [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3, -3, 2/3, -2/3, 3/2, -3/2, 1/4, -1/4]
            sage: [a.height() for a in lst]
            [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4]
        """

        #The previous version of this function, implemented by Nils
        #Bruin, used the sequence defined by $a_0=0$ and
        #$a_{n+1}=\frac{1}{2\lfloor a_n\rfloor+1-a_n}$ and generated the
        #sequence $$a_0,a_1,-a_1,a_2,-a_2,\ldots$$.  This is [A002487]
        #in Sloane's encyclopedia, attributed to [Stern].  It is not
        #monotone in height, but has other interesting properties
        #described in [CalkinWilf].
        #REFERENCES:
        #  [A002487] Sloane's OLEIS,
        #    http://www.research.att.com/~njas/sequences/A002487
        #  [CalkinWilf] N. Calkin and H.S. Wilf, Recounting the
        #    rationals, American Mathematical Monthly 107 (2000),
        #    360--363
        #  [Stern] M.A. Stern, Ueber eine zahlentheoretische Funktion,
        #    Journal fuer die reine und angewandte Mathematik 55
        #    (1858), 193--220
        #
        # [beginning of Nils' code]
        #from sage.rings.arith import integer_floor as floor
        #
        #n=self(0)
        #yield n
        #while True:
        #  n=1/(2*floor(n)+1-n)
        #  yield n
        #  yield -n
        # [end of Nils' code]

        yield self(0)
        yield self(1)
        yield self(-1)
        height = integer.Integer(1)
        while True:
            height = height + 1
            for other in range(1, height):
                if height.gcd(other) == 1:
                    yield self(other/height)
                    yield self(-other/height)
                    yield self(height/other)
                    yield self(-height/other)

    def range_by_height(self, start, end=None):
        r"""
        Range function for rational numbers, ordered by height.

        Returns a Python generator for the list of rational numbers with
        heights in ``range(start, end)``. Follows the same
        convention as Python range, see range? for details.

        See also QQ.__iter__?

        EXAMPLES:

        All rational numbers with height strictly less than 4::

            sage: list(QQ.range_by_height(4))
            [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3, -3, 2/3, -2/3, 3/2, -3/2]
            sage: [a.height() for a in QQ.range_by_height(4)]
            [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3]

        All rational numbers with height 2::

            sage: list(QQ.range_by_height(2, 3))
            [1/2, -1/2, 2, -2]

        Nonsensical integer arguments will return an empty generator::

            sage: list(QQ.range_by_height(3, 3))
            []
            sage: list(QQ.range_by_height(10, 1))
            []

        There are no rational numbers with height `\leq 0`::

            sage: list(QQ.range_by_height(-10, 1))
            []
        """
        if end is None:
            end = start
            start = 1
        if start < 1:
            start = 1
        for height in ZZ.range(start, end):
            if height == 1:
                yield self(0)
                yield self(1)
                yield self(-1)
            for other in ZZ.range(1, height):
                if height.gcd(other) == 1:
                    yield self(other/height)
                    yield self(-other/height)
                    yield self(height/other)
                    yield self(-height/other)

    def discriminant(self):
        """
        Return the discriminant of the field of rational numbers, which is
        1.

        EXAMPLES::

            sage: QQ.discriminant()
            1
        """
        return sage.rings.integer.Integer(1)

    def class_number(self):
        """
        Return the class number of the field of rational numbers, which is
        1.

        EXAMPLES::

            sage: QQ.class_number()
            1
        """
        return integer.Integer(1)

    def signature(self):
        """
        Return the signature of the rational field, which is (1,0), since
        there are 1 real and no complex embeddings.

        EXAMPLES::

            sage: QQ.signature()
            (1, 0)
        """
        return (integer.Integer(1), integer.Integer(0))


    def discriminant(self):
        """
        Return the discriminant of the field of rational numbers, which is
        1.

        EXAMPLES::

            sage: QQ.discriminant()
            1
        """
        return integer.Integer(1)

    def signature(self):
        """
        Return the signature of the rational field, which is (1,0), since
        there are 1 real and no complex embeddings.

        EXAMPLES::

            sage: QQ.signature()
            (1, 0)
        """
        return (integer.Integer(1), integer.Integer(0))


    def embeddings(self, K):
        """
        Return list of the one embedding of `\mathbb{Q}` into
        `K`, if it exists.

        EXAMPLES::

            sage: QQ.embeddings(QQ)
            [Ring Coercion endomorphism of Rational Field]
            sage: QQ.embeddings(CyclotomicField(5))
            [Ring Coercion morphism:
              From: Rational Field
              To:   Cyclotomic Field of order 5 and degree 4]

        `K` must have characteristic 0::

            sage: QQ.embeddings(GF(3))
            Traceback (most recent call last):
            ...
            ValueError: no embeddings of the rational field into K.
        """
        if K.characteristic() != 0:
            raise ValueError, "no embeddings of the rational field into K."
        return [self.hom(K)]

    def complex_embedding(self, prec=53):
        """
        Return embedding of the rational numbers into the complex numbers.

        EXAMPLES::

            sage: QQ.complex_embedding()
            Ring morphism:
              From: Rational Field
              To:   Complex Field with 53 bits of precision
              Defn: 1 |--> 1.00000000000000
            sage: QQ.complex_embedding(20)
            Ring morphism:
              From: Rational Field
              To:   Complex Field with 20 bits of precision
              Defn: 1 |--> 1.0000
        """
        import complex_field
        CC = complex_field.ComplexField(prec)
        return self.hom([CC(1)])

    def gens(self):
        """
        EXAMPLES::

            sage: QQ.gens()
            (1,)
        """
        return (self(1), )

    def gen(self, n=0):
        """
        EXAMPLES::

            sage: QQ.gen()
            1
        """
        if n == 0:
            return self(1)
        else:
            raise IndexError, "n must be 0"

    def degree(self):
        """
        EXAMPLES::

            sage: QQ.degree()
            1
        """
        return integer.Integer(1)

    def absolute_degree(self):
        """
        EXAMPLES::

            sage: QQ.absolute_degree()
            1
        """
        return integer.Integer(1)

    def ngens(self):
        """
        EXAMPLES::

            sage: QQ.ngens()
            1
        """
        return integer.Integer(1)

    def is_absolute(self):
        """
        `\mathbb{Q}` is an absolute extension of
        `\mathbb{Q}`.

        EXAMPLES::

            sage: QQ.is_absolute()
            True
        """
        return True

    def is_subring(self, K):
        """
        Return ``True`` if `\mathbb{Q}` is a subring of
        `K`.

        We are only able to determine this in some cases, e.g., when
        `K` is a field or of positive characteristic.

        EXAMPLES::

            sage: QQ.is_subring(QQ)
            True
            sage: QQ.is_subring(QQ['x'])
            True
            sage: QQ.is_subring(GF(7))
            False
            sage: QQ.is_subring(CyclotomicField(7))
            True
            sage: QQ.is_subring(ZZ)
            False
            sage: QQ.is_subring(Frac(ZZ))
            True
        """
        if K.is_field():
            return K.characteristic() == 0
        if K.characteristic() != 0:
            return False
        try:
            self.embeddings(K)
        except (TypeError, ValueError):
            return False
        return True

    def is_field(self):
        """
        Return ``True``, since the rational field is a field.

        EXAMPLES::

            sage: QQ.is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Return ``False``, since the rational field is not
        finite.

        EXAMPLES::

            sage: QQ.is_finite()
            False
        """
        return False

    def is_prime_field(self):
        """
        Return ``True``, since `\mathbb{Q}` is a prime
        field.

        EXAMPLES::

            sage: QQ.is_prime_field()
            True
        """
        return True

    def is_atomic_repr(self):
        return True

    def characteristic(self):
        """
        Return 0, since the rational field has characteristic 0.

        EXAMPLES::

            sage: c = QQ.characteristic(); c
            0
            sage: parent(c)
            Integer Ring
        """
        return integer.Integer(0)

    def maximal_order(self):
        """
        Return the maximal order of the rational numbers, i.e., the ring
        `\mathbb{Z}` of integers.

        EXAMPLES::

            sage: QQ.maximal_order()
            Integer Ring
            sage: QQ.ring_of_integers ()
            Integer Ring
        """
        from integer_ring import ZZ
        return ZZ

    def number_field(self):
        """
        Return the number field associated to `\mathbb{Q}`. Since
        `\mathbb{Q}` is a number field, this just returns
        `\mathbb{Q}` again.

        EXAMPLES::

            sage: QQ.number_field() is QQ
            True
        """
        return self

    def power_basis(self):
        r"""
        Return a power basis for this number field over its base field.

        The power basis is always [1] for the rational field. This method
        is defined to make the rational field behave more like a number
        field.

        EXAMPLES::

            sage: QQ.power_basis()
            [1]
        """
        return [ self.gen() ]

    def extension(self, poly, names, check=True, embedding=None):
        """
        EXAMPLES:

        We make a single absolute extension::

            sage: K.<a> = QQ.extension(x^3 + 5); K
            Number Field in a with defining polynomial x^3 + 5

        We make an extension generated by roots of two polynomials::

            sage: K.<a,b> = QQ.extension([x^3 + 5, x^2 + 3]); K
            Number Field in a with defining polynomial x^3 + 5 over its base field
            sage: b^2
            -3
            sage: a^3
            -5
        """
        from sage.rings.number_field.all import NumberField
        return NumberField(poly, names=names, check=check, embedding=embedding)

    def order(self):
        """
        EXAMPLES::

            sage: QQ.order()
            +Infinity
        """
        return infinity.infinity

    def _an_element_(self):
        return rational.Rational((1,2))

    def random_element(self, num_bound=None, den_bound=None, distribution=None):
        """
        EXAMPLES::

            sage: QQ.random_element(10,10) # random output
            -5/3
        """
        global ZZ
        if ZZ is None:
            import integer_ring
            ZZ = integer_ring.ZZ
        if num_bound == None:
            num = ZZ.random_element(distribution=distribution)
            den = ZZ.random_element(distribution=distribution)
            while den == 0: den = ZZ.random_element(distribution=distribution)
            return self((num, den))
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
        Return a root of unity in ``self``.

        INPUT:


        -  ``n`` - integer (default: 2) order of the root of
           unity


        EXAMPLES::

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
            return rational.Rational(1)
        elif n == 2:
            return rational.Rational(-1)
        else:
            raise ValueError, "no n-th root of unity in rational field"

    #################################
    ## Coercions to interfaces
    #################################
    def _gap_init_(self):
        """
        EXAMPLES::

            sage: gap(QQ)
            Rationals
        """
        return 'Rationals'

    def _magma_init_(self, magma):
        """
        EXAMPLES::

            sage: magma(QQ)                       # optional - magma
            Rational Field
        """
        return 'RationalField()'

    def _macaulay2_init_(self):
        """
        EXAMPLES::

            sage: macaulay2(QQ)                   # optional- macaulay2
            QQ
        """
        return "QQ"

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES:
            sage: sage_input(QQ, verify=True)
            # Verified
            QQ
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: QQ._sage_input_(SageInputBuilder(), False)
            {atomic:QQ}
        """
        return sib.name('QQ')


QQ = RationalField()
Q = QQ

def is_RationalField(x):
    return isinstance(x, RationalField)

def frac(n,d):
    return rational.Rational(n)/rational.Rational(d)
