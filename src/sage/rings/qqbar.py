"""
Field of Algebraic Numbers

AUTHOR:
    -- Carl Witty (2007-01-27): initial version
    -- Carl Witty (2007-10-29): massive rewrite to support complex
         as well as real numbers

This is an implementation of the algebraic numbers (the complex
numbers which are the zero of a polynomial in ZZ[x]; in other words,
the algebraic closure of QQ, with an embedding into CC).  All
computations are exact.  We also include an implementation of the
algebraic reals (the intersection of the algebraic numbers with RR).
The field of algebraic numbers is available with abbreviation QQbar;
the field of algebraic reals has abbreviation AA.

As with many other implementations of the algebraic numbers, we try
hard to avoid computing a number field and working in the number
field; instead, we use floating-point interval arithmetic whenever
possible (basically whenever we need to prove non-equalities), and
resort to symbolic computation only as needed (basically to prove
equalities).

Algebraic numbers exist in one of the following forms:

* a rational number
* the product of a rational number and an n'th root of unity
* the sum, difference, product, or quotient of algebraic numbers
* the negation, inverse, absolute value, norm, real part,
imaginary part, or complex conjugate of an algebraic number
* a particular root of a polynomial, given as a polynomial with
algebraic coefficients, an isolating interval (given as a
RealIntervalFieldElement) which encloses exactly one root, and
the multiplicity of the root
* a polynomial in one generator, where the generator is an algebraic
number given as the root of an irreducible polynomial with integral
coefficients and the polynomial is given as a NumberFieldElement

The multiplicative subgroup of the algebraic numbers generated
by the rational numbers and the roots of unity is handled particularly
efficiently, as long as these roots of unity come from the QQbar.zeta()
method.  Cyclotomic fields in general are fairly efficient, again
as long as they are derived from QQbar.zeta().

An algebraic number can be coerced into ComplexIntervalField (or
RealIntervalField, for algebraic reals); every algebraic number has a
cached interval of the highest precision yet calculated.

Everything is done with intervals except for comparisons.  By default,
comparisons compute the two algebraic numbers with 128-bit precision
intervals; if this does not suffice to prove that the numbers are different,
then we fall back on exact computation.

Note that division involves an implicit comparison of the divisor against
zero, and may thus trigger exact computation.

Also, using an algebraic number in the leading coefficient of
a polynomial also involves an implicit comparison against zero, which
again may trigger exact computation.

Note that we work fairly hard to avoid computing new number fields;
to help, we keep a lattice of already-computed number fields and
their inclusions.

EXAMPLES:
    sage: sqrt(AA(2)) > 0
    True
    sage: (sqrt(5 + 2*sqrt(QQbar(6))) - sqrt(QQbar(3)))^2 == 2
    True
    sage: AA((sqrt(5 + 2*sqrt(6)) - sqrt(3))^2) == 2
    True

For a monic cubic polynomial x^3 + b*x^2 + c*x + d with roots
s1,s2,s3, the discriminant is defined as (s1-s2)^2(s1-s3)^2(s2-s3)^2
and can be computed as b^2c^2 - 4*b^3d - 4*c^3 + 18*bcd - 27*d^2.
We can test that these definitions do give the same result.

    sage: def disc1(b, c, d):
    ...       return b^2*c^2 - 4*b^3*d - 4*c^3 + 18*b*c*d - 27*d^2
    sage: def disc2(s1, s2, s3):
    ...       return ((s1-s2)*(s1-s3)*(s2-s3))^2
    sage: x = polygen(AA)
    sage: p = x*(x-2)*(x-4)
    sage: cp = AA.common_polynomial(p)
    sage: d, c, b, _ = p.list()
    sage: s1 = AA.polynomial_root(cp, RIF(-1, 1))
    sage: s2 = AA.polynomial_root(cp, RIF(1, 3))
    sage: s3 = AA.polynomial_root(cp, RIF(3, 5))
    sage: disc1(b, c, d) == disc2(s1, s2, s3)
    True
    sage: p = p + 1
    sage: cp = AA.common_polynomial(p)
    sage: d, c, b, _ = p.list()
    sage: s1 = AA.polynomial_root(cp, RIF(-1, 1))
    sage: s2 = AA.polynomial_root(cp, RIF(1, 3))
    sage: s3 = AA.polynomial_root(cp, RIF(3, 5))
    sage: disc1(b, c, d) == disc2(s1, s2, s3)
    True
    sage: p = (x-sqrt(AA(2)))*(x-AA(2).nth_root(3))*(x-sqrt(AA(3)))
    sage: cp = AA.common_polynomial(p)
    sage: d, c, b, _ = p.list()
    sage: s1 = AA.polynomial_root(cp, RIF(1.4, 1.5))
    sage: s2 = AA.polynomial_root(cp, RIF(1.7, 1.8))
    sage: s3 = AA.polynomial_root(cp, RIF(1.2, 1.3))
    sage: disc1(b, c, d) == disc2(s1, s2, s3)
    True

We can coerce from symbolic expressions:

    sage: QQbar(sqrt(-5))
    [2.2360679774997893 .. 2.2360679774997899]*I
    sage: AA(sqrt(2) + sqrt(3))
    [3.1462643699419721 .. 3.1462643699419726]
    sage: AA(sqrt(2)) + sqrt(3)
    [3.1462643699419721 .. 3.1462643699419726]
    sage: sqrt(2) + QQbar(sqrt(3))
    [3.1462643699419721 .. 3.1462643699419726]
    sage: QQbar(I)
    1*I
    sage: AA(I)
    Traceback (most recent call last):
    ...
    TypeError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
    sage: QQbar(I * golden_ratio)
    [1.6180339887498946 .. 1.6180339887498950]*I
    sage: AA(golden_ratio)^2 - AA(golden_ratio)
    1
    sage: QQbar((-8)^(1/3))
    [0.99999999999999988 .. 1.0000000000000003] + [1.7320508075688771 .. 1.7320508075688775]*I
    sage: AA((-8)^(1/3))
    -2
    sage: QQbar((-4)^(1/4))
    [1.0000000000000000 .. 1.0000000000000000] + [1.0000000000000000 .. 1.0000000000000000]*I
    sage: AA((-4)^(1/4))
    Traceback (most recent call last):
    ...
    TypeError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real

Note the different behavior in taking roots: for AA we prefer real roots if they exist, but for QQbar we take the principal root:

    sage: AA(-1)^(1/3)
    -1
    sage: QQbar(-1)^(1/3)
    [0.49999999999999994 .. 0.50000000000000012] + [0.86602540378443859 .. 0.86602540378443871]*I

We can explicitly coerce from QQ[I].  (Technically, this is not quite
kosher, since QQ[I] doesn't come with an embedding; we don't know
whether the field generator is supposed to map to +I or -I.  We assume
that for any quadratic field with polynomial x^2+1, the generator maps
to +I.)

    sage: K.<im> = QQ[I]
    sage: pythag = QQbar(3/5 + 4*im/5); pythag
    4/5*I + 3/5
    sage: pythag.abs() == 1
    True

However, implicit coercion from QQ[I] is not allowed.

    sage: QQbar(1) + im
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '+': 'Algebraic Field' and 'Number Field in I with defining polynomial x^2 + 1'

We can implicitly coerce from algebraic reals to algebraic numbers:

    sage: a = QQbar(1); print a, a.parent()
    1 Algebraic Field
    sage: b = AA(1); print b, b.parent()
    1 Algebraic Real Field
    sage: c = a + b; print c, c.parent()
    2 Algebraic Field

Some computation with radicals:

    sage: phi = (1 + sqrt(AA(5))) / 2
    sage: phi^2 == phi + 1
    True
    sage: tau = (1 - sqrt(AA(5))) / 2
    sage: tau^2 == tau + 1
    True
    sage: phi + tau == 1
    True
    sage: tau < 0
    True

    sage: rt23 = sqrt(AA(2/3))
    sage: rt35 = sqrt(AA(3/5))
    sage: rt25 = sqrt(AA(2/5))
    sage: rt23 * rt35 == rt25
    True

The Sage rings AA and QQbar can decide equalities between radical
expressions (over the reals and complex numbers respectively):

    sage: a = AA((2/(3*sqrt(3)) + 10/27)^(1/3) - 2/(9*(2/(3*sqrt(3)) + 10/27)^(1/3)) + 1/3)
    sage: a
    [0.99999999999999988 .. 1.0000000000000003]
    sage: a == 1
    True

Algebraic numbers which are known to be rational print as rationals;
otherwise they print as intervals (with 53-bit precision).

    sage: AA(2)/3
    2/3
    sage: QQbar(5/7)
    5/7
    sage: QQbar(1/3 - 1/4*I)
    -1/4*I + 1/3
    sage: two = QQbar(4).nth_root(4)^2; two
    [1.9999999999999997 .. 2.0000000000000005]
    sage: two == 2; two
    True
    2
    sage: phi
    [1.6180339887498946 .. 1.6180339887498950]

We can find the real and imaginary parts of an algebraic number (exactly).

    sage: r = QQbar.polynomial_root(x^5 - x - 1, CIF(RIF(0.1, 0.2), RIF(1.0, 1.1))); r
    [0.18123244446987538 .. 0.18123244446987541] + [1.0839541013177105 .. 1.0839541013177108]*I
    sage: r.real()
    [0.18123244446987538 .. 0.18123244446987541]
    sage: r.imag()
    [1.0839541013177105 .. 1.0839541013177108]
    sage: r.minpoly()
    x^5 - x - 1
    sage: r.real().minpoly()
    x^10 + 3/16*x^6 + 11/32*x^5 - 1/64*x^2 + 1/128*x - 1/1024
    sage: r.imag().minpoly() # this takes a long time (143s on my laptop)
    x^20 - 5/8*x^16 - 95/256*x^12 - 625/1024*x^10 - 5/512*x^8 - 1875/8192*x^6 + 25/4096*x^4 - 625/32768*x^2 + 2869/1048576

We can find the absolute value and norm of an algebraic number exactly.
(Note that we define the norm as the product of a number and its
complex conjugate; this is the algebraic definition of norm, if we
view QQbar as AA[I].)

    sage: r = QQbar((-8)^(1/3)); r
    [0.99999999999999988 .. 1.0000000000000003] + [1.7320508075688771 .. 1.7320508075688775]*I
    sage: r.abs() == 2
    True
    sage: r.norm() == 4
    True
    sage: (r+I).norm().minpoly()
    x^2 - 10*x + 13
    sage: r = AA.polynomial_root(x^2 - x - 1, RIF(-1, 0)); r
    [-0.61803398874989491 .. -0.61803398874989479]
    sage: r.abs().minpoly()
    x^2 + x - 1

We can compute the multiplicative order of an algebraic number.

    sage: QQbar(-1/2 + I*sqrt(3)/2).multiplicative_order()
    3
    sage: QQbar(-sqrt(3)/2 + I/2).multiplicative_order()
    12
    sage: QQbar.zeta(12345).multiplicative_order()
    12345

Cyclotomic fields are very fast as long as we only multiply and divide:

    sage: z3_3 = QQbar.zeta(3) * 3
    sage: z4_4 = QQbar.zeta(4) * 4
    sage: z5_5 = QQbar.zeta(5) * 5
    sage: z6_6 = QQbar.zeta(6) * 6
    sage: z20_20 = QQbar.zeta(20) * 20
    sage: z3_3 * z4_4 * z5_5 * z6_6 * z20_20
    7200

And they're still pretty fast even if you add and subtract, and trigger
exact computation:

    sage: (z3_3 + z4_4 + z5_5 + z6_6 + z20_20)._exact_value()
    4*zeta60^15 + 5*zeta60^12 + 9*zeta60^10 + 20*zeta60^3 - 3 where a^16 + a^14 - a^10 - a^8 - a^6 + a^2 + 1 = 0 and a in [0.99452189536827328 .. 0.99452189536827341] + [0.10452846326765345 .. 0.10452846326765352]*I

The paper _ARPREC: An Arbitrary Precision Computation Package_ discusses
this result.  Evidently it is difficult to find, but we can easily
verify it.

    sage: alpha = QQbar.polynomial_root(x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1, RIF(1, 1.2))
    sage: lhs = alpha^630 - 1
    sage: rhs_num = (alpha^315 - 1) * (alpha^210 - 1) * (alpha^126 - 1)^2 * (alpha^90 - 1) * (alpha^3 - 1)^3 * (alpha^2 - 1)^5 * (alpha - 1)^3
    sage: rhs_den = (alpha^35 - 1) * (alpha^15 - 1)^2 * (alpha^14 - 1)^2 * (alpha^5 - 1)^6 * alpha^68
    sage: rhs = rhs_num / rhs_den
    sage: lhs
    [2.6420403358193507e44 .. 2.6420403358193520e44]
    sage: rhs
    [2.6420403358193503e44 .. 2.6420403358193520e44]
    sage: lhs - rhs
    [-7.9344219392947342e28 .. 8.1800756658404269e28]
    sage: lhs == rhs
    True
    sage: lhs - rhs
    0
    sage: lhs._exact_value()
    -242494609856316402264822833062350847769474540*a^9 + 862295472068289472491654837785947906234680703*a^8 - 829559238431038252116584538075753012193290520*a^7 - 125882239615006638366472766103700441555126185*a^6 + 1399067970863104691667276008776398309383579345*a^5 - 1561176687069361567616835847286958553574223422*a^4 + 761706318888840943058230840550737823821027895*a^3 + 580740464974951394762758666210754821723780266*a^2 - 954587496403409756503464154898858512440951323*a + 546081123623099782018260884934770383777092602 where a^10 - 4*a^9 + 5*a^8 - a^7 - 6*a^6 + 9*a^5 - 6*a^4 - a^3 + 5*a^2 - 4*a + 1 = 0 and a in [0.44406334400909258 .. 0.44406334400909265]
"""

import sage.rings.ring
from sage.structure.sage_object import SageObject
from sage.structure.parent_gens import ParentWithGens
from sage.rings.real_mpfr import RR
from sage.rings.real_mpfi import RealIntervalField, RIF, RealIntervalFieldElement, is_RealIntervalFieldElement
from sage.rings.complex_field import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField, is_ComplexIntervalField
from sage.rings.complex_interval import ComplexIntervalFieldElement, is_ComplexIntervalFieldElement
from sage.rings.polynomial.all import PolynomialRing, is_Polynomial
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field import NumberField, NumberField_quadratic, QuadraticField, CyclotomicField
from sage.rings.number_field.number_field_element import is_NumberFieldElement
from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic
from sage.rings.arith import factor
from sage.libs.pari.gen import pari
from sage.structure.element import generic_power
import infinity
from sage.misc.functional import cyclotomic_polynomial

CC = ComplexField()
CIF = ComplexIntervalField()

# Singleton object implementation copied from integer_ring.py
_obj = None
class _uniq_alg(object):
    def __new__(cls):
        global _obj
        if _obj is None:
            _obj = sage.rings.ring.Field.__new__(cls)
        return _obj

_obj_r = None
class _uniq_alg_r(object):
    def __new__(cls):
        global _obj_r
        if _obj_r is None:
            _obj_r = sage.rings.ring.Field.__new__(cls)
        return _obj_r

is_SymbolicExpressionRing = None

def late_import():
    global is_SymbolicExpressionRing
    if is_SymbolicExpressionRing is None:
        import sage.calculus.calculus
        is_SymbolicExpressionRing = sage.calculus.calculus.is_SymbolicExpressionRing

class AlgebraicField_common(sage.rings.ring.Field):
    def default_interval_prec(self):
        return 64

    def is_finite(self):
        return False

    def characteristic(self):
        return sage.rings.integer.Integer(0)

    def order(self):
        return infinity.infinity

    def common_polynomial(self, poly):
        """
        Given a polynomial with algebraic coefficients, returns a
        wrapper that caches high-precision calculations and
        factorizations.  This wrapper can be passed to polynomial_root
        in place of the polynomial.

        Using common_polynomial makes no semantic difference, but will
        improve efficiency if you are dealing with multiple roots
        of a single polynomial.

        EXAMPLES:
            sage: x = polygen(ZZ)
            sage: p = AA.common_polynomial(x^2 - x - 1)
            sage: phi = AA.polynomial_root(p, RIF(1, 2))
            sage: tau = AA.polynomial_root(p, RIF(-1, 0))
            sage: phi + tau == 1
            True
            sage: phi * tau == -1
            True

            sage: x = polygen(SR)
            sage: p = (x - sqrt(-5)) * (x - sqrt(3)); p
            x^2 + ((1 - sqrt(3))*(1 - sqrt(5)*I) - sqrt(3)*sqrt(5)*I - 1)*x + sqrt(3)*sqrt(5)*I
            sage: p = QQbar.common_polynomial(p)
            sage: a = QQbar.polynomial_root(p, CIF(RIF(-0.1, 0.1), RIF(2, 3))); a
            [-6.1814409042940318e-19 .. 6.0853195377215860e-19] + [2.2360679774997893 .. 2.2360679774997899]*I
            sage: b = QQbar.polynomial_root(p, RIF(1, 2)); b
            [1.7320508075688771 .. 1.7320508075688775]

        These "common polynomials" can be shared between real and
        complex roots:
             sage: p = AA.common_polynomial(x^3 - x - 1)
             sage: r1 = AA.polynomial_root(p, RIF(1.3, 1.4)); r1
             [1.3247179572447458 .. 1.3247179572447461]
             sage: r2 = QQbar.polynomial_root(p, CIF(RIF(-0.7, -0.6), RIF(0.5, 0.6))); r2
             [-0.66235897862237303 .. -0.66235897862237291] + [0.56227951206230120 .. 0.56227951206230132]*I
        """
        return AlgebraicPolynomialTracker(poly)

class AlgebraicRealField(_uniq_alg_r, AlgebraicField_common):
    r"""
    The field of algebraic reals.
    """

    def __init__(self):
        ParentWithGens.__init__(self, self, ('x',), normalize=False)

    def __call__(self, x):
        """
        Coerce x into the field of algebraic real numbers.

        """

        if isinstance(x, AlgebraicReal):
            return x
        elif isinstance(x, AlgebraicNumber):
            if x.imag().is_zero():
                return x.real()
            else:
                raise TypeError, "Cannot coerce algebraic number with non-zero imaginary part to algebraic real"
        elif hasattr(x, '_algebraic_'):
            return x._algebraic_(AA)
        return AlgebraicReal(x)

    def _repr_(self):
        return "Algebraic Real Field"

    # Is there a standard representation for this?
    def _latex_(self):
        return "\\mathbf{A}"

    def _coerce_impl(self, x):
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            return self(x)
        elif hasattr(x, '_algebraic_'):
            return x._algebraic_(AA)
        raise TypeError, 'no implicit coercion of element to the algebraic numbers'

    def has_coerce_map_from_impl(self, from_par):
        if from_par == ZZ or from_par == QQ or from_par == int or from_par == long:
            return True
        if from_par == AA:
            return True
        late_import()
        if is_SymbolicExpressionRing(from_par):
            return True
        return False

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def gens(self):
        return (self(1), )

    def gen(self, n=0):
        if n == 0:
            return self(1)
        else:
            raise IndexError, "n must be 0"

    def ngens(self):
        return 1

    def is_atomic_repr(self):
        return True

    def zeta(self, n=2):
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        else:
            raise ValueError, "no n-th root of unity in algebraic reals"

    def polynomial_root(self, poly, interval, multiplicity=1):
        """
        Given a polynomial with algebraic coefficients and an interval
        enclosing exactly one root of the polynomial, constructs
        an algebraic real representation of that root.

        The polynomial need not be irreducible, or even squarefree; but
        if the given root is a multiple root, its multiplicity must be
        specified.  (IMPORTANT NOTE: Currently, multiplicity-k roots
        are handled by taking the (k-1)th derivative of the polynomial.
        This means that the interval must enclose exactly one root
        of this derivative.)

        The conditions on the arguments (that the interval encloses exactly
        one root, and that multiple roots match the given multiplicity)
        are not checked; if they are not satisfied, an error may be
        thrown (possibly later, when the algebraic number is used),
        or wrong answers may result.

        Note that if you are constructing multiple roots of a single
        polynomial, it is better to use AA.common_polynomial (or
        QQbar.common_polynomial; the two are equivalent) to get a
        shared polynomial.

        EXAMPLES:
            sage: x = polygen(AA)
            sage: phi = AA.polynomial_root(x^2 - x - 1, RIF(1, 2)); phi
            [1.6180339887498946 .. 1.6180339887498950]
            sage: p = (x-1)^7 * (x-2)
            sage: r = AA.polynomial_root(p, RIF(9/10, 11/10), multiplicity=7)
            sage: r; r == 1
            [0.99999999999999988 .. 1.0000000000000003]
            True
            sage: p = (x-phi)*(x-sqrt(AA(2)))
            sage: r = AA.polynomial_root(p, RIF(1, 3/2))
            sage: r; r == sqrt(AA(2))
            [1.4142135623730949 .. 1.4142135623730952]
            True

        We allow complex polynomials, as long as the particular root
        in question is real.
            sage: K.<im> = QQ[I]
            sage: x = polygen(K)
            sage: p = (im + 1) * (x^3 - 2); p
            (I + 1)*x^3 - 2*I - 2
            sage: r = AA.polynomial_root(p, RIF(1, 2)); r^3
            [1.9999999999999997 .. 2.0000000000000005]
        """
        if not is_RealIntervalFieldElement(interval):
            raise ValueError, "interval argument of .polynomial_root on algebraic real field must be real"

        return AlgebraicReal(ANRoot(poly, interval, multiplicity))

def is_AlgebraicRealField(F):
    return isinstance(F, AlgebraicRealField)

AA = AlgebraicRealField()

class AlgebraicField(_uniq_alg, AlgebraicField_common):
    """
    The field of algebraic numbers.
    """

    def __init__(self):
        ParentWithGens.__init__(self, AA, ('I',), normalize=False)

    def __call__(self, x):
        """
        Coerce x into the field of algebraic numbers.

        """

        if isinstance(x, AlgebraicNumber):
            return x
        elif isinstance(x, AlgebraicReal):
            return AlgebraicNumber(x._descr)
        elif hasattr(x, '_algebraic_'):
            return x._algebraic_(QQbar)
        return AlgebraicNumber(x)

    def _repr_(self):
        return "Algebraic Field"

    def _coerce_impl(self, x):
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            return self(x)
        elif hasattr(x, '_algebraic_'):
            return x._algebraic_(QQbar)
        elif isinstance(x, AlgebraicReal):
            return AlgebraicNumber(x._descr)
        raise TypeError, 'no implicit coercion of element to the algebraic numbers'

    def has_coerce_map_from_impl(self, from_par):
        if from_par == ZZ or from_par == QQ or from_par == int or from_par == long:
            return True
        if from_par == AA or from_par == QQbar:
            return True
        late_import()
        if is_SymbolicExpressionRing(from_par):
            return True
        return False

    def gens(self):
        return(QQbar_I, )

    def gen(self, n=0):
        if n == 0:
            return QQbar_I
        else:
            raise IndexError, "n must be 0"

    def ngens(self):
        return 1

    def is_atomic_repr(self):
        return False

    def zeta(self, n=4):
        """
        Returns a primitive n'th root of unity.  (In fact, returns
        exp(2*pi*I/n).)

        EXAMPLES:
            sage: QQbar.zeta(1)
            1
            sage: QQbar.zeta(2)
            -1
            sage: QQbar.zeta(3)
            [-0.50000000000000012 .. -0.49999999999999994] + [0.86602540378443859 .. 0.86602540378443871]*I
            sage: QQbar.zeta(4)
            1*I
            sage: QQbar.zeta()
            1*I
            sage: QQbar.zeta(5)
            [0.30901699437494739 .. 0.30901699437494746] + [0.95105651629515353 .. 0.95105651629515365]*I
            sage: QQbar.zeta(314159)
            [0.99999999979999965 .. 0.99999999979999977] + [0.000020000016891958231 .. 0.000020000016891958236]*I
        """
        return AlgebraicNumber(ANRootOfUnity(QQ_1/n, QQ_1))

    def polynomial_root(self, poly, interval, multiplicity=1):
        """
        Given a polynomial with algebraic coefficients and an interval
        enclosing exactly one root of the polynomial, constructs
        an algebraic real representation of that root.

        The polynomial need not be irreducible, or even squarefree; but
        if the given root is a multiple root, its multiplicity must be
        specified.  (IMPORTANT NOTE: Currently, multiplicity-k roots
        are handled by taking the (k-1)th derivative of the polynomial.
        This means that the interval must enclose exactly one root
        of this derivative.)

        The conditions on the arguments (that the interval encloses exactly
        one root, and that multiple roots match the given multiplicity)
        are not checked; if they are not satisfied, an error may be
        thrown (possibly later, when the algebraic number is used),
        or wrong answers may result.

        Note that if you are constructing multiple roots of a single
        polynomial, it is better to use QQbar.common_polynomial
        to get a shared polynomial.

        EXAMPLES:
            sage: x = polygen(QQbar)
            sage: phi = QQbar.polynomial_root(x^2 - x - 1, RIF(0, 2)); phi
            [1.6180339887498946 .. 1.6180339887498950]
            sage: p = (x-1)^7 * (x-2)
            sage: r = QQbar.polynomial_root(p, RIF(9/10, 11/10), multiplicity=7)
            sage: r; r == 1
            [1.0000000000000000 .. 1.0000000000000000]
            True
            sage: p = (x-phi)*(x-sqrt(QQbar(2)))
            sage: r = QQbar.polynomial_root(p, RIF(1, 3/2))
            sage: r; r == sqrt(QQbar(2))
            [1.4142135623730949 .. 1.4142135623730952]
            True
        """
        return AlgebraicNumber(ANRoot(poly, interval, multiplicity))

def is_AlgebraicField(F):
    return isinstance(F, AlgebraicField)

QQbar = AlgebraicField()

def is_AlgebraicField_common(F):
    return isinstance(F, AlgebraicField_common)

def prec_seq():
    # XXX Should do some testing to see where the efficiency breaks are
    # in MPFR.  We could also test variants like "bits = bits + bits // 2"
    # (I think this is what MPFR uses internally).
    bits = 64
    while True:
        yield bits
        bits = bits * 2

_short_prec_seq = (64, 128, None)
def short_prec_seq(): return _short_prec_seq

def tail_prec_seq():
    bits = 256
    while True:
        yield bits
        bits = bits * 2

def rational_exact_root(r, d):
    """
    Checks whether the rational r is an exact d'th power.  If so, returns
    the d'th root of r; otherwise, returns None.

    EXAMPLES:
        sage: from sage.rings.qqbar import rational_exact_root
        sage: rational_exact_root(16/81, 4)
        2/3
        sage: rational_exact_root(8/81, 3) is None
        True
    """
    num = r.numerator()
    den = r.denominator()

    (num_rt, num_exact) = num.nth_root(d, report_exact=True)
    if not num_exact: return None
    (den_rt, den_exact) = den.nth_root(d, report_exact=True)
    if not den_exact: return None
    return (num_rt / den_rt)

def clear_denominators(poly):
    """
    Takes a monic polynomial and rescales the variable to get a monic
    polynomial with "integral" coefficients.  Works on any univariate
    polynomial whose base ring has a denominator() method that returns
    integers; for example, the base ring might be QQ or a number
    field.

    Returns the scale factor and the new polynomial.

    (Inspired by Pari's primitive_pol_to_monic().)

    We assume that coefficient denominators are "small"; the algorithm
    factors the denominators, to give the smallest possible scale factor.

    EXAMPLES:
        sage: from sage.rings.qqbar import clear_denominators

        sage: _.<x> = QQ['x']
        sage: clear_denominators(x + 3/2)
        (2, x + 3)
        sage: clear_denominators(x^2 + x/2 + 1/4)
        (2, x^2 + x + 1)

    """

    # This algorithm factors the polynomial denominators.
    # We should check the size of the denominators and switch to
    # an alternate, less precise algorithm if we decide factoring
    # would be too slow.

    d = poly.denominator()
    if d == 1:
        return d, poly
    deg = poly.degree()
    factors = {}
    for i in range(deg):
        d = poly[i].denominator()
        df = factor(d)
        for f, e in df:
            oe = 0
            if f in factors:
                oe = factors[f]
            min_e = (e + (deg-i) - 1) // (deg-i)
            factors[f] = max(oe, min_e)
    change = 1
    for f, e in factors.iteritems():
        change = change * f**e
    poly = poly * (change ** deg)
    poly = poly(poly.parent().gen() / change)
    return change, poly

def do_polred(poly):
    """
    Find the polynomial of lowest discriminant that generates
    the same field as poly, out of those returned by the Pari
    polred routine.  Returns a triple: (elt_fwd, elt_back, new_poly),
    where new_poly is the new polynomial, elt_fwd is a polynomial
    expression for a root of the new polynomial in terms of a root
    of the original polynomial, and elt_back is a polynomial expression
    for a root of the original polynomial in terms of a root of the
    new polynomial.

    EXAMPLES:
        sage: from sage.rings.qqbar import do_polred

        sage: _.<x> = QQ['x']
        sage: do_polred(x^2-5)
        (-1/2*x + 1/2, -2*x + 1, x^2 - x - 1)
        sage: do_polred(x^2-x-11)
        (-1/3*x + 2/3, -3*x + 2, x^2 - x - 1)
        sage: do_polred(x^3 + 123456)
        (-1/4*x, -4*x, x^3 - 1929)
    """
    degree = poly.degree()
    pari_poly = pari(poly)

    red_table = pari_poly.polred(3)

    best = None
    best_discr = None

    for i in range(red_table.nrows()):
        red_poly = red_table[i,1]
        if red_poly.poldegree() < degree:
            continue
        red_discr = red_poly.poldisc().abs()
        if best_discr is None or red_discr < best_discr:
            best = red_poly
            best_discr = red_discr
            best_elt = red_table[i,0]

    assert(best is not None)
    parent = poly.parent()
    rev = parent(best_elt.Mod(pari_poly).modreverse().lift())
    return parent(best_elt), rev, parent(best)

def isolating_interval(intv_fn, pol):
    """
    intv_fn is a function that takes a precision and returns an
    interval of that precision containing some particular root of pol.
    (It must return better approximations as the precision increases.)
    pol is an irreducible polynomial with rational coefficients.

    Returns an interval containing at most one root of pol.

    EXAMPLES:
        sage: from sage.rings.qqbar import isolating_interval

        sage: _.<x> = QQ['x']
        sage: isolating_interval(lambda prec: sqrt(RealIntervalField(prec)(2)), x^2 - 2)
        [1.41421356237309504876 .. 1.41421356237309504888]

    And an example that requires more precision:
        sage: delta = 10^(-70)
        sage: p = (x - 1) * (x - 1 - delta) * (x - 1 + delta)
        sage: isolating_interval(lambda prec: RealIntervalField(prec)(1 + delta), p)
        [1.00000000000000000000000000000000000000000000000000000000000000000000009999999999999999999999999999999999999999999999999999999999999999999999999999999999998 .. 1.00000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000014]

    The function also works with complex intervals and complex roots:
        sage: p = x^2 - x + 13/36
        sage: isolating_interval(lambda prec: ComplexIntervalField(prec)(1/2, 1/3), p)
        [0.500000000000000000000 .. 0.500000000000000000000] + [0.333333333333333333315 .. 0.333333333333333333343]*I
    """
    dpol = pol.derivative()

    for prec in prec_seq():
        intv = intv_fn(prec)
        ifld = intv.parent()

        # We need to verify that pol has exactly one root in the
        # interval intv.  We know (because it is a precondition of
        # calling this function) that it has at least one root in the
        # interval, so we only need to verify that it has at most one
        # root (that the interval is sufficiently narrow).

        # We do this by computing the derivative of the polynomial
        # over the interval.  If the derivative is bounded away from zero,
        # then we know there can be at most one root.

        if not dpol(intv).contains_zero():
            return intv

def find_zero_result(fn, l):
    """
    l is a list of some sort.  fn is a function which maps an element
    of l and a precision into an interval (either real or complex) of
    that precision, such that for sufficient precision, exactly one
    element of l results in an interval containing 0.  Returns that
    one element of l.

    EXAMPLES:
        sage: from sage.rings.qqbar import find_zero_result
        sage: _.<x> = QQ['x']
        sage: delta = 10^(-70)
        sage: p1 = x - 1
        sage: p2 = x - 1 - delta
        sage: p3 = x - 1 + delta
        sage: p2 == find_zero_result(lambda p, prec: p(RealIntervalField(prec)(1 + delta)), [p1, p2, p3])
        True
    """
    for prec in prec_seq():
        result = None
        ambig = False
        for v in l:
            intv = fn(v, prec)
            if intv.contains_zero():
                if result is not None:
                    ambig = True
                    break
                result = v
        if ambig:
            continue
        if result is None:
            raise ValueError, 'find_zero_result could not find any zeroes'
        return result

def conjugate_expand(v):
    """
    If the interval v (which may be real or complex) includes some
    purely real numbers, return v' containing v such that
    v' == v'.conjugate().  Otherwise return v unchanged.  (Note that if
    v' == v'.conjugate(), and v' includes one non-real root of a real
    polynomial, then v' also includes the conjugate of that root.
    Also note that the diameter of the return value is at most twice
    the diameter of the input.)

    EXAMPLES:
        sage: from sage.rings.qqbar import conjugate_expand
        sage: conjugate_expand(CIF(RIF(0, 1), RIF(1, 2)))
        [0.00000000000000000 .. 1.0000000000000000] + [1.0000000000000000 .. 2.0000000000000000]*I
        sage: conjugate_expand(CIF(RIF(0, 1), RIF(0, 1)))
        [0.00000000000000000 .. 1.0000000000000000] + [-1.0000000000000000 .. 1.0000000000000000]*I
        sage: conjugate_expand(CIF(RIF(0, 1), RIF(-2, 1)))
        [0.00000000000000000 .. 1.0000000000000000] + [-2.0000000000000000 .. 2.0000000000000000]*I
        sage: conjugate_expand(RIF(1, 2))
        [1.0000000000000000 .. 2.0000000000000000]
    """
    if is_RealIntervalFieldElement(v):
        return v
    im = v.imag()
    if not im.contains_zero():
        return v
    re = v.real()
    fld = ComplexIntervalField(v.prec())
    return fld(re, im.union(-im))

def conjugate_shrink(v):
    """
    If the interval v includes some purely real numbers, return
    a real interval containing only those real numbers.  Otherwise
    return v unchanged.

    If v includes exactly one root of a real polynomial, and v was
    returned by conjugate_expand(), then conjugate_shrink(v) still
    includes that root, and is a RealIntervalFieldElement iff the root
    in question is real.

    EXAMPLES:
        sage: from sage.rings.qqbar import conjugate_shrink
        sage: conjugate_shrink(RIF(3, 4))
        [3.0000000000000000 .. 4.0000000000000000]
        sage: conjugate_shrink(CIF(RIF(1, 2), RIF(1, 2)))
        [1.0000000000000000 .. 2.0000000000000000] + [1.0000000000000000 .. 2.0000000000000000]*I
        sage: conjugate_shrink(CIF(RIF(1, 2), RIF(0, 1)))
        [1.0000000000000000 .. 2.0000000000000000]
        sage: conjugate_shrink(CIF(RIF(1, 2), RIF(-1, 2)))
        [1.0000000000000000 .. 2.0000000000000000]
    """
    if is_RealIntervalFieldElement(v):
        return v
    im = v.imag()
    if im.contains_zero():
        return v.real()
    return v

# Cache some commonly-used polynomial rings
QQx = QQ['x']
QQx_x = QQx.gen()
QQy = QQ['y']
QQy_y = QQy.gen()
QQxy = QQ['x', 'y']
QQxy_x = QQxy.gen(0)
QQxy_y = QQxy.gen(1)

class AlgebraicGeneratorRelation(SageObject):
    """
    A simple class for maintaining relations in the lattice of algebraic
    extensions.
    """
    def __init__(self, child1, child1_poly, child2, child2_poly, parent):
        self.child1 = child1
        self.child1_poly = child1_poly
        self.child2 = child2
        self.child2_poly = child2_poly
        self.parent = parent

algebraic_generator_counter = 0

class AlgebraicGenerator(SageObject):
    """
    An AlgebraicGenerator represents both an algebraic number alpha and
    the number field QQ[alpha].  There is a single AlgebraicGenerator
    representing QQ (with alpha==0).

    The AlgebraicGenerator class is private, and should not be used
    directly.
    """

    def __init__(self, field, root):
        """
        Construct an AlgebraicGenerator object.

        sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
        sage: _.<y> = QQ['y']
        sage: x = polygen(QQbar)
        sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
        sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
        sage: gen = AlgebraicGenerator(nf, root)
        sage: gen
        Number Field in a with defining polynomial y^2 - y - 1 with a in [1.6180339887498946 .. 1.6180339887498950]
        sage: gen.field()
        Number Field in a with defining polynomial y^2 - y - 1
        sage: gen.is_trivial()
        False
        sage: gen.union(qq_generator) is gen
        True
        sage: qq_generator.union(gen) is gen
        True
        sage: nf = NumberField(y^2 + 1, name='a', check=False)
        sage: root = ANRoot(x^2 + 1, CIF(0, 1))
        sage: x = AlgebraicGenerator(nf, root); x
        Number Field in a with defining polynomial y^2 + 1 with a in [1.0000000000000000 .. 1.0000000000000000]*I
        """
        self._field = field
        self._pari_field = None
        self._trivial = (field is None)
        self._root = root
        self._unions = {}
        self._cyclotomic = False
        global algebraic_generator_counter
        self._index = algebraic_generator_counter
        algebraic_generator_counter += 1

    def __hash__(self):
        return self._index

    def __cmp__(self, other):
        return cmp(self._index, other._index)

    def set_cyclotomic(self, n):
        self._cyclotomic = True
        self._cyclotomic_order = ZZ(n)

    def is_complex(self):
        return self._root.is_complex()

    def _repr_(self):
        if self._trivial:
            return 'Trivial generator'
        else:
            if isinstance(self._root, ANRootOfUnity):
                return str(self._root)
            else:
                return '%s with a in %s'%(self._field, self._root._interval_fast(53))

    def is_trivial(self):
        """
        Returns true iff this is the trivial generator (alpha == 0), which
        does not actually extend the rationals.

        EXAMPLES:
            sage: from sage.rings.qqbar import qq_generator
            sage: qq_generator.is_trivial()
            True
        """
        return self._trivial

    def field(self):
        return self._field

    def pari_field(self):
        if self._pari_field is None:
            pari_pol = self._field.pari_polynomial().subst('x','y')
            self._pari_field = pari_pol.nfinit(1)
        return self._pari_field

    def conjugate(self):
        """
        If this generator is for the algebraic number alpha, return a
        generator for the complex conjugate of alpha.
        """
        try:
            return self._conjugate
        except AttributeError:
            if not self.is_complex():
                self._conjugate = self
            else:
                conj = AlgebraicGenerator(self._field, self._root.conjugate(None))
                self._conjugate = conj
                conj._conjugate = self
            if self._cyclotomic:
                conj_rel = QQx_x ** (self._cyclotomic_order - 1)
                rel = AlgebraicGeneratorRelation(self, conj_rel, conj, conj_rel, self)
                self._unions[conj] = rel
                conj._unions[self] = rel
            return self._conjugate

    def _interval_fast(self, prec):
        """
        Returns an interval containing this generator, to the specified
        precision.
        """
        return self._root._interval_fast(prec)

    def union(self, other):
        """
        Given generators alpha and beta, alpha.union(beta) gives a generator
        for the number field QQ[alpha][beta].

        EXAMPLES:
            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in [1.4142135623730949 .. 1.4142135623730952]
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = ANRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in [1.7320508075688771 .. 1.7320508075688775]
            sage: gen2.union(qq_generator) is gen2
            True
            sage: qq_generator.union(gen3) is gen3
            True
            sage: gen2.union(gen3)
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in [0.51763809020504147 .. 0.51763809020504159]
        """
        if self._trivial:
            return other
        if other._trivial:
            return self
        if self is other:
            return self
        if other in self._unions:
            return self._unions[other].parent
        if self._field.polynomial().degree() < other._field.polynomial().degree():
            self, other = other, self
        elif other._cyclotomic:
            self, other = other, self

        if self._cyclotomic and other._cyclotomic:
            parent_order = self._cyclotomic_order.lcm(other._cyclotomic_order)
            new_gen = cyclotomic_generator(parent_order)
            rel = AlgebraicGeneratorRelation(self, QQx_x ** (parent_order // self._cyclotomic_order),
                                             other, QQx_x ** (parent_order // other._cyclotomic_order),
                                             new_gen)
            self._unions[other] = rel
            other._unions[self] = rel
            return new_gen

        sp = self._field.polynomial()
        op = other._field.polynomial()
        op = QQx(op)
        # print sp
        # print op
        # print self._field.polynomial()
        # print self._field.polynomial().degree()
        # pari_nf = self._field.pari_nf()
        pari_nf = self.pari_field()
        # print pari_nf[0]
        factors = list(pari_nf.nffactor(op).lift())[0]
        # print factors
        x, y = QQxy.gens()
        factors_sage = [QQxy(p) for p in factors]
        # print factors_sage
        def find_fn(p, prec):
            ifield = RealIntervalField(prec)
            if_poly = ifield['x', 'y']
            ip = if_poly(p)
            return ip(other._root._interval_fast(prec), self._root._interval_fast(prec))
        my_factor = find_zero_result(find_fn, factors_sage)

        if my_factor.degree(x) == 1 and my_factor.coefficient(x) == 1:
            value = (-my_factor + x).univariate_polynomial(QQy)
            # print value
            rel = AlgebraicGeneratorRelation(self, QQy_y,
                                             other, value,
                                             self)
            self._unions[other] = rel
            other._unions[self] = rel
            return rel.parent

        # XXX need more caching here
        newpol, self_pol, k = pari_nf.rnfequation(my_factor, 1)
        k = int(k)
        # print newpol
        # print self_pol
        # print k

        newpol_sage = QQx(newpol)
        newpol_sage_y = QQy(newpol_sage)

        red_elt, red_back, red_pol = do_polred(newpol_sage_y)

        red_back_x = QQx(red_back)

        new_nf = NumberField(red_pol, name='a', check=False)

        self_pol_sage = QQx(self_pol.lift())

        new_nf_a = new_nf.gen()

        def intv_fn(prec):
            return conjugate_expand(red_elt(self._root._interval_fast(prec) * k + other._root._interval_fast(prec)))
        new_intv = conjugate_shrink(isolating_interval(intv_fn, red_pol))

        new_gen = AlgebraicGenerator(new_nf, ANRoot(QQx(red_pol), new_intv))
        rel = AlgebraicGeneratorRelation(self, self_pol_sage(red_back_x),
                                         other, (QQx_x - k*self_pol_sage)(red_back_x),
                                         new_gen)
        self._unions[other] = rel
        other._unions[self] = rel
        return new_gen

    def super_poly(self, super, checked=None):
        """
        Given a generator gen and another generator super, where super
        is the result of a tree of union() operations where one of the
        leaves is gen, gen.super_poly(super) returns a polynomial
        expressing the value of gen in terms of the value of super.
        (Except that if gen is qq_generator, super_poly() always
        returns None.)

        EXAMPLES:
            sage: from sage.rings.qqbar import AlgebraicGenerator, ANRoot, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in [1.4142135623730949 .. 1.4142135623730952]
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = ANRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in [1.7320508075688771 .. 1.7320508075688775]
            sage: gen2_3 = gen2.union(gen3)
            sage: gen2_3
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in [0.51763809020504147 .. 0.51763809020504159]
            sage: qq_generator.super_poly(gen2) is None
            True
            sage: gen2.super_poly(gen2_3)
            -a^3 + 3*a
            sage: gen3.super_poly(gen2_3)
            -a^2 + 2

        """
        if checked is None:
            checked = {}
        checked[self] = True
        if super is self:
            return self._field.gen()
        for u in self._unions.values():
            # print self, u.parent
            # print checked
            if u.parent in checked:
                continue
            poly = u.parent.super_poly(super, checked)
            # print poly, u.child1, u.child2
            if poly is None:
                continue
            if self is u.child1:
                return u.child1_poly(poly)
            assert(self is u.child2)
            return u.child2_poly(poly)
        return None

    def __call__(self, elt):
        """
        Takes an algebraic number which is represented as either a rational
        or a number field element, and which is in a subfield of the
        field generated by this generator.  Lifts the number into the
        field of this generator, and returns either a Rational or a
        NumberFieldElement depending on whether this is the trivial generator.

        EXAMPLES:
            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, ANExtensionElement, ANRational
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in [1.4142135623730949 .. 1.4142135623730952]
            sage: sqrt2 = ANExtensionElement(gen2, nf2.gen())
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = ANRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in [1.7320508075688771 .. 1.7320508075688775]
            sage: sqrt3 = ANExtensionElement(gen3, nf3.gen())
            sage: gen2_3 = gen2.union(gen3)
            sage: gen2_3
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in [0.51763809020504147 .. 0.51763809020504159]
            sage: gen2_3(sqrt2)
            -a^3 + 3*a
            sage: gen2_3(ANRational(1/7))
            1/7
            sage: gen2_3(sqrt3)
            -a^2 + 2
        """
        if self._trivial:
            return elt.rational_value()
        if elt.is_rational():
            return self._field(elt.rational_value())
        if elt.generator() is self:
            return elt.field_element_value()
        gen = elt.generator()
        sp = gen.super_poly(self)
        # print gen
        # print self
        assert(not(sp is None))
        # print sp
        # print self._field
        return self._field(elt.field_element_value().polynomial()(sp))

# These are the functions used to add, subtract, multiply, and divide
# algebraic numbers.  Basically, we try to compute exactly if the
# result would be a Gaussian rational, or a rational times a root
# of unity; or if both arguments are already known to be in the same
# number field.  Otherwise we fall back to floating-point computation,
# to be backed up by exact symbolic computation only as required.

# These choices are motivated partly by efficiency considerations
# (not backed up by benchmarks, so other possibilities might be more
# efficient), and partly by concerns for the prettiness of output:
# we want algebraic numbers to print as Gaussian rationals, rather
# than as intervals, as often as possible.

def an_addsub_rational(a, b, sub):
    va = a._descr._value
    vb = b._descr._value
    if sub:
        v = va - vb
    else:
        v = va + vb
    return ANRational(v)

def an_muldiv_rational(a, b, div):
    va = a._descr._value
    vb = b._descr._value
    if div:
        v = va / vb
    else:
        v = va * vb
    return ANRational(v)

def an_addsub_zero(a, b, sub):
    if b._descr.is_rational() and b._descr.rational_value().is_zero():
        return a._descr
    # we know a is 0
    if sub:
        return b._descr.neg(b)
    else:
        return b._descr

def an_muldiv_zero(a, b, div):
    if b._descr.is_rational() and b._descr.rational_value().is_zero():
        if div:
            raise ValueError, "algebraic number division by zero"
        else:
            return ANRational(0)
    # we know a is 0
    return ANRational(0)

def an_addsub_gaussian(a, b, sub):
    va = a._descr.gaussian_value()
    vb = b._descr.gaussian_value()
    if sub:
        v = va - vb
    else:
        v = va + vb
    return ANExtensionElement(QQbar_I_generator, v)

def an_muldiv_gaussian(a, b, div):
    va = a._descr.gaussian_value()
    vb = b._descr.gaussian_value()
    if div:
        v = va / vb
    else:
        v = va * vb
    return ANExtensionElement(QQbar_I_generator, v)

def an_addsub_expr(a, b, sub):
    return ANBinaryExpr(a, b, ('-' if sub else '+'))

def an_muldiv_expr(a, b, div):
    return ANBinaryExpr(a, b, ('/' if div else '*'))

def an_muldiv_rootunity(a, b, div):
    ad = a._descr
    bd = b._descr
    if div:
        return ANRootOfUnity(ad.angle() - bd.angle(), ad.scale() / bd.scale())
    else:
        return ANRootOfUnity(ad.angle() + bd.angle(), ad.scale() * bd.scale())

def an_addsub_rootunity(a, b, sub):
    ad = a._descr
    bd = b._descr
    if ad._angle == bd._angle:
        if sub:
            return ANRootOfUnity(ad.angle(), ad.scale() - bd.scale())
        else:
            return ANRootOfUnity(ad.angle(), ad.scale() + bd.scale())
    else:
        return an_addsub_expr(a, b, sub)

def an_muldiv_element(a, b, div):
    ad = a._descr
    bd = b._descr
    adg = ad.generator()
    bdg = bd.generator()
    if adg == qq_generator or adg == bdg:
        if div:
            return ANExtensionElement(bdg, ad._value / bd._value)
        else:
            return ANExtensionElement(bdg, ad._value * bd._value)
    if bdg == qq_generator:
        if div:
            return ANExtensionElement(adg, ad._value / bd._value)
        else:
            return ANExtensionElement(adg, ad._value * bd._value)
    return ANBinaryExpr(a, b, ('/' if div else '*'))

def an_addsub_element(a, b, sub):
    ad = a._descr
    bd = b._descr
    adg = ad.generator()
    bdg = bd.generator()
    if adg == qq_generator or adg == bdg:
        if sub:
            return ANExtensionElement(bdg, ad._value - bd._value)
        else:
            return ANExtensionElement(bdg, ad._value + bd._value)
    if bdg == qq_generator:
        if sub:
            return ANExtensionElement(adg, ad._value - bd._value)
        else:
            return ANExtensionElement(adg, ad._value + bd._value)
    return ANBinaryExpr(a, b, ('-' if sub else '+'))

# Here we hand-craft a simple multimethod dispatch.
_mul_algo = {}
_add_algo = {}
_descriptors = ('zero', 'rational', 'imaginary', 'gaussian', 'rootunity', 'element', 'other')
for a in _descriptors:
    for b in _descriptors:
        key = (a, b)
        if a == 'zero' or b == 'zero':
            _mul_algo[key] = an_muldiv_zero
            _add_algo[key] = an_addsub_zero
            continue
        if a == 'rational' and b == 'rational':
            _mul_algo[key] = an_muldiv_rational
            _add_algo[key] = an_addsub_rational
            continue
        if b == 'rational':
            a1, b1 = b, a
        else:
            a1, b1 = a, b
        if a1 == 'rational':
            if b1 == 'imaginary':
                _mul_algo[key] = an_muldiv_rootunity
                _add_algo[key] = an_addsub_gaussian
                continue
            if b1 == 'gaussian':
                _mul_algo[key] = an_muldiv_gaussian
                _add_algo[key] = an_addsub_gaussian
                continue
            if b1 == 'rootunity':
                _mul_algo[key] = an_muldiv_rootunity
                _add_algo[key] = an_addsub_expr
                continue
            if b1 == 'element':
                _mul_algo[key] = an_muldiv_element
                _add_algo[key] = an_addsub_element
                continue
        if b1 == 'imaginary':
            a1, b1 = b1, a1
        if a1 == 'imaginary':
            if b1 == 'imaginary' or b1 == 'rootunity':
                _mul_algo[key] = an_muldiv_rootunity
                _add_algo[key] = an_addsub_rootunity
                continue
            if b1 == 'gaussian':
                _mul_algo[key] = an_muldiv_gaussian
                _add_algo[key] = an_addsub_gaussian
                continue
        if a1 == 'gaussian' and b1 == 'gaussian':
            _mul_algo[key] = an_muldiv_gaussian
            _add_algo[key] = an_addsub_gaussian
            continue
        if a1 == 'rootunity' and b1 == 'rootunity':
            _mul_algo[key] = an_muldiv_rootunity
            _add_algo[key] = an_addsub_rootunity
            continue
        if a1 == 'element' and b1 == 'element':
            _mul_algo[key] = an_muldiv_element
            _add_algo[key] = an_addsub_element
            continue
        _mul_algo[key] = an_muldiv_expr
        _add_algo[key] = an_addsub_expr

class ANDescr(SageObject):
    """
    An AlgebraicNumber or AlgebraicReal is a wrapper around an ANDescr object.
    ANDescr is an abstract base class, which should never
    be directly instantiated; its concrete subclasses are
    ANRational, ANBinaryExpr, ANUnaryExpr, ANRootOfUnity,
    ANRoot, and ANExtensionElement.
    ANDescr and all of its subclasses are private, and
    should not be used directly.
    """
    def is_exact(self):
        """
        Returns True if self is an ANRational, ANRootOfUnity, or
        ANExtensionElement.

        EXAMPLES:
            sage: from sage.rings.qqbar import ANRational
            sage: ANRational(1/2).is_exact()
            True
            sage: QQbar(3+I)._descr.is_exact()
            True
            sage: QQbar.zeta(17)._descr.is_exact()
            True
        """
        return False

    def is_rational(self):
        """
        Returns True if self is an ANRational object.  (Note that
        the constructors for ANExtensionElement and ANRootOfUnity
        will actually return ANRational objects for rational numbers.)

        EXAMPLES:
            sage: from sage.rings.qqbar import ANRational
            sage: ANRational(3/7).is_rational()
            True
        """
        return False

    def is_field_element(self):
        """
        Returns True if self is an ANExtensionElement.

            sage: from sage.rings.qqbar import ANExtensionElement, ANRoot, AlgebraicGenerator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: sqrt2 = ANExtensionElement(gen2, nf2.gen())
            sage: sqrt2.is_field_element()
            True
        """
        return False

    def neg(self, n):
        return ANUnaryExpr(n, '-')

    def invert(self, n):
        return ANUnaryExpr(n, '~')

    def abs(self, n):
        return ANUnaryExpr(n, 'abs')

    def real(self, n):
        if self.is_complex():
            return ANUnaryExpr(n, 'real')
        else:
            return self

    def imag(self, n):
        if self.is_complex():
            return ANUnaryExpr(n, 'imag')
        else:
            return ANRational(0)

    def conjugate(self, n):
        if self.is_complex():
            return ANUnaryExpr(n, 'conjugate')
        else:
            return self

    def norm(self, n):
        if self.is_complex():
            return ANUnaryExpr(n, 'norm')
        else:
            return (n*n)._descr

class AlgebraicNumber_base(sage.structure.element.FieldElement):
    """
    This is the common base class for algebraic numbers (complex
    numbers which are the zero of a polynomial in ZZ[x]) and algebraic
    reals (algebraic numbers which happen to be real).

    AlgebraicNumber objects can be created using QQbar (==
    AlgebraicNumberField()), and AlgebraicReal objects can be created
    using AA (== AlgebraicRealField()).  They can be created either by
    coercing a rational or a symbolic expression, or by using the
    QQbar.polynomial_root() or AA.polynomial_root() method to
    construct a particular root of a polynomial with algebraic
    coefficients.  Also, AlgebraicNumber and AlgebraicReal are closed
    under addition, subtraction, multiplication, division (except by
    0), and rational powers (including roots), except that for a
    negative AlgebraicReal, taking a power with an even denominator returns
    an AlgebraicNumber instead of an AlgebraicReal.

    AlgebraicNumber and AlgebraicReal objects can be approximated to
    any desired precision.  They can be compared exactly; if the two
    numbers are very close, or are equal, this may require exact
    computation, which can be extremely slow.

    As long as exact computation is not triggered, computation with
    algebraic numbers should not be too much slower than computation with
    intervals.  As mentioned above, exact computation is triggered
    when comparing two algebraic numbers which are very close together.
    This can be an explicit comparison in user code, but the following
    list of actions (not necessarily complete) can also trigger exact
    computation:
    Dividing by an algebraic number which is very close to 0.
    Using an algebraic number which is very close to 0 as the leading
    coefficient in a polynomial.
    Taking a root of an alebraic number which is very close to 0.

    The exact definition of "very close" is subject to change; currently,
    we compute our best approximation of the two numbers using 128-bit
    arithmetic, and see if that's sufficient to decide the comparison.
    Note that comparing two algebraic numbers which are actually equal will
    always trigger exact computation, unless they are actually the same object.

    EXAMPLES:
        sage: sqrt(QQbar(2))
        [1.4142135623730949 .. 1.4142135623730952]
        sage: sqrt(QQbar(2))^2 == 2
        True
        sage: x = polygen(QQbar)
        sage: phi = QQbar.polynomial_root(x^2 - x - 1, RIF(1, 2))
        sage: phi
        [1.6180339887498946 .. 1.6180339887498950]
        sage: phi^2 == phi+1
        True
        sage: AA(sqrt(65537))
        [256.00195311754941 .. 256.00195311754948]
    """

    def __init__(self, parent, x):
        """
        Initialize an algebraic number.  The argument must be either
        a rational number, a Gaussian rational, or a subclass of ANDescr.

        sage: from sage.rings.qqbar import ANRootOfUnity
        sage: AlgebraicReal(22/7)
        22/7
        sage: AlgebraicNumber(ANRootOfUnity(2/5, 1))
        [-0.80901699437494746 .. -0.80901699437494734] + [0.58778525229247302 .. 0.58778525229247314]*I
        """
        sage.structure.element.FieldElement.__init__(self, parent)
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            self._descr = ANRational(x)
        elif isinstance(x, (ANDescr)):
            self._descr = x
        elif parent is QQbar and \
                 isinstance(x, NumberFieldElement_quadratic) and \
                 list(x.parent().polynomial()) == [1, 0, 1]:
            self._descr = ANExtensionElement(QQbar_I_generator, QQbar_I_nf(x.list()))
        else:
            raise TypeError, "Illegal initializer for algebraic number"

        self._value = self._descr._interval_fast(64)

    def _repr_(self):
        """
        Returns the print representation of this number.

        EXAMPLES:
            sage: AA(22/7)
            22/7
            sage: QQbar(1/3 + 2/7*I)
            2/7*I + 1/3
            sage: QQbar.zeta(4) + 5
            I + 5
            sage: QQbar.zeta(17)
            [0.93247222940435570 .. 0.93247222940435582] + [0.36124166618715292 .. 0.36124166618715298]*I
            sage: AA(19).sqrt()
            [4.3588989435406730 .. 4.3588989435406740]
        """
        if self._descr.is_rational():
            return repr(self._descr)
        if isinstance(self._descr, ANRootOfUnity) and self._descr._angle == QQ_1_4:
            return '%s*I'%self._descr._scale
        if isinstance(self._descr, ANExtensionElement) and self._descr._generator is QQbar_I_generator:
            return repr(self._descr._value)
        if self.parent() is QQbar:
            return repr(CIF(self._value))
        else:
            return repr(RIF(self._value))

    def _mul_(self, other):
        """
        TESTS:
            sage: AA(sqrt(2)) * AA(sqrt(8))
            [3.9999999999999995 .. 4.0000000000000009]
        """
        sd = self._descr
        od = other._descr
        sdk = sd.kind()
        odk = od.kind()
        return type(self)(_mul_algo[sdk, odk](self, other, False))

    def _div_(self, other):
        """
        TESTS:
            sage: AA(sqrt(2)) / AA(sqrt(8))
            [0.49999999999999994 .. 0.50000000000000012]
        """
        sd = self._descr
        od = other._descr
        sdk = sd.kind()
        odk = od.kind()
        return type(self)(_mul_algo[sdk, odk](self, other, True))

    def __invert__(self):
        """
        TESTS:
            sage: ~AA(sqrt(~2))
            [1.4142135623730949 .. 1.4142135623730952]
        """
        sd = self._descr
        return type(self)(self._descr.invert(self))

    def _add_(self, other):
        """
        TESTS:
            sage: x = polygen(ZZ)
            sage: rt1, rt2 = (x^2 - x - 1).roots(ring=AA, multiplicities=False)
            sage: rt1 + rt2
            [0.99999999999999988 .. 1.0000000000000003]
        """
        sd = self._descr
        od = other._descr
        sdk = sd.kind()
        odk = od.kind()
        return type(self)(_add_algo[sdk, odk](self, other, False))

    def _sub_(self, other):
        """
        TESTS:
            sage: AA(golden_ratio) * 2 - AA(5).sqrt()
            [0.99999999999999988 .. 1.0000000000000003]
        """
        sd = self._descr
        od = other._descr
        sdk = sd.kind()
        odk = od.kind()
        return type(self)(_add_algo[sdk, odk](self, other, True))

    def _neg_(self):
        """
        TESTS:
            sage: -QQbar(I)
            -1*I
        """
        return type(self)(self._descr.neg(self))

    def __abs__(self):
        """
        TESTS:
            sage: abs(AA(sqrt(2) - sqrt(3)))
            [0.31783724519578221 .. 0.31783724519578227]
            sage: abs(QQbar(3+4*I))
            5
            sage: v = QQbar.zeta(3) + 1
            sage: v.exactify()
            sage: v.abs().minpoly()
            x - 1
        """
        return AlgebraicReal(self._descr.abs(self))

    def __hash__(self):
        """
        Compute a hash code for this number (equal algebraic numbers will
        have the same hash code, different algebraic numbers are likely
        to have different hash codes).

        This may trigger exact computation, but that is very unlikely.

        TESTS:
        The hash code is stable, even when the representation changes.
            sage: two = QQbar(4).nth_root(4)^2
            sage: two
            [1.9999999999999997 .. 2.0000000000000005]
            sage: h1 = hash(two)
            sage: two == 2
            True
            sage: two
            2
            sage: h2 = hash(two)
            sage: h1 == h2
            True

            sage: h1 = hash(QQbar.zeta(6))
            sage: h2 = hash(QQbar(1/2 + I*sqrt(3)/2))
            sage: h1 == h2
            True

        Unfortunately, the hash code for algebraic numbers which are close
        enough to each other are the same.  (This is inevitable, if
        equal algebraic reals give the same hash code and hashing does
        not always trigger exact computation.)
            sage: h1 = hash(QQbar(0))
            sage: h2 = hash(QQbar(1/2^100))
            sage: hash(h1) == hash(h2)
            True

        """

        # The only way I can think of to hash algebraic numbers without
        # always triggering exact computation is to use interval_exact().
        # However, interval_exact() always triggers exact computation
        # if the number is exactly representable in floating point, which
        # is presumably not too unlikely (algebraic reals like 0, 1/2,
        # 1, or 2 are presumably not uncommon).

        # So I modify the algebraic real by adding 1/123456789 to it before
        # calling interval_exact().  Then, exact computation will be triggered
        # by algebraic reals which are sufficiently close to
        # (some floating point number minus 1/123456789).  Hopefully,
        # -1/123456789 comes up in algebraic real computations far less
        # often than 0 does.  Algebraic numbers have a similar offset added,
        # with an additional complex component of 1/987654321*I.

        # All of this effort to avoid exact computation is probably wasted,
        # anyway... in almost all uses of hash codes, if the hash codes
        # match, the next step is to compare for equality; and comparing
        # for equality often requires exact computation.  (If a==b,
        # then checking a==b requires exact computation unless (a is b).)

        if self.parent() is AA:
            return hash((self + AA_hash_offset).interval_exact(RIF))
        else:
            return hash((self + QQbar_hash_offset).interval_exact(CIF))

    def sqrt(self):
        """
        Return the square root of this number.

        EXAMPLES:
            sage: AA(2).sqrt()
            [1.4142135623730949 .. 1.4142135623730952]
            sage: QQbar(I).sqrt()
            [0.70710678118654746 .. 0.70710678118654758] + [0.70710678118654746 .. 0.70710678118654758]*I
        """
        return self.__pow__(~ZZ(2))

    def nth_root(self, n):
        """
        Return the nth root of this number.

        Note that for odd n and negative real numbers, AlgebraicReal
        and AlgebraicNumber values give different answers: AlgebraicReal
        values prefer real results, and AlgebraicNumber values
        return the principal root.

        EXAMPLES:
            sage: AA(-8).nth_root(3)
            -2
            sage: QQbar(-8).nth_root(3)
            [0.99999999999999988 .. 1.0000000000000003] + [1.7320508075688771 .. 1.7320508075688775]*I
            sage: QQbar.zeta(12).nth_root(15)
            [0.99939082701909565 .. 0.99939082701909577] + [0.034899496702500969 .. 0.034899496702500977]*I
        """
        return self.__pow__(~ZZ(n))

    def exactify(self):
        """
        Compute an exact representation for this number.

        EXAMPLES:
            sage: two = QQbar(4).nth_root(4)^2
            sage: two
            [1.9999999999999997 .. 2.0000000000000005]
            sage: two.exactify()
            sage: two
            2
        """
        od = self._descr
        if od.is_exact(): return
        self._descr = self._descr.exactify()
        new_val = self._descr._interval_fast(self.parent().default_interval_prec())
        if is_RealIntervalFieldElement(new_val) and is_ComplexIntervalFieldElement(self._value):
            self._value = self._value.real().intersection(new_val)
        elif is_RealIntervalFieldElement(self._value) and is_ComplexIntervalFieldElement(new_val):
            self._value = self._value.intersection(new_val.real())
        else:
            self._value = self._value.intersection(new_val)

    def _exact_field(self):
        """
        Returns a generator for a number field that includes this number
        (not necessarily the smallest such number field).

        EXAMPLES:
            sage: QQbar(2)._exact_field()
            Trivial generator
            sage: (sqrt(QQbar(2)) + sqrt(QQbar(19)))._exact_field()
            Number Field in a with defining polynomial y^4 - 20*y^2 + 81 with a in [3.7893137826710354 .. 3.7893137826710360]
            sage: (QQbar(7)^(3/5))._exact_field()
            Number Field in a with defining polynomial y^5 - 7 with a in [1.4757731615945519 .. 1.4757731615945522]
        """

        sd = self._descr
        if sd.is_exact():
            return sd.generator()
        self.exactify()
        return self._exact_field()

    def _exact_value(self):
        """
        Returns an ANRational, an ANRootOfUnity, or an
        ANExtensionElement representing this value.

        EXAMPLES:
            sage: QQbar(2)._exact_value()
            2
            sage: (sqrt(QQbar(2)) + sqrt(QQbar(19)))._exact_value()
            1/9*a^3 + a^2 - 11/9*a - 10 where a^4 - 20*a^2 + 81 = 0 and a in [3.7893137826710354 .. 3.7893137826710360]
            sage: (QQbar(7)^(3/5))._exact_value()
            a^3 where a^5 - 7 = 0 and a in [1.4757731615945519 .. 1.4757731615945522]
        """
        sd = self._descr
        if sd.is_exact():
            return sd
        self.exactify()
        return self._descr

    def _more_precision(self):
        """
        Recompute the interval bounding this number with higher-precision
        interval arithmetic.

        EXAMPLES:
            sage: rt2 = sqrt(QQbar(2))
            sage: rt2._value
            [1.41421356237309504876 .. 1.41421356237309504888]
            sage: rt2._more_precision()
            sage: rt2._value
            [1.414213562373095048801688724209698078568 .. 1.414213562373095048801688724209698078575]
            sage: rt2._more_precision()
            sage: rt2._value
            [1.414213562373095048801688724209698078569671875376948073176679737990732478462101 .. 1.414213562373095048801688724209698078569671875376948073176679737990732478462120]
        """
        prec = self._value.prec()
        self._value = self._descr._interval_fast(prec*2)

    def minpoly(self):
        """
        Compute the minimal polynomial of this algebraic number.
        The minimal polynomial is the monic polynomial of least degree
        having this number as a root; it is unique.

        EXAMPLES:
            sage: QQbar(4).sqrt().minpoly()
            x - 2
            sage: ((QQbar(2).nth_root(4))^2).minpoly()
            x^2 - 2
            sage: v = sqrt(QQbar(2)) + sqrt(QQbar(3)); v
            [3.1462643699419721 .. 3.1462643699419726]
            sage: p = v.minpoly(); p
            x^4 - 10*x^2 + 1
            sage: p(RR(v.real()))
            1.31006316905768e-14
        """
        try:
            return self._minimal_polynomial
        except AttributeError:
            self.exactify()
            self._minimal_polynomial = self._descr.minpoly()
            return self._minimal_polynomial

    def degree(self):
        """
        Return the degree of this algebraic number (the degree of its
        minimal polynomial, or equivalently, the degree of the smallest
        algebraic extension of the rationals containing this number).

        EXAMPLES:
            sage: QQbar(5/3).degree()
            1
            sage: sqrt(QQbar(2)).degree()
            2
            sage: QQbar(17).nth_root(5).degree()
            5
            sage: sqrt(3+sqrt(QQbar(8))).degree()
            2
        """
        return self.minpoly().degree()

    def interval_fast(self, field):
        """
        Given a RealIntervalField, compute the value of this number
        using interval arithmetic of at least the precision of the field,
        and return the value in that field.  (More precision may be used
        in the computation.)  The returned interval may be arbitrarily
        imprecise, if this number is the result of a sufficiently long
        computation chain.

        EXAMPLES:
            sage: x = AA(2).sqrt()
            sage: x.interval_fast(RIF)
            [1.4142135623730949 .. 1.4142135623730952]
            sage: x.interval_fast(RealIntervalField(200))
            [1.4142135623730950488016887242096980785696718753769480731766796 .. 1.4142135623730950488016887242096980785696718753769480731766809]
            sage: x = QQbar(I).sqrt()
            sage: x.interval_fast(CIF)
            [0.70710678118654746 .. 0.70710678118654758] + [0.70710678118654746 .. 0.70710678118654758]*I
            sage: x.interval_fast(RIF)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert number to real interval.
        """
        if field.prec() == self._value.prec():
            return field(self._value)
        elif field.prec() > self._value.prec():
            self._more_precision()
            return self.interval_fast(field)
        else:
            return field(self._value)

    def interval_diameter(self, diam):
        """
        Compute an interval representation of self with diameter() at
        most diam.  The precision of the returned value is unpredictable.

        EXAMPLES:
            sage: AA(2).sqrt().interval_diameter(1e-10)
            [1.41421356237309504876 .. 1.41421356237309504888]
            sage: AA(2).sqrt().interval_diameter(1e-30)
            [1.414213562373095048801688724209698078568 .. 1.414213562373095048801688724209698078575]
            sage: QQbar(2).sqrt().interval_diameter(1e-10)
            [1.41421356237309504876 .. 1.41421356237309504888]
            sage: QQbar(2).sqrt().interval_diameter(1e-30)
            [1.414213562373095048801688724209698078568 .. 1.414213562373095048801688724209698078575]
        """
        if diam <= 0:
            raise ValueError, 'diameter must be positive in interval_diameter'

        while self._value.diameter() > diam:
            self._more_precision()

        return self._value

    def interval(self, field):
        """
        Given an interval field (real or complex, as appropriate) of
        precision p, compute an interval representation of self with
        diameter() at most 2^-p; then round that representation into
        the given field.  Here diameter() is relative diameter for
        intervals not containing 0, and absolute diameter for
        intervals that do contain 0; thus, if the returned interval
        does not contain 0, it has at least p-1 good bits.

        EXAMPLES:
            sage: RIF64 = RealIntervalField(64)
            sage: x = AA(2).sqrt()
            sage: y = x*x
            sage: y = 1000 * y - 999 * y
            sage: y.interval_fast(RIF64)
            [1.99999999999999966693 .. 2.00000000000000033307]
            sage: y.interval(RIF64)
            [1.99999999999999999989 .. 2.00000000000000000022]
            sage: CIF64 = ComplexIntervalField(64)
            sage: x = QQbar.zeta(11)
            sage: x.interval_fast(CIF64)
            [0.841253532831181168808 .. 0.841253532831181168918] + [0.540640817455597581992 .. 0.540640817455597582210]*I
            sage: x.interval(CIF64)
            [0.841253532831181168808 .. 0.841253532831181168864] + [0.540640817455597582101 .. 0.540640817455597582156]*I
        """
        target = RR(1.0) >> field.prec()
        val = self.interval_diameter(target)
        return field(val)

class AlgebraicNumber(AlgebraicNumber_base):
    """
    The class for algebraic numbers (complex numbers which are the roots
    of a polynomial with integer coefficients).  Much of its functionality
    is inherited from AlgebraicNumber_base.
    """
    def __init__(self, x):
        AlgebraicNumber_base.__init__(self, QQbar, x)

    def __cmp__(self, other):
        """
        Compare two algebraic numbers, lexicographically.  (That is,
        first compare the real components; if the real componetns are
        equal, compare the imaginary components.)

        EXAMPLES:
            sage: x = QQbar.zeta(3); x
            [-0.50000000000000012 .. -0.49999999999999994] + [0.86602540378443859 .. 0.86602540378443871]*I
            sage: cmp(QQbar(-1), x)
            -1
            sage: cmp(QQbar(-1/2), x)
            -1
            sage: cmp(QQbar(0), x)
            1
        """
        if self is other: return 0
        rcmp = cmp(self.real(), other.real())
        if rcmp != 0:
            return rcmp
        return cmp(self.imag(), other.imag())

    def __eq__(self, other):
        """
        Test two algebraic numbers for equality.

        EXAMPLES:
            sage: QQbar.zeta(6) == QQbar(1/2 + I*sqrt(3)/2)
            True
            sage: QQbar(I) == QQbar(I * (2^100+1)/(2^100))
            False
        """
        if not isinstance(other, AlgebraicNumber):
            other = self.parent()(other)
        if self is other: return True
        if other._descr.is_rational() and other._descr.rational_value() == 0:
            return not self.__nonzero__()
        if self._descr.is_rational() and self._descr.rational_value() == 0:
            return not other.__nonzero__()
        return not self._sub_(other).__nonzero__()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __nonzero__(self):
        """
        Check whether self is equal is nonzero.  This is fast if
        interval arithmetic proves that self is nonzero, but may be
        slow if the number actually is very close to zero.

        EXAMPLES:
            sage: (QQbar.zeta(2) + 1).__nonzero__()
            False
            sage: (QQbar.zeta(7) / (2^500)).__nonzero__()
            True
        """
        val = self._value
        if not val.contains_zero():
            return True
        if self._descr.is_field_element():
            # The ANExtensionElement constructor returns an ANRational
            # instead, if the number is zero.
            return True
        if self._descr.is_rational():
            return self._descr.rational_value().__nonzero__()
        if self._value.prec() < 128:
            self._more_precision()
            return self.__nonzero__()

        # Sigh...
        self.exactify()
        return self.__nonzero__()

    def __pow__(self, e):
        """
        self^p returns the p'th power of self (where p can be an arbitrary
        rational).  If p is (a/b), takes the principal b'th root of self,
        then takes that to the a'th power.  (Note that this differs
        from __pow__ on algebraic reals, where real roots are preferred
        over principal roots if they exist.)

        EXAMPLES:
            sage: QQbar(2)^(1/2)
            [1.4142135623730949 .. 1.4142135623730952]
            sage: QQbar(8)^(2/3)
            4
            sage: QQbar(8)^(2/3) == 4
            True
            sage: x = polygen(QQbar)
            sage: phi = QQbar.polynomial_root(x^2 - x - 1, RIF(1, 2))
            sage: tau = QQbar.polynomial_root(x^2 - x - 1, RIF(-1, 0))
            sage: rt5 = QQbar(5)^(1/2)
            sage: phi^10 / rt5
            [55.003636123247410 .. 55.003636123247418]
            sage: tau^10 / rt5
            [0.0036361232474132654 .. 0.0036361232474132659]
            sage: (phi^10 - tau^10) / rt5
            [54.999999999999992 .. 55.000000000000008]
            sage: (phi^10 - tau^10) / rt5 == fibonacci(10)
            True
            sage: (phi^50 - tau^50) / rt5 == fibonacci(50)
            True
            sage: QQbar(-8)^(1/3)
            [0.99999999999999988 .. 1.0000000000000003] + [1.7320508075688771 .. 1.7320508075688775]*I
            sage: (QQbar(-8)^(1/3))^3
            -8
            sage: QQbar(32)^(1/5)
            2
            sage: a = QQbar.zeta(7)^(1/3); a
            [0.95557280578614067 .. 0.95557280578614079] + [0.29475517441090420 .. 0.29475517441090427]*I
            sage: a == QQbar.zeta(21)
            True
            sage: QQbar.zeta(7)^6
            [0.62348980185873348 .. 0.62348980185873360] - [0.78183148246802980 .. 0.78183148246802992]*I
            sage: (QQbar.zeta(7)^6)^(1/3) * QQbar.zeta(21)
            1
        """
        e = QQ._coerce_(e)
        n = e.numerator()
        d = e.denominator()
        if d == 1:
            return generic_power(self, n)

        # First, check for exact roots.
        if isinstance(self._descr, ANRational):
            rt = rational_exact_root(abs(self._descr._value), d)
            if rt is not None:
                if self._descr._value < 0:
                    return AlgebraicNumber(ANRootOfUnity(~(2*d), rt))**n
                else:
                    return AlgebraicNumber(ANRational(rt))**n
        elif isinstance(self._descr, ANRootOfUnity):
            rt = rational_exact_root(abs(self._descr._scale), d)
            if rt is not None:
                if self._descr._scale < 0:
                    return AlgebraicNumber(ANRootOfUnity((self._descr._angle - QQ_1_2)/d, rt))**n
                else:
                    return AlgebraicNumber(ANRootOfUnity(self._descr._angle/d, rt))**n

        # Without this special case, we don't know the multiplicity
        # of the desired root
        if self.is_zero():
            return AlgebriacNumber(0)
        argument_is_pi = False
        for prec in short_prec_seq():
            if prec is None:
                # We know that self.real() < 0, since self._value
                # crosses the negative real line and self._value
                # is known to be non-zero.
                isgn = self.imag().sign()
                val = self._value
                argument = val.argument()
                if isgn == 0:
                    argument = argument.parent().pi()
                    argument_is_pi = True
                elif isgn > 0:
                    if argument < 0:
                        argument = argument + 2 * argument.parent().pi()
                else:
                    if argument > 0:
                        argument = argument - 2 * argument.parent().pi()
            else:
                val = self._interval_fast(prec)
                if not val.crosses_log_branch_cut():
                    argument = val.argument()
                    if val.imag().is_zero() and val.real() < 0:
                        argument_is_pi = True
                    break

        target_abs = abs(val) ** e
        target_arg = argument * e

        for prec in tail_prec_seq():
            if target_abs.relative_diameter() < RR_1_10 and (target_arg * d).absolute_diameter() < RR_1_10:
                break

            val = self._interval_fast(prec)

            target_abs = abs(val) ** e
            argument = val.argument()
            if argument_is_pi:
                argument = argument.parent().pi()
            target_arg = argument * e

        pow_n = self**n
        poly = QQbarPoly.gen()**d - pow_n

        prec = target_abs.prec()
        if argument_is_pi and d == 2:
            target_real = 0
        else:
            target_real = target_arg.cos() * target_abs
        target = ComplexIntervalField(prec)(target_real,
                                            target_arg.sin() * target_abs)

        return AlgebraicNumber(ANRoot(poly, target))

    def _interval_fast(self, prec):
        return self.interval_fast(ComplexIntervalField(prec))

    def real(self):
        return AlgebraicReal(self._descr.real(self))

    def imag(self):
        return AlgebraicReal(self._descr.imag(self))

    def conjugate(self):
        """
        Returns the complex conjugate of self.

        EXAMPLES:
            sage: QQbar(3 + 4*I).conjugate()
            [3.0000000000000000 .. 3.0000000000000000] - [4.0000000000000000 .. 4.0000000000000000]*I
            sage: QQbar.zeta(7).conjugate()
            [0.62348980185873348 .. 0.62348980185873360] - [0.78183148246802980 .. 0.78183148246802992]*I
            sage: QQbar.zeta(7) + QQbar.zeta(7).conjugate()
            [1.2469796037174669 .. 1.2469796037174672] + [-3.2526065174565134e-19 .. 3.7947076036992656e-19]*I
        """
        return AlgebraicNumber(self._descr.conjugate(self))

    def norm(self):
        """
        Returns self * self.conjugate().  This is the algebraic
        definition of norm, if we view QQbar as AA[I].

        EXAMPLES:
            sage: QQbar(3 + 4*I).norm()
            25
            sage: type(QQbar(I).norm())
            <class 'sage.rings.qqbar.AlgebraicReal'>
            sage: QQbar.zeta(1007).norm()
            1
        """
        return AlgebraicReal(self._descr.norm(self))

    def interval_exact(self, field):
        """
        Given a ComplexIntervalField, compute the best possible
        approximation of this number in that field.  Note that if
        either the real or imaginary parts of this number are
        sufficiently close to some floating-point number (and, in
        particular, if either is exactly representable in floating-point),
        then this will trigger exact computation, which may be very slow.

        EXAMPLES:
            sage: a = QQbar(I).sqrt(); a
            [0.70710678118654746 .. 0.70710678118654758] + [0.70710678118654746 .. 0.70710678118654758]*I
            sage: a.interval_exact(CIF)
            [0.70710678118654746 .. 0.70710678118654758] + [0.70710678118654746 .. 0.70710678118654758]*I
            sage: b = QQbar((1+I)*sqrt(2)/2)
            sage: (a - b).interval(CIF)
            [-5.4210108624275222e-20 .. 5.4210108624275222e-20] + [-1.0842021724855045e-19 .. 5.4210108624275222e-20]*I
            sage: (a - b).interval_exact(CIF)
            0
        """
        if not is_ComplexIntervalField(field):
            raise ValueError, "AlgebraicNumber interval_exact requires a ComplexIntervalField"
        rfld = field._real_field()
        re = self.real().interval_exact(rfld)
        im = self.imag().interval_exact(rfld)
        return field(re, im)

    def _complex_mpfr_field_(self, field):
        if is_ComplexIntervalField(field):
            return self.interval(field)
        else:
            return self.complex_number(field)


    def complex_number(self, field):
        """
        Given a ComplexField, compute a good approximation to self in that
        field.  The approximation will be off by at most two ulp's in
        each component, except for components which are very close to
        zero, which will have an abolute error at most 2^(-(field.prec()-1)).

        EXAMPLES:
            sage: a = QQbar.zeta(5)
            sage: a.complex_number(CIF)
            0.309016994374947 + 0.951056516295154*I
            sage: (a + a.conjugate()).complex_number(CIF)
            0.618033988749895 - 5.42101086242752e-20*I
        """
        v = self.interval(ComplexIntervalField(field.prec()))
        return v.center()

    def complex_exact(self, field):
        """
        Given a ComplexField, return the best possible approximation of
        this number in that field.  Note that if either component is
        sufficiently close to the halfway point between two floating-point
        numbers in the corresponding RealField, then this will trigger
        exact computation, which may be very slow.

        EXAMPLES:
            sage: a = QQbar.zeta(9) + I + QQbar.zeta(9).conjugate(); a
            [1.5320888862379560 .. 1.5320888862379563] + [0.99999999999999988 .. 1.0000000000000003]*I
            sage: a.complex_exact(CIF)
            [1.5320888862379560 .. 1.5320888862379563] + [1.0000000000000000 .. 1.0000000000000000]*I
        """
        rfld = field._real_field()
        re = self.real().real_exact(rfld)
        im = self.imag().real_exact(rfld)
        return field(re, im)

    def multiplicative_order(self):
        """
        Compute the multiplicative order of this algebraic real
        number.  That is, find the smallest positive integer n such
        that x^n == 1.  If there is no such n, returns +Infinity.

        We first check that abs(x) is very close to 1.  If so, we compute
        x exactly and examine its argument.

        EXAMPLES:
            sage: QQbar(-sqrt(3)/2 - I/2).multiplicative_order()
            12
            sage: QQbar(1).multiplicative_order()
            1
            sage: QQbar(-I).multiplicative_order()
            4
            sage: QQbar(707/1000 + 707/1000*I).multiplicative_order()
            +Infinity
            sage: QQbar(3/5 + 4/5*I).multiplicative_order()
            +Infinity
        """
        if not (1 in CIF(self).norm()):
            return infinity.infinity
        if self.norm() != 1:
            return infinity.infinity
        ra = self.rational_argument()
        if ra is None:
            return infinity.infinity
        return ra.denominator()

    def rational_argument(self):
        """
        Returns the argument of self, divided by 2*pi, as long as this
        result is rational.  Otherwise returns None.  Always triggers
        exact computation.

        EXAMPLES:
            sage: QQbar((1+I)*(sqrt(2)+sqrt(5))).rational_argument()
            1/8
            sage: QQbar(-1 + I*sqrt(3)).rational_argument()
            1/3
            sage: QQbar(-1 - I*sqrt(3)).rational_argument()
            -1/3
            sage: QQbar(3+4*I).rational_argument() is None
            True
            sage: (QQbar.zeta(7654321)^65536).rational_argument()
            65536/7654321
            sage: (QQbar.zeta(3)^65536).rational_argument()
            1/3
        """
        # This always triggers exact computation.  An alternate method
        # could almost always avoid exact computation when the result
        # is None: if we can compute an upper bound on the degree of
        # this algebraic number without exact computation, we can use
        # the method of ANExtensionElement.rational_argument().

        # Even a very loose upper bound would suffice; for instance,
        # an upper bound of 2^100, when the true degree was 8, would
        # still be efficient.

        self.exactify()
        return self._descr.rational_argument(self)

class AlgebraicReal(AlgebraicNumber_base):
    def __init__(self, x):
        AlgebraicNumber_base.__init__(self, AA, x)

    def __cmp__(self, other):
        """
        Compare two algebraic reals.

        EXAMPLES:
            sage: cmp(AA(golden_ratio), AA(sqrt(5)))
            -1
            sage: cmp(AA(golden_ratio), AA((sqrt(5)+1)/2))
            0
            sage: cmp(AA(7), AA(50/7))
            -1
        """
        if self is other: return 0
        if other._descr.is_rational() and other._descr.rational_value() == 0:
            return self.sign()
        elif self._descr.is_rational() and self._descr.rational_value() == 0:
            return -other.sign()
        else:
            return self._sub_(other).sign()

    def __pow__(self, e):
        """
        self^p returns the p'th power of self (where p can be an
        arbitrary rational).  If p is (a/b), takes the b'th root of
        self, then takes that to the a'th power.  If self is negative
        and b is odd, it takes the real b'th root; if self is odd and
        b is even, this takes a complex root.  Note that the behavior
        when self is negative and b is odd differs from the complex
        case; algebraic numbers select the principal complex b'th
        root, but algebraic reals select the real root.

        EXAMPLES:
            sage: AA(2)^(1/2)
            [1.4142135623730949 .. 1.4142135623730952]
            sage: AA(8)^(2/3)
            4
            sage: AA(8)^(2/3) == 4
            True
            sage: x = polygen(AA)
            sage: phi = AA.polynomial_root(x^2 - x - 1, RIF(0, 2))
            sage: tau = AA.polynomial_root(x^2 - x - 1, RIF(-2, 0))
            sage: rt5 = AA(5)^(1/2)
            sage: phi^10 / rt5
            [55.003636123247410 .. 55.003636123247418]
            sage: tau^10 / rt5
            [0.0036361232474132654 .. 0.0036361232474132659]
            sage: (phi^10 - tau^10) / rt5
            [54.999999999999992 .. 55.000000000000008]
            sage: (phi^10 - tau^10) / rt5 == fibonacci(10)
            True
            sage: (phi^50 - tau^50) / rt5 == fibonacci(50)
            True

        TESTS:
            sage: AA(-8)^(1/3)
            -2
            sage: AA(-8)^(2/3)
            4
            sage: AA(32)^(3/5)
            8
            sage: AA(-16)^(1/2)
            4*I
            sage: AA(-16)^(1/4)
            [1.4142135623730949 .. 1.4142135623730952] + [1.4142135623730949 .. 1.4142135623730952]*I
            sage: AA(-16)^(1/4)/QQbar.zeta(8)
            2
        """
        e = QQ._coerce_(e)
        n = e.numerator()
        d = e.denominator()
        if d == 1:
            return generic_power(self, n)

        # First, check for exact roots.
        if isinstance(self._descr, ANRational):
            rt = rational_exact_root(abs(self._descr._value), d)
            if rt is not None:
                if self._descr._value < 0:
                    if d % 2 == 0:
                        return AlgebraicNumber(ANRootOfUnity(~(2*d), rt))**n
                    else:
                        return AlgebraicReal(ANRational(-rt))**n
                else:
                    return AlgebraicReal(ANRational(rt))**n

        # Without this special case, we don't know the multiplicity
        # of the desired root
        if self.sign() == 0:
            return AlgebriacNumber(0)
        if d % 2 == 0:
            if self.sign() < 0:
                return QQbar(self).__pow__(e)
        pow_n = self**n
        poly = AAPoly.gen()**d - pow_n
        range = pow_n.interval_fast(RIF)
        if d % 2 == 0:
            result_min = 0
        else:
            result_min = min(range.lower(), -1)
        result_max = max(range.upper(), 1)
        return AlgebraicReal(ANRoot(poly, RIF(result_min, result_max)))

    def sign(self):
        """
        Compute the sign of this algebraic number (return -1 if negative,
        0 if zero, or 1 if positive).

        Computes an interval enclosing this number using 128-bit interval
        arithmetic; if this interval includes 0, then fall back to
        exact computation (which can be very slow).

        EXAMPLES:
            sage: AA(-5).nth_root(7).sign()
            -1
            sage: (AA(2).sqrt() - AA(2).sqrt()).sign()
            0
        """
        if self._value.lower() > 0:
            return 1
        elif self._value.upper() < 0:
            return -1
        elif self._descr.is_rational():
            val = self._descr.rational_value()
            if val > 0:
                return 1
            elif val < 0:
                return -1
            else:
                return 0
        elif self._descr.is_field_element():
            # All field elements are irrational by construction
            # (the ANExtensionElement constructor will return an ANRational
            # instead, if the number is actually rational).
            # An irrational number must eventually be different from 0
            self._more_precision()
            return self.sign()
        elif self._value.prec() < 128:
            # OK, we'll try adding precision one more time
            # print self._value
            self._more_precision()
            return self.sign()
        else:
            # Sigh...
            self.exactify()
            return self.sign()

    def _interval_fast(self, prec):
        return self.interval_fast(RealIntervalField(prec))

    def interval_exact(self, field):
        """
        Given a RealIntervalField, compute the best possible
        approximation of this number in that field.  Note that if this
        number is sufficiently close to some floating-point number
        (and, in particular, if this number is exactly representable in
        floating-point), then this will trigger exact computation, which
        may be very slow.

        EXAMPLES:
            sage: x = AA(2).sqrt()
            sage: y = x*x
            sage: x.interval(RIF)
            [1.4142135623730949 .. 1.4142135623730952]
            sage: x.interval_exact(RIF)
            [1.4142135623730949 .. 1.4142135623730952]
            sage: y.interval(RIF)
            [1.9999999999999997 .. 2.0000000000000005]
            sage: y.interval_exact(RIF)
            [2.0000000000000000 .. 2.0000000000000000]
            sage: z = 1 + AA(2).sqrt() / 2^200
            sage: z.interval(RIF)
            [1.0000000000000000 .. 1.0000000000000003]
            sage: z.interval_exact(RIF)
            [1.0000000000000000 .. 1.0000000000000003]
        """
        for extra in (0, 40):
            target = RR(1.0) >> field.prec()
            # p==precise; pr==precise rounded
            pval = self.interval_diameter(target)
            pbot = pval.lower()
            ptop = pval.upper()
            val = field(pval)
            bot = val.lower()
            top = val.upper()
            prbot = pbot.parent()(bot)
            prtop = ptop.parent()(top)
            if bot == top or (bot.nextabove() == top and
                              prbot < pbot and ptop < prtop):
                return val

        # Even 40 extra bits of precision aren't enough to prove that
        # self is not an exactly representable float.
        self.exactify()
        while True:
            # p==precise; pr==precise rounded
            pval = self._value
            pbot = pval.lower()
            ptop = pval.upper()
            val = field(pval)
            bot = val.lower()
            top = val.upper()
            prbot = pbot.parent()(bot)
            prtop = ptop.parent()(top)
            if bot == top or (bot.nextabove() == top and
                              prbot < pbot and ptop < prtop):
                return val

            self._more_precision()

    def real_number(self, field):
        """
        Given a RealField, compute a good approximation to self in that field.
        The approximation will be off by at most two ulp's, except for
        numbers which are very close to 0, which will have an absolute
        error at most 2^(-(field.prec()-1)).  Also, the rounding mode of the
        field is respected.

        EXAMPLES:
            sage: x = AA(2).sqrt()^2
            sage: x.real_number(RR)
            2.00000000000000
            sage: x.real_number(RealField(53, rnd='RNDD'))
            1.99999999999999
            sage: x.real_number(RealField(53, rnd='RNDU'))
            2.00000000000001
            sage: x.real_number(RealField(53, rnd='RNDZ'))
            1.99999999999999
            sage: (-x).real_number(RR)
            -2.00000000000000
            sage: (-x).real_number(RealField(53, rnd='RNDD'))
            -2.00000000000001
            sage: (-x).real_number(RealField(53, rnd='RNDU'))
            -1.99999999999999
            sage: (-x).real_number(RealField(53, rnd='RNDZ'))
            -1.99999999999999
            sage: (x-2).real_number(RR)
            5.42101086242752e-20
            sage: (x-2).real_number(RealField(53, rnd='RNDD'))
            -1.08420217248551e-19
            sage: (x-2).real_number(RealField(53, rnd='RNDU'))
            2.16840434497101e-19
            sage: (x-2).real_number(RealField(53, rnd='RNDZ'))
            0.000000000000000
            sage: y = AA(2).sqrt()
            sage: y.real_number(RR)
            1.41421356237309
            sage: y.real_number(RealField(53, rnd='RNDD'))
            1.41421356237309
            sage: y.real_number(RealField(53, rnd='RNDU'))
            1.41421356237310
            sage: y.real_number(RealField(53, rnd='RNDZ'))
            1.41421356237309
        """
        v = self.interval(RealIntervalField(field.prec()))

        mode = field.rounding_mode()
        if mode == 'RNDN':
            return v.center()
        if mode == 'RNDD':
            return v.lower()
        if mode == 'RNDU':
            return v.upper()
        if mode == 'RNDZ':
            if v > 0:
                return field(v.lower())
            elif v < 0:
                return field(v.upper())
            else:
                return field(0)

    _mpfr_ = real_number

    def _complex_mpfr_field_(self, field):
        if is_ComplexIntervalField(field):
            return field(self.interval(field._real_field()))
        else:
            return field(self.real_number(field._real_field()))

    def real_exact(self, field):
        """
        Given a RealField, compute the best possible approximation of
        this number in that field.  Note that if this number is
        sufficiently close to the halfway point between two
        floating-point numbers in the field (for the default
        round-to-nearest mode) or if the number is sufficiently close
        to a floating-point number in the field (for directed rounding
        modes), then this will trigger exact computation, which may be
        very slow.

        The rounding mode of the field is respected.

        EXAMPLES:
            sage: x = AA(2).sqrt()^2
            sage: x.real_exact(RR)
            2.00000000000000
            sage: x.real_exact(RealField(53, rnd='RNDD'))
            2.00000000000000
            sage: x.real_exact(RealField(53, rnd='RNDU'))
            2.00000000000000
            sage: x.real_exact(RealField(53, rnd='RNDZ'))
            2.00000000000000
            sage: (-x).real_exact(RR)
            -2.00000000000000
            sage: (-x).real_exact(RealField(53, rnd='RNDD'))
            -2.00000000000000
            sage: (-x).real_exact(RealField(53, rnd='RNDU'))
            -2.00000000000000
            sage: (-x).real_exact(RealField(53, rnd='RNDZ'))
            -2.00000000000000
            sage: (x-2).real_exact(RR)
            0.000000000000000
            sage: (x-2).real_exact(RealField(53, rnd='RNDD'))
            0.000000000000000
            sage: (x-2).real_exact(RealField(53, rnd='RNDU'))
            0.000000000000000
            sage: (x-2).real_exact(RealField(53, rnd='RNDZ'))
            0.000000000000000
            sage: y = AA(2).sqrt()
            sage: y.real_exact(RR)
            1.41421356237310
            sage: y.real_exact(RealField(53, rnd='RNDD'))
            1.41421356237309
            sage: y.real_exact(RealField(53, rnd='RNDU'))
            1.41421356237310
            sage: y.real_exact(RealField(53, rnd='RNDZ'))
            1.41421356237309
        """
        for extra in (0, 40):
            target = RR(1.0) >> field.prec()
            val = self.interval_diameter(target)
            fbot = field(val.lower())
            ftop = field(val.upper())
            if fbot == ftop:
                return ftop

        # Even 40 extra bits of precision aren't enough to determine the
        # answer.
        rifp1 = RealIntervalField(field.prec() + 1)
        rifp2 = RealIntervalField(field.prec() + 2)

        val = self.interval_exact(rifp1)

        # Call the largest floating-point number <= self 'x'.  Then
        # val may be [x .. x], [x .. x + 1/2 ulp],
        # [x + 1/2 ulp .. x + 1/2 ulp], or [x + 1/2 ulp .. x + 1 ulp];
        # in the second and fourth cases, the true value is not equal
        # to either of the interval endpoints.

        mid = rifp2(val).center()

        # Now mid may be x, x + 1/4 ulp, x + 1/2 ulp, or x + 3/4 ulp; in
        # the first and third cases, mid is the exact, true value of self;
        # in the second and fourth cases, self is close to mid, and is
        # neither x, x + 1/2 ulp, nor x + 1 ulp.

        # In all of these cases, in all rounding modes, the rounded value
        # of mid is the same as the rounded value of self.

        return field(mid)

class ANRational(ANDescr):
    """
    The subclass of ANDescr that represents an arbitrary
    rational.  This class is private, and should not be used directly.
    """

    def __init__(self, x):
        """
        TESTS:
            sage: polygen(QQbar) / int(3)
            1/3*x
            sage: QQbar(int(7)) / QQbar(long(2))
            7/2
        """
        if isinstance(x, (sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            self._value = x
        elif isinstance(x, (int, long)):
            self._value = ZZ(x)
        else:
            raise TypeError, "Illegal initializer for algebraic number rational"

    def _repr_(self):
        return repr(self._value)

    def kind(self):
        if self._value.is_zero():
            return 'zero'
        else:
            return 'rational'

    def _interval_fast(self, prec):
        return RealIntervalField(prec)(self._value)

    def generator(self):
        return qq_generator

    def is_complex(self):
        return False

    def is_rational(self):
        return True

    def rational_value(self):
        return self._value

    def exactify(self):
        return self

    def is_exact(self):
        return True

    def minpoly(self):
        return QQx_x - self._value

    def neg(self, n):
        return ANRational(-self._value)

    def invert(self, n):
        return ANRational(~self._value)

    def abs(self, n):
        return ANRational(abs(self._value))

    def rational_argument(self, n):
        if self._value > 0:
            return QQ(0)
        if self._value < 0:
            return QQ(1)/2
        return None

    def gaussian_value(self):
        return QQbar_I_nf(self._value)

    def angle(self):
        return QQ_0

    def scale(self):
        return self._value

class ANRootOfUnity(ANDescr):
    """
    The subclass of ANDescr that represents a rational multiplied
    by a root of unity.  This class is private, and should not be
    used directly.

    Such numbers are represented by a "rational angle" and a rational
    scale.  The "rational angle" is the argument of the number, divided by
    2*pi; so given angle and scale, the number is:
    scale*(cos(2*pi*angle) + sin(2*pi*angle)*I); or equivalently
    scale*(e^(2*pi*angle*I)).

    We normalize so that 0<angle<1/2; this requires allowing both positive
    and negative scales.  (Attempts to create an ANRootOfUnity with an
    angle which is a multiple of 1/2 end up creating an ANRational instead.)
    """

    def __new__(self, angle, scale):
        """
        Construct an ANRootOfUnity from a rational angle and a rational
        scale.  If the number is actually a real rational, returns an
        ANRational instead.
        """
        if scale.is_zero():
            return ANRational(0)
        try:
            int_angle = ZZ(angle*2)
        except TypeError:
            return ANDescr.__new__(self)

        if int_angle & 1:
            # int_angle is odd
            return ANRational(-scale)
        else:
            # int_angle is even
            return ANRational(scale)

    def __init__(self, angle, scale):
        """
        Construct an ANRootOfUnity from a rational angle and a rational
        scale.
        """

        angle2 = angle * 2
        fl2 = angle2.floor()
        angle2 = angle2 - fl2
        angle = angle2 / 2
        if fl2 & 1:
            scale = -scale
        self._angle = angle
        self._scale = scale

    def _repr_(self):
        return "%s*e^(2*pi*I*%s)"%(self._scale, self._angle)

    def kind(self):
        if self._angle == QQ_1_4:
            return 'imaginary'
        else:
            return 'rootunity'

    def _interval_fast(self, prec):
        argument = self._angle * RealIntervalField(prec).pi() * 2
        if self._angle == QQ_1_4:
            return ComplexIntervalField(prec)(0, self._scale)
        else:
            return ComplexIntervalField(prec)(argument.cos(), argument.sin()) * self._scale

    def generator(self):
        return cyclotomic_generator(self._angle.denominator())

    def field_element_value(self):
        gen = self.generator()
        f = gen._field
        a = f.gen()
        return self._scale * a ** self._angle.numerator()

    def is_complex(self):
        return True

    def exactify(self):
        return self

    def is_exact(self):
        return True

    def minpoly(self):
        """
        EXAMPLES:
            sage: a = QQbar.zeta(7) * 2; a
            [1.2469796037174669 .. 1.2469796037174672] + [1.5636629649360596 .. 1.5636629649360599]*I
            sage: a.minpoly()
            x^6 + 2*x^5 + 4*x^4 + 8*x^3 + 16*x^2 + 32*x + 64
            sage: a.minpoly()(a)
            [-3.1918911957973251e-16 .. 3.4694469519536142e-16] + [-3.3133218391157016e-16 .. 3.2786273695961655e-16]*I
            sage: a.minpoly()(a) == 0
            True
        """
        # This could be more efficient...
        p = cyclotomic_polynomial(self._angle.denominator())
        p = p(p.parent().gen() / self._scale)
        p = p / p.leading_coefficient()
        return p

    def neg(self, n):
        return ANRootOfUnity(self._angle, -self._scale)

    def invert(self, n):
        # We want ANRootOfUnity(-self._angle, ~self._scale);
        # but that's not normalized, so we pre-normalize it to:
        return ANRootOfUnity(QQ_1_2 - self._angle, -~self._scale)

    def conjugate(self, n):
        # We want ANRootOfUnity(-self._angle, self._scale);
        # but that's not normalized, so we pre-normalize it to:
        return ANRootOfUnity(QQ_1_2 - self._angle, -self._scale)

    def abs(self, n):
        return ANRational(abs(self._scale))

    def norm(self, n):
        return ANRational(self._scale * self._scale)

    def rational_argument(self, n):
        if self._scale > 0:
            return self._angle
        else:
            return self._angle - QQ_1_2

    def gaussian_value(self):
        assert(self._angle == QQ_1_4)
        return QQbar_I_nf(self._scale * QQbar_I_nf.gen())

    def angle(self):
        return self._angle

    def scale(self):
        return self._scale

def is_AlgebraicReal(x):
    return isinstance(x, AlgebraicReal)

def is_AlgebraicNumber(x):
    return isinstance(x, AlgebraicNumber)

QQbarPoly = PolynomialRing(QQbar, 'x')
AAPoly = PolynomialRing(AA, 'x')

class AlgebraicPolynomialTracker(SageObject):
    """
    Keeps track of a polynomial used for algebraic numbers.

    If multiple algebraic numbers are created as roots of a single
    polynomial, this allows the polynomial and information about
    the polynomial to be shared.  This reduces work if the polynomial
    must be recomputed at higher precision, or if it must be factored.

    This class is private, and should only be constructed by
    AA.common_polynomial() or QQbar.common_polynomial(), and should
    only be used as an argument to AA.polynomial_root() or
    QQbar.polynomial_root().  (It doesn't matter whether you create
    the common polynomial with AA.common_polynomial() or
    QQbar.common_polynomial().)

    EXAMPLES:
        sage: x = polygen(QQbar)
        sage: P = QQbar.common_polynomial(x^2 - x - 1)
        sage: P
        x^2 + (-1)*x - 1
        sage: QQbar.polynomial_root(P, RIF(1, 2))
        [1.6180339887498946 .. 1.6180339887498950]
    """

    def __init__(self, poly):
        if not is_Polynomial(poly):
            raise ValueError, "Trying to create AlgebraicPolynomialTracker on non-Polynomial"
        if isinstance(poly.base_ring(), AlgebraicField_common):
            complex = is_AlgebraicField(poly.base_ring())
        else:
            try:
                poly = poly.change_ring(AA)
                complex = False
            except TypeError:
                poly = poly.change_ring(QQbar)
                complex = True
        self._poly = poly
        self._complex = complex
        self._exact = False
        self._roots_cache = {}

    def _repr_(self):
        return repr(self._poly)

    def poly(self):
        return self._poly

    def is_complex(self):
        return self._complex

    def complex_roots(self, prec, multiplicity):
        """
        EXAMPLES:
            sage: x = polygen(ZZ)
            sage: cp = AA.common_polynomial(x^4 - 2)
            sage: cp.complex_roots(30, 1)
            [[1.1892071150027208 .. 1.1892071150027213], [-1.1892071150027213 .. -1.18920711500272...], [1.18920711500272... .. 1.1892071150027213]*I, [-1.1892071150027213 .. -1.1892071150027208]*I]
        """
        if self._roots_cache.has_key(multiplicity):
            roots = self._roots_cache[multiplicity]
            if roots[0] >= prec:
                return roots[1]

        p = self._poly
        for i in range(multiplicity - 1):
            p = p.derivative()

        from sage.rings.polynomial.complex_roots import complex_roots
        roots_mult = complex_roots(p, min_prec=prec)
        roots = [rt for (rt, mult) in roots_mult if mult == 1]
        self._roots_cache[multiplicity] = (prec, roots)
        return roots

    def exactify(self):
        """
        Compute a common field that holds all of the algebraic coefficients
        of this polynomial, then factor the polynomial over that field.
        Store the factors for later use (ignoring multiplicity).
        """
        if self._exact:
            return

        self._exact = True

        gen = qq_generator

        for c in self._poly.list():
            c.exactify()
            gen = gen.union(c._exact_field())

        self._gen = gen

        coeffs = [gen(c._exact_value()) for c in self._poly.list()]

        if gen.is_trivial():
            qp = QQy(coeffs)
            self._factors = [fac_exp[0] for fac_exp in qp.factor()]
        else:
            fld = gen.field()
            fld_poly = fld['x']

            fp = fld_poly(coeffs)

            self._factors = [fac_exp[0] for fac_exp in fp.factor()]

    def factors(self):
        self.exactify()
        return self._factors

    def generator(self):
        self.exactify()
        return self._gen

class ANRoot(ANDescr):
    """
    The subclass of ANDescr that represents a particular
    root of a polynomial with algebraic coefficients.
    This class is private, and should not be used directly.
    """
    def __init__(self, poly, interval, multiplicity=1):
        if not isinstance(poly, AlgebraicPolynomialTracker):
            poly = AlgebraicPolynomialTracker(poly)
        self._poly = poly
        self._multiplicity = multiplicity
        self._complex = is_ComplexIntervalFieldElement(interval)
        self._complex_poly = poly.is_complex()
        self._interval = self.refine_interval(interval, 64)

    def _repr_(self):
        return 'Root %s of %s'%(self._interval, self._poly)

    def kind(self):
        return 'other'

    def is_complex(self):
        return self._complex

    def conjugate(self, n):
        if not self._complex:
            return self
        if not self._complex_poly:
            return ANRoot(self._poly, self._interval.conjugate(), self._multiplicity)

        raise NotImplementedError

    def refine_interval(self, interval, prec):
        if self._complex or self._complex_poly:
            v = self._complex_refine_interval(interval, prec)
            if self._complex:
                return v
            else:
                return v.real()
        else:
            return self._real_refine_interval(interval, prec)

    def _real_refine_interval(self, interval, prec):
        """
        Takes an interval which is assumed to enclose exactly one root
        of the polynomial (or, with multiplicity=k, exactly one root
        of the k-1'th derivative); and a precision, in bits.

        Tries to find a narrow interval enclosing the root using
        interval arithmetic of the given precision.  (No particular
        number of resulting bits of precision is guaranteed.)

        Uses a combination of Newton's method (adapted for interval
        arithmetic) and bisection.  The algorithm will converge very
        quickly if started with a sufficiently narrow interval.

        EXAMPLES:
            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(AA)
            sage: rt2 = ANRoot(x^2 - 2, RIF(0, 2))
            sage: rt2.refine_interval(RIF(0, 2), 75)
            [1.41421356237309504880163 .. 1.41421356237309504880175]
        """
        # Don't throw away bits in the original interval; doing so might
        # invalidate it (include an extra root)

        field = RealIntervalField(max(prec, interval.prec()))
        interval = field(interval)
        if interval.is_exact():
            return interval

        p = self._poly.poly()
        dp = p.derivative()
        for i in xrange(0, self._multiplicity - 1):
            p = dp
            dp = p.derivative()

        zero = field(0)

        poly_ring = field['x']

        # interval_p = poly_ring(p)
        coeffs = [c._interval_fast(prec) for c in p.list()]
        interval_p = poly_ring(coeffs)

        # This special case is important: this is the only way we could
        # refine "infinitely deep" (we could get an interval of diameter
        # about 2^{-2^31}, and then hit floating-point underflow); avoiding
        # this case here means we don't have to worry about iterating too
        # many times later
        if coeffs[0].is_zero() and interval.contains_zero():
            return zero

        # interval_dp = poly_ring(dp)
        dcoeffs = [c.interval_fast(field) for c in dp.list()]
        interval_dp = poly_ring(dcoeffs)

        linfo = {}
        uinfo = {}

        def update_info(info, x):
            info['endpoint'] = x
            val = interval_p(field(x))
            info['value'] = val
            # sign == 1 if val is bounded greater than 0
            # sign == -1 if val is bounded less than 0
            # sign == 0 if val might be 0
            if val > zero:
                info['sign'] = 1
            elif val < zero:
                info['sign'] = -1
            else:
                info['sign'] = 0

        update_info(linfo, interval.lower())
        update_info(uinfo, interval.upper())

        newton_lower = True

        while True:
            if linfo['sign'] == 0 and uinfo['sign'] == 0:
                # We take it on faith that there is a root in interval,
                # even though we can't prove it at the current precision.
                # We can't do any more refining...
                return interval

            if linfo['sign'] == uinfo['sign']:
                # Oops...
                # print self._poly.poly()
                # print interval_p
                # print linfo['endpoint'], linfo['value'], linfo['sign']
                # print uinfo['endpoint'], uinfo['value'], uinfo['sign']
                raise ValueError, "Refining interval that does not bound unique root!"

            # Use a simple algorithm:
            # Try an interval Newton-Raphson step.  If this does not add at
            # least one bit of information, or if it fails (because the
            # slope is not bounded away from zero), then try bisection.
            # If this fails because the value at the midpoint is not
            # bounded away from zero, then also try the 1/4 and 3/4 points.
            # If all of these fail, then return the current interval.

            slope = interval_dp(interval)

            newton_success = False
            diam = interval.diameter()

            if not (zero in slope):
                # OK, we try Newton-Raphson.
                # I have no idea if it helps, but each time through the loop,
                # we either do Newton-Raphson from the left endpoint or
                # the right endpoint, alternating.

                newton_lower = not newton_lower
                if newton_lower:
                    new_range = linfo['endpoint'] - linfo['value'] / slope
                else:
                    new_range = uinfo['endpoint'] - uinfo['value'] / slope

                if new_range.lower() in interval:
                    interval = field(new_range.lower(), interval.upper())
                    update_info(linfo, interval.lower())
                if new_range.upper() in interval:
                    interval = field(interval.lower(), new_range.upper())
                    update_info(uinfo, interval.upper())

                new_diam = interval.diameter()

                if new_diam == 0:
                    # Wow, we managed to find the answer exactly.
                    # (I think this can only happen with a linear polynomial,
                    # in which case we shouldn't have been in this
                    # function at all; but oh well.)
                    return interval

                if (new_diam << 1) <= diam:
                    # We got at least one bit.
                    newton_success = True

            if not newton_success:
                center = interval.center()

                def try_bisection(mid):
                    minfo = {}
                    update_info(minfo, mid)
                    if minfo['sign'] == 0:
                        return interval, False
                    # We check to make sure the new interval is actually
                    # narrower; this might not be true if the interval
                    # is less than 4 ulp's wide
                    if minfo['sign'] == linfo['sign'] and mid > interval.lower():
                        linfo['endpoint'] = minfo['endpoint']
                        linfo['value'] = minfo['value']
                        linfo['sign'] = minfo['sign']
                        return field(mid, interval.upper()), True
                    if minfo['sign'] == uinfo['sign'] and mid < interval.upper():
                        uinfo['endpoint'] = minfo['endpoint']
                        uinfo['value'] = minfo['value']
                        uinfo['sign'] = minfo['sign']
                        return field(interval.lower(), mid), True
                    return interval, False

                interval, bisect_success = try_bisection(center)

                if not bisect_success:
                    uq = (center + interval.upper()) / 2
                    interval, bisect_success = try_bisection(uq)
                if not bisect_success:
                    lq = (center + interval.lower()) / 2
                    interval, bisect_success = try_bisection(lq)

                if not bisect_success:
                    # OK, we've refined about as much as we can.
                    # (We might be able to trim a little more off the edges,
                    # but the interval is no more than twice as wide as the
                    # narrowest possible.)
                    return interval

    def _complex_refine_interval(self, interval, prec):
        """
        Takes an interval which is assumed to enclose exactly one root
        of the polynomial (or, with multiplicity=k, exactly one root
        of the k-1'th derivative); and a precision, in bits.

        Tries to find a narrow interval enclosing the root using
        interval arithmetic of the given precision.  (No particular
        number of resulting bits of precision is guaranteed.)

        Uses Newton's method (adapted for interval arithmetic).  The
        algorithm will converge very quickly if started with a
        sufficiently narrow interval.  If Newton's method fails, then
        we falls back on computing all the roots of the polynomial
        numerically, and select the appropriate root.

        EXAMPLES:
            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(QQbar)
            sage: intv = CIF(RIF(0, 1), RIF(0.1, 1))
            sage: rt = ANRoot(x^5 - 1, intv)
            sage: new_intv = rt.refine_interval(intv, 53); new_intv
            [0.309016994374947424022 .. 0.309016994374947424185] + [0.951056516295153572002 .. 0.951056516295153572219]*I
            sage: rt.refine_interval(new_intv, 70)
            [0.30901699437494742410148 .. 0.30901699437494742410319] + [0.95105651629515357211565 .. 0.95105651629515357211735]*I
        """
        # Don't throw away bits in the original interval; doing so might
        # invalidate it (include an extra root)

        field = ComplexIntervalField(max(prec, interval.prec()))
        interval = field(interval)
        if interval.is_exact():
            return interval

        p = self._poly.poly()
        dp = p.derivative()
        for i in xrange(0, self._multiplicity - 1):
            p = dp
            dp = p.derivative()

        zero = field(0)

        poly_ring = field['x']

        # interval_p = poly_ring(p)
        coeffs = [c.interval_fast(field) for c in p.list()]
        interval_p = poly_ring(coeffs)

        # This special case is important: this is the only way we could
        # refine "infinitely deep" (we could get an interval of diameter
        # about 2^{-2^31}, and then hit floating-point underflow); avoiding
        # this case here means we don't have to worry about iterating too
        # many times later
        if coeffs[0].is_zero() and zero in interval:
            return zero

        # interval_dp = poly_ring(dp)
        dcoeffs = [c.interval_fast(field) for c in dp.list()]
        interval_dp = poly_ring(dcoeffs)

        while True:
            center = field(interval.center())
            val = interval_p(center)

            slope = interval_dp(interval)

            diam = interval.diameter()

            if zero in slope:
                # Give up and fall back on root isolation.
                return self._complex_isolate_interval(interval, prec)

            if not (zero in slope):
                new_range = center - val / slope
                interval = interval.intersection(new_range)

                new_diam = interval.diameter()

                if new_diam == 0:
                    # Wow; we nailed it exactly.  (This may happen
                    # whenever the root is exactly equal to some
                    # floating-point number, and cannot happen
                    # if the root is not equal to a floating-point
                    # number.)  We just return the perfect answer.
                    return interval

                if new_diam == diam:
                    # We're not getting any better.  There are two
                    # possible reasons for this.  Either we have
                    # refined as much as possible given the imprecision
                    # of our interval polynomial, and we have the best
                    # answer we're going to get at this precision;
                    # or we started with a poor approximation to the
                    # root, resulting in a broad range of possible
                    # slopes in this interval, and Newton-Raphson refining
                    # is not going to help.

                    # I don't have a formal proof, but I believe the
                    # following test differentiates between these two
                    # behaviors.  (If I'm wrong, we might get bad behavior
                    # like infinite loops, but we still won't actually
                    # return a wrong answer.)

                    if val.contains_zero():
                        # OK, center must be a good approximation
                        # to the root (in the current precision, anyway).
                        # And the expression "center - val / slope"
                        # above means that we have a pretty good interval,
                        # even if slope is a poor estimate.
                        return interval

                    # The center of the current interval is known
                    # not to be a root.  This should let us divide
                    # the interval in half, and improve on our previous
                    # estimates.  I can only think of two reasons why
                    # it might not:
                    # 1) the "center" of the interval is actually
                    # on one of the edges of the interval (because the
                    # interval is only one ulp wide), or
                    # 2) the slope estimate is so bad that
                    # "center - val / slope" doesn't give us information.

                    # With complex intervals implemented as
                    # rectangular regions of the complex plane, it's
                    # possible that "val / slope" includes zero even
                    # if both "val" and "slope" are bounded away from
                    # zero, if the diameter of the (interval) argument
                    # of val or slope is large enough.

                    # So we test the diameter of the argument of
                    # slope; if it's small, we decide that we must
                    # have a good interval, but if it's big, we decide
                    # that we probably can't make progress with
                    # Newton-Raphson.

                    # I think the relevant measure of "small" is
                    # whether the diameter is less than pi/2; in that
                    # case, no matter the value of "val" (as long as
                    # "val" is fairly precise), "val / slope" should
                    # be bounded away from zero.  But we compare
                    # against 1 instead, in the hopes that this might
                    # be slightly faster.

                    if slope.argument().absolute_diameter() < 1:
                        return interval

                    # And now it's time to give up.
                    return self._complex_isolate_interval(interval, prec)

    def _complex_isolate_interval(self, interval, prec):
        """
        Find a precise approximation to the unique root in interval,
        by finding a precise approximation to all roots of the
        polynomial, and checking which one is in interval.  Slow but sure.

        EXAMPLES:
            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(QQbar)
            sage: intv = CIF(RIF(0, 1), RIF(0.1, 1))
            sage: rt = ANRoot(x^5 - 1, intv)
            sage: rt._complex_isolate_interval(intv, 53)
            [0.309016994374947424022 .. 0.309016994374947424185] + [0.951056516295153572002 .. 0.951056516295153572219]*I
        """
        rts = self._poly.complex_roots(prec, self._multiplicity)

        # Find all the roots that overlap interval.
        our_root = [rt for rt in rts if rt.overlaps(interval)]

        if len(our_root) == 1:
            return our_root[0]

        if len(our_root) == 0:
            raise ValueError, "Complex root interval does not include any roots"

        # We have more than one root that overlap the current interval.
        # Technically, this might not be an error; perhaps the actual
        # root is just outside our interval, even though the (presumably
        # tight) interval containing that root touches our interval.

        # But it seems far more likely that the provided interval is
        # just too big.

        raise ValueError, "Complex root interval probably includes multiple roots"

    def exactify(self):
        """
        Returns either an ANRational or an
        ANExtensionElement with the same value as this number.

        EXAMPLES:
            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(QQbar)
            sage: two = ANRoot((x-2)*(x-sqrt(QQbar(2))), RIF(1.9, 2.1))
            sage: two.exactify()
            2
            sage: two.exactify().rational_value()
            2
            sage: strange = ANRoot(x^2 + sqrt(QQbar(3))*x - sqrt(QQbar(2)), RIF(-0, 1))
            sage: strange.exactify()
            a where a^8 - 6*a^6 + 5*a^4 - 12*a^2 + 4 = 0 and a in [0.60510122651395104 .. 0.60510122651395116]
        """
        gen = self._poly.generator()

        if gen.is_trivial():
            qpf = self._poly.factors()
            def find_fn(factor, prec):
                return factor(self._interval_fast(prec))
            my_factor = find_zero_result(find_fn, qpf)

            # Factoring always returns monic polynomials over the rationals
            assert(my_factor.is_monic())

            if my_factor.degree() == 1:
                return ANRational(-my_factor[0])

            den, my_factor = clear_denominators(my_factor)

            red_elt, red_back, red_pol = do_polred(my_factor)

            field = NumberField(red_pol, 'a', check=False)

            def intv_fn(rif):
                return conjugate_expand(red_elt(self._interval_fast(rif) * den))
            new_intv = conjugate_shrink(isolating_interval(intv_fn, red_pol))
            root = ANRoot(QQx(red_pol), new_intv)
            new_gen = AlgebraicGenerator(field, root)

            return ANExtensionElement(new_gen, red_back(field.gen())/den)
        else:
            fld = gen.field()

            fpf = self._poly.factors()
            # print fpf
            def find_fn(factor, prec):
                # XXX
                ifield = (ComplexIntervalField if self.is_complex() else RealIntervalField)(prec)
                if_poly = ifield['x']
                gen_val = gen._interval_fast(prec)
                self_val = self._interval_fast(prec)
                ip = if_poly([c.polynomial()(gen_val) for c in factor])
                return ip(self_val)
            my_factor = find_zero_result(find_fn, fpf)

            # print my_factor
            assert(my_factor.is_monic())

            if my_factor.degree() == 1:
                return ANExtensionElement(gen, -my_factor[0])

            # rnfequation needs a monic polynomial with integral coefficients.
            # We achieve this with a change of variables.

            den, my_factor = clear_denominators(my_factor)


            pari_nf = gen.pari_field()

            # print pari_nf[0]
            x, y = QQxy.gens()
            my_factor = QQxy['z']([c.polynomial()(y) for c in my_factor])(x)
            # print my_factor

            # XXX much duplicate code with AlgebraicGenerator.union()

            # XXX need more caching here
            newpol, self_pol, k = pari_nf.rnfequation(my_factor, 1)
            k = int(k)
            # print newpol
            # print self_pol
            # print k

            newpol_sage = QQx(newpol)
            newpol_sage_y = QQy(newpol_sage)

            red_elt, red_back, red_pol = do_polred(newpol_sage_y)

            new_nf = NumberField(red_pol, name='a', check=False)

            self_pol_sage = QQx(self_pol.lift())

            new_nf_a = new_nf.gen()

            def intv_fn(prec):
                return conjugate_expand(red_elt(gen._interval_fast(prec) * k + self._interval_fast(prec) * den))
            new_intv = conjugate_shrink(isolating_interval(intv_fn, red_pol))

            root = ANRoot(QQx(red_pol), new_intv)
            new_gen = AlgebraicGenerator(new_nf, root)
            red_back_a = red_back(new_nf.gen())
            new_poly = ((QQx_x - k * self_pol_sage)(red_back_a)/den)
            return ANExtensionElement(new_gen, new_poly)

    def _more_precision(self):
        """
        Recompute the interval enclosing this root at higher
        precision.
        """
        prec = self._interval.prec()
        self._interval = self.refine_interval(self._interval, prec*2)

    def _interval_fast(self, prec):
        """
        Given a RealIntervalField, compute the value of this number
        using interval arithmetic of at least the precision of the field,
        and return the value in that field.  (More precision may be used
        in the computation.)
        """
        if prec == self._interval.prec():
            return self._interval
        if prec < self._interval.prec():
            return type(self._interval.parent())(prec)(self._interval)
        self._more_precision()
        return self._interval_fast(prec)

qq_generator = AlgebraicGenerator(None, ANRoot(AAPoly.gen(), RIF(0)))

_cyclotomic_gen_cache = {}
def cyclotomic_generator(n):
    try:
        return _cyclotomic_gen_cache[n]
    except KeyError:
        assert(n > 2 and n != 4)
        n = ZZ(n)
        f = CyclotomicField(n)
        v = ANRootOfUnity(~n, QQ_1)
        g = AlgebraicGenerator(f, v)
        g.set_cyclotomic(n)
        _cyclotomic_gen_cache[n] = g
        return g

class ANExtensionElement(ANDescr):
    """
    The subclass of ANDescr that represents a number field
    element in terms of a specific generator.  Consists of a polynomial
    with rational coefficients in terms of the generator, and the
    generator itself, an AlgebraicGenerator.
    """

    def __new__(self, generator, value):
        try:
            return ANRational(value._rational_())
        except TypeError:
            if generator is QQbar_I_generator and value[0].is_zero():
                return ANRootOfUnity(QQ_1_4, value[1])
            return ANDescr.__new__(self)

    def __init__(self, generator, value):
        self._generator = generator
        self._value = value
        self._exactly_real = not generator.is_complex()

    def _repr_(self):
        return '%s where %s = 0 and a in %s'%(self._value,
                                              self._generator.field().polynomial()._repr(name='a'),
                                              self._generator._interval_fast(53))

    def kind(self):
        if self._generator is QQbar_I_generator:
            return 'gaussian'
        else:
            return 'element'

    def is_complex(self):
        return not self._exactly_real

    def is_exact(self):
        return True

    def is_field_element(self):
        return True

    def generator(self):
        return self._generator

    def exactify(self):
        return self

    def field_element_value(self):
        return self._value

    def minpoly(self):
        return self._value.minpoly()

    def _interval_fast(self, prec):
        gen_val = self._generator._interval_fast(prec)
        v = self._value.polynomial()(gen_val)
        if self._exactly_real and is_ComplexIntervalFieldElement(v):
            return v.real()
        return v

    def neg(self, n):
        return ANExtensionElement(self._generator, -self._value)

    def invert(self, n):
        return ANExtensionElement(self._generator, ~self._value)

    def conjugate(self, n):
        if self._exactly_real:
            return self
        else:
            return ANExtensionElement(self._generator.conjugate(), self._value)

    def norm(self, n):
        if self._exactly_real:
            return (n*n)._descr
        elif self._generator is QQbar_I_generator:
            return ANRational(self._value.norm())
        else:
            return ANUnaryExpr(n, 'norm')

    def abs(self, n):
        return AlgebraicReal(self.norm(n)).sqrt()._descr

    def rational_argument(self, n):
        """
        If the argument of self is 2*pi times some rational number,
        return that rational; otherwise, return None.
        """

        # If the argument of self is 2*pi times some rational number a/b,
        # then self/abs(self) is a root of the b'th cyclotomic polynomial.
        # This implies that the algebraic degree of self is at least
        # phi(b).  Working backward, we know that the algebraic degree
        # of self is at most the degree of the generator, so that gives
        # an upper bound on phi(b).  According to
        # http://mathworld.wolfram.com/TotientFunction.html,
        # phi(b) >= sqrt(b) for b > 6; this gives us an upper bound on b.
        # We then check to see if this is possible; if so, we test
        # to see if it actually holds.

        if self._exactly_real:
            if n > 0:
                return 0
            else:
                return QQ(1)/2

        gen_degree = self._generator._field.degree()
        if gen_degree <= 2:
            max_b = 6
        else:
            max_b = gen_degree*gen_degree
        rat_arg_fl = ComplexIntervalField(100)(n).argument() / RealIntervalField(100).pi() / 2
        rat_arg = rat_arg_fl.simplest_rational()
        if rat_arg.denominator() > max_b:
            return None
        n_exp = n ** rat_arg.denominator()
        if n_exp.real() > 0 and n_exp.imag() == 0:
            return rat_arg
        # Strictly speaking, we need to look for the second-simplest
        # rational in rat_arg_fl and make sure its denominator is > max_b.
        # For now, we just punt.
        raise NotImplementedError

    def gaussian_value(self):
        assert(self._generator is QQbar_I_generator)
        return self._value

class ANUnaryExpr(ANDescr):
    def __init__(self, arg, op):
        self._arg = arg
        self._op = op

    def kind(self):
        return 'other'

    def _interval_fast(self, prec):
        op = self._op

        if op == '-':
            return -self._arg._interval_fast(prec)

        if op == '~':
            return ~self._arg._interval_fast(prec)

        if op == 'real':
            v = self._arg._interval_fast(prec)
            if is_ComplexIntervalFieldElement(v):
                return v.real()
            else:
                return v

        if op == 'imag':
            v = self._arg._interval_fast(prec)
            if is_ComplexIntervalFieldElement(v):
                return v.imag()
            else:
                return RealIntervalField(prec)(0)

        if op == 'abs':
            return abs(self._arg._interval_fast(prec))

        if op == 'norm':
            v = self._arg._interval_fast(prec)
            if is_ComplexIntervalFieldElement(v):
                return v.norm()
            else:
                return v.square()

        if op == 'conjugate':
            v = self._arg._interval_fast(prec)
            if is_ComplexIntervalFieldElement(v):
                return v.conjugate()
            else:
                return v

        raise NotImplementedError

    def exactify(self):
        op = self._op
        arg = self._arg

        if op == '-':
            arg.exactify()
            return arg._descr.neg(None)

        if op == '~':
            arg.exactify()
            return arg._descr.invert(None)

        if op == 'real':
            arg.exactify()
            rv = (arg + arg.conjugate()) / 2
            rv.exactify()
            rvd = rv._descr
            rvd._exactly_real = True
            return rvd

        if op == 'imag':
            arg.exactify()
            iv = QQbar_I * (arg.conjugate() - arg) / 2
            iv.exactify()
            ivd = iv._descr
            ivd._exactly_real = True
            return ivd

        if op == 'abs':
            arg.exactify()
            if arg.parent() is AA:
                if arg.sign() > 0:
                    return arg._descr
                else:
                    return arg._descr.neg(None)

            v = (arg * arg.conjugate()).sqrt()
            v.exactify()
            vd = v._descr
            vd._exactly_real = True
            return vd

        if op == 'norm':
            arg.exactify()
            v = arg * arg.conjugate()
            v.exactify()
            vd = v._descr
            vd._exactly_real = True
            return vd

        if op == 'conjugate':
            arg.exactify()
            return arg._descr.conjugate(None)

class ANBinaryExpr(ANDescr):
    def __init__(self, left, right, op):
        self._left = left
        self._right = right
        self._op = op
        self._complex = True

    def kind(self):
        return 'other'

    def is_complex(self):
        return self._complex

    def _interval_fast(self, prec):
        op = self._op

        lv = self._left._interval_fast(prec)
        rv = self._right._interval_fast(prec)

        if not (is_ComplexIntervalFieldElement(lv) or is_ComplexIntervalFieldElement(rv)):
            self._complex = False

        if op == '-':
            return lv - rv

        if op == '+':
            return lv + rv

        if op == '/':
            return lv / rv

        if op == '*':
            return lv * rv

        raise NotImplementedError

    def exactify(self):
        left = self._left
        right = self._right
        left.exactify()
        right.exactify()
        gen = left._exact_field().union(right._exact_field())
        left_value = gen(left._exact_value())
        right_value = gen(right._exact_value())

        op = self._op

        if op == '+':
            value = left_value + right_value
        if op == '-':
            value = left_value - right_value
        if op == '*':
            value = left_value * right_value
        if op == '/':
            value = left_value / right_value

        if gen.is_trivial():
            return ANRational(value)
        else:
            return ANExtensionElement(gen, value)

ax = QQbarPoly.gen()
# def heptadecagon():
#     # Compute the exact (x,y) coordinates of the vertices of a 34-gon.
#     # (Take every other coordinate to get the vertices of a
#     # heptadecagon.)
#     # Formulas from:
#     # Weisstein, Eric W. "Trigonometry Angles--Pi/17." From
#     # MathWorld--A Wolfram Web Resource.
#     # http://mathworld.wolfram.com/TrigonometryAnglesPi17.html

#     rt17 = AA(17).sqrt()
#     rt2 = AA(2).sqrt()
#     eps = (17 + rt17).sqrt()
#     epss = (17 - rt17).sqrt()
#     delta = rt17 - 1
#     alpha = (34 + 6*rt17 + rt2*delta*epss - 8*rt2*eps).sqrt()
#     beta = 2*(17 + 3*rt17 - 2*rt2*eps - rt2*epss).sqrt()
#     x = rt2*(15 + rt17 + rt2*(alpha + epss)).sqrt()/8
#     y = rt2*(epss**2 - rt2*(alpha + epss)).sqrt()/8

#     cx, cy = 1, 0
#     for i in range(34):
#         cx, cy = x*cx-y*cy, x*cy+y*cx
#     print cx, cy
#     print cx.sign(), cy.sign()
#     print "Yo!"
#     print (cx-1).sign()
#     print "OK"
# # heptadecagon()

# def heptadecagon2():
#     # Compute the exact (x,y) coordinates of the vertices of a 34-gon.
#     # (Take every other coordinate to get the vertices of a
#     # heptadecagon.)
#     # Formulas from:
#     # Weisstein, Eric W. "Heptadecagon." From MathWorld--A Wolfram
#     # Web Resource. http://mathworld.wolfram.com/Heptadecagon.html

#     x = AA.polynomial_root(256*ax**8 - 128*ax**7 - 448*ax**6 + 192*ax**5 + 240*ax**4 - 80*ax**3 - 40*ax**2 + 8*ax + 1, RIF(0.9829, 0.983))
#     y = (1-x**2).sqrt()

#     cx, cy = 1, 0
#     for i in range(34):
#         cx, cy = x*cx-y*cy, x*cy+y*cx
#     print cx, cy
#     print cx.sign(), cy.sign()
#     print (cx-1).sign()
#     return x, y
# # heptadecagon2()


RR_1_10 = RR(ZZ(1)/10)
QQ_0 = QQ(0)
QQ_1 = QQ(1)
QQ_1_2 = QQ(1)/2
QQ_1_4 = QQ(1)/4

QQbar_I_nf = QuadraticField(-1, 'I')
# XXX change ANRoot to ANRootOfUnity below
QQbar_I_generator = AlgebraicGenerator(QQbar_I_nf, ANRoot(AAPoly.gen()**2 + 1, CIF(0, 1)))
QQbar_I = AlgebraicNumber(ANExtensionElement(QQbar_I_generator, QQbar_I_nf.gen()))
_cyclotomic_gen_cache[4] = QQbar_I_generator
QQbar_I_generator.set_cyclotomic(4)

AA_hash_offset = AA(~ZZ(123456789))

QQbar_hash_offset = AlgebraicNumber(ANExtensionElement(QQbar_I_generator, ~ZZ(123456789) + QQbar_I_nf.gen()/ZZ(987654321)))

ZZX_x = ZZ['x'].gen()

# This is used in the _algebraic_ method of the golden_ratio constant,
# in sage/functions/constants.py
AA_golden_ratio = None

def get_AA_golden_ratio():
    global AA_golden_ratio
    if AA_golden_ratio is None:
        AA_golden_ratio_nf = NumberField(ZZX_x**2 - ZZX_x - 1, 'phi')
        AA_golden_ratio_generator = AlgebraicGenerator(AA_golden_ratio_nf, ANRoot(AAPoly.gen()**2 - AAPoly.gen() - 1, RIF(1.618, 1.6181)))
        AA_golden_ratio = AlgebraicReal(ANExtensionElement(AA_golden_ratio_generator, AA_golden_ratio_nf.gen()))
    return AA_golden_ratio
