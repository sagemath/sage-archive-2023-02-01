"""
Field of Algebraic Reals

AUTHOR:
    -- Carl Witty (2007-01-27): initial version

This is an implementation of the algebraic reals (the real numbers which
are the zero of a polynomial in ZZ[x]).  All computations are exact.

As with many other implementations of the algebraic numbers, we avoid
computing a number field and working in the number field whenever possible.

Algebraic numbers exist in one of the following forms:

* a rational number
* the sum, difference, product, or quotient of algebraic numbers
* a particular root of a polynomial, given as a polynomial with
algebraic coefficients, an isolating interval (given as a
RealIntervalFieldElement) which encloses exactly one root, and
the multiplicity of the root
* a polynomial in a generator, where the generator is an algebraic
number given as the root of an irreducible polynomial with integral
coefficients and the polynomial is given as a NumberFieldElement

An algebraic number can be coerced into RealIntervalField; every algebraic
number has a cached interval of the highest precision yet calculated.

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
to halp, we keep a lattice of already-computed number fields and
their inclusions.

EXAMPLES:
    sage: sqrt(Alg(2)) > 0
    True
    sage: (sqrt(5 + 2*sqrt(Alg(6))) - sqrt(Alg(3)))**2 == 2
    True

For a monic cubic polynomial x^3 + b*x^2 + c*x + d with roots
s1,s2,s3, the discriminant is defined as (s1-s2)^2(s1-s3)^2(s2-s3)^2
and can be computed as b^2c^2 - 4*b^3d - 4*c^3 + 18*bcd - 27*d^2.
We can test that these definitions do give the same result.

    sage: def disc1(b, c, d):
    ...       return b^2*c^2 - 4*b^3*d - 4*c^3 + 18*b*c*d - 27*d^2
    sage: def disc2(s1, s2, s3):
    ...       return ((s1-s2)*(s1-s3)*(s2-s3))^2
    sage: x = polygen(Alg)
    sage: p = x*(x-2)*(x-4)
    sage: cp = Alg.common_polynomial(p)
    sage: d, c, b, _ = p.list()
    sage: s1 = Alg.polynomial_root(cp, RIF(-1, 1))
    sage: s2 = Alg.polynomial_root(cp, RIF(1, 3))
    sage: s3 = Alg.polynomial_root(cp, RIF(3, 5))
    sage: disc1(b, c, d) == disc2(s1, s2, s3)
    True
    sage: p = p + 1
    sage: cp = Alg.common_polynomial(p)
    sage: d, c, b, _ = p.list()
    sage: s1 = Alg.polynomial_root(cp, RIF(-1, 1))
    sage: s2 = Alg.polynomial_root(cp, RIF(1, 3))
    sage: s3 = Alg.polynomial_root(cp, RIF(3, 5))
    sage: disc1(b, c, d) == disc2(s1, s2, s3)
    True
    sage: p = (x-sqrt(Alg(2)))*(x-Alg(2).nth_root(3))*(x-sqrt(Alg(3)))
    sage: cp = Alg.common_polynomial(p)
    sage: d, c, b, _ = p.list()
    sage: s1 = Alg.polynomial_root(cp, RIF(1.4, 1.5))
    sage: s2 = Alg.polynomial_root(cp, RIF(1.7, 1.8))
    sage: s3 = Alg.polynomial_root(cp, RIF(1.2, 1.3))
    sage: disc1(b, c, d) == disc2(s1, s2, s3)
    True

Some computation with radicals:

    sage: phi = (1 + sqrt(Alg(5))) / 2
    sage: phi^2 == phi + 1
    True
    sage: tau = (1 - sqrt(Alg(5))) / 2
    sage: tau^2 == tau + 1
    True
    sage: phi + tau == 1
    True
    sage: tau < 0
    True

    sage: rt23 = sqrt(Alg(2/3))
    sage: rt35 = sqrt(Alg(3/5))
    sage: rt25 = sqrt(Alg(2/5))
    sage: rt23 * rt35 == rt25
    True

Algebraic reals which are known to be rational print as rationals; otherwise
they print as intervals (with 53-bit precision).

    sage: Alg(2)/3
    2/3
    sage: Alg(5/7)
    5/7
    sage: two = sqrt(Alg(4)); two
    [1.9999999999999997 ... 2.0000000000000005]
    sage: two == 2; two
    True
    2
    sage: phi
    [1.6180339887498946 ... 1.6180339887498950]

The paper _ARPREC: An Arbitrary Precision Computation Package_ discusses
this result.  Evidently it is difficult to find, but we can easily
verify it.

    sage: alpha = Alg.polynomial_root(x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1, RIF(1, 1.2))
    sage: lhs = alpha^630 - 1
    sage: rhs_num = (alpha^315 - 1) * (alpha^210 - 1) * (alpha^126 - 1)^2 * (alpha^90 - 1) * (alpha^3 - 1)^3 * (alpha^2 - 1)^5 * (alpha - 1)^3
    sage: rhs_den = (alpha^35 - 1) * (alpha^15 - 1)^2 * (alpha^14 - 1)^2 * (alpha^5 - 1)^6 * alpha^68
    sage: rhs = rhs_num / rhs_den
    sage: lhs
    [2.6420403358193507e44 ... 2.6420403358193520e44]
    sage: rhs
    [2.6420403358193507e44 ... 2.6420403358193520e44]
    sage: lhs - rhs
    [-62883485433074552000000000000 ... 63753912023197085000000000000]
    sage: lhs == rhs
    True
    sage: lhs - rhs
    [0.00000000000000000 ... 0.00000000000000000]
    sage: lhs._exact_value()
    -242494609856316402264822833062350847769474540*a^9 + 862295472068289472491654837785947906234680703*a^8 - 829559238431038252116584538075753012193290520*a^7 - 125882239615006638366472766103700441555126185*a^6 + 1399067970863104691667276008776398309383579345*a^5 - 1561176687069361567616835847286958553574223422*a^4 + 761706318888840943058230840550737823821027895*a^3 + 580740464974951394762758666210754821723780266*a^2 - 954587496403409756503464154898858512440951323*a + 546081123623099782018260884934770383777092602 where a^10 - 4*a^9 + 5*a^8 - a^7 - 6*a^6 + 9*a^5 - 6*a^4 - a^3 + 5*a^2 - 4*a + 1 = 0 and a in [0.44406334400909258 ... 0.44406334400909265]
"""

import sage.rings.ring
from sage.structure.sage_object import SageObject
from sage.structure.parent_gens import ParentWithGens
from sage.rings.real_mpfr import RR
from sage.rings.real_mpfi import RealIntervalField, RIF, RealIntervalFieldElement
from sage.rings.polynomial_ring import PolynomialRing
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field import NumberField
from sage.rings.number_field.number_field_element import is_NumberFieldElement
from sage.rings.arith import factor
from sage.libs.pari.gen import pari

# Singleton object implementation copied from integer_ring.py
_obj = None
class _uniq_alg(object):
    def __new__(cls):
        global _obj
        if _obj is None:
            _obj = sage.rings.ring.Field.__new__(cls)
        return _obj

class AlgebraicField(_uniq_alg, sage.rings.ring.Field):
    r"""
    The field of algebraic reals.

    """

    def __init__(self):
        ParentWithGens.__init__(self, self, ('x',), normalize=False)
        self._default_interval_field = RealIntervalField(64)

    def _repr_(self):
        return "Algebraic Field"

    # Is there a standard representation for this?
    def _latex_(self):
        return "\\mathbf{A}"

    def __call__(self, x):
        """
        Coerce x into the field of algebraic numbers.

        """

        if isinstance(x, AlgebraicNumber):
            return x
        return AlgebraicNumber(x)

    def _coerce_impl(self, x):
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            return self(x)
        raise TypeError, 'no implicit coercion of element to the algebraic numbers'

    def _is_valid_homomorphism_(self, codomain, im_gens):
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def default_interval_field(self):
        return self._default_interval_field

    def gens(self):
        return (self(1), )

    def gen(self, n=0):
        if n == 0:
            return self(1)
        else:
            raise IndexError, "n must be 0"

    def ngens(self):
        return 1

    def is_finite(self):
        return False

    def is_atomic_repr(self):
        return True

    def characteristic(self):
        return sage.rings.integer.Integer(0)

    def order(self):
        return infinity.infinity

    def zeta(self, n=2):
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        else:
            raise ValueError, "no n-th root of unity in algebraic reals"

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
            sage: x = polygen(Alg)
            sage: p = Alg.common_polynomial(x^2 - x - 1)
            sage: phi = Alg.polynomial_root(p, RIF(0, 2))
            sage: tau = Alg.polynomial_root(p, RIF(-2, 0))
            sage: phi + tau == 1
            True
            sage: phi * tau == -1
            True
        """
        return AlgebraicPolynomialTracker(poly)

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
        polynomial, it is better to use Alg.common_polynomial
        to get a shared polynomial.

        EXAMPLES:
            sage: x = polygen(Alg)
            sage: phi = Alg.polynomial_root(x^2 - x - 1, RIF(0, 2)); phi
            [1.6180339887498946 ... 1.6180339887498950]
            sage: p = (x-1)^7 * (x-2)
            sage: r = Alg.polynomial_root(p, RIF(9/10, 11/10), multiplicity=7)
            sage: r; r == 1
            [0.99999999999999988 ... 1.0000000000000003]
            True
            sage: p = (x-phi)*(x-sqrt(Alg(2)))
            sage: r = Alg.polynomial_root(p, RIF(1, 3/2))
            sage: r; r == sqrt(Alg(2))
            [1.4142135623730949 ... 1.4142135623730952]
            True
        """
        return AlgebraicNumber(AlgebraicNumberRoot(poly, interval, multiplicity))

def is_AlgebraicField(F):
    return isinstance(F, AlgebraicField)

Alg = AlgebraicField()

def rif_seq():
    # XXX Should do some testing to see where the efficiency breaks are
    # in MPFR.
    bits = 64
    while True:
        yield RealIntervalField(bits)
        bits = bits * 2
#         # XXX temporary debugging:
#         if bits > 1024:
#             break

def clear_denominators(poly):
    """
    Takes a monic polynomial and rescales the variable to get a monic
    polynomial with "integral" coefficients.  Works on any univariate
    polynomial whose base ring has a denominator() method that returns
    integers; for example, the base ring might be QQ or a number
    field.

    Returns the scale factor and the new polynomial.

    (Inspired by Pari's primitive_pol_to_monic().)

    EXAMPLES:
        sage: from sage.rings.algebraic import clear_denominators

        sage: _.<x> = QQ['x']
        sage: clear_denominators(x + 3/2)
        (2, x + 3)
        sage: clear_denominators(x^2 + x/2 + 1/4)
        (2, x^2 + x + 1)

    """

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
        sage: from sage.rings.algebraic import do_polred

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
    return parent(best_elt), parent(best_elt.Mod(pari_poly).modreverse().lift()), parent(best)

def isolating_interval(intv_fn, pol):
    """
    intv_fn is a function that takes a RealIntervalField and returns an
    interval containing some particular root of pol.  (It must return
    better approximations as the RealIntervalField precision increases.)
    pol is an irreducible polynomial with rational coefficients.

    Returns an interval containing at most one root of pol.

    EXAMPLES:
        sage: from sage.rings.algebraic import isolating_interval

        sage: _.<x> = QQ['x']
        sage: isolating_interval(lambda rif: sqrt(rif(2)), x^2 - 2)
        [1.41421356237309504876 ... 1.41421356237309504888]

    And an example that requires more precision:
        sage: delta = 10^(-70)
        sage: p = (x - 1) * (x - 1 - delta) * (x - 1 + delta)
        sage: isolating_interval(lambda rif: rif(1 + delta), p)
        [1.00000000000000000000000000000000000000000000000000000000000000000000009999999999999999999999999999999999999999999999999999999999999999999999999999999999998 ... 1.00000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000014]
    """
    for rif in rif_seq():
        intv = intv_fn(rif)

        # We need to verify that pol has exactly one root in the
        # interval intv.  We know (because it is a precondition of
        # calling this function) that it has at least one root in the
        # interval, so we only need to verify that it has at most one
        # root (that the interval is sufficiently narrow).

        # Call the root we want alpha, and the bounding interval
        # [bot ... top].  By Budan's theorem, if pol(x+bot)
        # has at most one more sign change than pol(x+top),
        # then pol has at most one root in this interval.
        # This is true if all coefficients in pol(x+bot)
        # except for the constant coefficient have the same
        # sign as the corresponding coefficients in pol(x+top).

        # We can check this by computing pol(x+new_intv).
        # This gives an "interval polynomial" which represents
        # a set of polynomials, including pol(x+bot) and
        # pol(x+top).  If all coefficients except the constant
        # coefficient of this interval polynomial are bounded away from
        # zero, then pol has at most one root in intv.

        # This is a sufficient condition, but not a necessary one.
        # It is not obvious that there always exists a bounding
        # interval for alpha with this property.  Fortunately,
        # in the current case with pol irreducible, such
        # a bounding interval does exist.

        # Consider the polynomial pol(x+alpha), where
        # pol has degree n.  The coefficient of x^k in this
        # polynomial has an exact representation as a polynomial
        # in QQ[alpha] of degree (n-k).  If this coefficient is
        # equal to zero, then that means that alpha is a zero
        # of this degree (n-k) polynomial.  But if k>0, then this
        # is impossible, since alpha is a root of an irreducible
        # degree n polynomial.

        # Thus, all coefficients of pol(x+alpha) are nonzero,
        # except for the constant coefficient, which is zero.  We
        # are guaranteed that if we select a sufficiently tight bounding
        # interval around alpha, then pol(x+[bot...top])
        # will have coefficients bounded away from zero, which proves
        # that there is only one root within that bounding interval.

        zero = rif(0)
        rif_poly_ring = rif['x']

        ip = pol(rif_poly_ring.gen() + intv)
        coeffs = ip.list()
        assert(zero in coeffs[0])

        ok = True

        for c in coeffs[1:]:
            if zero in c:
                ok = False
                break
        if ok:
            return intv

def find_zero_result(fn, l):
    """
    l is a list of some sort.
    fn is a function which maps an element of l and a RealIntervalField
    into an interval, such that for a sufficiently precise RealIntervalField,
    exactly one element of l results in an interval containing 0.
    Returns that one element of l.

    EXAMPLES:
        sage: from sage.rings.algebraic import find_zero_result
        sage: _.<x> = QQ['x']
        sage: delta = 10^(-70)
        sage: p1 = x - 1
        sage: p2 = x - 1 - delta
        sage: p3 = x - 1 + delta
        sage: p2 == find_zero_result(lambda p, rif: p(rif(1 + delta)), [p1, p2, p3])
        True
    """
    for rif in rif_seq():
        result = None
        ambig = False
        for v in l:
            intv = fn(v, rif)
            if 0 in intv:
                if result is not None:
                    ambig = True
                    break
                result = v
        if ambig:
            continue
        if result is None:
            raise ValueError, 'find_zero_result could not find any zeroes'
        return result

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
    An AlgebraicGenerator represents both an algebraic real alpha and
    the number field QQ[alpha].  There is a single AlgebraicGenerator
    representing QQ (with alpha==1).

    The AlgebraicGenerator class is private, and should not be used
    directly.
    """

    def __init__(self, field, root):
        """
        Construct an AlgebraicGenerator object.

        sage: from sage.rings.algebraic import AlgebraicNumberRoot, AlgebraicGenerator, unit_generator
        sage: _.<y> = QQ['y']
        sage: x = polygen(Alg)
        sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
        sage: root = AlgebraicNumberRoot(x^2 - x - 1, RIF(1, 2))
        sage: x = AlgebraicGenerator(nf, root)
        sage: x
        Number Field in a with defining polynomial y^2 - y - 1 with a in [1.6180339887498946 ... 1.6180339887498950]
        sage: x.field()
        Number Field in a with defining polynomial y^2 - y - 1
        sage: x.is_unit()
        False
        sage: x.union(unit_generator) is x
        True
        sage: unit_generator.union(x) is x
        True
        """
        self._field = field
        self._unit = (field is None)
        self._root = root
        self._unions = {}
        global algebraic_generator_counter
        self._index = algebraic_generator_counter
        algebraic_generator_counter += 1

    def __hash__(self):
        return self._index

    def __cmp__(self, other):
        return cmp(self._index, other._index)

    def _repr_(self):
        if self._unit:
            return 'Unit generator'
        else:
            return '%s with a in %s'%(self._field, self._root.interval_fast(RIF))

    def is_unit(self):
        """
        Returns true iff this is the unit generator (alpha == 1), which
        does not actually extend the rationals.

        EXAMPLES:
            sage: from sage.rings.algebraic import unit_generator
            sage: unit_generator.is_unit()
            True
        """
        return self._unit

    def field(self):
        return self._field

    def interval_fast(self, field):
        """
        Returns an interval containing this generator, to the specified
        precision.
        """
        return self._root.interval_fast(field)

    def union(self, other):
        """
        Given generators alpha and beta, alpha.union(beta) gives a generator
        for the number field QQ[alpha][beta].

        EXAMPLES:
            sage: from sage.rings.algebraic import AlgebraicNumberRoot, AlgebraicGenerator, unit_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(Alg)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = AlgebraicNumberRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in [1.4142135623730949 ... 1.4142135623730952]
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = AlgebraicNumberRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in [1.7320508075688771 ... 1.7320508075688775]
            sage: gen2.union(unit_generator) is gen2
            True
            sage: unit_generator.union(gen3) is gen3
            True
            sage: gen2.union(gen3)
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in [0.51763809020504147 ... 0.51763809020504159]
        """
        if self._unit:
            return other
        if other._unit:
            return self
        if self is other:
            return self
        if other in self._unions:
            return self._unions[other].parent
        if self._field.polynomial().degree() < other._field.polynomial().degree():
            self, other = other, self
        sp = self._field.polynomial()
        op = other._field.polynomial()
        op = QQx(op)
        # print sp
        # print op
        # print self._field.polynomial()
        # print self._field.polynomial().degree()
        pari_nf = self._field.pari_nf()
        # print pari_nf[0]
        factors = list(pari_nf.nffactor(op).lift())[0]
        # print factors
        x, y = QQxy.gens()
        # XXX Go through strings because multivariate polynomial coercion from
        # Pari is not implemented
        factors_sage = [QQxy(str(p)) for p in factors]
        # print factors_sage
        def find_fn(p, rif):
            rif_poly = rif['x', 'y']
            ip = rif_poly(p)
            return ip(other._root.interval_fast(rif), self._root.interval_fast(rif))
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

        def intv_fn(rif):
            return red_elt(self._root.interval_fast(rif) * k + other._root.interval_fast(rif))
        new_intv = isolating_interval(intv_fn, red_pol)

        new_gen = AlgebraicGenerator(new_nf, AlgebraicNumberRoot(QQx(red_pol), new_intv))
        rel = AlgebraicGeneratorRelation(self, self_pol_sage(red_back_x),
                                         other, (QQx_x - k*self_pol_sage)(red_back_x),
                                         new_gen)
        self._unions[other] = rel
        other._unions[self] = rel
        return rel.parent

    def super_poly(self, super, checked=None):
        """
        Given a generator gen and another generator super, where super
        is the result of a tree of union() operations where one of the
        leaves is gen, gen.super_poly(super) returns a polynomial
        expressing the value of gen in terms of the value of super.
        (Except that if gen is unit_generator, super_poly() always
        returns None.)

        EXAMPLES:
            sage: from sage.rings.algebraic import AlgebraicGenerator, AlgebraicNumberRoot, unit_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(Alg)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = AlgebraicNumberRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in [1.4142135623730949 ... 1.4142135623730952]
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = AlgebraicNumberRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in [1.7320508075688771 ... 1.7320508075688775]
            sage: gen2_3 = gen2.union(gen3)
            sage: gen2_3
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in [0.51763809020504147 ... 0.51763809020504159]
            sage: unit_generator.super_poly(gen2) is None
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
        Takes an AlgebraicNumber which is represented as either a rational
        or a number field element, and which is in a subfield of the
        field generated by this generator.  Lifts the number into the
        field of this generator, and returns either a Rational or a
        NumberFieldElement depending on whether this is the unit generator.

        EXAMPLES:
            sage: from sage.rings.algebraic import AlgebraicNumberRoot, AlgebraicGenerator, AlgebraicNumberExtensionElement, AlgebraicNumberRational
            sage: _.<y> = QQ['y']
            sage: x = polygen(Alg)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = AlgebraicNumberRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in [1.4142135623730949 ... 1.4142135623730952]
            sage: sqrt2 = AlgebraicNumberExtensionElement(gen2, nf2.gen())
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = AlgebraicNumberRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in [1.7320508075688771 ... 1.7320508075688775]
            sage: sqrt3 = AlgebraicNumberExtensionElement(gen3, nf3.gen())
            sage: gen2_3 = gen2.union(gen3)
            sage: gen2_3
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in [0.51763809020504147 ... 0.51763809020504159]
            sage: gen2_3(sqrt2)
            -a^3 + 3*a
            sage: gen2_3(AlgebraicNumberRational(1/7))
            1/7
            sage: gen2_3(sqrt3)
            -a^2 + 2
        """
        if self._unit:
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

class AlgebraicNumberDescr(SageObject):
    """
    An AlgebraicNumber is a wrapper around an AlgebraicNumberDescr object.
    AlgebraicNumberDescr is an abstract base class, which should never
    be directly instantiated; its concrete subclasses are
    AlgebraicNumberRational, AlgebraicNumberExpression,
    AlgebraicNumberRoot, and AlgebraicNumberExtensionElement.
    AlgebraicNumberDescr and all of its subclasses are private, and
    should not be used directly.
    """
    def is_exact(self):
        """
        Returns True if self is an AlgebraicNumberRational or an
        AlgebraicNumberExtensionElement.

        EXAMPLES:
            sage: from sage.rings.algebraic import AlgebraicNumberRational
            sage: AlgebraicNumberRational(1/2).is_exact()
            True
        """
        return False

    def is_rational(self):
        """
        Returns True if self is an AlgebraicNumberRational or
        an AlgebraicNumberExtensionElement which is actually rational.

        EXAMPLES:
            sage: from sage.rings.algebraic import AlgebraicNumberRational
            sage: AlgebraicNumberRational(3/7).is_rational()
            True
        """
        return False

    def is_field_element(self):
        """
        Returns True if self is an AlgebraicNumberExtensionElement.

            sage: from sage.rings.algebraic import AlgebraicNumberExtensionElement, AlgebraicNumberRoot, AlgebraicGenerator
            sage: _.<y> = QQ['y']
            sage: x = polygen(Alg)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = AlgebraicNumberRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: sqrt2 = AlgebraicNumberExtensionElement(gen2, nf2.gen())
            sage: sqrt2.is_field_element()
            True
        """
        return False

class AlgebraicNumberRational(AlgebraicNumberDescr):
    """
    The subclass of AlgebraicNumberDescr that represents an arbitrary
    rational.  This class is private, and should not be used directly.
    """

    def __init__(self, x):
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            self._value = x
        else:
            raise TypeError, "Illegal initializer for algebraic number rational"

    def _repr_(self):
        return repr(self._value)

    def interval_fast(self, rif):
        return rif(self._value)

    def is_rational(self):
        return True

    def rational_value(self):
        return self._value

    def exactify(self):
        return self

    def is_exact(self):
        return True

class AlgebraicNumberExpression(AlgebraicNumberDescr):
    """
    The subclass of AlgebraicNumberDescr that represents the sum,
    difference, product, or quotient of two algebraic numbers.
    This class is private, and should not be used directly.
    """

    def __init__(self, left, right, op):
        """
        op can be '+', '-', '*', or '/'.
        left and right are algebraic numbers.
        If op is '-', then left can be None (in which case it is taken to
        mean 0).
        if op is '/', then left can be None (in which case it is taken to
        mean 1).

        EXAMPLES:
            sage: from sage.rings.algebraic import AlgebraicNumberExpression
            sage: AlgebraicNumberExpression(Alg(1/3), Alg(1/2), '+')
            [0.83333333333333325 ... 0.83333333333333338] (1/3 + 1/2)
            sage: AlgebraicNumberExpression(Alg(1/3), Alg(1/2), '-')
            [-0.16666666666666669 ... -0.16666666666666662] (1/3 - 1/2)
            sage: AlgebraicNumberExpression(Alg(1/3), Alg(1/2), '*')
            [0.16666666666666665 ... 0.16666666666666669] (1/3 * 1/2)
            sage: AlgebraicNumberExpression(Alg(1/3), Alg(1/2), '/')
            [0.66666666666666662 ... 0.66666666666666675] (1/3 / 1/2)
            sage: AlgebraicNumberExpression(None, Alg(1/2), '-')
            [-0.50000000000000000 ... -0.50000000000000000] (None - 1/2)
            sage: AlgebraicNumberExpression(None, Alg(1/2), '/')
            [2.0000000000000000 ... 2.0000000000000000] (None / 1/2)
        """
        self._left = left
        self._right = right
        self._op = op

    def _repr_(self):
        return '%s (%s %s %s)'%(self.interval_fast(RIF), self._left, self._op, self._right)

    def interval_fast(self, field):
        """
        Returns an interval containing self with precision given by the
        field argument (which must be a RealIntervalField).  Note that
        the result may not be a minimal-width interval.

        EXAMPLES:
            sage: from sage.rings.algebraic import AlgebraicNumberExpression
            sage: five_sixths = AlgebraicNumberExpression(Alg(1/2), Alg(1/3), '+')
            sage: five_sixths.interval_fast(RealIntervalField(4))
            [0.812 ... 0.875]
            sage: five_sixths.interval_fast(RealIntervalField(70))
            [0.83333333333333333333305 ... 0.83333333333333333333390]
        """
        op = self._op
        if op == '+':
            return self._left.interval_fast(field) + self._right.interval_fast(field)
        if op == '-':
            if self._left is None:
                return -self._right.interval_fast(field)
            else:
                return self._left.interval_fast(field) - self._right.interval_fast(field)
        if op == '*':
            return self._left.interval_fast(field) * self._right.interval_fast(field)
        if op == '/':
            # Note: this may trigger exact computation
            if self._right.sign() == 0:
                raise ZeroDivisionError, 'Division by zero in algebraic number field'
            if self._left is None:
                return ~self._right.interval_fast(field)
            else:
                return self._left.interval_fast(field) / self._right.interval_fast(field)
        raise ValueError, 'Illegal operation for AlgebraicNumberExpression'

    def exactify(self):
        """
        Return a new exact AlgebraicNumberDescr with the same value as self.

        EXAMPLES:
            sage: from sage.rings.algebraic import AlgebraicNumberExpression
            sage: five_sixths = AlgebraicNumberExpression(Alg(1/2), Alg(1/3), '+')
            sage: five_sixths.exactify()
            5/6
        """
        left = self._left
        if left is None:
            if op == '-':
                left = Alg(0)
            else:
                left = Alg(1)
        left.exactify()

        right = self._right
        right.exactify()

        gen = left._exact_field().union(right._exact_field())

        # print gen
        # print left
        # print right

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

        # print gen
        # print left, left_value
        # print right, right_value

        if gen.is_unit():
            return AlgebraicNumberRational(value)
        else:
            return AlgebraicNumberExtensionElement(gen, value)


class AlgebraicNumber(sage.structure.element.FieldElement):
    """
    An algebraic real (a real number which is the zero of a polynomial
    in ZZ[x]).

    AlgebraicNumber objects can be created using Alg (== AlgebraicNumberField);
    either by coercing a rational, or by using the Alg.polynomial_root()
    method to construct a particular root of a polynomial with algebraic
    coefficients.  Also, AlgebraicNumber is closed under addition,
    subtraction, multiplication, division (except by 0), and rational
    powers (including roots), except for negative numbers and powers
    with an even denominator.

    AlgebraicNumber objects can be approximated to any desired precision.
    They can be compared exactly; if the two numbers are very close,
    this may require exact computation, which can be extremely slow.

    As long as exact computation is not triggered, computation with
    algebraic reals should not be too much slower than computation with
    intervals.  As mentioned above, exact computation is triggered
    when comparing two algebraic reals which are very close together.
    This can be an explicit comparison in user code, but the following
    list of actions (not necessarily complete) can also trigger exact
    computation:
    Dividing by an algebraic real which is very close to 0.
    Using an algebraic real which is very close to 0 as the leading coefficient
    in a polynomial.
    Taking a root of an alebraic real which is very close to 0.

    The exact definition of "very close" is subject to change; currently,
    we compute our best approximation of the two numbers using 128-bit
    arithmetic, and see if that's sufficient to decide the comparison.
    Note that comparing two algebraic reals which are actually equal will
    always trigger exact computation.

    EXAMPLES:
        sage: sqrt(Alg(2))
        [1.4142135623730949 ... 1.4142135623730952]
        sage: sqrt(Alg(2))^2 == 2
        True
        sage: x = polygen(Alg)
        sage: phi = Alg.polynomial_root(x^2 - x - 1, RIF(0, 2))
        sage: phi
        [1.6180339887498946 ... 1.6180339887498950]
        sage: phi^2 == phi+1
        True
    """

    def __init__(self, x):
        """
        Initialize an algebraic number.  The argument must be either
        a rational number or a subclass of AlgebraicNumberDescr.

        sage: from sage.rings.algebraic import AlgebraicNumberExpression
        sage: AlgebraicNumber(22/7)
        22/7
        sage: AlgebraicNumber(AlgebraicNumberExpression(Alg(1/2), Alg(1/5), '+'))
        [0.69999999999999995 ... 0.70000000000000007]
        """
        sage.structure.element.FieldElement.__init__(self, Alg)
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            self._descr = AlgebraicNumberRational(x)
        elif isinstance(x, (AlgebraicNumberDescr)):
            self._descr = x
        else:
            raise TypeError, "Illegal initializer for algebraic number"

        self._value = self._descr.interval_fast(Alg.default_interval_field())

    def _repr_(self):
        if self._descr.is_rational():
            return repr(self._descr)
        return repr(RIF(self._value))

    def _mul_(self, other):
        sd = self._descr
        od = other._descr
        if sd.is_rational() and od.is_rational():
            value = sd.rational_value() * od.rational_value()
            return AlgebraicNumber(AlgebraicNumberRational(value))
        elif sd.is_field_element() and \
                od.is_field_element() and \
                sd.field_parent() == od.field_parent():
            value = sd.field_element_value() * od.field_element_value()
            return AlgebraicNumber(AlgebraicNumberExtensionElement(sd.field_parent(), value))
        else:
            value = AlgebraicNumberExpression(self, other, '*')
            return AlgebraicNumber(value)

    def _div_(self, other):
        sd = self._descr
        od = other._descr
        if sd.is_rational() and od.is_rational():
            value = sd.rational_value() / od.rational_value()
            return AlgebraicNumber(AlgebraicNumberRational(value))
        elif sd.is_field_element() and \
                od.is_field_element() and \
                sd.field_parent() == od.field_parent():
            value = sd.field_element_value() / od.field_element_value()
            return AlgebraicNumber(AlgebraicNumberExtensionElement(sd.field_parent(), value))
        else:
            value = AlgebraicNumberExpression(self, other, '/')
            return AlgebraicNumber(value)

    def __invert__(self):
        sd = self._descr
        if sd.is_rational:
            value = ~sd.rational_value()
            return AlgebraicNumber(AlgebraicNumberRational(value))
        elif sd.is_field_element():
            value = ~sd.field_element_value()
            return AlgebraicNumber(AlgebraicNumberExtensionElement(sd.field_parent(), value))
        else:
            value = AlgebraicNumberExpression(None, self, '/')
            return AlgebraicNumber(value)

    def _add_(self, other):
        sd = self._descr
        od = other._descr
        if sd.is_rational() and od.is_rational():
            value = sd.rational_value() + od.rational_value()
            return AlgebraicNumber(AlgebraicNumberRational(value))
        elif sd.is_field_element() and \
                od.is_field_element() and \
                sd.field_parent() == od.field_parent():
            value = sd.field_element_value() + od.field_element_value()
            return AlgebraicNumber(AlgebraicNumberExtensionElement(sd.field_parent(), value))
        else:
            value = AlgebraicNumberExpression(self, other, '+')
            return AlgebraicNumber(value)

    def _sub_(self, other):
        sd = self._descr
        od = other._descr
        if sd.is_rational() and od.is_rational():
            value = sd.rational_value() - od.rational_value()
            return AlgebraicNumber(AlgebraicNumberRational(value))
        elif sd.is_field_element() and \
                od.is_field_element() and \
                sd.field_parent() == od.field_parent():
            value = sd.field_element_value() - od.field_element_value()
            return AlgebraicNumber(AlgebraicNumberExtensionElement(sd.field_parent(), value))
        else:
            value = AlgebraicNumberExpression(self, other, '-')
            return AlgebraicNumber(value)

    def _neg_(self):
        sd = self._descr
        if sd.is_rational():
            value = -sd.rational_value()
            return AlgebraicNumber(AlgebraicNumberRational(value))
        elif sd.is_field_element():
            value = -sd.field_element_value()
            return AlgebraicNumber(AlgebraicNumberExtensionElement(sd.field_parent(), value))
        else:
            value = AlgebraicNumberExpression(None, self, '-')
            return AlgebraicNumber(value)

    def __cmp__(self, other):
        if other._descr.is_rational() and other._descr.rational_value() == 0:
            return self.sign()
        elif self._descr.is_rational() and self._descr.rational_value() == 0:
            return -other.sign()
        else:
            return self._sub_(other).sign()

    def __pow__(self, e):
        """
        self^p returns the p'th power of self (where p can be an arbitrary
        rational).

        EXAMPLES:
            sage: Alg(2)^(1/2)
            [1.4142135623730949 ... 1.4142135623730952]
            sage: Alg(8)^(2/3)
            [3.9999999999999995 ... 4.0000000000000009]
            sage: Alg(8)^(2/3) == 4
            True
            sage: x = polygen(Alg)
            sage: phi = Alg.polynomial_root(x^2 - x - 1, RIF(0, 2))
            sage: tau = Alg.polynomial_root(x^2 - x - 1, RIF(-2, 0))
            sage: rt5 = Alg(5)^(1/2)
            sage: phi^10 / rt5
            [55.003636123247410 ... 55.003636123247418]
            sage: tau^10 / rt5
            [0.0036361232474132654 ... 0.0036361232474132659]
            sage: (phi^10 - tau^10) / rt5
            [54.999999999999992 ... 55.000000000000008]
            sage: (phi^10 - tau^10) / rt5 == fibonacci(10)
            True
            sage: (phi^50 - tau^50) / rt5 == fibonacci(50)
            True
        """
        e = QQ._coerce_(e)
        n = e.numerator()
        d = e.denominator()
        if d == 1:
            if n == 0:
                # implements 0^0 == 1
                return AlgebraicNumber(1)
            elif n < 0:
                return (~self).__pow__(-n)
            elif n == 1:
                return self
            else:
                pow_n2 = self.__pow__(n//2)
                if n % 2 == 1:
                    return pow_n2 * pow_n2 * self
                else:
                    return pow_n2 * pow_n2
        # Without this special case, we don't know the multiplicity
        # of the desired root
        if self.sign() == 0:
            return AlgebriacNumber(0)
        if d % 2 == 0:
            if self.sign() < 0:
                raise ValueError, 'complex algebraics not implemented'
        pow_n = self**n
        poly = AlgPoly.gen()**d - pow_n
        range = pow_n.interval_fast(RIF)
        if d % 2 == 0:
            result_min = 0
        else:
            result_min = min(range.lower(), -1)
        result_max = max(range.upper(), 1)
        return AlgebraicNumber(AlgebraicNumberRoot(poly, RIF(result_min, result_max)))

    def sqrt(self):
        return self.__pow__(~ZZ(2))

    def nth_root(self, n):
        return self.__pow__(~ZZ(n))

    def exactify(self):
        """
        Compute an exact representation for this number.

        EXAMPLES:
            sage: two = sqrt(Alg(4))
            sage: two
            [1.9999999999999997 ... 2.0000000000000005]
            sage: two.exactify()
            sage: two
            2
        """
        od = self._descr
        self._descr = self._descr.exactify()
        new_val = self._descr.interval_fast(self.parent().default_interval_field())
        self._value = self._value.intersection(new_val)
#         if self._value != self._descr.interval_fast(RIF):
#             v1 = self._value
#             v2 = self._descr.interval_fast(RIF)
#             # print (v1 != v2)
#             # print v1 - v2
#             # print v1, v2
#             # print 'try1', self._descr.interval_fast(RIF), 'done'
#             global od2
#             od2 = od
#             # print 'a', od.interval_fast(RIF), self._value
#             # print 'b', self._descr.interval_fast(RIF), v2
#             # print 'try2', self._descr.interval_fast(RIF), 'done'
#             raise ValueError

    def _exact_field(self):
        """
        Returns a generator for a number field that includes this number
        (not necessarily the smallest such number field).

        EXAMPLES:
            sage: Alg(2)._exact_field()
            Unit generator
            sage: (sqrt(Alg(2)) + sqrt(Alg(19)))._exact_field()
            Number Field in a with defining polynomial y^4 - 20*y^2 + 81 with a in [3.7893137826710354 ... 3.7893137826710360]
            sage: (Alg(7)^(3/5))._exact_field()
            Number Field in a with defining polynomial y^5 - 7 with a in [1.4757731615945519 ... 1.4757731615945522]
        """

        sd = self._descr
        if sd.is_rational():
            return unit_generator
        if sd.is_field_element():
            return sd.field_parent()
        self.exactify()
        return self._exact_field()

    def _exact_value(self):
        """
        Returns either an AlgebraicNumberRational or an
        AlgebraicNumberExtensionElement representing this value.

        EXAMPLES:
            sage: Alg(2)._exact_value()
            2
            sage: (sqrt(Alg(2)) + sqrt(Alg(19)))._exact_value()
            1/9*a^3 + a^2 - 11/9*a - 10 where a^4 - 20*a^2 + 81 = 0 and a in [3.7893137826710354 ... 3.7893137826710360]
            sage: (Alg(7)^(3/5))._exact_value()
            a^3 where a^5 - 7 = 0 and a in [1.4757731615945519 ... 1.4757731615945522]
        """
        sd = self._descr
        if sd.is_exact():
            return sd
        self.exactify()
        return self._exact_value()

    def _more_precision(self):
        """
        Recompute the interval bounding this number with higher-precision
        interval arithmetic.

        EXAMPLES:
            sage: rt2 = sqrt(Alg(2))
            sage: rt2._value
            [1.41421356237309504876 ... 1.41421356237309504888]
            sage: rt2._more_precision()
            sage: rt2._value
            [1.414213562373095048801688724209698078568 ... 1.414213562373095048801688724209698078575]
            sage: rt2._more_precision()
            sage: rt2._value
            [1.414213562373095048801688724209698078569671875376948073176679737990732478462101 ... 1.414213562373095048801688724209698078569671875376948073176679737990732478462120]
        """
        prec = self._value.prec()
        field = RealIntervalField(prec * 2)
        self._value = self._descr.interval_fast(field)

    def sign(self):
        """
        Compute the sign of this algebraic number (return -1 if negative,
        0 if zero, or 1 if positive).

        Computes an interval enclosing this number using 128-bit interval
        arithmetic; if this interval includes 0, then fall back to
        exact computation (which can be very slow).

        EXAMPLES:
            sage: Alg(-5).nth_root(7).sign()
            -1
            sage: (Alg(2).sqrt() - Alg(2).sqrt()).sign()
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
            if not self._descr.is_irrational():
                self._descr = AlgebraicNumberRational(self._descr.rational_value())
                return self.sign()
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

    def interval_fast(self, field):
        """
        Given a RealIntervalField, compute the value of this number
        using interval arithmetic of at least the precision of the field,
        and return the value in that field.  (More precision may be used
        in the computation.)  The returned interval may be arbitrarily
        imprecise, if this number is the result of a sufficiently long
        computation chain.

        EXAMPLES:
            sage: x = Alg(2).sqrt()
            sage: x.interval_fast(RIF)
            [1.4142135623730949 ... 1.4142135623730952]
            sage: x.interval_fast(RealIntervalField(200))
            [1.4142135623730950488016887242096980785696718753769480731766796 ... 1.4142135623730950488016887242096980785696718753769480731766809]
            sage: x = Alg(4).sqrt()
            sage: (x-2).interval_fast(RIF)
            [-1.0842021724855045e-19 ... 2.1684043449710089e-19]
        """
        if field.prec() == self._value.prec():
            return self._value
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
            sage: Alg(2).sqrt().interval_diameter(1e-10)
            [1.41421356237309504876 ... 1.41421356237309504888]
            sage: Alg(2).sqrt().interval_diameter(1e-30)
            [1.414213562373095048801688724209698078568 ... 1.414213562373095048801688724209698078575]
        """
        if diam <= 0:
            raise ValueError, 'diameter must be positive in interval_diameter'

        while self._value.diameter() > diam:
            self._more_precision()

        return self._value

    def interval(self, field):
        """
        Given a RealIntervalField of precision p, compute an interval
        representation of self with diameter() at most 2^-p; then round
        that representation into the given field.  Here diameter() is
        relative diameter for intervals not containing 0, and absolute
        diameter for intervals that do contain 0; thus, if the returned
        interval does not contain 0, it has at least p-1 good bits.

        EXAMPLES:
            sage: RIF64 = RealIntervalField(64)
            sage: x = Alg(2).sqrt()
            sage: y = x*x
            sage: y = 1000 * y - 999 * y
            sage: y.interval_fast(RIF64)
            [1.99999999999999966693 ... 2.00000000000000033307]
            sage: y.interval(RIF64)
            [1.99999999999999999989 ... 2.00000000000000000022]
        """
        target = RR(1.0) >> field.prec()
        val = self.interval_diameter(target)
        return field(val)

    def interval_exact(self, field):
        """
        Given a RealIntervalField, compute the best possible
        approximation of this number in that field.  Note that if this
        number is sufficiently close to some floating-point number
        (and, in particular, if this number is exactly representable in
        floating-point), then this will trigger exact computation, which
        may be very slow.

        EXAMPLES:
            sage: x = Alg(2).sqrt()
            sage: y = x*x
            sage: x.interval(RIF)
            [1.4142135623730949 ... 1.4142135623730952]
            sage: x.interval_exact(RIF)
            [1.4142135623730949 ... 1.4142135623730952]
            sage: y.interval(RIF)
            [1.9999999999999997 ... 2.0000000000000005]
            sage: y.interval_exact(RIF)
            [2.0000000000000000 ... 2.0000000000000000]
            sage: z = 1 + Alg(2).sqrt() / 2^200
            sage: z.interval(RIF)
            [1.0000000000000000 ... 1.0000000000000003]
            sage: z.interval_exact(RIF)
            [1.0000000000000000 ... 1.0000000000000003]
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

    def real(self, field):
        """
        Given a RealField, compute a good approximation to self in that field.
        The approximation will be off by at most two ulp's, except for
        numbers which are very close to 0, which will have an absolute
        error at most 2^(-(field.prec()-1)).  Also, the rounding mode of the
        field is respected.

        EXAMPLES:
            sage: x = Alg(2).sqrt()^2
            sage: x.real(RR)
            2.00000000000000
            sage: x.real(RealField(53, rnd='RNDD'))
            1.99999999999999
            sage: x.real(RealField(53, rnd='RNDU'))
            2.00000000000001
            sage: x.real(RealField(53, rnd='RNDZ'))
            1.99999999999999
            sage: (-x).real(RR)
            -2.00000000000000
            sage: (-x).real(RealField(53, rnd='RNDD'))
            -2.00000000000001
            sage: (-x).real(RealField(53, rnd='RNDU'))
            -1.99999999999999
            sage: (-x).real(RealField(53, rnd='RNDZ'))
            -1.99999999999999
            sage: (x-2).real(RR)
            5.42101086242752e-20
            sage: (x-2).real(RealField(53, rnd='RNDD'))
            -1.08420217248551e-19
            sage: (x-2).real(RealField(53, rnd='RNDU'))
            2.16840434497101e-19
            sage: (x-2).real(RealField(53, rnd='RNDZ'))
            0.000000000000000
            sage: y = Alg(2).sqrt()
            sage: y.real(RR)
            1.41421356237309
            sage: y.real(RealField(53, rnd='RNDD'))
            1.41421356237309
            sage: y.real(RealField(53, rnd='RNDU'))
            1.41421356237310
            sage: y.real(RealField(53, rnd='RNDZ'))
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

    def real_exact(self, field):
        """
        Given a RealField, compute the best possible approximation of
        this number in that field.  Note that if this number is sufficiently
        close to some floating-point number in the field (and, in particular,
        if this number is exactly representable in the field), then
        this will trigger exact computation, which may be very slow.

        The rounding mode of the field is respected.

        EXAMPLES:
            sage: x = Alg(2).sqrt()^2
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
            sage: y = Alg(2).sqrt()
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
        # val may be [x ... x], [x ... x + 1/2 ulp],
        # [x + 1/2 ulp ... x + 1/2 ulp], or [x + 1/2 ulp ... x + 1 ulp];
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

def is_AlgebraicNumber(x):
    return isinstance(x, AlgebraicNumber)

AlgPoly = PolynomialRing(Alg, 'x')

class AlgebraicPolynomialTracker(SageObject):
    """
    Keeps track of a polynomial used for algebraic numbers.

    If multiple algebraic numbers are created as roots of a single
    polynomial, this allows the polynomial and information about
    the polynomial to be shared.  This reduces work if the polynomial
    must be recomputed at higher precision, or if it must be made
    exact.

    This class is private, and should only be constructed by
    AlgebraicField.common_polynomial(), and should only be used as
    an argument to AlgebraicField.polynomial_root().

    EXAMPLES:
        sage: x = polygen(Alg)
        sage: P = Alg.common_polynomial(x^2 - x - 1)
        sage: P
        x^2 - x - 1
        sage: Alg.polynomial_root(P, RIF(0, 2))
        [1.6180339887498946 ... 1.6180339887498950]
    """

    def __init__(self, poly):
        poly = AlgPoly._coerce_(poly)
        self._poly = poly
        self._exact = False

    def _repr_(self):
        return repr(self._poly)

    def poly(self):
        return self._poly

    def exactify(self):
        """
        Compute a common field that holds all of the algebraic coefficients
        of this polynomial, then factor the polynomial over that field.
        Store the factors for later use (ignoring multiplicity).
        """
        if self._exact:
            return

        self._exact = True

        gen = unit_generator

        for c in self._poly.list():
            c.exactify()
            gen = gen.union(c._exact_field())

        self._gen = gen

        coeffs = [gen(c._exact_value()) for c in self._poly.list()]

        if gen.is_unit():
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

class AlgebraicNumberRoot(AlgebraicNumberDescr):
    """
    The subclass of AlgebraicNumberDescr that represents a particular
    root of a polynomial with algebraic coefficients.
    This class is private, and should not be used directly.
    """
    def __init__(self, poly, interval, multiplicity=1):
        if not isinstance(poly, AlgebraicPolynomialTracker):
            poly = AlgebraicPolynomialTracker(poly)
        self._poly = poly
        # Rethink precision control...integer bitcounts vs.
        # RealIntervalField values
        self._multiplicity = multiplicity
        self._interval = self.refine_interval(interval, 64)

    def _repr_(self):
        return 'Root %s of %s'%(self._interval, self._poly)

    def refine_interval(self, interval, precision):
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
            sage: from sage.rings.algebraic import AlgebraicNumberRoot
            sage: x = polygen(Alg)
            sage: rt2 = AlgebraicNumberRoot(x^2 - 2, RIF(0, 2))
            sage: rt2.refine_interval(RIF(0, 2), 75)
            [1.41421356237309504880163 ... 1.41421356237309504880175]
        """
        p = self._poly.poly()
        dp = p.derivative()
        for i in xrange(0, self._multiplicity - 1):
            p = dp
            dp = p.derivative()

        # Don't throw away bits in the original interval; doing so might
        # invalidate it (include an extra root)
        field = RealIntervalField(max(precision, interval.prec()))
        zero = field(0)
        interval = field(interval)
        poly_ring = field['x']

        # XXX Once this is properly integrated into SAGE,
        # then the following should be replaced with:
        # interval_p = poly_ring(p)
        # (and similarly for interval_dp)
        coeffs = [c.interval_fast(field) for c in p.list()]
        interval_p = poly_ring(coeffs)

        # This special case is important: this is the only way we could
        # refine "infinitely deep" (we could get an interval of diameter
        # about 2^{-2^31}, and then hit floating-point underflow); avoiding
        # this case here means we don't have to worry about iterating too
        # many times later
        if coeffs[0] == zero and zero in interval:
            return zero

        dcoeffs = [c.interval_fast(field) for c in dp.list()]
        interval_dp = poly_ring(dcoeffs)

        linfo = {}
        uinfo = {}

        # print coeffs
        # print p

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

            # print interval

            if linfo['sign'] == uinfo['sign']:
                # Oops...
                print self._poly.poly()
                print interval_p
                print linfo['endpoint'], linfo['value'], linfo['sign']
                print uinfo['endpoint'], uinfo['value'], uinfo['sign']
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


    def exactify(self):
        """
        Returns either an AlgebraicNumberRational or an
        AlgebraicNumberExtensionElement with the same value as this number.

        EXAMPLES:
            sage: from sage.rings.algebraic import AlgebraicNumberRoot
            sage: x = polygen(Alg)
            sage: two = AlgebraicNumberRoot((x-2)*(x-sqrt(Alg(2))), RIF(1.5, 3))
            sage: two.exactify()
            2 where a^2 - 2 = 0 and a in [1.4142135623730949 ... 1.4142135623730952]
            sage: two.exactify().rational_value()
            2
            sage: strange = AlgebraicNumberRoot(x^2 + sqrt(Alg(3))*x - sqrt(Alg(2)), RIF(-1, 3))
            sage: strange.exactify()
            a where a^8 - 6*a^6 + 5*a^4 - 12*a^2 + 4 = 0 and a in [0.60510122651395104 ... 0.60510122651395116]
        """
        gen = self._poly.generator()

        if gen.is_unit():
            qpf = self._poly.factors()
            def find_fn(factor, rif):
                return factor(self.interval_fast(rif))
            my_factor = find_zero_result(find_fn, qpf)

            # Factoring always returns monic polynomials over the rationals
            assert(my_factor.is_monic())

            if my_factor.degree() == 1:
                return AlgebraicNumberRational(-my_factor[0])

            den, my_factor = clear_denominators(my_factor)

            red_elt, red_back, red_pol = do_polred(my_factor)

            field = NumberField(red_pol, 'a', check=False)

            def intv_fn(rif):
                return red_elt(self.interval_fast(rif) * den)
            new_intv = isolating_interval(intv_fn, red_pol)
            root = AlgebraicNumberRoot(AlgPoly(red_pol), new_intv)
            new_gen = AlgebraicGenerator(field, root)

            return AlgebraicNumberExtensionElement(new_gen, red_back(field.gen())/den)
        else:
            fld = gen.field()

            fpf = self._poly.factors()
            # print fpf
            def find_fn(factor, rif):
                rif_poly = rif['x']
                gen_val = gen.interval_fast(rif)
                self_val = self.interval_fast(rif)
                ip = rif_poly([c.polynomial()(gen_val) for c in factor])
                return ip(self_val)
            my_factor = find_zero_result(find_fn, fpf)

            # print my_factor
            assert(my_factor.is_monic())

            if my_factor.degree() == 1:
                return AlgebraicNumberExtensionElement(gen, -my_factor[0])

            # rnfequation needs a monic polynomial with integral coefficients.
            # We achieve this with a change of variables.

            den, my_factor = clear_denominators(my_factor)

            pari_nf = gen.field().pari_nf()
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

            def intv_fn(rif):
                return red_elt(gen.interval_fast(rif) * k + self.interval_fast(rif) * den)
            new_intv = isolating_interval(intv_fn, red_pol)

            root = AlgebraicNumberRoot(QQx(red_pol), new_intv)
            new_gen = AlgebraicGenerator(new_nf, root)
            red_back_a = red_back(new_nf.gen())
            new_poly = ((QQx_x - k * self_pol_sage)(red_back_a)/den)
            return AlgebraicNumberExtensionElement(new_gen, new_poly)

    def _more_precision(self):
        """
        Recompute the interval enclosing this root at higher
        precision.
        """
        prec = self._interval.prec()
        self._interval = self.refine_interval(self._interval, prec*2)

    def interval_fast(self, field):
        """
        Given a RealIntervalField, compute the value of this number
        using interval arithmetic of at least the precision of the field,
        and return the value in that field.  (More precision may be used
        in the computation.)
        """
        if field.prec() == self._interval.prec():
            return self._interval
        if field.prec() < self._interval.prec():
            return field(self._interval)
        self._more_precision()
        return self.interval_fast(field)

unit_generator = AlgebraicGenerator(None, AlgebraicNumberRoot(AlgPoly.gen() - 1, RIF(1)))

class AlgebraicNumberExtensionElement(AlgebraicNumberDescr):
    """
    The subclass of AlgebraicNumberDescr that represents a number field
    element in terms of a specific generator.  Consists of a polynomial
    with rational coefficients in terms of the generator, and the
    generator itself, an AlgebraicGenerator.
    """

    # XXX Should override __new__, and return an AlgebraicNumberRational
    # if the value is rational.

    def __init__(self, generator, value):
        self._generator = generator
        self._value = value

    def _repr_(self):
        return '%s where %s = 0 and a in %s'%(self._value,
                                              self._generator.field().polynomial()._repr(name='a'),
                                              self._generator.interval_fast(RIF))

    def generator(self):
        return self._generator

    def is_exact(self):
        return True

    def is_field_element(self):
        return True

    def field_parent(self):
        return self._generator

    def exactify(self):
        return self

    def rational_value(self):
        poly = self._value.polynomial()
        assert(poly.is_constant())
        return poly[0]

    def is_irrational(self):
        return self._value.polynomial().degree() >= 1

    def field_element_value(self):
        return self._value

    def interval_fast(self, field):
        gen_val = self._generator.interval_fast(field)
        # XXX Coercion to field() below is necessary in case this is
        # a constant polynomial (that is, this is a rational number).
        # If we maintain the invariant that AlgebraicNumberExtensionElement
        # values are never rational, then the coercion is redundant.
        return field(self._value.polynomial()(gen_val))

ax = AlgPoly.gen()
# def heptadecagon():
#     # Compute the exact (x,y) coordinates of the vertices of a 34-gon.
#     # (Take every other coordinate to get the vertices of a
#     # heptadecagon.)
#     # Formulas from:
#     # Weisstein, Eric W. "Trigonometry Angles--Pi/17." From
#     # MathWorld--A Wolfram Web Resource.
#     # http://mathworld.wolfram.com/TrigonometryAnglesPi17.html

#     rt17 = Alg(17).sqrt()
#     rt2 = Alg(2).sqrt()
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

#     x = Alg.polynomial_root(256*ax**8 - 128*ax**7 - 448*ax**6 + 192*ax**5 + 240*ax**4 - 80*ax**3 - 40*ax**2 + 8*ax + 1, RIF(0.9829, 0.983))
#     y = (1-x**2).sqrt()

#     cx, cy = 1, 0
#     for i in range(34):
#         cx, cy = x*cx-y*cy, x*cy+y*cx
#     print cx, cy
#     print cx.sign(), cy.sign()
#     print (cx-1).sign()
#     return x, y
# # heptadecagon2()
