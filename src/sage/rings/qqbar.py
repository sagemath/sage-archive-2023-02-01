"""
Field of Algebraic Numbers

AUTHOR:

- Carl Witty (2007-01-27): initial version
- Carl Witty (2007-10-29): massive rewrite to support complex as well as real numbers

This is an implementation of the algebraic numbers (the complex
numbers which are the zero of a polynomial in `\ZZ[x]`; in other
words, the algebraic closure of `\QQ`, with an embedding into `\CC`).
All computations are exact. We also include an implementation of the
algebraic reals (the intersection of the algebraic numbers with
`\RR`). The field of algebraic numbers `\QQbar` is available with
abbreviation ``QQbar``; the field of algebraic reals has abbreviation
``AA``.

As with many other implementations of the algebraic numbers, we try
hard to avoid computing a number field and working in the number
field; instead, we use floating-point interval arithmetic whenever
possible (basically whenever we need to prove non-equalities), and
resort to symbolic computation only as needed (basically to prove
equalities).

Algebraic numbers exist in one of the following forms:

- a rational number

- the sum, difference, product, or quotient of algebraic numbers

- the negation, inverse, absolute value, norm, real part,
  imaginary part, or complex conjugate of an algebraic number

- a particular root of a polynomial, given as a polynomial with
  algebraic coefficients together with an isolating interval (given as
  a ``RealIntervalFieldElement``) which encloses exactly one root, and
  the multiplicity of the root

- a polynomial in one generator, where the generator is an algebraic
  number given as the root of an irreducible polynomial with integral
  coefficients and the polynomial is given as a
  ``NumberFieldElement``.

An algebraic number can be coerced into ``ComplexIntervalField`` (or
``RealIntervalField``, for algebraic reals); every algebraic number has a
cached interval of the highest precision yet calculated.

In most cases, computations that need to compare two algebraic numbers
compute them with 128-bit precision intervals; if this does not suffice to
prove that the numbers are different, then we fall back on exact
computation.

Note that division involves an implicit comparison of the divisor against
zero, and may thus trigger exact computation.

Also, using an algebraic number in the leading coefficient of
a polynomial also involves an implicit comparison against zero, which
again may trigger exact computation.

Note that we work fairly hard to avoid computing new number fields;
to help, we keep a lattice of already-computed number fields and
their inclusions.

EXAMPLES::

    sage: sqrt(AA(2)) > 0
    True
    sage: (sqrt(5 + 2*sqrt(QQbar(6))) - sqrt(QQbar(3)))^2 == 2
    True
    sage: AA((sqrt(5 + 2*sqrt(6)) - sqrt(3))^2) == 2
    True

For a monic cubic polynomial `x^3 + bx^2 + cx + d` with roots `s1`,
`s2`, `s3`, the discriminant is defined as
`(s1-s2)^2(s1-s3)^2(s2-s3)^2` and can be computed as `b^2c^2 - 4b^3d -
4c^3 + 18bcd - 27d^2`. We can test that these definitions do give the
same result::

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

We can convert from symbolic expressions::

    sage: QQbar(sqrt(-5))
    2.236067977499790?*I
    sage: AA(sqrt(2) + sqrt(3))
    3.146264369941973?
    sage: QQbar(I)
    I
    sage: AA(I)
    Traceback (most recent call last):
    ...
    ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
    sage: QQbar(I * golden_ratio)
    1.618033988749895?*I
    sage: AA(golden_ratio)^2 - AA(golden_ratio)
    1
    sage: QQbar((-8)^(1/3))
    1.000000000000000? + 1.732050807568878?*I
    sage: AA((-8)^(1/3))
    -2
    sage: QQbar((-4)^(1/4))
    1 + 1*I
    sage: AA((-4)^(1/4))
    Traceback (most recent call last):
    ...
    ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real

The coercion, however, goes in the other direction, since not all
symbolic expressions are algebraic numbers::

    sage: QQbar(sqrt(2)) + sqrt(3)
    sqrt(3) + 1.414213562373095?
    sage: QQbar(sqrt(2) + QQbar(sqrt(3)))
    3.146264369941973?

Note the different behavior in taking roots: for ``AA`` we prefer real
roots if they exist, but for ``QQbar`` we take the principal root::

    sage: AA(-1)^(1/3)
    -1
    sage: QQbar(-1)^(1/3)
    0.500000000000000? + 0.866025403784439?*I

We can explicitly coerce from `\QQ[I]`. (Technically, this is not quite
kosher, since `\QQ[I]` doesn't come with an embedding; we do not know
whether the field generator is supposed to map to `+I` or `-I`. We assume
that for any quadratic field with polynomial `x^2+1`, the generator maps
to `+I`.)::

    sage: K.<im> = QQ[I]
    sage: pythag = QQbar(3/5 + 4*im/5); pythag
    4/5*I + 3/5
    sage: pythag.abs() == 1
    True

However, implicit coercion from `\QQ[I]` is not allowed::

    sage: QQbar(1) + im
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '+': 'Algebraic Field' and 'Number Field in I with defining polynomial x^2 + 1'

We can implicitly coerce from algebraic reals to algebraic numbers::

    sage: a = QQbar(1); print a, a.parent()
    1 Algebraic Field
    sage: b = AA(1); print b, b.parent()
    1 Algebraic Real Field
    sage: c = a + b; print c, c.parent()
    2 Algebraic Field

Some computation with radicals::

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

The Sage rings ``AA`` and ``QQbar`` can decide equalities between radical
expressions (over the reals and complex numbers respectively)::

    sage: a = AA((2/(3*sqrt(3)) + 10/27)^(1/3) - 2/(9*(2/(3*sqrt(3)) + 10/27)^(1/3)) + 1/3)
    sage: a
    1.000000000000000?
    sage: a == 1
    True

Algebraic numbers which are known to be rational print as rationals;
otherwise they print as intervals (with 53-bit precision)::

    sage: AA(2)/3
    2/3
    sage: QQbar(5/7)
    5/7
    sage: QQbar(1/3 - 1/4*I)
    -1/4*I + 1/3
    sage: two = QQbar(4).nth_root(4)^2; two
    2.000000000000000?
    sage: two == 2; two
    True
    2
    sage: phi
    1.618033988749895?

We can find the real and imaginary parts of an algebraic number (exactly)::

    sage: r = QQbar.polynomial_root(x^5 - x - 1, CIF(RIF(0.1, 0.2), RIF(1.0, 1.1))); r
    0.1812324444698754? + 1.083954101317711?*I
    sage: r.real()
    0.1812324444698754?
    sage: r.imag()
    1.083954101317711?
    sage: r.minpoly()
    x^5 - x - 1
    sage: r.real().minpoly()
    x^10 + 3/16*x^6 + 11/32*x^5 - 1/64*x^2 + 1/128*x - 1/1024
    sage: r.imag().minpoly()  # long time (10s on sage.math, 2013)
    x^20 - 5/8*x^16 - 95/256*x^12 - 625/1024*x^10 - 5/512*x^8 - 1875/8192*x^6 + 25/4096*x^4 - 625/32768*x^2 + 2869/1048576

We can find the absolute value and norm of an algebraic number exactly.
(Note that we define the norm as the product of a number and its
complex conjugate; this is the algebraic definition of norm, if we
view ``QQbar`` as ``AA[I]``.)::

    sage: R.<x> = QQ[]
    sage: r = (x^3 + 8).roots(QQbar, multiplicities=False)[2]; r
    1.000000000000000? + 1.732050807568878?*I
    sage: r.abs() == 2
    True
    sage: r.norm() == 4
    True
    sage: (r+QQbar(I)).norm().minpoly()
    x^2 - 10*x + 13
    sage: r = AA.polynomial_root(x^2 - x - 1, RIF(-1, 0)); r
    -0.618033988749895?
    sage: r.abs().minpoly()
    x^2 + x - 1

We can compute the multiplicative order of an algebraic number::

    sage: QQbar(-1/2 + I*sqrt(3)/2).multiplicative_order()
    3
    sage: QQbar(-sqrt(3)/2 + I/2).multiplicative_order()
    12
    sage: (QQbar.zeta(23)**5).multiplicative_order()
    23

The paper "ARPREC: An Arbitrary Precision Computation Package" by
Bailey, Yozo, Li and Thompson discusses this result. Evidently it is
difficult to find, but we can easily verify it. ::

    sage: alpha = QQbar.polynomial_root(x^10 + x^9 - x^7 - x^6 - x^5 - x^4 - x^3 + x + 1, RIF(1, 1.2))
    sage: lhs = alpha^630 - 1
    sage: rhs_num = (alpha^315 - 1) * (alpha^210 - 1) * (alpha^126 - 1)^2 * (alpha^90 - 1) * (alpha^3 - 1)^3 * (alpha^2 - 1)^5 * (alpha - 1)^3
    sage: rhs_den = (alpha^35 - 1) * (alpha^15 - 1)^2 * (alpha^14 - 1)^2 * (alpha^5 - 1)^6 * alpha^68
    sage: rhs = rhs_num / rhs_den
    sage: lhs
    2.642040335819351?e44
    sage: rhs
    2.642040335819351?e44
    sage: lhs - rhs
    0.?e29
    sage: lhs == rhs
    True
    sage: lhs - rhs
    0
    sage: lhs._exact_value()
    10648699402510886229334132989629606002223831*a^9 + 23174560249100286133718183712802529035435800*a^8 + 27259790692625442252605558473646959458901265*a^7 + 21416469499004652376912957054411004410158065*a^6 + 14543082864016871805545108986578337637140321*a^5 + 6458050008796664339372667222902512216589785*a^4 - 3052219053800078449122081871454923124998263*a^3 - 14238966128623353681821644902045640915516176*a^2 - 16749022728952328254673732618939204392161001*a - 9052854758155114957837247156588012516273410 where a^10 + a^9 - a^7 - a^6 - a^5 - a^4 - a^3 + a + 1 = 0 and a in 1.176280818259918?

Given an algebraic number, we can produce a string that will reproduce
that algebraic number if you type the string into Sage. We can see
that until exact computation is triggered, an algebraic number keeps
track of the computation steps used to produce that number::

    sage: rt2 = AA(sqrt(2))
    sage: rt3 = AA(sqrt(3))
    sage: n = (rt2 + rt3)^5; n
    308.3018001722975?
    sage: sage_input(n)
    R.<x> = AA[]
    v1 = AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951))) + AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))
    v2 = v1*v1
    v2*v2*v1

But once exact computation is triggered, the computation tree is discarded,
and we get a way to produce the number directly::

    sage: n == 109*rt2 + 89*rt3
    True
    sage: sage_input(n)
    R.<x> = AA[]
    v = AA.polynomial_root(AA.common_polynomial(x^4 - 4*x^2 + 1), RIF(RR(0.51763809020504148), RR(0.51763809020504159)))
    -109*v^3 - 89*v^2 + 327*v + 178

We can also see that some computations (basically, those which are
easy to perform exactly) are performed directly, instead of storing
the computation tree::

    sage: z3_3 = QQbar.zeta(3) * 3
    sage: z4_4 = QQbar.zeta(4) * 4
    sage: z5_5 = QQbar.zeta(5) * 5
    sage: sage_input(z3_3 * z4_4 * z5_5)
    R.<x> = AA[]
    3*QQbar.polynomial_root(AA.common_polynomial(x^2 + x + 1), CIF(RIF(-RR(0.50000000000000011), -RR(0.49999999999999994)), RIF(RR(0.8660254037844386), RR(0.86602540378443871))))*QQbar(4*I)*(5*QQbar.polynomial_root(AA.common_polynomial(x^4 + x^3 + x^2 + x + 1), CIF(RIF(RR(0.3090169943749474), RR(0.30901699437494745)), RIF(RR(0.95105651629515353), RR(0.95105651629515364)))))

Note that the ``verify=True`` argument to ``sage_input`` will always trigger
exact computation, so running ``sage_input`` twice in a row on the same number
will actually give different answers. In the following, running ``sage_input``
on ``n`` will also trigger exact computation on ``rt2``, as you can see by the
fact that the third output is different than the first::

    sage: rt2 = AA(sqrt(2))
    sage: n = rt2^2
    sage: sage_input(n, verify=True)
    # Verified
    R.<x> = AA[]
    v = AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
    v*v
    sage: sage_input(n, verify=True)
    # Verified
    AA(2)
    sage: n = rt2^2
    sage: sage_input(n, verify=True)
    # Verified
    AA(2)

Just for fun, let's try ``sage_input`` on a very complicated expression. The
output of this example changed with the rewriting of polynomial multiplication
algorithms in :trac:`10255`::

    sage: rt2 = sqrt(AA(2))
    sage: rt3 = sqrt(QQbar(3))
    sage: x = polygen(QQbar)
    sage: nrt3 = AA.polynomial_root((x-rt2)*(x+rt3), RIF(-2, -1))
    sage: one = AA.polynomial_root((x-rt2)*(x-rt3)*(x-nrt3)*(x-1-rt3-nrt3), RIF(0.9, 1.1))
    sage: one
    1.000000000000000?
    sage: sage_input(one, verify=True)
    # Verified
    R1 = QQbar['x']
    x1 = R1.gen()
    R2 = AA['x']
    x2 = R2.gen()
    cp1 = AA.common_polynomial(x2^2 - 2)
    v1 = QQbar.polynomial_root(cp1, RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
    v2 = QQbar.polynomial_root(AA.common_polynomial(x1^2 - 3), CIF(RIF(RR(1.7320508075688772), RR(1.7320508075688774)), RIF(RR(0))))
    v3 = -v1 - v2
    v4 = QQbar.polynomial_root(cp1, RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
    cp2 = AA.common_polynomial(x1^2 + (-v4 + v2)*x1 - v4*v2)
    v5 = QQbar.polynomial_root(cp2, RIF(-RR(1.7320508075688774), -RR(1.7320508075688772)))
    v6 = v3 - v5
    v7 = -1 - v2 - QQbar.polynomial_root(cp2, RIF(-RR(1.7320508075688774), -RR(1.7320508075688772)))
    v8 = v1*v2
    v9 = v8 - v3*v5
    si = v8*v5
    AA.polynomial_root(AA.common_polynomial(x1^4 + (v6 + v7)*x1^3 + (v9 + v6*v7)*x1^2 + (-si + v9*v7)*x1 - si*v7), RIF(RR(0.99999999999999989), RR(1.0000000000000002)))
    sage: one
    1

We can pickle and unpickle algebraic fields (and they are globally unique)::

    sage: loads(dumps(AlgebraicField())) is AlgebraicField()
    True
    sage: loads(dumps(AlgebraicRealField())) is AlgebraicRealField()
    True

We can pickle and unpickle algebraic numbers::

    sage: loads(dumps(QQbar(10))) == QQbar(10)
    True
    sage: loads(dumps(QQbar(5/2))) == QQbar(5/2)
    True
    sage: loads(dumps(QQbar.zeta(5))) == QQbar.zeta(5)
    True

    sage: t = QQbar(sqrt(2)); type(t._descr)
    <class 'sage.rings.qqbar.ANRoot'>
    sage: loads(dumps(t)) == QQbar(sqrt(2))
    True

    sage: t.exactify(); type(t._descr)
    <class 'sage.rings.qqbar.ANExtensionElement'>
    sage: loads(dumps(t)) == QQbar(sqrt(2))
    True

    sage: t = ~QQbar(sqrt(2)); type(t._descr)
    <class 'sage.rings.qqbar.ANUnaryExpr'>
    sage: loads(dumps(t)) == 1/QQbar(sqrt(2))
    True

    sage: t = QQbar(sqrt(2)) + QQbar(sqrt(3)); type(t._descr)
    <class 'sage.rings.qqbar.ANBinaryExpr'>
    sage: loads(dumps(t)) == QQbar(sqrt(2)) + QQbar(sqrt(3))
    True

We can convert elements of ``QQbar`` and ``AA`` into the following
types: ``float``, ``complex``, ``RDF``, ``CDF``, ``RR``, ``CC``,
``RIF``, ``CIF``, ``ZZ``, and ``QQ``, with a few exceptions. (For the
arbitrary-precision types, ``RR``, ``CC``, ``RIF``, and ``CIF``, it
can convert into a field of arbitrary precision.)

Converting from ``QQbar`` to a real type (``float``, ``RDF``, ``RR``,
``RIF``, ``ZZ``, or ``QQ``) succeeds only if the ``QQbar`` is actually
real (has an imaginary component of exactly zero). Converting from
either ``AA`` or ``QQbar`` to ``ZZ`` or ``QQ`` succeeds only if the
number actually is an integer or rational. If conversion fails, a
ValueError will be raised.

Here are examples of all of these conversions::

    sage: all_vals = [AA(42), AA(22/7), AA(golden_ratio), QQbar(-13), QQbar(89/55), QQbar(-sqrt(7)), QQbar.zeta(5)]
    sage: def convert_test_all(ty):
    ....:     def convert_test(v):
    ....:         try:
    ....:             return ty(v)
    ....:         except (TypeError, ValueError):
    ....:             return None
    ....:     return [convert_test(_) for _ in all_vals]
    sage: convert_test_all(float)
    [42.0, 3.1428571428571432, 1.618033988749895, -13.0, 1.6181818181818182, -2.6457513110645907, None]
    sage: convert_test_all(complex)
    [(42+0j), (3.1428571428571432+0j), (1.618033988749895+0j), (-13+0j), (1.6181818181818182+0j), (-2.6457513110645907+0j), (0.30901699437494745+0.9510565162951536j)]
    sage: convert_test_all(RDF)
    [42.0, 3.1428571428571432, 1.618033988749895, -13.0, 1.6181818181818182, -2.6457513110645907, None]
    sage: convert_test_all(CDF)
    [42.0, 3.1428571428571432, 1.618033988749895, -13.0, 1.6181818181818182, -2.6457513110645907, 0.30901699437494745 + 0.9510565162951536*I]
    sage: convert_test_all(RR)
    [42.0000000000000, 3.14285714285714, 1.61803398874989, -13.0000000000000, 1.61818181818182, -2.64575131106459, None]
    sage: convert_test_all(CC)
    [42.0000000000000, 3.14285714285714, 1.61803398874989, -13.0000000000000, 1.61818181818182, -2.64575131106459, 0.309016994374947 + 0.951056516295154*I]
    sage: convert_test_all(RIF)
    [42, 3.142857142857143?, 1.618033988749895?, -13, 1.618181818181819?, -2.645751311064591?, None]
    sage: convert_test_all(CIF)
    [42, 3.142857142857143?, 1.618033988749895?, -13, 1.618181818181819?, -2.645751311064591?, 0.3090169943749474? + 0.9510565162951536?*I]
    sage: convert_test_all(ZZ)
    [42, None, None, -13, None, None, None]
    sage: convert_test_all(QQ)
    [42, 22/7, None, -13, 89/55, None, None]

Compute the exact coordinates of a 34-gon (the formulas used are from
Weisstein, Eric W. "Trigonometry Angles--Pi/17." and can be found at
http://mathworld.wolfram.com/TrigonometryAnglesPi17.html)::

    sage: rt17 = AA(17).sqrt()
    sage: rt2 = AA(2).sqrt()
    sage: eps = (17 + rt17).sqrt()
    sage: epss = (17 - rt17).sqrt()
    sage: delta = rt17 - 1
    sage: alpha = (34 + 6*rt17 + rt2*delta*epss - 8*rt2*eps).sqrt()
    sage: beta = 2*(17 + 3*rt17 - 2*rt2*eps - rt2*epss).sqrt()
    sage: x = rt2*(15 + rt17 + rt2*(alpha + epss)).sqrt()/8
    sage: y = rt2*(epss**2 - rt2*(alpha + epss)).sqrt()/8

    sage: cx, cy = 1, 0
    sage: for i in range(34):
    ....:    cx, cy = x*cx-y*cy, x*cy+y*cx
    sage: cx
    1.000000000000000?
    sage: cy
    0.?e-15

    sage: ax = polygen(AA)
    sage: x2 = AA.polynomial_root(256*ax**8 - 128*ax**7 - 448*ax**6 + 192*ax**5 + 240*ax**4 - 80*ax**3 - 40*ax**2 + 8*ax + 1, RIF(0.9829, 0.983))
    sage: y2 = (1-x2**2).sqrt()
    sage: x - x2
    0.?e-18
    sage: y - y2
    0.?e-17

Ideally, in the above example we should be able to test ``x == x2`` and ``y ==
y2`` but this is currently infinitely long.

TESTS:

Verify that :trac:`10981` is fixed::

    sage: x = AA['x'].gen()
    sage: P = 1/(1+x^4)
    sage: P.partial_fraction_decomposition()
    (0, [(-0.3535533905932738?*x + 1/2)/(x^2 - 1.414213562373095?*x + 1), (0.3535533905932738?*x + 1/2)/(x^2 + 1.414213562373095?*x + 1)])
"""
import itertools
import operator

import sage.rings.ring
from sage.misc.fast_methods import Singleton
from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject
from sage.rings.real_mpfr import RR
from sage.rings.real_mpfi import RealIntervalField, RIF, is_RealIntervalFieldElement, RealIntervalField_class
from sage.rings.complex_field import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField, is_ComplexIntervalField
from sage.rings.complex_interval import is_ComplexIntervalFieldElement
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.number_field.number_field import NumberField, QuadraticField, CyclotomicField
from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_quadratic
from sage.arith.all import factor
from sage.structure.element import generic_power, canonical_coercion
import infinity
from sage.misc.functional import cyclotomic_polynomial

CC = ComplexField()
CIF = ComplexIntervalField()

is_SymbolicExpressionRing = None

class AlgebraicField_common(sage.rings.ring.Field):
    r"""
    Common base class for the classes :class:`~AlgebraicRealField` and
    :class:`~AlgebraicField`.
    """

    def default_interval_prec(self):
        r"""
        Return the default interval precision used for root isolation.

        EXAMPLES::

            sage: AA.default_interval_prec()
            64
        """

        return 64

    def is_finite(self):
        r"""
        Check whether this field is finite. Since this class is only used for
        fields of characteristic 0, always returns False.

        EXAMPLE::

            sage: QQbar.is_finite()
            False
        """
        return False

    def characteristic(self):
        r"""
        Return the characteristic of this field. Since this class is only used
        for fields of characteristic 0, always returns 0.

        EXAMPLES::

            sage: AA.characteristic()
            0
        """
        return sage.rings.integer.Integer(0)

    def order(self):
        r"""
        Return the cardinality of self. Since this class is only used for
        fields of characteristic 0, always returns Infinity.

        EXAMPLE::

            sage: QQbar.order()
            +Infinity
        """
        return infinity.infinity

    def common_polynomial(self, poly):
        """
        Given a polynomial with algebraic coefficients, returns a
        wrapper that caches high-precision calculations and
        factorizations. This wrapper can be passed to polynomial_root
        in place of the polynomial.

        Using ``common_polynomial`` makes no semantic difference, but will
        improve efficiency if you are dealing with multiple roots
        of a single polynomial.

        EXAMPLES::

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
            x^2 + (-sqrt(3) - sqrt(-5))*x + sqrt(3)*sqrt(-5)
            sage: p = QQbar.common_polynomial(p)
            sage: a = QQbar.polynomial_root(p, CIF(RIF(-0.1, 0.1), RIF(2, 3))); a
            0.?e-18 + 2.236067977499790?*I
            sage: b = QQbar.polynomial_root(p, RIF(1, 2)); b
            1.732050807568878?

        These "common polynomials" can be shared between real and
        complex roots::

             sage: p = AA.common_polynomial(x^3 - x - 1)
             sage: r1 = AA.polynomial_root(p, RIF(1.3, 1.4)); r1
             1.324717957244746?
             sage: r2 = QQbar.polynomial_root(p, CIF(RIF(-0.7, -0.6), RIF(0.5, 0.6))); r2
             -0.6623589786223730? + 0.5622795120623013?*I
        """
        return AlgebraicPolynomialTracker(poly)

class AlgebraicRealField(Singleton, AlgebraicField_common):
    r"""
    The field of algebraic reals.

    TESTS::

        sage: AA == loads(dumps(AA))
        True
    """

    def __new__(cls):
        r"""
        This method is there to ensure that pickles created before this class
        was made a :class:`~sage.misc.fast_methods.Singleton` still load.

        TESTS::

            sage: s = loads('x\x9cmQ\xcbR\x141\x14\xad\x11A\x083\xe2\x03T|'
            ....: '\x82l`\xd3\xff\xe0\x86\x8de/\xba*\xcb\xa9[\xe9\xf4'
            ....: '\xa5;e:=\'I+,\xa6J\x17B\xf9\xd7f\x08\xe2s\x95\xa4\xee9\xf7<'
            ....: '\xf2\xe5\x8e\x0e\xaa\xe5"D?\xea8z.\x9a\x0b\xa7z\xa3I[\x15'
            ....: '\x82\xf8\xf3\x85\xc9\xb1<xg[\xae\xbd2\xbabeO\r\xdb\x86>\x9b'
            ....: '\xd8\x91V\x91\xdb\xc1_\xe0f\xa57\xae\r\x05P+/\xfe\xe5\x08'
            ....: '\xaci\xa2z46\x1aG$Z\x8e*F/p\xf7oC\xa33\x18\x99</<\x07v\tf'
            ....: '\x06\'F\xe7\xb9\x195\x0b\xacg\xc2\x8d\xbc\xe1P\x9c\xad\x04'
            ....: '\x828\xcd\x076N\x96W\xb8WaSN\x17\xca\xa7\r9\r\xb6.+\x88Kl'
            ....: '\x97e\xb7\x16+LO\xbeb\xb6\xc4\xfdc)\x88\xfb\x9a\x9b&\x05'
            ....: '\xc0N)wI\x0f\xee\x13\xfbH=\xc7nh(U\xc2xP\xca\r\xd2\x8d'
            ....: '\x8a\n\x0fK\xb9\xf5+\xfe\xa3n3MV\x98\x80\xc7rr\xfe\r\xbbr'
            ....: '\x9bZv\xecU\x1c|\xc0\xde\x12O\xe4:\xd5*0\x9ev3\xb9C\x0b'
            ....: '\xa3?Z\xa6\xa4\x11R6<{?I\xa2l\xb9\xbf6;\xb8\\\xc6\xe0\xb1'
            ....: '\x9f\xb3\xf6&\xe8\xe2,\xb3R\x13\xf9\xf2\xe1\xda\x9c\xc0s'
            ....: '\xb9\xf7?.\xe1E7\xeb\xa6W\x15^&\x80q&\x1aeo\x93Y\x13"^\xcd'
            ....: '\xf1Z\xee\xdf\x92W\x18Z\xa4\xa6(\xd7\x867\xdf\x93\xad\x9fL'
            ....: '\xa5W\xff\x90\x89\x07s\x1c\xfe6\xd2\x03{\xcdy\xf4v\x8e\xa3'
            ....: '\xb1.~\x000\xc2\xe0\xa1')
            sage: s is AA
            True

        """
        try: return AA
        except BaseException: return AlgebraicField_common.__new__(cls)

    def __init__(self):
        r"""
        Standard initialization function.

        EXAMPLE:

        This function calls functions in superclasses which set the category, so we check that.

            sage: QQbar.category() # indirect doctest
            Category of fields
        """
        AlgebraicField_common.__init__(self, self, ('x',), normalize=False)

    def _element_constructor_(self, x):
        r"""
        Construct an element of the field of algebraic real numbers from ``x``.

        EXAMPLE::

            sage: QQbar(sqrt(2)) in AA # indirect doctest
            True
            sage: QQbar(I) in AA
            False
            sage: AA in AA
            False

        The following should both return True (this is a bug). ::

            sage: sqrt(2) in AA # not tested
            False
            sage: K.<z> = CyclotomicField(5); z + 1/z in AA # not tested
            False
        """
        if isinstance(x, AlgebraicReal):
            return x
        elif isinstance(x, AlgebraicNumber):
            if x.imag().is_zero():
                return x.real()
            else:
                raise ValueError("Cannot coerce algebraic number with non-zero imaginary part to algebraic real")
        elif hasattr(x, '_algebraic_'):
            return x._algebraic_(AA)
        return AlgebraicReal(x)

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: AA._repr_()
            'Algebraic Real Field'
        """
        return "Algebraic Real Field"

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: AA._repr_option('element_is_atomic')
            True
        """
        if key == 'element_is_atomic':
            return True
        return super(AlgebraicRealField, self)._repr_option(key)

    # Is there a standard representation for this?
    def _latex_(self):
        r"""
        Latex representation of self.

        EXAMPLE::

            sage: AA._latex_()
            '\\mathbf{A}'
        """
        return "\\mathbf{A}"

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(AA, verify=True)
            # Verified
            AA
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: AA._sage_input_(SageInputBuilder(), False)
            {atomic:AA}
        """
        return sib.name('AA')

    def _coerce_map_from_(self, from_par):
        r"""
        Set up the coercion model.

        TESTS::

            sage: AA.has_coerce_map_from(ZZ) # indirect doctest
            True
            sage: K.<a> = QuadraticField(7, embedding=AA(7).sqrt()); AA.has_coerce_map_from(K)
            True
            sage: a in AA
            True
            sage: a + AA(3)
            5.645751311064590?
            sage: AA.has_coerce_map_from(SR)
            False
        """
        return (from_par is ZZ or from_par is QQ
                or from_par is AA)

    def completion(self, p, prec, extras = {}):
        r"""
        Return the completion of self at the place `p`. Only implemented for `p
        = \infty` at present.

        INPUT:

        - ``p`` -- either a prime (not implemented at present) or Infinity
        - ``prec`` -- precision of approximate field to return
        - ``extras`` -- a dict of extra keyword arguments for the ``RealField``
          constructor

        EXAMPLES::

            sage: AA.completion(infinity, 500)
            Real Field with 500 bits of precision
            sage: AA.completion(infinity, prec=53, extras={'type':'RDF'})
            Real Double Field
            sage: AA.completion(infinity, 53) is RR
            True
            sage: AA.completion(7, 10)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if p == infinity.Infinity:
            from sage.rings.real_mpfr import create_RealField
            return create_RealField(prec, **extras)
        else:
            raise NotImplementedError

    def algebraic_closure(self):
        """
        Return the algebraic closure of this field, which is the field
        `\overline{\QQ}` of algebraic numbers.

        EXAMPLES::

            sage: AA.algebraic_closure()
            Algebraic Field
        """
        return QQbar

    def _is_valid_homomorphism_(self, codomain, im_gens):
        r"""
        Attempt to construct a homomorphism from self to codomain sending the
        generators to ``im_gens``. Since this field is not finitely generated,
        this cannot be implemented in a mathematically sensible way, and we
        just test that there exists a canonical coercion.

        EXAMPLE::

            sage: AA._is_valid_homomorphism_(QQbar, [QQbar(1)])
            True
            sage: AA._is_valid_homomorphism_(QQ, [QQ(1)])
            False
        """
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

    def gens(self):
        r"""
        Return a set of generators for this field. As this field is not
        finitely generated, we opt for just returning 1.

        EXAMPLE::

            sage: AA.gens()
            (1,)
        """
        return (self(1), )

    def gen(self, n=0):
        r"""
        Return the `n`-th element of the tuple returned by :meth:`gens`.

        EXAMPLE::

            sage: AA.gen(0)
            1
            sage: AA.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: n must be 0
        """
        if n == 0:
            return self(1)
        else:
            raise IndexError("n must be 0")

    def ngens(self):
        r"""
        Return the size of the tuple returned by :meth:`gens`.

        EXAMPLE::

            sage: AA.ngens()
            1
        """
        return 1

    def zeta(self, n=2):
        r"""
        Return an `n`-th root of unity in this field. This will raise a
        ``ValueError`` if `n \ne \{1, 2\}` since no such root exists.

        INPUT:

        - ``n`` (integer) -- default 2

        EXAMPLE::

            sage: AA.zeta(1)
            1
            sage: AA.zeta(2)
            -1
            sage: AA.zeta()
            -1
            sage: AA.zeta(3)
            Traceback (most recent call last):
            ...
            ValueError: no n-th root of unity in algebraic reals

        Some silly inputs::

            sage: AA.zeta(Mod(-5, 7))
            -1
            sage: AA.zeta(0)
            Traceback (most recent call last):
            ...
            ValueError: no n-th root of unity in algebraic reals
        """
        if n == 1:
            return self(1)
        elif n == 2:
            return self(-1)
        else:
            raise ValueError("no n-th root of unity in algebraic reals")

    def polynomial_root(self, poly, interval, multiplicity=1):
        r"""
        Given a polynomial with algebraic coefficients and an interval
        enclosing exactly one root of the polynomial, constructs
        an algebraic real representation of that root.

        The polynomial need not be irreducible, or even squarefree; but
        if the given root is a multiple root, its multiplicity must be
        specified. (IMPORTANT NOTE: Currently, multiplicity-`k` roots
        are handled by taking the `(k-1)`-st derivative of the polynomial.
        This means that the interval must enclose exactly one root
        of this derivative.)

        The conditions on the arguments (that the interval encloses exactly
        one root, and that multiple roots match the given multiplicity)
        are not checked; if they are not satisfied, an error may be
        thrown (possibly later, when the algebraic number is used),
        or wrong answers may result.

        Note that if you are constructing multiple roots of a single
        polynomial, it is better to use ``AA.common_polynomial`` (or
        ``QQbar.common_polynomial``; the two are equivalent) to get a
        shared polynomial.

        EXAMPLES::

            sage: x = polygen(AA)
            sage: phi = AA.polynomial_root(x^2 - x - 1, RIF(1, 2)); phi
            1.618033988749895?
            sage: p = (x-1)^7 * (x-2)
            sage: r = AA.polynomial_root(p, RIF(9/10, 11/10), multiplicity=7)
            sage: r; r == 1
            1.000000000000000?
            True
            sage: p = (x-phi)*(x-sqrt(AA(2)))
            sage: r = AA.polynomial_root(p, RIF(1, 3/2))
            sage: r; r == sqrt(AA(2))
            1.414213562373095?
            True

        We allow complex polynomials, as long as the particular root
        in question is real. ::

            sage: K.<im> = QQ[I]
            sage: x = polygen(K)
            sage: p = (im + 1) * (x^3 - 2); p
            (I + 1)*x^3 - 2*I - 2
            sage: r = AA.polynomial_root(p, RIF(1, 2)); r^3
            2.000000000000000?
        """
        if not is_RealIntervalFieldElement(interval):
            raise ValueError("interval argument of .polynomial_root on algebraic real field must be real")

        return AlgebraicReal(ANRoot(poly, interval, multiplicity))

    def _factor_univariate_polynomial(self, f):
        """
        Factor the univariate polynomial ``f``.

        INPUT:

        - ``f`` -- a univariate polynomial defined over the real algebraic field

        OUTPUT:

        - A factorization of ``f`` over the real algebraic numbers into a unit
          and monic irreducible factors

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.factor`.

        TESTS::

            sage: R.<x> = AA[]
            sage: AA._factor_univariate_polynomial(x)
            x
            sage: AA._factor_univariate_polynomial(2*x)
            (2) * x
            sage: AA._factor_univariate_polynomial((x^2 + 1)^2)
            (x^2 + 1)^2
            sage: AA._factor_univariate_polynomial(x^8 + 1)
            (x^2 - 1.847759065022574?*x + 1.000000000000000?) * (x^2 - 0.7653668647301795?*x + 1.000000000000000?) * (x^2 + 0.7653668647301795?*x + 1.000000000000000?) * (x^2 + 1.847759065022574?*x + 1.000000000000000?)
            sage: AA._factor_univariate_polynomial(R(3))
            3
            sage: AA._factor_univariate_polynomial(12*x^2 - 4)
            (12) * (x - 0.5773502691896258?) * (x + 0.5773502691896258?)
            sage: AA._factor_univariate_polynomial(12*x^2 + 4)
            (12) * (x^2 + 0.3333333333333334?)
            sage: AA._factor_univariate_polynomial(EllipticCurve('11a1').change_ring(AA).division_polynomial(5))
            (5) * (x - 16.00000000000000?) * (x - 5.000000000000000?) * (x - 1.959674775249769?) * (x + 2.959674775249769?) * (x^2 - 2.854101966249685?*x + 15.47213595499958?) * (x^2 + 1.909830056250526?*x + 1.660606461254312?) * (x^2 + 3.854101966249685?*x + 6.527864045000421?) * (x^2 + 13.09016994374948?*x + 93.33939353874569?)

        """
        rr = f.roots()
        cr = [(r,e) for r,e in f.roots(QQbar) if r.imag()>0]

        from sage.structure.factorization import Factorization
        return Factorization(
            [(f.parent()([-r,1]),e) for r,e in rr] +
            [(f.parent()([r.norm(),-2*r.real(),1]),e) for r,e in cr],
            unit=f.leading_coefficient())

def is_AlgebraicRealField(F):
    r"""
    Check whether ``F`` is an :class:`~AlgebraicRealField` instance. For internal use.

    EXAMPLE::

        sage: from sage.rings.qqbar import is_AlgebraicRealField
        sage: [is_AlgebraicRealField(x) for x in [AA, QQbar, None, 0, "spam"]]
        [True, False, False, False, False]
    """
    return isinstance(F, AlgebraicRealField)

# Create the globally unique AlgebraicRealField object.
AA = AlgebraicRealField()

class AlgebraicField(Singleton, AlgebraicField_common):
    """
    The field of all algebraic complex numbers.
    """

    def __new__(cls):
        r"""
        This method is there to ensure that pickles created before this class
        was made a :class:`~sage.misc.fast_methods.Singleton` still load.

        TESTS::

            sage: s = loads('x\x9c}RMo\x131\x10U(-\xad\x9b\x92\x16ZJh\x80~'
            ....: '\x00MZX~\x03\x97J\x08\xb1\x87H>F\x96\xd7;\xdd\xb1\xd8x3\xb6'
            ....: '\x17\xe8!\x12\x1c\xda\xaa\xff\x9aI\xb7\x04\x8a*N\xb65\xef'
            ....: '\xcd\xbc\xf7\xc6?\xee\x99\xa0\x0bHB\xf4\xb5\x89\xb5'
            ....: '\x87$?szl\x8d2\xa5\x0eA\xdc~Q\xab/{\x1f\xca\x022\xaf\xad9'
            ....: '\xb1P\xe6\xea\x9b\x8d\xa8\x8c\x8ePT\xfe\x8cn\xday\xeb\x8a'
            ....: '\x90\x10e\xda\x8b\xdbxA\x0bF\xa9\xac\xb6e\xb4N)Q@\xd41zA'
            ....: '\xf7\xff\x15R;K5(\x0f\x13\x0f\x01\x1c\xc3l\xe5D\xed<\xe4'
            ....: '\xb5\x01A\x8b\r\xe1f\xb4\x85\x90\x9c\xce\x06\x04q\xd2\x1c'
            ....: '\xb44\x98^\xd2\x83!-\xcb\xf6D{\xee\xd0\xb8\xa0\x95\x8b!\x89'
            ....: '\x0bZMS\\\x88Cj\x0f~\xd2\xda\x94\x1e\xf6\xa5P0\xce \xcfY<uR'
            ....: '\xb9\xa9L\xe5\xbe\x82\x8fj\x0c\x11\xab\\q\x14@\xeb\xa9\\R&'
            ....: '\xd7Q\xd3F*W\xfeX\x7f\x84\xcb\\\x99a\x02=\x96\xad\x8f\xe7'
            ....: '\xb4)WU\x01\x0e\xbc\x8e\x95\x0f\xb45\xa5\'rQe:\x00m#G\xb9;'
            ....: '\x8ff\x08\xba\xbc+\xce\xa7\xff\x89s\xce\x11\xd4E\xf6\xf3'
            ....: '\x8c\xfdt\xd9\xcf\x0e\xfb\xe9M\xe9y\x1f;)\xae\xa7\xb8'
            ....: '\x91"KC\x96\xf4\xfd\x9c^ \xabx\x89\xdb\xd8\x93\x1d5\xb1'
            ....: '\xe6K\t\x8a-\x06\x8e\x96v?\xb5\xd83\x940\xbe\xce\xaar'
            ....: '\xcd.*O{\x8d\x8c\xb1\r&9mX\xbc\x88\xe6\xf2\xf9:\x1bA\xfbr'
            ....: '\xeb.\xae\xa2\x03\xec\xe1\xce\xe5\x90^1\xc0:\x1b\xad.\xe7'
            ....: '\xc1\x966Dz=\xa27\xb2;\'\xcf0j\xc2\x8bR\xcd\xd6\xe8\xf0'
            ....: '\x8ae\xfdfj3\xfb\x06\r\xb1?\xa2\xc1_%S\x817\xd0\x94'
            ....: '\x8eFt\\g\xc8\x96p\x0f\xf7\xf1\x00\xd7\xb0\xcd\x1a\xde"'
            ....: '\x0f{\x87\x87W\xc8\xdc\x04\x19\xf5\xbe\xce\x92_p\'\x13\xc5')
            sage: s is QQbar
            True
        """
        try: return QQbar
        except BaseException: return AlgebraicField_common.__new__(cls)

    def __init__(self):
        r"""
        Standard init function.

        We test by setting the category::

            sage: QQbar.category() # indirect doctest
            Category of fields
            sage: QQbar.base_ring()
            Algebraic Real Field

        TESTS::

            sage: QQbar._repr_option('element_is_atomic')
            False
        """
        AlgebraicField_common.__init__(self, AA, ('I',), normalize=False)

    def _element_constructor_(self, x):
        """
        Try to construct an element of the field of algebraic numbers from `x`.

        EXAMPLES::

            sage: sqrt(2) in QQbar # indirect doctest
            True
            sage: 22/7 in QQbar
            True
            sage: pi in QQbar
            False
        """
        if isinstance(x, AlgebraicNumber):
            return x
        elif isinstance(x, AlgebraicReal):
            return AlgebraicNumber(x._descr)
        elif hasattr(x, '_algebraic_'):
            return x._algebraic_(QQbar)
        return AlgebraicNumber(x)

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: QQbar._repr_()
            'Algebraic Field'
        """
        return "Algebraic Field"

    def _latex_(self):
        r"""
        Latex representation of self.

        EXAMPLE::

            sage: QQbar._latex_()
            '\\overline{\\QQ}'
        """
        return "\\overline{\\QQ}"

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(QQbar, verify=True)
            # Verified
            QQbar
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: QQbar._sage_input_(SageInputBuilder(), False)
            {atomic:QQbar}
        """
        return sib.name('QQbar')

    def _coerce_map_from_(self, from_par):
        r"""
        Set up the coercion model.

        TESTS::

            sage: QQbar.has_coerce_map_from(ZZ) # indirect doctest
            True
            sage: QQbar.has_coerce_map_from(AA)
            True
            sage: QQbar.has_coerce_map_from(CC)
            False
            sage: QQbar.has_coerce_map_from(SR)
            False
        """
        return (from_par is ZZ or from_par is QQ
                or from_par is AA or from_par is QQbar)

    def completion(self, p, prec, extras = {}):
        r"""
        Return the completion of self at the place `p`. Only implemented for `p
        = \infty` at present.

        INPUT:

        - ``p`` -- either a prime (not implemented at present) or Infinity
        - ``prec`` -- precision of approximate field to return
        - ``extras`` -- a dict of extra keyword arguments for the ``RealField``
          constructor

        EXAMPLES::

            sage: QQbar.completion(infinity, 500)
            Complex Field with 500 bits of precision
            sage: QQbar.completion(infinity, prec=53, extras={'type':'RDF'})
            Complex Double Field
            sage: QQbar.completion(infinity, 53) is CC
            True
            sage: QQbar.completion(3, 20)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if p == infinity.Infinity:
            from sage.rings.real_mpfr import create_RealField
            return create_RealField(prec, **extras).complex_field()
        else:
            raise NotImplementedError

    def algebraic_closure(self):
        """
        Return the algebraic closure of this field. As this field is already
        algebraically closed, just returns self.

        EXAMPLES::

            sage: QQbar.algebraic_closure()
            Algebraic Field
        """
        return self

    def construction(self):
        """
        Return a functor that constructs self (used by the coercion machinery).

        EXAMPLE::

            sage: QQbar.construction()
            (AlgebraicClosureFunctor, Rational Field)
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        from sage.all import QQ
        return (AlgebraicClosureFunctor(), QQ)

    def gens(self):
        r"""
        Return a set of generators for this field. As this field is not
        finitely generated over its prime field, we opt for just returning I.

        EXAMPLE::

            sage: QQbar.gens()
            (I,)
        """
        return(QQbar_I, )

    def gen(self, n=0):
        r"""
        Return the `n`-th element of the tuple returned by :meth:`gens`.

        EXAMPLE::

            sage: QQbar.gen(0)
            I
            sage: QQbar.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: n must be 0
        """
        if n == 0:
            return QQbar_I
        else:
            raise IndexError("n must be 0")

    def ngens(self):
        r"""
        Return the size of the tuple returned by :meth:`gens`.

        EXAMPLE::

            sage: QQbar.ngens()
            1
        """
        return 1

    @cached_method
    def zeta(self, n=4):
        r"""
        Returns a primitive `n`'th root of unity, specifically `\exp(2*\pi*i/n)`.

        INPUT:

        - ``n`` (integer) -- default 4

        EXAMPLES::

            sage: QQbar.zeta(1)
            1
            sage: QQbar.zeta(2)
            -1
            sage: QQbar.zeta(3)
            -0.500000000000000? + 0.866025403784439?*I
            sage: QQbar.zeta(4)
            I
            sage: QQbar.zeta()
            I
            sage: QQbar.zeta(5)
            0.3090169943749474? + 0.9510565162951536?*I
            sage: QQbar.zeta(3000)
            0.999997806755380? + 0.002094393571219374?*I
        """
        if n == 1:
            return self.one()
        elif n == 2:
            return -self.one()
        elif n == 4:
            return self.gen()
        else:
            nf = CyclotomicField(n, embedding=CC.zeta(n))
            p = nf.polynomial()
            root = ANRoot(p, ComplexIntervalField(64).zeta(n))
            gen = AlgebraicGenerator(nf, root)
            return AlgebraicNumber(ANExtensionElement(gen, nf.gen()))

    def polynomial_root(self, poly, interval, multiplicity=1):
        r"""
        Given a polynomial with algebraic coefficients and an interval
        enclosing exactly one root of the polynomial, constructs
        an algebraic real representation of that root.

        The polynomial need not be irreducible, or even squarefree; but
        if the given root is a multiple root, its multiplicity must be
        specified. (IMPORTANT NOTE: Currently, multiplicity-`k` roots
        are handled by taking the `(k-1)`-st derivative of the polynomial.
        This means that the interval must enclose exactly one root
        of this derivative.)

        The conditions on the arguments (that the interval encloses exactly
        one root, and that multiple roots match the given multiplicity)
        are not checked; if they are not satisfied, an error may be
        thrown (possibly later, when the algebraic number is used),
        or wrong answers may result.

        Note that if you are constructing multiple roots of a single
        polynomial, it is better to use ``QQbar.common_polynomial``
        to get a shared polynomial.

        EXAMPLES::

            sage: x = polygen(QQbar)
            sage: phi = QQbar.polynomial_root(x^2 - x - 1, RIF(0, 2)); phi
            1.618033988749895?
            sage: p = (x-1)^7 * (x-2)
            sage: r = QQbar.polynomial_root(p, RIF(9/10, 11/10), multiplicity=7)
            sage: r; r == 1
            1
            True
            sage: p = (x-phi)*(x-sqrt(QQbar(2)))
            sage: r = QQbar.polynomial_root(p, RIF(1, 3/2))
            sage: r; r == sqrt(QQbar(2))
            1.414213562373095?
            True
        """
        return AlgebraicNumber(ANRoot(poly, interval, multiplicity))

    def random_element(self, poly_degree=2, *args, **kwds):
        r"""
        Returns a random algebraic number.

        INPUT:

        - ``poly_degree`` - default: 2 - degree of the random polynomial over
          the integers of which the returned algebraic number is a root. This
          is not necessarily the degree of the minimal polynomial of the
          number. Increase this parameter to achieve a greater diversity of
          algebraic numbers, at a cost of greater computation time. You can
          also vary the distribution of the coefficients but that will not vary
          the degree of the extension containing the element.

        - ``args``, ``kwds`` - arguments and keywords passed to the random
          number generator for elements of ``ZZ``, the integers. See
          :meth:`~sage.rings.integer_ring.IntegerRing_class.random_element` for
          details, or see example below.

        OUTPUT:

        An element of ``QQbar``, the field of algebraic numbers (see
        :mod:`sage.rings.qqbar`).

        ALGORITHM:

        A polynomial with degree between 1 and ``poly_degree``,
        with random integer coefficients is created. A root of this
        polynomial is chosen at random. The default degree is
        2 and the integer coefficients come from a distribution
        heavily weighted towards `0, \pm 1, \pm 2`.

        EXAMPLES::

            sage: a = QQbar.random_element()
            sage: a                         # random
            0.2626138748742799? + 0.8769062830975992?*I
            sage: a in QQbar
            True

            sage: b = QQbar.random_element(poly_degree=20)
            sage: b                         # random
            -0.8642649077479498? - 0.5995098147478391?*I
            sage: b in QQbar
            True

        Parameters for the distribution of the integer coefficients
        of the polynomials can be passed on to the random element method
        for integers. For example, current default behavior of this method
        returns zero about 15% of the time; if we do not include zero as a
        possible coefficient, there will never be a zero constant term, and
        thus never a zero root. ::

            sage: z = [QQbar.random_element(x=1, y=10) for _ in range(20)]
            sage: QQbar(0) in z
            False

         If you just want real algebraic numbers you can filter them out.
         Using an odd degree for the polynomials will insure some degree of
         success. ::

            sage: r = []
            sage: while len(r) < 3:
            ...     x = QQbar.random_element(poly_degree=3)
            ...     if x in AA:
            ...       r.append(x)
            sage: (len(r) == 3) and all([z in AA for z in r])
            True

        TESTS:

            sage: QQbar.random_element('junk')
            Traceback (most recent call last):
            ...
            TypeError: polynomial degree must be an integer, not junk
            sage: QQbar.random_element(poly_degree=0)
            Traceback (most recent call last):
            ...
            ValueError: polynomial degree must be greater than zero, not 0

        Random vectors already have a 'degree' keyword, so
        we cannot use that for the polynomial's degree. ::

            sage: v = random_vector(QQbar, degree=2, poly_degree=3)
            sage: v                                 # random
            (0.4694381338921299?, -0.500000000000000? + 0.866025403784439?*I)
        """
        import sage.rings.all
        import sage.misc.prandom
        try:
            poly_degree = sage.rings.all.ZZ(poly_degree)
        except TypeError:
            msg = "polynomial degree must be an integer, not {0}"
            raise TypeError(msg.format(poly_degree))
        if poly_degree < 1:
            msg = "polynomial degree must be greater than zero, not {0}"
            raise ValueError(msg.format(poly_degree))
        R = PolynomialRing(sage.rings.all.ZZ, 'x')
        p = R.random_element(degree=poly_degree, *args, **kwds)
        # degree zero polynomials have no roots
        # totally zero poly has degree -1
        # add a random leading term
        if p.degree() < 1:
            g = R.gen(0)
            m = sage.misc.prandom.randint(1, poly_degree)
            p = p + g**m
        roots = p.roots(ring=QQbar, multiplicities=False)

        # p will have at least one root; pick one at random
        # could we instead just compute one root "randomly"?
        m = sage.misc.prandom.randint(0, len(roots)-1)
        return roots[m]

    def _factor_univariate_polynomial(self, f):
        """
        Factor the univariate polynomial ``f``.

        INPUT:

        - ``f`` -- a univariate polynomial defined over the algebraic field

        OUTPUT:

        - A factorization of ``f`` over the algebraic numbers into a unit and
          monic irreducible factors

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.factor`.

        TESTS::

            sage: R.<x> = QQbar[]
            sage: QQbar._factor_univariate_polynomial(x)
            x
            sage: QQbar._factor_univariate_polynomial(2*x)
            (2) * x
            sage: QQbar._factor_univariate_polynomial((x^2 + 1)^2)
            (x - I)^2 * (x + I)^2
            sage: QQbar._factor_univariate_polynomial(x^8 - 1)
            (x - 1) * (x - 0.7071067811865475? - 0.7071067811865475?*I) * (x - 0.7071067811865475? + 0.7071067811865475?*I) * (x - I) * (x + I) * (x + 0.7071067811865475? - 0.7071067811865475?*I) * (x + 0.7071067811865475? + 0.7071067811865475?*I) * (x + 1)
            sage: QQbar._factor_univariate_polynomial(12*x^2 - 4)
            (12) * (x - 0.5773502691896258?) * (x + 0.5773502691896258?)
            sage: QQbar._factor_univariate_polynomial(R(-1))
            -1
            sage: QQbar._factor_univariate_polynomial(EllipticCurve('11a1').change_ring(QQbar).division_polynomial(5))
            (5) * (x - 16) * (x - 5) * (x - 1.959674775249769?) * (x - 1.427050983124843? - 3.665468789467727?*I) * (x - 1.427050983124843? + 3.665468789467727?*I) * (x + 0.9549150281252629? - 0.8652998037182486?*I) * (x + 0.9549150281252629? + 0.8652998037182486?*I) * (x + 1.927050983124843? - 1.677599044300515?*I) * (x + 1.927050983124843? + 1.677599044300515?*I) * (x + 2.959674775249769?) * (x + 6.545084971874737? - 7.106423590645660?*I) * (x + 6.545084971874737? + 7.106423590645660?*I)

        """
        from sage.structure.factorization import Factorization
        return Factorization([(f.parent()([-r,1]),e) for r,e in f.roots()], unit=f.leading_coefficient())

def is_AlgebraicField(F):
    r"""
    Check whether ``F`` is an :class:`~AlgebraicField` instance.

    EXAMPLE::

        sage: from sage.rings.qqbar import is_AlgebraicField
        sage: [is_AlgebraicField(x) for x in [AA, QQbar, None, 0, "spam"]]
        [False, True, False, False, False]
    """
    return isinstance(F, AlgebraicField)

# Create the globally unique AlgebraicField object.
QQbar = AlgebraicField()

def is_AlgebraicField_common(F):
    r"""
    Check whether ``F`` is an :class:`~AlgebraicField_common` instance.

    EXAMPLE::

        sage: from sage.rings.qqbar import is_AlgebraicField_common
        sage: [is_AlgebraicField_common(x) for x in [AA, QQbar, None, 0, "spam"]]
        [True, True, False, False, False]
    """
    return isinstance(F, AlgebraicField_common)

def prec_seq():
    r"""
    Return a generator object which iterates over an infinite increasing
    sequence of precisions to be tried in various numerical computations.

    Currently just returns powers of 2 starting at 64.

    EXAMPLE::

        sage: g = sage.rings.qqbar.prec_seq()
        sage: [next(g), next(g), next(g)]
        [64, 128, 256]
    """
    # XXX Should do some testing to see where the efficiency breaks are
    # in MPFR. We could also test variants like "bits = bits + bits // 2"
    # (I think this is what MPFR uses internally).
    bits = 64
    while True:
        yield bits
        bits = bits * 2

_short_prec_seq = (64, 128, None)
def short_prec_seq():
    r"""
    Return a sequence of precisions to try in cases when an infinite-precision
    computation is possible: returns a couple of small powers of 2 and then
    ``None``.

    EXAMPLE::

        sage: from sage.rings.qqbar import short_prec_seq
        sage: short_prec_seq()
        (64, 128, None)
    """
    return _short_prec_seq

def tail_prec_seq():
    r"""
    A generator over precisions larger than those in :func:`~short_prec_seq`.

    EXAMPLE::

        sage: from sage.rings.qqbar import tail_prec_seq
        sage: g = tail_prec_seq()
        sage: [next(g), next(g), next(g)]
        [256, 512, 1024]
    """
    bits = 256
    while True:
        yield bits
        bits = bits * 2

def rational_exact_root(r, d):
    r"""
    Checks whether the rational `r` is an exact `d`'th power. If so, returns
    the `d`'th root of `r`; otherwise, returns None.

    EXAMPLES::

        sage: from sage.rings.qqbar import rational_exact_root
        sage: rational_exact_root(16/81, 4)
        2/3
        sage: rational_exact_root(8/81, 3) is None
        True
    """
    num = r.numerator()
    den = r.denominator()

    (num_rt, num_exact) = num.nth_root(d, truncate_mode=1)
    if not num_exact: return None
    (den_rt, den_exact) = den.nth_root(d, truncate_mode=1)
    if not den_exact: return None
    return (num_rt / den_rt)

def clear_denominators(poly):
    """
    Takes a monic polynomial and rescales the variable to get a monic
    polynomial with "integral" coefficients. Works on any univariate
    polynomial whose base ring has a ``denominator()`` method that returns
    integers; for example, the base ring might be `\QQ` or a number
    field.

    Returns the scale factor and the new polynomial.

    (Inspired by Pari's ``primitive_pol_to_monic()``.)

    We assume that coefficient denominators are "small"; the algorithm
    factors the denominators, to give the smallest possible scale factor.

    EXAMPLES::

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
    r"""
    Find a polynomial of reasonably small discriminant that generates
    the same number field as ``poly``, using the PARI ``polredbest``
    function.

    INPUT:

    - ``poly`` - a monic irreducible polynomial with integer coefficients.

    OUTPUT:

    A triple (``elt_fwd``, ``elt_back``, ``new_poly``), where:

    - ``new_poly`` is the new polynomial defining the same number field,
    - ``elt_fwd`` is a polynomial expression for a root of the new
      polynomial in terms of a root of the original polynomial,
    - ``elt_back`` is a polynomial expression for a root of the original
      polynomial in terms of a root of the new polynomial.

    EXAMPLES::

        sage: from sage.rings.qqbar import do_polred
        sage: R.<x> = QQ['x']
        sage: oldpol = x^2 - 5
        sage: fwd, back, newpol = do_polred(oldpol)
        sage: newpol
        x^2 - x - 1
        sage: Kold.<a> = NumberField(oldpol)
        sage: Knew.<b> = NumberField(newpol)
        sage: newpol(fwd(a))
        0
        sage: oldpol(back(b))
        0
        sage: do_polred(x^2 - x - 11)
        (1/3*x + 1/3, 3*x - 1, x^2 - x - 1)
        sage: do_polred(x^3 + 123456)
        (1/4*x, 4*x, x^3 + 1929)

    This shows that :trac:`13054` has been fixed::

        sage: do_polred(x^4 - 4294967296*x^2 + 54265257667816538374400)
        (1/4*x, 4*x, x^4 - 268435456*x^2 + 211973662764908353025)
    """
    new_poly, elt_back = poly._pari_().polredbest(flag=1)

    parent = poly.parent()
    elt_fwd = elt_back.modreverse()
    return parent(elt_fwd.lift()), parent(elt_back.lift()), parent(new_poly)

def isolating_interval(intv_fn, pol):
    """
    ``intv_fn`` is a function that takes a precision and returns an
    interval of that precision containing some particular root of pol.
    (It must return better approximations as the precision increases.)
    pol is an irreducible polynomial with rational coefficients.

    Returns an interval containing at most one root of pol.

    EXAMPLES::

        sage: from sage.rings.qqbar import isolating_interval

        sage: _.<x> = QQ['x']
        sage: isolating_interval(lambda prec: sqrt(RealIntervalField(prec)(2)), x^2 - 2)
        1.4142135623730950488?

    And an example that requires more precision::

        sage: delta = 10^(-70)
        sage: p = (x - 1) * (x - 1 - delta) * (x - 1 + delta)
        sage: isolating_interval(lambda prec: RealIntervalField(prec)(1 + delta), p)
        1.000000000000000000000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000000000000000000000000000000000?

    The function also works with complex intervals and complex roots::

        sage: p = x^2 - x + 13/36
        sage: isolating_interval(lambda prec: ComplexIntervalField(prec)(1/2, 1/3), p)
        0.500000000000000000000? + 0.3333333333333333334?*I
    """
    dpol = pol.derivative()

    for prec in prec_seq():
        intv = intv_fn(prec)
        ifld = intv.parent()

        # We need to verify that pol has exactly one root in the
        # interval intv. We know (because it is a precondition of
        # calling this function) that it has at least one root in the
        # interval, so we only need to verify that it has at most one
        # root (that the interval is sufficiently narrow).

        # We do this by computing the derivative of the polynomial
        # over the interval. If the derivative is bounded away from zero,
        # then we know there can be at most one root.

        if not dpol(intv).contains_zero():
            return intv

def find_zero_result(fn, l):
    """
    ``l`` is a list of some sort. ``fn`` is a function which maps an element of
    ``l`` and a precision into an interval (either real or complex) of that
    precision, such that for sufficient precision, exactly one element of ``l``
    results in an interval containing 0. Returns that one element of ``l``.

    EXAMPLES::

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
            raise ValueError('find_zero_result could not find any zeroes')
        return result

def conjugate_expand(v):
    r"""
    If the interval ``v`` (which may be real or complex) includes some
    purely real numbers, return ``v'`` containing ``v`` such that
    ``v' == v'.conjugate()``. Otherwise return ``v`` unchanged. (Note that if
    ``v' == v'.conjugate()``, and ``v'`` includes one non-real root of a real
    polynomial, then ``v'`` also includes the conjugate of that root.
    Also note that the diameter of the return value is at most twice
    the diameter of the input.)

    EXAMPLES::

        sage: from sage.rings.qqbar import conjugate_expand
        sage: conjugate_expand(CIF(RIF(0, 1), RIF(1, 2))).str(style='brackets')
        '[0.00000000000000000 .. 1.0000000000000000] + [1.0000000000000000 .. 2.0000000000000000]*I'
        sage: conjugate_expand(CIF(RIF(0, 1), RIF(0, 1))).str(style='brackets')
        '[0.00000000000000000 .. 1.0000000000000000] + [-1.0000000000000000 .. 1.0000000000000000]*I'
        sage: conjugate_expand(CIF(RIF(0, 1), RIF(-2, 1))).str(style='brackets')
        '[0.00000000000000000 .. 1.0000000000000000] + [-2.0000000000000000 .. 2.0000000000000000]*I'
        sage: conjugate_expand(RIF(1, 2)).str(style='brackets')
        '[1.0000000000000000 .. 2.0000000000000000]'
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
    r"""
    If the interval ``v`` includes some purely real numbers, return
    a real interval containing only those real numbers. Otherwise
    return ``v`` unchanged.

    If ``v`` includes exactly one root of a real polynomial, and ``v`` was
    returned by ``conjugate_expand()``, then ``conjugate_shrink(v)`` still
    includes that root, and is a ``RealIntervalFieldElement`` iff the root
    in question is real.

    EXAMPLES::

        sage: from sage.rings.qqbar import conjugate_shrink
        sage: conjugate_shrink(RIF(3, 4)).str(style='brackets')
        '[3.0000000000000000 .. 4.0000000000000000]'
        sage: conjugate_shrink(CIF(RIF(1, 2), RIF(1, 2))).str(style='brackets')
        '[1.0000000000000000 .. 2.0000000000000000] + [1.0000000000000000 .. 2.0000000000000000]*I'
        sage: conjugate_shrink(CIF(RIF(1, 2), RIF(0, 1))).str(style='brackets')
        '[1.0000000000000000 .. 2.0000000000000000]'
        sage: conjugate_shrink(CIF(RIF(1, 2), RIF(-1, 2))).str(style='brackets')
        '[1.0000000000000000 .. 2.0000000000000000]'
    """
    if is_RealIntervalFieldElement(v):
        return v
    im = v.imag()
    if im.contains_zero():
        return v.real()
    return v

def number_field_elements_from_algebraics(numbers, minimal=False):
    r"""
    Given a sequence of elements of either ``AA`` or ``QQbar``
    (or a mixture), computes a number field containing all of these
    elements, these elements as members of that number field, and a
    homomorphism from the number field back to ``AA`` or
    ``QQbar``.

    This may not return the smallest such number field, unless
    ``minimal=True`` is specified.

    Also, a single number can be passed, rather than a sequence; and
    any values which are not elements of ``AA`` or ``QQbar``
    will automatically be coerced to ``QQbar``

    This function may be useful for efficiency reasons: doing exact
    computations in the corresponding number field will be faster
    than doing exact computations directly in ``AA`` or ``QQbar``.

    EXAMPLES:

    We can use this to compute the splitting field of a polynomial.
    (Unfortunately this takes an unreasonably long time for non-toy
    examples.)::

        sage: x = polygen(QQ)
        sage: p = x^3 + x^2 + x + 17
        sage: rts = p.roots(ring=QQbar, multiplicities=False)
        sage: splitting = number_field_elements_from_algebraics(rts)[0]; splitting
        Number Field in a with defining polynomial y^6 - 40*y^4 - 22*y^3 + 873*y^2 + 1386*y + 594
        sage: p.roots(ring=splitting)
        [(361/29286*a^5 - 19/3254*a^4 - 14359/29286*a^3 + 401/29286*a^2 + 18183/1627*a + 15930/1627, 1), (49/117144*a^5 - 179/39048*a^4 - 3247/117144*a^3 + 22553/117144*a^2 + 1744/4881*a - 17195/6508, 1), (-1493/117144*a^5 + 407/39048*a^4 + 60683/117144*a^3 - 24157/117144*a^2 - 56293/4881*a - 53033/6508, 1)]
        sage: rt2 = AA(sqrt(2)); rt2
        1.414213562373095?
        sage: rt3 = AA(sqrt(3)); rt3
        1.732050807568878?
        sage: qqI = QQbar.zeta(4); qqI
        I
        sage: z3 = QQbar.zeta(3); z3
        -0.500000000000000? + 0.866025403784439?*I
        sage: rt2b = rt3 + rt2 - rt3; rt2b
        1.414213562373095?
        sage: rt2c = z3 + rt2 - z3; rt2c
        1.414213562373095? + 0.?e-19*I

        sage: number_field_elements_from_algebraics(rt2)
        (Number Field in a with defining polynomial y^2 - 2, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 2
            To:   Algebraic Real Field
            Defn: a |--> 1.414213562373095?)

        sage: number_field_elements_from_algebraics((rt2,rt3))
        (Number Field in a with defining polynomial y^4 - 4*y^2 + 1, [-a^3 + 3*a, -a^2 + 2], Ring morphism:
            From: Number Field in a with defining polynomial y^4 - 4*y^2 + 1
            To:   Algebraic Real Field
            Defn: a |--> 0.5176380902050415?)

    We've created ``rt2b`` in such a way that \sage doesn't initially know
    that it's in a degree-2 extension of `\QQ`::

        sage: number_field_elements_from_algebraics(rt2b)
        (Number Field in a with defining polynomial y^4 - 4*y^2 + 1, -a^3 + 3*a, Ring morphism:
            From: Number Field in a with defining polynomial y^4 - 4*y^2 + 1
            To:   Algebraic Real Field
            Defn: a |--> 0.5176380902050415?)

    We can specify ``minimal=True`` if we want the smallest number field::

        sage: number_field_elements_from_algebraics(rt2b, minimal=True)
        (Number Field in a with defining polynomial y^2 - 2, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 2
            To:   Algebraic Real Field
            Defn: a |--> 1.414213562373095?)

    Things work fine with rational numbers, too::

        sage: number_field_elements_from_algebraics((QQbar(1/2), AA(17)))
        (Rational Field, [1/2, 17], Ring morphism:
            From: Rational Field
            To:   Algebraic Real Field
            Defn: 1 |--> 1)

    Or we can just pass in symbolic expressions, as long as they can be
    coerced into ``QQbar``::

        sage: number_field_elements_from_algebraics((sqrt(7), sqrt(9), sqrt(11)))
        (Number Field in a with defining polynomial y^4 - 9*y^2 + 1, [-a^3 + 8*a, 3, -a^3 + 10*a], Ring morphism:
            From: Number Field in a with defining polynomial y^4 - 9*y^2 + 1
            To:   Algebraic Real Field
            Defn: a |--> 0.3354367396454047?)

    Here we see an example of doing some computations with number field
    elements, and then mapping them back into ``QQbar``::

        sage: (fld,nums,hom) = number_field_elements_from_algebraics((rt2, rt3, qqI, z3))
        sage: fld,nums,hom  # random
        (Number Field in a with defining polynomial y^8 - y^4 + 1, [-a^5 + a^3 + a, a^6 - 2*a^2, a^6, -a^4], Ring morphism:
          From: Number Field in a with defining polynomial y^8 - y^4 + 1
          To:   Algebraic Field
          Defn: a |--> -0.2588190451025208? - 0.9659258262890683?*I)
        sage: (nfrt2, nfrt3, nfI, nfz3) = nums
        sage: hom(nfrt2)
        1.414213562373095? + 0.?e-18*I
        sage: nfrt2^2
        2
        sage: nfrt3^2
        3
        sage: nfz3 + nfz3^2
        -1
        sage: nfI^2
        -1
        sage: sum = nfrt2 + nfrt3 + nfI + nfz3; sum
        2*a^6 + a^5 - a^4 - a^3 - 2*a^2 - a
        sage: hom(sum)
        2.646264369941973? + 1.866025403784439?*I
        sage: hom(sum) == rt2 + rt3 + qqI + z3
        True
        sage: [hom(n) for n in nums] == [rt2, rt3, qqI, z3]
        True

    TESTS::

        sage: number_field_elements_from_algebraics(rt3)
        (Number Field in a with defining polynomial y^2 - 3, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 3
            To:   Algebraic Real Field
            Defn: a |--> 1.732050807568878?)
        sage: number_field_elements_from_algebraics((rt2,qqI))
        (Number Field in a with defining polynomial y^4 + 1, [-a^3 + a, a^2], Ring morphism:
            From: Number Field in a with defining polynomial y^4 + 1
            To:   Algebraic Field
            Defn: a |--> 0.7071067811865475? + 0.7071067811865475?*I)

    Note that for the first example, where \sage doesn't realize that
    the number is real, we get a homomorphism to ``QQbar``; but with
    ``minimal=True``, we get a homomorphism to ``AA``. Also note
    that the exact answer depends on a Pari function that gives
    different answers for 32-bit and 64-bit machines::

        sage: number_field_elements_from_algebraics(rt2c)
        (Number Field in a with defining polynomial y^4 + 2*y^2 + 4, 1/2*a^3, Ring morphism:
            From: Number Field in a with defining polynomial y^4 + 2*y^2 + 4
            To:   Algebraic Field
            Defn: a |--> -0.7071067811865475? - 1.224744871391589?*I)
        sage: number_field_elements_from_algebraics(rt2c, minimal=True)
        (Number Field in a with defining polynomial y^2 - 2, a, Ring morphism:
            From: Number Field in a with defining polynomial y^2 - 2
            To:   Algebraic Real Field
            Defn: a |--> 1.414213562373095?)

    """
    gen = qq_generator

    # Keep track of whether we were given a single value or a list.
    single_number = False
    try:
        len(numbers)
    except TypeError:
        numbers = [numbers]
        single_number = True

    def mk_algebraic(x):
        if isinstance(x, AlgebraicNumber_base):
            return x
        return QQbar(x)

    numbers = [mk_algebraic(_) for _ in numbers]

    for v in numbers:
        if minimal:
            v.simplify()
        gen = gen.union(v._exact_field())

    fld = gen._field

    nums = [gen(v._exact_value()) for v in numbers]

    if single_number:
        nums = nums[0]

    hom = fld.hom([gen.root_as_algebraic()])

    return (fld, nums, hom)

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
        r"""
        EXAMPLE::

            sage: from sage.rings.qqbar import AlgebraicGeneratorRelation
            sage: c = AlgebraicGeneratorRelation(None, None, None, None, None)
            sage: c
            <class 'sage.rings.qqbar.AlgebraicGeneratorRelation'>
        """
        self.child1 = child1
        self.child1_poly = child1_poly
        self.child2 = child2
        self.child2_poly = child2_poly
        self.parent = parent

algebraic_generator_counter = 0

class AlgebraicGenerator(SageObject):
    r"""
    An ``AlgebraicGenerator`` represents both an algebraic number `\alpha` and
    the number field `\QQ[\alpha]`. There is a single ``AlgebraicGenerator``
    representing `\QQ` (with `\alpha=0`).

    The ``AlgebraicGenerator`` class is private, and should not be used
    directly.
    """

    def __init__(self, field, root):
        """
        Construct an ``AlgebraicGenerator`` object.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: y = polygen(QQ, 'y')
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen
            Number Field in a with defining polynomial y^2 - y - 1 with a in 1.618033988749895?
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
            Number Field in a with defining polynomial y^2 + 1 with a in 1*I
        """
        self._field = field
        self._pari_field = None
        self._trivial = (field is QQ)
        self._root = root
        self._root_as_algebraic = (QQbar if root.is_complex() else AA)(root)
        self._unions = {}
        self._cyclotomic = False
        global algebraic_generator_counter
        self._index = algebraic_generator_counter
        algebraic_generator_counter += 1

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = QQbar(sqrt(2)) + QQbar(sqrt(3))
            sage: t.exactify()
            sage: type(t._descr._generator)
            <class 'sage.rings.qqbar.AlgebraicGenerator'>
            sage: loads(dumps(t)) == t
            True
        """
        return (AlgebraicGenerator, (self._field, self._root))

    def __hash__(self):
        r"""
        Return a hash value for self. This will depend on the order that
        commands get executed at load time, so we do not test the value that is
        returned, just that it doesn't raise an error.

        EXAMPLE::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: hash(gen) # random
        """
        return self._index

    def __cmp__(self, other):
        r"""
        Compare self with another AlgebraicGenerator object.

        EXAMPLE::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen.__cmp__(qq_generator)
            1
        """
        return cmp(self._index, other._index)

    def is_complex(self):
        r"""
        Return True if this is a generator for a non-real number field.

        EXAMPLE::

            sage: z7 = QQbar.zeta(7)
            sage: g = z7._descr._generator
            sage: g.is_complex()
            True

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator
            sage: y = polygen(QQ, 'y')
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen.is_complex()
            False
        """
        return self._root.is_complex()

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: from sage.rings.qqbar import qq_generator
            sage: qq_generator._repr_()
            'Trivial generator'

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: y = polygen(QQ)
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen._repr_()
            'Number Field in a with defining polynomial x^2 - x - 1 with a in 1.618033988749895?'
        """
        if self._trivial:
            return 'Trivial generator'
        else:
            return '%s with a in %s'%(self._field, self._root._interval_fast(53))

    def root_as_algebraic(self):
        r"""
        Return the root attached to self as an algebraic number.

        EXAMPLE::

            sage: t = sage.rings.qqbar.qq_generator.root_as_algebraic(); t
            1
            sage: t.parent()
            Algebraic Real Field
        """
        return self._root_as_algebraic

    def is_trivial(self):
        """
        Returns true iff this is the trivial generator (alpha == 1), which
        does not actually extend the rationals.

        EXAMPLES::

            sage: from sage.rings.qqbar import qq_generator
            sage: qq_generator.is_trivial()
            True
        """
        return self._trivial

    def field(self):
        r"""
        Return the number field attached to self.

        EXAMPLE::

            sage: from sage.rings.qqbar import qq_generator
            sage: qq_generator.field()
            Rational Field
        """
        return self._field

    def pari_field(self):
        r"""
        Return the PARI field attached to this generator.

        EXAMPLE::


            sage: from sage.rings.qqbar import qq_generator
            sage: qq_generator.pari_field()
            Traceback (most recent call last):
            ...
            ValueError: No PARI field attached to trivial generator

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: y = polygen(QQ)
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen.pari_field()
             [y^2 - y - 1, [2, 0], ...]
        """
        if self.is_trivial(): raise ValueError("No PARI field attached to trivial generator")
        if self._pari_field is None:
            pari_pol = self._field.pari_polynomial("y")
            self._pari_field = pari_pol.nfinit(1)
        return self._pari_field

    def conjugate(self):
        r"""
        If this generator is for the algebraic number `\alpha`, return a
        generator for the complex conjugate of `\alpha`.

        EXAMPLE::

            sage: from sage.rings.qqbar import AlgebraicGenerator
            sage: x = polygen(QQ); f = x^4 + x + 17
            sage: nf = NumberField(f,name='a')
            sage: b = f.roots(QQbar)[0][0]
            sage: root = b._descr
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen.conjugate()
            Number Field in a with defining polynomial x^4 + x + 17 with a in -1.436449997483091? + 1.374535713065812?*I
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

        EXAMPLE::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: y = polygen(QQ, 'y')
            sage: x = polygen(QQbar)
            sage: nf = NumberField(y^2 - y - 1, name='a', check=False)
            sage: root = ANRoot(x^2 - x - 1, RIF(1, 2))
            sage: gen = AlgebraicGenerator(nf, root)
            sage: gen._interval_fast(128)
            1.61803398874989484820458683436563811773?
        """
        return self._root._interval_fast(prec)

    def union(self, other):
        r""" Given generators ``alpha`` and ``beta``,
        ``alpha.union(beta)`` gives a generator for the number field
        `\QQ[\alpha][\beta]`.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in 1.414213562373095?
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = ANRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in 1.732050807568878?
            sage: gen2.union(qq_generator) is gen2
            True
            sage: qq_generator.union(gen3) is gen3
            True
            sage: gen2.union(gen3)
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in 0.5176380902050415?
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
        r"""
        Given a generator ``gen`` and another generator ``super``, where ``super``
        is the result of a tree of ``union()`` operations where one of the
        leaves is ``gen``, ``gen.super_poly(super)`` returns a polynomial
        expressing the value of ``gen`` in terms of the value of ``super``
        (except that if ``gen`` is ``qq_generator``, ``super_poly()`` always
        returns None.)

        EXAMPLES::

            sage: from sage.rings.qqbar import AlgebraicGenerator, ANRoot, qq_generator
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in 1.414213562373095?
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = ANRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in 1.732050807568878?
            sage: gen2_3 = gen2.union(gen3)
            sage: gen2_3
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in 0.5176380902050415?
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
        Takes an algebraic number which is represented as either a
        rational or a number field element, and which is in a subfield
        of the field generated by this generator. Lifts the number
        into the field of this generator, and returns either a
        ``Rational`` or a ``NumberFieldElement`` depending on whether
        this is the trivial generator.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot, AlgebraicGenerator, ANExtensionElement, ANRational
            sage: _.<y> = QQ['y']
            sage: x = polygen(QQbar)
            sage: nf2 = NumberField(y^2 - 2, name='a', check=False)
            sage: root2 = ANRoot(x^2 - 2, RIF(1, 2))
            sage: gen2 = AlgebraicGenerator(nf2, root2)
            sage: gen2
            Number Field in a with defining polynomial y^2 - 2 with a in 1.414213562373095?
            sage: sqrt2 = ANExtensionElement(gen2, nf2.gen())
            sage: nf3 = NumberField(y^2 - 3, name='a', check=False)
            sage: root3 = ANRoot(x^2 - 3, RIF(1, 2))
            sage: gen3 = AlgebraicGenerator(nf3, root3)
            sage: gen3
            Number Field in a with defining polynomial y^2 - 3 with a in 1.732050807568878?
            sage: sqrt3 = ANExtensionElement(gen3, nf3.gen())
            sage: gen2_3 = gen2.union(gen3)
            sage: gen2_3
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in 0.5176380902050415?
            sage: gen2_3(sqrt2)
            -a^3 + 3*a
            sage: gen2_3(ANRational(1/7))
            1/7
            sage: gen2_3(sqrt3)
            -a^2 + 2
        """
        if self._trivial:
            return elt._value
        if isinstance(elt, ANRational):
            return self._field(elt._value)
        if elt.generator() is self:
            return elt.field_element_value()
        gen = elt.generator()
        sp = gen.super_poly(self)
        assert(not(sp is None))
        return self._field(elt.field_element_value().polynomial()(sp))


# dictionary for multimethod dispatch
_binop_algo = {}

class ANDescr(SageObject):
    r"""
    An ``AlgebraicNumber`` or ``AlgebraicReal`` is a wrapper around an
    ``ANDescr`` object. ``ANDescr`` is an abstract base class, which should
    never be directly instantiated; its concrete subclasses are ``ANRational``,
    ``ANBinaryExpr``, ``ANUnaryExpr``, ``ANRoot``, and ``ANExtensionElement``.
    ``ANDescr`` and all of its subclasses are for internal use, and should not
    be used directly.
    """
    def is_simple(self):
        r"""
        Checks whether this descriptor represents a value with the same
        algebraic degree as the number field associated with the descriptor.

        Returns ``True`` if self is an ``ANRational``, or a minimal
        ``ANExtensionElement``.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRational
            sage: ANRational(1/2).is_simple()
            True
            sage: rt2 = AA(sqrt(2))
            sage: rt3 = AA(sqrt(3))
            sage: rt2b = rt3 + rt2 - rt3
            sage: rt2.exactify()
            sage: rt2._descr.is_simple()
            True
            sage: rt2b.exactify()
            sage: rt2b._descr.is_simple()
            False
            sage: rt2b.simplify()
            sage: rt2b._descr.is_simple()
            True
        """
        return False

    # Unitary operators: the second argument "n" is an AlgebraicNumber_base
    # wrapper around self.

    def neg(self, n):
        r"""
        Negation of self.

        EXAMPLE::

            sage: a = QQbar(sqrt(2))
            sage: b = a._descr
            sage: b.neg(a)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        return ANUnaryExpr(n, '-')

    def invert(self, n):
        r"""
        1/self.

        EXAMPLE::

            sage: a = QQbar(sqrt(2))
            sage: b = a._descr
            sage: b.invert(a)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        return ANUnaryExpr(n, '~')

    def abs(self, n):
        r"""
        Absolute value of self.

        EXAMPLE::

            sage: a = QQbar(sqrt(2))
            sage: b = a._descr
            sage: b.abs(a)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        return ANUnaryExpr(n, 'abs')

    def real(self, n):
        r"""
        Real part of self.

        EXAMPLE::

            sage: a = QQbar(sqrt(-7))
            sage: b = a._descr
            sage: b.real(a)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        if self.is_complex():
            return ANUnaryExpr(n, 'real')
        else:
            return self

    def imag(self, n):
        r"""
        Imaginary part of self.

        EXAMPLE::

            sage: a = QQbar(sqrt(-7))
            sage: b = a._descr
            sage: b.imag(a)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        if self.is_complex():
            return ANUnaryExpr(n, 'imag')
        else:
            return ANRational(0)

    def conjugate(self, n):
        r"""
        Complex conjugate of self.

        EXAMPLE::

            sage: a = QQbar(sqrt(-7))
            sage: b = a._descr
            sage: b.conjugate(a)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        if self.is_complex():
            return ANUnaryExpr(n, 'conjugate')
        else:
            return self

    def norm(self, n):
        r"""
        Field norm of self from `\overline{\QQ}` to its real subfield
        `\mathbf{A}`, i.e.~the square of the usual complex absolute value.

        EXAMPLE::

            sage: a = QQbar(sqrt(-7))
            sage: b = a._descr
            sage: b.norm(a)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        if self.is_complex():
            return ANUnaryExpr(n, 'norm')
        else:
            return (n*n)._descr


class AlgebraicNumber_base(sage.structure.element.FieldElement):
    r"""
    This is the common base class for algebraic numbers (complex
    numbers which are the zero of a polynomial in `\ZZ[x]`) and algebraic
    reals (algebraic numbers which happen to be real).

    ``AlgebraicNumber`` objects can be created using ``QQbar`` (==
    ``AlgebraicNumberField()``), and ``AlgebraicReal`` objects can be created
    using ``AA`` (== ``AlgebraicRealField()``). They can be created either by
    coercing a rational or a symbolic expression, or by using the
    ``QQbar.polynomial_root()`` or ``AA.polynomial_root()`` method to
    construct a particular root of a polynomial with algebraic
    coefficients. Also, ``AlgebraicNumber`` and ``AlgebraicReal`` are closed
    under addition, subtraction, multiplication, division (except by
    0), and rational powers (including roots), except that for a
    negative ``AlgebraicReal``, taking a power with an even denominator returns
    an ``AlgebraicNumber`` instead of an ``AlgebraicReal``.

    ``AlgebraicNumber``   and   ``AlgebraicReal``   objects   can   be
    approximated  to  any desired  precision. They  can be  compared
    exactly; if the two numbers are very close, or are equal, this may
    require exact computation, which can be extremely slow.

    As long as exact computation is not triggered, computation with
    algebraic numbers should not be too much slower than computation with
    intervals. As mentioned above, exact computation is triggered
    when comparing two algebraic numbers which are very close together.
    This can be an explicit comparison in user code, but the following
    list of actions (not necessarily complete) can also trigger exact
    computation:

    - Dividing by an algebraic number which is very close to 0.

    - Using an algebraic number which is very close to 0 as the leading
      coefficient in a polynomial.

    - Taking a root of an algebraic number which is very close to 0.

    The exact definition of "very close" is subject to change; currently,
    we compute our best approximation of the two numbers using 128-bit
    arithmetic, and see if that's sufficient to decide the comparison.
    Note that comparing two algebraic numbers which are actually equal will
    always trigger exact computation, unless they are actually the same object.

    EXAMPLES::

        sage: sqrt(QQbar(2))
        1.414213562373095?
        sage: sqrt(QQbar(2))^2 == 2
        True
        sage: x = polygen(QQbar)
        sage: phi = QQbar.polynomial_root(x^2 - x - 1, RIF(1, 2))
        sage: phi
        1.618033988749895?
        sage: phi^2 == phi+1
        True
        sage: AA(sqrt(65537))
        256.0019531175495?
    """

    def __init__(self, parent, x):
        r"""
        Initialize an algebraic number. The argument must be either
        a rational number, a Gaussian rational, or a subclass of ``ANDescr``.

        EXAMPLES::

            sage: AlgebraicReal(22/7)
            22/7
        """
        sage.structure.element.FieldElement.__init__(self, parent)
        if isinstance(x, (int, long, sage.rings.integer.Integer,
                          sage.rings.rational.Rational)):
            self._descr = ANRational(x)
        elif isinstance(x, ANDescr):
            self._descr = x
        elif parent is QQbar and \
                 isinstance(x, NumberFieldElement_quadratic) and \
                 list(x.parent().polynomial()) == [1, 0, 1]:
            self._descr = ANExtensionElement(QQbar_I_generator, QQbar_I_nf(x.list()))
        else:
            raise TypeError("Illegal initializer for algebraic number")

        self._value = self._descr._interval_fast(64)

    def _repr_(self):
        """
        Returns the print representation of this number.

        EXAMPLES::

            sage: AA(22/7) # indirect doctest
            22/7
            sage: QQbar(1/3 + 2/7*I)
            2/7*I + 1/3
            sage: QQbar.zeta(4) + 5
            I + 5
            sage: QQbar.zeta(4)
            I
            sage: 3*QQbar.zeta(4)
            3*I
            sage: QQbar.zeta(17)
            0.9324722294043558? + 0.3612416661871530?*I
            sage: AA(19).sqrt()
            4.358898943540674?
        """
        if isinstance(self._descr, ANRational):
            return repr(self._descr)
        if isinstance(self._descr, ANExtensionElement) and self._descr._generator is QQbar_I_generator:
            return repr(self._descr._value)
        if self.parent() is QQbar:
            return repr(CIF(self._value))
        else:
            return repr(RIF(self._value))

    def _latex_(self):
        r"""
        Returns the latex representation of this number.

        EXAMPLES::

            sage: latex(AA(22/7))
            \frac{22}{7}
            sage: latex(QQbar(1/3 + 2/7*I))
            \frac{2}{7} \sqrt{-1} + \frac{1}{3}
            sage: latex(QQbar.zeta(4) + 5)
            \sqrt{-1} + 5
            sage: latex(QQbar.zeta(4))
            \sqrt{-1}
            sage: latex(3*QQbar.zeta(4))
            3 \sqrt{-1}
            sage: latex(QQbar.zeta(17))
            0.9324722294043558? + 0.3612416661871530? \sqrt{-1}
            sage: latex(AA(19).sqrt())
            4.358898943540674?
        """
        from sage.misc.latex import latex
        if isinstance(self._descr, ANRational):
            return latex(self._descr._value)
        if isinstance(self._descr, ANExtensionElement) and self._descr._generator is QQbar_I_generator:
            return latex(self._descr._value)
        return repr(self).replace('*I', r' \sqrt{-1}')

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES:

        These examples are mostly copied from the doctests of
        the ``handle_sage_input`` functions; see those for more examples::

            sage: sage_input(QQbar(3))
            QQbar(3)
            sage: sage_input(AA(22/7))
            AA(22/7)
            sage: sage_input(22/7*QQbar.zeta(4))
            QQbar(22/7*I)
            sage: sage_input(QQbar.zeta(5)^3)
            R.<x> = AA[]
            QQbar.polynomial_root(AA.common_polynomial(x^4 + x^3 + x^2 + x + 1), CIF(RIF(RR(0.3090169943749474), RR(0.30901699437494745)), RIF(RR(0.95105651629515353), RR(0.95105651629515364))))^3
            sage: sage_input((AA(3)^(1/2))^(1/3))
            R.<x> = AA[]
            AA.polynomial_root(AA.common_polynomial(x^3 - AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))), RIF(RR(1.2009369551760025), RR(1.2009369551760027)))
            sage: sage_input(QQbar(3+4*I))
            QQbar(3 + 4*I)
            sage: sage_input(-sqrt(AA(2)))
            R.<x> = AA[]
            -AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(2 + sqrt(AA(2)))
            R.<x> = AA[]
            2 + AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))

        And a nice big example::

            sage: K.<x> = QQ[]
            sage: p = K.random_element(4); p
            1/2*x^4 - 1/95*x^3 - 1/2*x^2 - 4
            sage: rts = p.roots(ring=QQbar, multiplicities=False); rts
            [-1.830225346898784?, 1.842584249981426?, 0.004346864248152390? - 1.540200655088741?*I, 0.004346864248152390? + 1.540200655088741?*I]
            sage: sage_input(rts, verify=True)  # long time (2s on sage.math, 2013)
            # Verified
            R.<x> = AA[]
            cp = AA.common_polynomial(1/2*x^4 - 1/95*x^3 - 1/2*x^2 - 4)
            [QQbar.polynomial_root(cp, CIF(RIF(-RR(1.8302253468987832), -RR(1.830225346898783)), RIF(RR(0)))), QQbar.polynomial_root(cp, CIF(RIF(RR(1.8425842499814258), RR(1.842584249981426)), RIF(RR(0)))), QQbar.polynomial_root(cp, CIF(RIF(RR(0.0043468642481523899), RR(0.0043468642481523908)), RIF(-RR(1.5402006550887404), -RR(1.5402006550887402)))), QQbar.polynomial_root(cp, CIF(RIF(RR(0.0043468642481523899), RR(0.0043468642481523908)), RIF(RR(1.5402006550887402), RR(1.5402006550887404))))]

            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sqrt(QQbar(7))._sage_input_(sib, False)
            {call: {getattr: {atomic:QQbar}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:QQbar}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:7}})}, {call: {atomic:CIF}({call: {atomic:RIF}({call: {atomic:RR}({atomic:2.6457513110645903})}, {call: {atomic:RR}({atomic:2.6457513110645907})})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:0})})})})}
        """
        (v, complicated) = \
            self._descr.handle_sage_input(sib, coerce, self.parent() is QQbar)
        if complicated or True:
            sib.id_cache(self, v, 'v')
        return v

    def _mul_(self, other):
        """
        TESTS::

            sage: AA(sqrt(2)) * AA(sqrt(8)) # indirect doctest
            4.000000000000000?
        """
        sk = type(self._descr)
        ok = type(other._descr)
        return type(self)(_binop_algo[sk,ok](self, other, operator.mul))

    def _div_(self, other):
        """
        TESTS::

            sage: AA(sqrt(2)) / AA(sqrt(8)) # indirect doctest
            0.500000000000000?
        """
        sk = type(self._descr)
        ok = type(other._descr)
        return type(self)(_binop_algo[sk,ok](self, other, operator.div))

    def __invert__(self):
        """
        TESTS::

            sage: ~AA(sqrt(~2))
            1.414213562373095?
        """
        sd = self._descr
        return type(self)(self._descr.invert(self))

    def _add_(self, other):
        """
        TESTS::

            sage: x = polygen(ZZ)
            sage: rt1, rt2 = (x^2 - x - 1).roots(ring=AA, multiplicities=False)
            sage: rt1 + rt2 # indirect doctest
            1.000000000000000?
        """
        sk = type(self._descr)
        ok = type(other._descr)
        return type(self)(_binop_algo[sk,ok](self, other, operator.add))

    def _sub_(self, other):
        """
        TESTS::

            sage: AA(golden_ratio) * 2 - AA(5).sqrt() # indirect doctest
            1.000000000000000?
        """
        sk = type(self._descr)
        ok = type(other._descr)
        return type(self)(_binop_algo[sk,ok](self, other, operator.sub))

    def _neg_(self):
        """
        TESTS::

            sage: -QQbar(I) # indirect doctest
            -I
        """
        return type(self)(self._descr.neg(self))

    def __abs__(self):
        """
        TESTS::

            sage: abs(AA(sqrt(2) - sqrt(3)))
            0.3178372451957823?
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

        The hash code is stable, even when the representation changes::

            sage: two = QQbar(4).nth_root(4)^2
            sage: two
            2.000000000000000?
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
        enough to each other are the same. (This is inevitable, if
        equal algebraic reals give the same hash code and hashing does
        not always trigger exact computation.)::

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
        # calling interval_exact(). Then, exact computation will be triggered
        # by algebraic reals which are sufficiently close to
        # (some floating point number minus 1/123456789). Hopefully,
        # -1/123456789 comes up in algebraic real computations far less
        # often than 0 does. Algebraic numbers have a similar offset added,
        # with an additional complex component of 1/987654321*I.

        # All of this effort to avoid exact computation is probably wasted,
        # anyway... in almost all uses of hash codes, if the hash codes
        # match, the next step is to compare for equality; and comparing
        # for equality often requires exact computation. (If a==b,
        # then checking a==b requires exact computation unless (a is b).)

        if self.parent() is AA:
            return hash((self + AA_hash_offset).interval_exact(RIF))
        else:
            return hash((self + QQbar_hash_offset).interval_exact(CIF))

    def is_square(self):
        """
        Return whether or not this number is square.

        OUTPUT:

        (boolean) True in all cases for elements of QQbar; True for
        non-negative elements of AA, otherwise False.

        EXAMPLES::

            sage: AA(2).is_square()
            True
            sage: AA(-2).is_square()
            False
            sage: QQbar(-2).is_square()
            True
            sage: QQbar(I).is_square()
            True
        """
        if self.parent() is AA:
            return bool(self >= 0)
        else:
            return True

    def is_integer(self):
        """
        Return True if this number is a integer

        EXAMPLES::

            sage: QQbar(2).is_integer()
            True
            sage: QQbar(1/2).is_integer()
            False
        """
        return self in ZZ

    def sqrt(self, all=False, extend=True):
        """
        Return the square root(s) of this number.

        INPUT:

        - ``extend`` - bool (default: True); ignored if self is in QQbar, or
          positive in AA. If self is negative in AA, do the following: if True,
          return a square root of self in QQbar, otherwise raise a ValueError.

        - ``all`` - bool (default: False); if True, return a list of all square
          roots. If False, return just one square root, or raise an ValueError
          if self is a negative element of AA and extend=False.

        OUTPUT:

        Either the principal square root of self, or a list of its
        square roots (with the principal one first).

        EXAMPLES::

            sage: AA(2).sqrt()
            1.414213562373095?

            sage: QQbar(I).sqrt()
            0.7071067811865475? + 0.7071067811865475?*I
            sage: QQbar(I).sqrt(all=True)
            [0.7071067811865475? + 0.7071067811865475?*I, -0.7071067811865475? - 0.7071067811865475?*I]

            sage: a = QQbar(0)
            sage: a.sqrt()
            0
            sage: a.sqrt(all=True)
            [0]

            sage: a = AA(0)
            sage: a.sqrt()
            0
            sage: a.sqrt(all=True)
            [0]

        This second example just shows that the program doesn't care where 0
        is defined, it gives the same answer regardless. After all, how many
        ways can you square-root zero?

        ::

            sage: AA(-2).sqrt()
            1.414213562373095?*I

            sage: AA(-2).sqrt(all=True)
            [1.414213562373095?*I, -1.414213562373095?*I]

            sage: AA(-2).sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: -2 is not a square in AA, being negative. Use extend = True for a square root in QQbar.


        """
        # deal with 0 first:

        if self.is_zero():
            if all:
                return [self]
            else:
                return self

        # raise an error if appropriate:

        if self.parent() is AA and self<0 and not extend:
            if not all:
                raise ValueError("%s is not a square in AA, being negative. Use extend = True for a square root in QQbar."%self)
            else:
                return []

        root = self ** ~ZZ(2)

        if all:
            return [root, -root]
        else:
           return root

    def nth_root(self, n, all=False):
        r"""
        Return the ``n``-th root of this number.

        INPUT:

        -  ``all`` - bool (default: ``False``). If ``True``, return a list of
           all `n`-th roots as complex algebraic numbers.

        .. WARNING::

            Note that for odd `n`, all=`False` and negative real numbers,
            ``AlgebraicReal`` and ``AlgebraicNumber`` values give different
            answers: ``AlgebraicReal`` values prefer real results, and
            ``AlgebraicNumber`` values return the principal root.

        EXAMPLES::

            sage: AA(-8).nth_root(3)
            -2
            sage: QQbar(-8).nth_root(3)
            1.000000000000000? + 1.732050807568878?*I
            sage: QQbar.zeta(12).nth_root(15)
            0.9993908270190957? + 0.03489949670250097?*I

        You can get all ``n``-th roots of algebraic numbers::

            sage: AA(-8).nth_root(3, all=True)
            [1.000000000000000? + 1.732050807568878?*I,
            -2.000000000000000? + 0.?e-18*I,
            1.000000000000000? - 1.732050807568878?*I]

            sage: QQbar(1+I).nth_root(4, all=True)
            [1.069553932363986? + 0.2127475047267431?*I,
             -0.2127475047267431? + 1.069553932363986?*I,
             -1.069553932363986? - 0.2127475047267431?*I,
             0.2127475047267431? - 1.069553932363986?*I]

        TESTS::

            sage: AA(-8).nth_root(3, all=True)[1]
            -2.000000000000000? + 0.?e-18*I
            sage: _.parent()
            Algebraic Field

            sage: AA(-2).nth_root(5, all=True) == QQbar(-2).nth_root(5, all=True)   # long time
            True
        """
        if not all:
            return self ** ~ZZ(n)
        else:
            root = QQbar(self) ** ~ZZ(n)
            zlist = [root]
            zeta = QQbar.zeta(n)
            for k in range(1, n):
                root *= zeta
                zlist.append(root)
            return zlist

    def as_number_field_element(self, minimal=False):
        r"""
        Returns a number field containing this value, a representation of
        this value as an element of that number field, and a homomorphism
        from the number field back to ``AA`` or ``QQbar``.

        This may not return the smallest such number field, unless
        ``minimal=True`` is specified.

        To compute a single number field containing multiple algebraic
        numbers, use the function
        ``number_field_elements_from_algebraics`` instead.

        EXAMPLES::

            sage: QQbar(sqrt(8)).as_number_field_element()
            (Number Field in a with defining polynomial y^2 - 2, 2*a, Ring morphism:
                From: Number Field in a with defining polynomial y^2 - 2
                To:   Algebraic Real Field
                Defn: a |--> 1.414213562373095?)
            sage: x = polygen(ZZ)
            sage: p = x^3 + x^2 + x + 17
            sage: (rt,) = p.roots(ring=AA, multiplicities=False); rt
            -2.804642726932742?
            sage: (nf, elt, hom) = rt.as_number_field_element()
            sage: nf, elt, hom
            (Number Field in a with defining polynomial y^3 - 2*y^2 - 31*y - 50, a^2 - 5*a - 19, Ring morphism:
              From: Number Field in a with defining polynomial y^3 - 2*y^2 - 31*y - 50
              To:   Algebraic Real Field
              Defn: a |--> 7.237653139801104?)
            sage: hom(elt) == rt
            True

        We see an example where we do not get the minimal number field unless
        we specify ``minimal=True``::

            sage: rt2 = AA(sqrt(2))
            sage: rt3 = AA(sqrt(3))
            sage: rt3b = rt2 + rt3 - rt2
            sage: rt3b.as_number_field_element()
            (Number Field in a with defining polynomial y^4 - 4*y^2 + 1, -a^2 + 2, Ring morphism:
                From: Number Field in a with defining polynomial y^4 - 4*y^2 + 1
                To:   Algebraic Real Field
                Defn: a |--> 0.5176380902050415?)
            sage: rt3b.as_number_field_element(minimal=True)
            (Number Field in a with defining polynomial y^2 - 3, a, Ring morphism:
                From: Number Field in a with defining polynomial y^2 - 3
                To:   Algebraic Real Field
                Defn: a |--> 1.732050807568878?)
        """
        return number_field_elements_from_algebraics(self, minimal=minimal)

    def exactify(self):
        """
        Compute an exact representation for this number.

        EXAMPLES::

            sage: two = QQbar(4).nth_root(4)^2
            sage: two
            2.000000000000000?
            sage: two.exactify()
            sage: two
            2
        """
        od = self._descr
        if isinstance(od, (ANRational, ANExtensionElement)): return
        self._set_descr(self._descr.exactify())

    def _set_descr(self, new_descr):
        """
        Set ``self._descr`` to ``new_descr``, and update
        ``self._value`` accordingly.

        EXAMPLES::

            sage: c = QQbar(-1)**(1/3) - QQbar(3)**(1/2)/2*QQbar.gen()
            sage: c._value
            0.5000000000000000000? + 0.?e-19*I
            sage: c.exactify()   # indirect doctest
            sage: c._value
            0.500000000000000000000?
        """
        self._descr = new_descr
        new_val = self._descr._interval_fast(self.parent().default_interval_prec())
        if is_RealIntervalFieldElement(new_val) and is_ComplexIntervalFieldElement(self._value):
            self._value = self._value.real().intersection(new_val)
        elif is_RealIntervalFieldElement(self._value) and is_ComplexIntervalFieldElement(new_val):
            self._value = self._value.intersection(new_val.real())
        else:
            self._value = self._value.intersection(new_val)

    def simplify(self):
        """
        Compute an exact representation for this number, in the
        smallest possible number field.

        EXAMPLES::

            sage: rt2 = AA(sqrt(2))
            sage: rt3 = AA(sqrt(3))
            sage: rt2b = rt3 + rt2 - rt3
            sage: rt2b.exactify()
            sage: rt2b._exact_value()
            a^3 - 3*a where a^4 - 4*a^2 + 1 = 0 and a in 1.931851652578137?
            sage: rt2b.simplify()
            sage: rt2b._exact_value()
            a where a^2 - 2 = 0 and a in 1.414213562373095?
        """
        self.exactify()
        od = self._descr
        if od.is_simple(): return
        self._set_descr(od.simplify(self))

    def _exact_field(self):
        """
        Returns a generator for a number field that includes this number
        (not necessarily the smallest such number field).

        EXAMPLES::

            sage: QQbar(2)._exact_field()
            Trivial generator
            sage: (sqrt(QQbar(2)) + sqrt(QQbar(19)))._exact_field()
            Number Field in a with defining polynomial y^4 - 20*y^2 + 81 with a in 2.375100220297941?
            sage: (QQbar(7)^(3/5))._exact_field()
            Number Field in a with defining polynomial y^5 - 2*y^4 - 18*y^3 + 38*y^2 + 82*y - 181 with a in 2.554256611698490?
        """

        sd = self._descr
        if isinstance(sd, (ANRational, ANExtensionElement)):
            return sd.generator()
        self.exactify()
        return self._exact_field()

    def _exact_value(self):
        r"""
        Returns an ``ANRational`` or an ``ANExtensionElement`` representing this
        value.

        EXAMPLES::

            sage: QQbar(2)._exact_value()
            2
            sage: (sqrt(QQbar(2)) + sqrt(QQbar(19)))._exact_value()
            -1/9*a^3 - a^2 + 11/9*a + 10 where a^4 - 20*a^2 + 81 = 0 and a in 2.375100220297941?
            sage: (QQbar(7)^(3/5))._exact_value()
            2*a^4 + 2*a^3 - 34*a^2 - 17*a + 150 where a^5 - 2*a^4 - 18*a^3 + 38*a^2 + 82*a - 181 = 0 and a in 2.554256611698490?
        """
        sd = self._descr
        if isinstance(sd, (ANRational, ANExtensionElement)):
            return sd
        self.exactify()
        return self._descr

    def _more_precision(self):
        """
        Recompute the interval bounding this number with higher-precision
        interval arithmetic.

        EXAMPLES::

            sage: rt2 = sqrt(QQbar(2))
            sage: rt2._value
            1.4142135623730950488?
            sage: rt2._more_precision()
            sage: rt2._value
            1.41421356237309504880168872420969807857?
            sage: rt2._more_precision()
            sage: rt2._value
            1.41421356237309504880168872420969807856967187537694807317667973799073247846211?
        """
        prec = self._value.prec()
        self._value = self._descr._interval_fast(prec*2)

    def minpoly(self):
        """
        Compute the minimal polynomial of this algebraic number.
        The minimal polynomial is the monic polynomial of least degree
        having this number as a root; it is unique.

        EXAMPLES::

            sage: QQbar(4).sqrt().minpoly()
            x - 2
            sage: ((QQbar(2).nth_root(4))^2).minpoly()
            x^2 - 2
            sage: v = sqrt(QQbar(2)) + sqrt(QQbar(3)); v
            3.146264369941973?
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

        EXAMPLES::

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
        r"""
        Given a :class:`RealIntervalField` or
        :class:`ComplexIntervalField`, compute the value of this number
        using interval arithmetic of at least the precision of the field,
        and return the value in that field. (More precision may be used
        in the computation.)  The returned interval may be arbitrarily
        imprecise, if this number is the result of a sufficiently long
        computation chain.

        EXAMPLES::

            sage: x = AA(2).sqrt()
            sage: x.interval_fast(RIF)
            1.414213562373095?
            sage: x.interval_fast(RealIntervalField(200))
            1.414213562373095048801688724209698078569671875376948073176680?
            sage: x = QQbar(I).sqrt()
            sage: x.interval_fast(CIF)
            0.7071067811865475? + 0.7071067811865475?*I
            sage: x.interval_fast(RIF)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 0.7071067811865475244? + 0.7071067811865475244?*I to real interval
        """
        while self._value.prec() < field.prec():
            self._more_precision()
        return field(self._value)

    def interval_diameter(self, diam):
        """
        Compute an interval representation of self with ``diameter()`` at
        most ``diam``. The precision of the returned value is unpredictable.

        EXAMPLES::

            sage: AA(2).sqrt().interval_diameter(1e-10)
            1.4142135623730950488?
            sage: AA(2).sqrt().interval_diameter(1e-30)
            1.41421356237309504880168872420969807857?
            sage: QQbar(2).sqrt().interval_diameter(1e-10)
            1.4142135623730950488?
            sage: QQbar(2).sqrt().interval_diameter(1e-30)
            1.41421356237309504880168872420969807857?
        """
        if diam <= 0:
            raise ValueError('diameter must be positive in interval_diameter')

        while self._value.diameter() > diam:
            self._more_precision()

        return self._value

    def interval(self, field):
        r"""
        Given an interval field (real or complex, as appropriate) of
        precision `p`, compute an interval representation of self with
        ``diameter()`` at most `2^{-p}`; then round that representation into
        the given field. Here ``diameter()`` is relative diameter for
        intervals not containing 0, and absolute diameter for
        intervals that do contain 0; thus, if the returned interval
        does not contain 0, it has at least `p-1` good bits.

        EXAMPLES::

            sage: RIF64 = RealIntervalField(64)
            sage: x = AA(2).sqrt()
            sage: y = x*x
            sage: y = 1000 * y - 999 * y
            sage: y.interval_fast(RIF64)
            2.000000000000000?
            sage: y.interval(RIF64)
            2.000000000000000000?
            sage: CIF64 = ComplexIntervalField(64)
            sage: x = QQbar.zeta(11)
            sage: x.interval_fast(CIF64)
            0.8412535328311811689? + 0.5406408174555975821?*I
            sage: x.interval(CIF64)
            0.8412535328311811689? + 0.5406408174555975822?*I

        The following implicitly use this method::

            sage: RIF(AA(5).sqrt())
            2.236067977499790?
            sage: AA(-5).sqrt().interval(RIF)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 2.236067977499790?*I to real interval

        TESTS:

        Check that :trac:`20209` is fixed::

            sage: RIF(QQbar(2).sqrt())
            1.414213562373095?
            sage: RIF(QQbar.gen() + QQbar(2).sqrt() - QQbar.gen())
            1.414213562373095?
            sage: RIF((QQbar.gen() + QQbar(2).sqrt() - QQbar.gen()).sqrt())
            1.189207115002722?

            sage: RealIntervalField(129)(QQbar(3).sqrt())
            1.73205080756887729352744634150587236695?
            sage: RIF(QQbar.gen())
            Traceback (most recent call last):
            ...
            TypeError: unable to convert I to real interval
        """
        target = RR(1.0) >> field.prec()
        val = self.interval_diameter(target)
        if isinstance(field, RealIntervalField_class) and is_ComplexIntervalFieldElement(val):
            if val.imag().is_zero():
                return field(val.real())
            elif self.imag().is_zero():
                return field(self.real())
            else:
                raise TypeError("unable to convert {} to real interval".format(self))
        else:
            return field(val)

    _real_mpfi_ = interval

    def radical_expression(self):
        r"""
        Attempt to obtain a symbolic expression using radicals. If no
        exact symbolic expression can be found, the algebraic number
        will be returned without modification.

        EXAMPLES::

            sage: AA(1/sqrt(5)).radical_expression()
            1/5*sqrt(5)
            sage: AA(sqrt(5 + sqrt(5))).radical_expression()
            sqrt(sqrt(5) + 5)
            sage: QQbar.zeta(5).radical_expression()
            1/4*sqrt(5) + 1/2*sqrt(-1/2*sqrt(5) - 5/2) - 1/4
            sage: a = QQ[x](x^7 - x - 1).roots(AA, False)[0]
            sage: a.radical_expression()
            1.112775684278706?
            sage: a.radical_expression().parent() == SR
            False
            sage: a = sorted(QQ[x](x^7-x-1).roots(QQbar, False), key=imag)[0]
            sage: a.radical_expression()
            -0.3636235193291805? - 0.9525611952610331?*I
            sage: QQbar.zeta(5).imag().radical_expression()
            1/2*sqrt(1/2*sqrt(5) + 5/2)
            sage: AA(5/3).radical_expression()
            5/3
            sage: AA(5/3).radical_expression().parent() == SR
            True
            sage: QQbar(0).radical_expression()
            0

        TESTS:

        In this example we find the correct answer despite the fact that
        multiple roots overlap with the current value. As a consequence,
        the precision of the evaluation will have to be increased.

        ::

            sage: a = AA(sqrt(2) + 10^25)
            sage: p = a.minpoly()
            sage: v = a._value
            sage: f = ComplexIntervalField(v.prec())
            sage: [f(b.rhs()).overlaps(f(v)) for b in SR(p).solve(x)]
            [True, True]
            sage: a.radical_expression()
            sqrt(2) + 10000000000000000000000000
        """
        from sage.symbolic.ring import SR # Lazy to avoid cyclic dependency

        # Adapted from NumberFieldElement._symbolic_()
        poly = self.minpoly()
        var = SR(poly.variable_name())
        if is_ComplexIntervalFieldElement(self._value):
            interval_field = self._value.parent()
        else:
            interval_field = ComplexIntervalField(self._value.prec())
        roots = poly.roots(SR, multiplicities=False)
        if len(roots) != poly.degree():
            return self
        while True:
            candidates = []
            for root in roots:
                if interval_field(root).overlaps(interval_field(self._value)):
                    candidates.append(root)
            if len(candidates) == 1:
                return candidates[0]
            roots = candidates
            interval_field = interval_field.to_prec(interval_field.prec()*2)

class AlgebraicNumber(AlgebraicNumber_base):
    r"""
    The class for algebraic numbers (complex numbers which are the roots
    of a polynomial with integer coefficients). Much of its functionality
    is inherited from ``AlgebraicNumber_base``.
    """
    def __init__(self, x):
        r"""
        Initialize this AlgebraicNumber object.

        EXAMPLE::

            sage: t = QQbar.zeta(5)
            sage: type(t)
            <class 'sage.rings.qqbar.AlgebraicNumber'>
        """
        AlgebraicNumber_base.__init__(self, QQbar, x)

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = QQbar.zeta(5)
            sage: loads(dumps(t)) == t
            True
        """
        return (AlgebraicNumber, (self._descr, ))

    def __cmp__(self, other):
        """
        Compare two algebraic numbers, lexicographically. (That is,
        first compare the real components; if the real components are
        equal, compare the imaginary components.)

        EXAMPLES::

            sage: x = QQbar.zeta(3); x
            -0.500000000000000? + 0.866025403784439?*I
            sage: cmp(QQbar(-1), x)
            -1
            sage: cmp(QQbar(-1/2), x)
            -1
            sage: cmp(QQbar(0), x)
            1

        One problem with this lexicographic ordering is the fact that if
        two algebraic numbers have the same real component, that real
        component has to be compared for exact equality, which can be
        a costly operation.  For the special case where both numbers
        have the same minimal polynomial, that cost can be avoided,
        though (see :trac:`16964`)::

            sage: x = polygen(ZZ)
            sage: p = 69721504*x^8 + 251777664*x^6 + 329532012*x^4 + 184429548*x^2 + 37344321
            sage: sorted(p.roots(QQbar,False))
            [-0.0221204634374361? - 1.090991904211621?*I,
             -0.0221204634374361? + 1.090991904211621?*I,
             -0.8088604911480535?*I,
             0.?e-215 - 0.7598602580415435?*I,
             0.?e-229 + 0.7598602580415435?*I,
             0.8088604911480535?*I,
             0.0221204634374361? - 1.090991904211621?*I,
             0.0221204634374361? + 1.090991904211621?*I]

        It also works for comparison of conjugate roots even in a degenerate
        situation where many roots have the same real part. In the following
        example, the polynomial ``p2`` is irreducible and all its roots have
        real part equal to `1`::

            sage: p1 = x^8 + 74*x^7 + 2300*x^6 + 38928*x^5 + \
            ....: 388193*x^4 + 2295312*x^3 + 7613898*x^2 + \
            ....: 12066806*x + 5477001
            sage: p2 = p1((x-1)^2)
            sage: sum(1 for r in p2.roots(CC,False) if abs(r.real() - 1) < 0.0001)
            16
            sage: r1 = QQbar.polynomial_root(p2, CIF(1, (-4.1,-4.0)))
            sage: r2 = QQbar.polynomial_root(p2, CIF(1, (4.0, 4.1)))
            sage: cmp(r1,r2), cmp(r1,r1), cmp(r2,r2), cmp(r2,r1)
            (-1, 0, 0, 1)

        Though, comparing roots which are not equal or conjugate is much
        slower because the algorithm needs to check the equality of the real
        parts::

            sage: sorted(p2.roots(QQbar,False))   # long time - 3 secs
            [1.000000000000000? - 4.016778562562223?*I,
             1.000000000000000? - 3.850538755978243?*I,
             1.000000000000000? - 3.390564396412898?*I,
             ...
             1.000000000000000? + 3.390564396412898?*I,
             1.000000000000000? + 3.850538755978243?*I,
             1.000000000000000? + 4.016778562562223?*I]
        """
        # case 0: same object
        if self is other: return 0

        # case 1: real parts are clearly distinct
        ri1 = self._value.real()
        ri2 = other._value.real()
        if not ri1.overlaps(ri2):
            return cmp(ri1, ri2)

        # case 2: possibly equal or conjugate values
        # (this case happen a lot when sorting the roots of a real polynomial)
        if is_RealIntervalFieldElement(self._value):
            ci1 = ri1.parent().zero()
        else:
            ci1 = self._value.imag().abs()
        if is_RealIntervalFieldElement(other._value):
            ci2 = ri2.parent().zero()
        else:
            ci2 = other._value.imag().abs()
        if ci1.overlaps(ci2) and self.minpoly() == other.minpoly():
            ri = ri1.union(ri2)
            ci = ci1.union(ci2)
            roots = self.minpoly().roots(QQbar, False)
            roots = [r for r in roots if r._value.real().overlaps(ri)
                     and r._value.imag().abs().overlaps(ci)]
            if len(roots) == 1:
                # There is only a single (real) root matching both descriptors
                # so they both must be that root and therefore equal.
                return 0
            if (len(roots) == 2 and
                not roots[0]._value.imag().contains_zero()):
                # There is a complex conjugate pair of roots matching both
                # descriptors, so compare by imaginary value.
                ii1 = self._value.imag()
                while ii1.contains_zero():
                    self._more_precision()
                    ii1 = self._value.imag()
                ii2 = other._value.imag()
                while ii2.contains_zero():
                    other._more_precision()
                    ii2 = other._value.imag()
                if ii1.overlaps(ii2):
                    return 0
                return cmp(ii1, ii2)

        # case 3: try hard to compare real parts and imaginary parts
        rcmp = cmp(self.real(), other.real())
        if rcmp != 0:
            return rcmp
        return cmp(self.imag(), other.imag())

    def __eq__(self, other):
        """
        Test two algebraic numbers for equality.

        EXAMPLES::

            sage: QQbar.zeta(6) == QQbar(1/2 + I*sqrt(3)/2)
            True
            sage: QQbar(I) == QQbar(I * (2^100+1)/(2^100))
            False
            sage: QQbar(2) == 2
            True
            sage: QQbar(2) == GF(7)(2)
            False
            sage: GF(7)(2) in QQbar
            False
        """
        if not isinstance(other, AlgebraicNumber):
            try:
                self, other = canonical_coercion(self, other)
                return self == other
            except TypeError:
                return False
        if self is other: return True
        if isinstance(other._descr, ANRational) and other._descr._value.is_zero():
            return not self
        if isinstance(self._descr, ANRational) and self._descr._value.is_zero():
            return not other
        return not self._sub_(other)

    def __ne__(self, other):
        r"""
        Test two algebraic numbers for inequality.

        EXAMPLES::

            sage: QQbar.zeta(6) != QQbar(1/2 + I*sqrt(3)/2)
            False
            sage: QQbar(I) != QQbar(I * (2^100+1)/(2^100))
            True
            sage: QQbar(2) != 2
            False
            sage: QQbar(2) != GF(7)(2)
            True
        """
        return not self == other

    def __nonzero__(self):
        """
        Check whether self is equal is nonzero. This is fast if
        interval arithmetic proves that self is nonzero, but may be
        slow if the number actually is very close to zero.

        EXAMPLES::

            sage: (QQbar.zeta(2) + 1).__nonzero__()
            False
            sage: (QQbar.zeta(7) / (2^500)).__nonzero__()
            True
        """
        val = self._value
        d = self._descr
        if not val.contains_zero() or isinstance(d, ANExtensionElement):
            return True
        elif isinstance(d, ANRational):
            return bool(d._value)

        while self._value.prec() < 128:
            self._more_precision()
            if not self._value.contains_zero():
                return True

        # Sigh...
        self.exactify()
        return self.__nonzero__()

    def __pow__(self, e):
        r""" ``self**p`` returns the `p`'th power of self (where `p` can
        be an arbitrary rational). If `p` is `(a/b)`, takes the principal
        `b`'th root of self, then takes that to the `a`'th power. (Note
        that this differs from ``__pow__`` on algebraic reals, where real
        roots are preferred over principal roots if they exist.)

        EXAMPLES::

            sage: QQbar(2)^(1/2)
            1.414213562373095?
            sage: QQbar(8)^(2/3)
            4
            sage: QQbar(8)^(2/3) == 4
            True
            sage: x = polygen(QQbar)
            sage: phi = QQbar.polynomial_root(x^2 - x - 1, RIF(1, 2))
            sage: tau = QQbar.polynomial_root(x^2 - x - 1, RIF(-1, 0))
            sage: rt5 = QQbar(5)^(1/2)
            sage: phi^10 / rt5
            55.00363612324742?
            sage: tau^10 / rt5
            0.003636123247413266?
            sage: (phi^10 - tau^10) / rt5
            55.00000000000000?
            sage: (phi^10 - tau^10) / rt5 == fibonacci(10)
            True
            sage: (phi^50 - tau^50) / rt5 == fibonacci(50)
            True
            sage: QQbar(-8)^(1/3)
            1.000000000000000? + 1.732050807568878?*I
            sage: (QQbar(-8)^(1/3))^3
            -8
            sage: QQbar(32)^(1/5)
            2
            sage: a = QQbar.zeta(7)^(1/3); a
            0.9555728057861407? + 0.2947551744109043?*I
            sage: a == QQbar.zeta(21)
            True
            sage: QQbar.zeta(7)^6
            0.6234898018587335? - 0.7818314824680299?*I
            sage: (QQbar.zeta(7)^6)^(1/3) * QQbar.zeta(21)
            1.000000000000000? + 0.?e-18*I
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
                    z = QQbar.zeta(2*d)**n
                    return z * AlgebraicNumber(ANRational(rt**n))
                else:
                    return AlgebraicNumber(ANRational(rt**n))

        # Without this special case, we do not know the multiplicity
        # of the desired root
        if self.is_zero():
            return QQbar.zero()
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

    def _mpfr_(self, field):
        r"""
        Given a ``RealField``, compute a good approximation to self in
        that field. Works only if the imaginary component of self is
        exactly zero; otherwise it raises a ``ValueError``.

        EXAMPLES::

            sage: QQbar(sqrt(2))._mpfr_(RR)
            1.41421356237309
            sage: QQbar(-22/7)._mpfr_(RR)
            -3.14285714285714
            sage: QQbar.zeta(3)._mpfr_(RR)
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
        """
        return AA(self)._mpfr_(field)

    def __float__(self):
        r"""
        Compute a good float approximation to self. Works only if the
        imaginary component of self is exactly zero; otherwise it
        raises a ``ValueError``.

        EXAMPLES::

            sage: QQbar(sqrt(2)).__float__()
            1.414213562373095
            sage: float(QQbar(-22/7))
            -3.1428571428571432
            sage: float(QQbar.zeta(3))
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
        """
        return AA(self).__float__()

    def __complex__(self):
        r"""
        Compute a good complex approximation to self.

        EXAMPLES::

            sage: QQbar(sqrt(2)).__complex__()
            (1.414213562373095+0j)
            sage: complex(QQbar.zeta(3))
            (-0.5+0.8660254037844386j)
        """
        return CC(self).__complex__()

    def _complex_double_(self, cdf):
        r"""
        Compute a good approximation to self in CDF.

        EXAMPLES::

            sage: QQbar(sqrt(-5))._complex_double_(CDF)
            2.23606797749979*I
            sage: CDF(QQbar.zeta(12))
            0.8660254037844386 + 0.5*I
        """
        return cdf(CC(self))

    def _interval_fast(self, prec):
        r"""
        Shortcut for :meth:`AlgebraicNumber_base.interval_fast` which uses the complex interval field.

        EXAMPLE::

            sage: QQbar(sqrt(-5))._interval_fast(100)
            2.236067977499789696409173...?*I
        """
        return self.interval_fast(ComplexIntervalField(prec))

    def _integer_(self, ZZ=None):
        """
        Return self as an Integer.

        EXAMPLES::

            sage: QQbar(0)._integer_()
            0
            sage: QQbar(0)._integer_().parent()
            Integer Ring
            sage: QQbar.zeta(6)._integer_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
            sage: QQbar(sqrt(17))._integer_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce non-integral Algebraic Real 4.123105625617660? to Integer
            sage: QQbar(sqrt(16))._integer_()
            4
            sage: v = QQbar(1 + I*sqrt(3))^5 + QQbar(16*sqrt(3)*I); v
            16.00000000000000? + 0.?e-17*I
            sage: v._integer_()
            16
        """
        return AA(self)._integer_(ZZ)

    def _rational_(self):
        """
        Return self as a Rational.

        EXAMPLES::

            sage: QQbar(-22/7)._rational_()
            -22/7
            sage: QQbar(3)._rational_().parent()
            Rational Field
            sage: (QQbar.zeta(7)^3)._rational_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary part to algebraic real
            sage: QQbar(sqrt(2))._rational_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce irrational Algebraic Real 1.414213562373095? to Rational
            sage: v1 = QQbar(1/3 + I*sqrt(5))^7
            sage: v2 = QQbar((100336/729*golden_ratio - 50168/729)*I)
            sage: v = v1 + v2; v
            -259.6909007773206? + 0.?e-15*I
            sage: v._rational_()
            -567944/2187
        """
        return AA(self)._rational_()

    def real(self):
        r"""
        Return the real part of self.

        EXAMPLE::

            sage: QQbar.zeta(5).real()
            0.3090169943749474?
        """
        return AlgebraicReal(self._descr.real(self))

    def imag(self):
        r"""
        Return the imaginary part of self.

        EXAMPLE::

            sage: QQbar.zeta(7).imag()
            0.7818314824680299?
        """
        return AlgebraicReal(self._descr.imag(self))

    def conjugate(self):
        """
        Returns the complex conjugate of self.

        EXAMPLES::

            sage: QQbar(3 + 4*I).conjugate()
            3 - 4*I
            sage: QQbar.zeta(7).conjugate()
            0.6234898018587335? - 0.7818314824680299?*I
            sage: QQbar.zeta(7) + QQbar.zeta(7).conjugate()
            1.246979603717467? + 0.?e-18*I
        """
        return AlgebraicNumber(self._descr.conjugate(self))

    def norm(self):
        r"""
        Returns ``self * self.conjugate()``. This is the algebraic
        definition of norm, if we view ``QQbar`` as ``AA[I]``.

        EXAMPLES::

            sage: QQbar(3 + 4*I).norm()
            25
            sage: type(QQbar(I).norm())
            <class 'sage.rings.qqbar.AlgebraicReal'>
            sage: QQbar.zeta(1007).norm()
            1.000000000000000?
        """
        return AlgebraicReal(self._descr.norm(self))

    def interval_exact(self, field):
        r"""
        Given a ``ComplexIntervalField``, compute the best possible
        approximation of this number in that field. Note that if
        either the real or imaginary parts of this number are
        sufficiently close to some floating-point number (and, in
        particular, if either is exactly representable in floating-point),
        then this will trigger exact computation, which may be very slow.

        EXAMPLES::

            sage: a = QQbar(I).sqrt(); a
            0.7071067811865475? + 0.7071067811865475?*I
            sage: a.interval_exact(CIF)
            0.7071067811865475? + 0.7071067811865475?*I
            sage: b = QQbar((1+I)*sqrt(2)/2)
            sage: (a - b).interval(CIF)
            0.?e-19 + 0.?e-18*I
            sage: (a - b).interval_exact(CIF)
            0
        """
        if not is_ComplexIntervalField(field):
            raise ValueError("AlgebraicNumber interval_exact requires a ComplexIntervalField")
        rfld = field._real_field()
        re = self.real().interval_exact(rfld)
        im = self.imag().interval_exact(rfld)
        return field(re, im)

    def _complex_mpfr_field_(self, field):
        r"""
        Compute an approximation to self in the given field, which may be
        either an interval field (in which case ``self.interval()`` is called)
        or any other complex field (in which case ``self.complex_number()`` is
        called).

        EXAMPLE::

            sage: a = QQbar(1 + I).sqrt()
            sage: t = a._complex_mpfr_field_(CIF); t
            1.098684113467810? + 0.4550898605622274?*I
            sage: parent(t)
            Complex Interval Field with 53 bits of precision
            sage: t = a._complex_mpfr_field_(ComplexField(100)); t
            1.0986841134678099660398011952 + 0.45508986056222734130435775782*I
            sage: parent(t)
            Complex Field with 100 bits of precision
        """
        if is_ComplexIntervalField(field):
            return self.interval(field)
        else:
            return self.complex_number(field)

    def complex_number(self, field):
        r"""
        Given the complex field ``field`` compute an accurate approximation of
        this element in that field.

        The approximation will be off by at most two ulp's in each component,
        except for components which are very close to zero, which will have an
        absolute error at most `2^{-prec+1}` where `prec` is the precision of
        the field.

        EXAMPLES::

            sage: a = QQbar.zeta(5)
            sage: a.complex_number(CC)
            0.309016994374947 + 0.951056516295154*I

            sage: b = QQbar(2).sqrt() + QQbar(3).sqrt() * QQbar.gen()
            sage: b.complex_number(ComplexField(128))
            1.4142135623730950488016887242096980786 + 1.7320508075688772935274463415058723669*I
        """
        v = self.interval(ComplexIntervalField(field.prec()))
        return field(v)

    def complex_exact(self, field):
        r"""
        Given a ``ComplexField``, return the best possible approximation of
        this number in that field. Note that if either component is
        sufficiently close to the halfway point between two floating-point
        numbers in the corresponding ``RealField``, then this will trigger
        exact computation, which may be very slow.

        EXAMPLES::

            sage: a = QQbar.zeta(9) + QQbar(I) + QQbar.zeta(9).conjugate(); a
            1.532088886237957? + 1.000000000000000?*I
            sage: a.complex_exact(CIF)
            1.532088886237957? + 1*I
        """
        rfld = field._real_field()
        re = self.real().real_exact(rfld)
        im = self.imag().real_exact(rfld)
        return field(re, im)

    def multiplicative_order(self):
        r"""
        Compute the multiplicative order of this algebraic real
        number. That is, find the smallest positive integer `n` such
        that `x^n = 1`. If there is no such `n`, returns ``+Infinity``.

        We first check that ``abs(x)`` is very close to 1. If so, we compute
        `x` exactly and examine its argument.

        EXAMPLES::

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
        d = self.minpoly().is_cyclotomic(True)
        return d if d else infinity.infinity

    def rational_argument(self):
        r"""
        Returns the argument of self, divided by `2\pi`, as long as this
        result is rational. Otherwise returns None. Always triggers
        exact computation.

        EXAMPLES::

            sage: QQbar((1+I)*(sqrt(2)+sqrt(5))).rational_argument()
            1/8
            sage: QQbar(-1 + I*sqrt(3)).rational_argument()
            1/3
            sage: QQbar(-1 - I*sqrt(3)).rational_argument()
            -1/3
            sage: QQbar(3+4*I).rational_argument() is None
            True
            sage: (QQbar(2)**(1/5) * QQbar.zeta(7)**2).rational_argument()
            2/7
            sage: (QQbar.zeta(73)**5).rational_argument()
            5/73
            sage: (QQbar.zeta(3)^65536).rational_argument()
            1/3
        """
        # This always triggers exact computation. An alternate method
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
        """
        Create an algebraic real from x, possibly taking the real part of x.

        TESTS:

        Both of the following examples, from :trac:`11728`, trigger
        taking the real part below. This is necessary because
        sometimes a very small (e.g., 1e-17) complex part appears in a
        complex interval used to create an AlgebraicReal.::

            sage: a = QQbar((-1)^(1/4)); b = AA(a^3-a); t = b.as_number_field_element()
            sage: b*1
            -1.414213562373095?
        """
        AlgebraicNumber_base.__init__(self, AA, x)
        self._ensure_real()

    def _ensure_real(self):
        """
        This is used internally by some methods to check if
        self._value is a complex interval, and if so, take the real
        part.

        EXAMPLES::

            sage: a = QQbar((-1)^(1/4)); b = AA(a^3-a); b._value
            -1.4142135623730950488?
            sage: b._value = a._value; b._value
            0.7071067811865475244? + 0.7071067811865475244?*I
            sage: b._ensure_real()
            sage: b._value
            0.7071067811865475244?
            sage: type(b._value)
            <type 'sage.rings.real_mpfi.RealIntervalFieldElement'>
        """
        if is_ComplexIntervalFieldElement(self._value):
            self._value = self._value.real()

    def _more_precision(self):
        """
        Recompute the interval bounding this number with higher-precision
        interval arithmetic.

        EXAMPLE::

            sage: a = QQbar(sqrt(2))
            sage: a._more_precision()

        TESTS:

        We have to ensure after doing this that self._value is still
        real which isn't the case without calling _ensure_real (see
        :trac:`11728`)::

            sage: P = AA[x](1+x^4); a1,a2 = P.factor()[0][0],P.factor()[1][0]; a1*a2
            x^4 + 1.000000000000000?
            sage: a1,a2
            (x^2 - 1.414213562373095?*x + 1, x^2 + 1.414213562373095?*x + 1)
            sage: a1*a2
            x^4 + 1
        """
        AlgebraicNumber_base._more_precision(self)
        self._ensure_real()

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = AA(sqrt(2))
            sage: loads(dumps(t)) == t
            True
        """
        return (AlgebraicReal, (self._descr, ))

    def __cmp__(self, other):
        """
        Compare two algebraic reals.

        EXAMPLES::

            sage: cmp(AA(golden_ratio), AA(sqrt(5)))
            -1
            sage: cmp(AA(golden_ratio), AA((sqrt(5)+1)/2))
            0
            sage: cmp(AA(7), AA(50/7))
            -1
        """
        if self is other: return 0
        if isinstance(other._descr, ANRational) and other._descr._value.is_zero():
            return self.sign()
        elif isinstance(self._descr, ANRational) and self._descr._value.is_zero():
            return -other.sign()
        else:
            return self._sub_(other).sign()

    def __pow__(self, e):
        """
        ``self**p`` returns the `p`'th power of self (where `p` can be an
        arbitrary rational). If `p` is `(a/b)`, takes the `b`'th root of
        self, then takes that to the `a`'th power. If self is negative
        and `b` is odd, it takes the real `b`'th root; if self is odd and
        `b` is even, this takes a complex root. Note that the behavior
        when self is negative and `b` is odd differs from the complex
        case; algebraic numbers select the principal complex `b`'th
        root, but algebraic reals select the real root.

        EXAMPLES::

            sage: AA(2)^(1/2)
            1.414213562373095?
            sage: AA(8)^(2/3)
            4
            sage: AA(8)^(2/3) == 4
            True
            sage: x = polygen(AA)
            sage: phi = AA.polynomial_root(x^2 - x - 1, RIF(0, 2))
            sage: tau = AA.polynomial_root(x^2 - x - 1, RIF(-2, 0))
            sage: rt5 = AA(5)^(1/2)
            sage: phi^10 / rt5
            55.00363612324742?
            sage: tau^10 / rt5
            0.003636123247413266?
            sage: (phi^10 - tau^10) / rt5
            55.00000000000000?
            sage: (phi^10 - tau^10) / rt5 == fibonacci(10)
            True
            sage: (phi^50 - tau^50) / rt5 == fibonacci(50)
            True

        TESTS::

            sage: AA(-8)^(1/3)
            -2
            sage: AA(-8)^(2/3)
            4
            sage: AA(32)^(3/5)
            8
            sage: AA(-16)^(1/2)
            4*I
            sage: AA(-16)^(1/4)
            1.414213562373095? + 1.414213562373095?*I
            sage: AA(-16)^(1/4)/QQbar.zeta(8)
            2

        We check that :trac:`7859` is fixed::

            sage: (AA(2)^(1/2)-AA(2)^(1/2))^(1/2)
            0
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
                        z = QQbar.zeta(2*d)**n
                        return z * AlgebraicNumber(ANRational(rt**n))
                    else:
                        return AlgebraicReal(ANRational((-rt)**n))
                else:
                    return AlgebraicReal(ANRational(rt**n))

        # Without this special case, we do not know the multiplicity
        # of the desired root
        if self.sign() == 0:
            return AA.zero()
        if d % 2 == 0:
            if self.sign() < 0:
                return QQbar(self) ** e
        pow_n = self**n
        poly = AAPoly.gen()**d - pow_n
        range = pow_n.interval_fast(RIF)
        if d % 2 == 0:
            result_min = 0
        else:
            result_min = min(range.lower(), -1)
        result_max = max(range.upper(), 1)
        return AlgebraicReal(ANRoot(poly, RIF(result_min, result_max)))

    def _integer_(self, Z=None):
        """
        Return self as an Integer.

        EXAMPLES::

            sage: AA(42)._integer_()
            42
            sage: AA(42)._integer_().parent()
            Integer Ring
            sage: AA(golden_ratio)._integer_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce non-integral Algebraic Real 1.618033988749895? to Integer
            sage: (AA(golden_ratio)^10 + AA(1-golden_ratio)^10)._integer_()
            123
            sage: AA(-22/7)._integer_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce non-integral Algebraic Real -22/7 to Integer
        """
        if self._value.lower().ceiling() > self._value.upper().floor():
            # The value is known to be non-integral.
            raise ValueError("Cannot coerce non-integral Algebraic Real %s to Integer" % self)

        self.exactify()
        if not isinstance(self._descr, ANRational):
            raise ValueError("Cannot coerce irrational Algebraic Real %s to Integer" % self)

        return ZZ(self._descr._value)

    def _floor_ceil(self, method):
        r"""
        Helper method used by :meth:`floor()`, :meth:`ceil()`,
        :meth:`round()`, and :meth:`trunc()`.

        TESTS::

            sage: x = polygen(QQ)
            sage: a = AA.polynomial_root(x^5 - (1-2^(-80)), RIF((0,2)))
            sage: b = AA.polynomial_root(x^5 - (1+2^(-80)), RIF((0,2)))
            sage: two = (a+b)^5 - 5*(a^4*b+a*b^4) - 10*(a^3*b^2+a^2*b^3)
            sage: one_half = 1/two
            sage: [[z.floor(), z.ceil(), z.round(), z.trunc()] # indirect doctest
            ....:  for z in [a, -a, b, -b, 6*(a+two),
            ....:            AA(0), AA(1), AA(-1), AA(1/2), AA(-1/2)]]
            [[0, 1, 1, 0], [-1, 0, -1, 0], [1, 2, 1, 1], [-2, -1, -1, -1],
            [17, 18, 18, 17], [0, 0, 0, 0], [1, 1, 1, 1], [-1, -1, -1, -1],
            [0, 1, 1, 0], [-1, 0, -1, 0]]
            sage: [[z.floor(), z.ceil(), z.trunc()] for z in [two, a*b]] # long time
            [[2, 2, 2], [0, 1, 0]]
            sage: [one_half.round(), (-one_half).round()] # long time
            [1, -1]
        """
        for i in itertools.count():
            candidate = method(self._value.lower())
            if candidate == method(self._value.upper()):
                return candidate
            self._more_precision()
            # field elements are irrational by construction
            if i == 2 and not isinstance(self._descr, ANExtensionElement):
                try:
                    return method(self._rational_())
                except (ValueError, TypeError):
                    pass

    def floor(self):
        r"""
        Return the largest integer not greater than ``self``.

        EXAMPLES::

            sage: AA(sqrt(2)).floor()
            1
            sage: AA(-sqrt(2)).floor()
            -2
            sage: AA(42).floor()
            42

        TESTS:

        Check that :trac:`15501` is fixed::

            sage: a = QQbar((-1)^(1/4)).real()
            sage: (floor(a-a) + a).parent()
            Algebraic Real Field
        """
        return self._floor_ceil(lambda x: x.floor())

    def ceil(self):
        r"""
        Return the smallest integer not smaller than ``self``.

        EXAMPLES::

            sage: AA(sqrt(2)).ceil()
            2
            sage: AA(-sqrt(2)).ceil()
            -1
            sage: AA(42).ceil()
            42
        """
        return self._floor_ceil(lambda x: x.ceil())

    def round(self):
        r"""
        Round ``self`` to the nearest integer.

        EXAMPLES::

            sage: AA(sqrt(2)).round()
            1
            sage: AA(1/2).round()
            1
            sage: AA(-1/2).round()
            -1
        """
        return self._floor_ceil(lambda x: x.round())

    def trunc(self):
        r"""
        Round ``self`` to the nearest integer toward zero.

        EXAMPLES::

            sage: AA(sqrt(2)).trunc()
            1
            sage: AA(-sqrt(2)).trunc()
            -1
            sage: AA(1).trunc()
            1
            sage: AA(-1).trunc()
            -1
        """
        return self._floor_ceil(lambda x: x.trunc())

    def _rational_(self):
        """
        Return self as a Rational.

        EXAMPLES::

            sage: AA(42)._rational_().parent()
            Rational Field
            sage: AA(-22/7)._rational_()
            -22/7
            sage: AA(sqrt(7))._rational_()
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce irrational Algebraic Real 2.645751311064591? to Rational
            sage: v = AA(1/2 + sqrt(2))^3 - AA(11/4*sqrt(2)); v
            3.125000000000000?
            sage: v._rational_()
            25/8
        """
        self.exactify()
        if not isinstance(self._descr, ANRational):
            raise ValueError("Cannot coerce irrational Algebraic Real %s to Rational" % self)

        return QQ(self._descr._value)

    def real(self):
        """
        Returns the real part of this algebraic real (so it always returns
        self).

        EXAMPLES::

            sage: a = AA(sqrt(2) + sqrt(3))
            sage: a.real()
            3.146264369941973?
            sage: a.real() is a
            True
        """
        return self

    def imag(self):
        """
        Returns the imaginary part of this algebraic real (so it always
        returns 0).

        EXAMPLES::

            sage: a = AA(sqrt(2) + sqrt(3))
            sage: a.imag()
            0
            sage: parent(a.imag())
            Algebraic Real Field
        """
        return AA_0

    def conjugate(self):
        """
        Returns the complex conjugate of self, i.e. returns itself.

        EXAMPLES::

            sage: a = AA(sqrt(2) + sqrt(3))
            sage: a.conjugate()
            3.146264369941973?
            sage: a.conjugate() is a
            True
        """
        return self

    def sign(self):
        """
        Compute the sign of this algebraic number (return -1 if negative,
        0 if zero, or 1 if positive).

        Computes an interval enclosing this number using 128-bit interval
        arithmetic; if this interval includes 0, then fall back to
        exact computation (which can be very slow).

        EXAMPLES::

            sage: AA(-5).nth_root(7).sign()
            -1
            sage: (AA(2).sqrt() - AA(2).sqrt()).sign()
            0
        """
        if self._value.lower() > 0:
            return 1
        elif self._value.upper() < 0:
            return -1
        elif isinstance(self._descr, ANRational):
            return self._descr._value.sign()
        elif isinstance(self._descr, ANExtensionElement):
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
        r"""
        Compute an approximation to this ``AlgebraicReal`` object in a real interval field of precision prec.

        EXAMPLE::

            sage: t = AA(sqrt(7))
            sage: t._interval_fast(100)
            2.64575131106459059050161575364?
        """
        return self.interval_fast(RealIntervalField(prec))

    def interval_exact(self, field):
        """
        Given a ``RealIntervalField``, compute the best possible
        approximation of this number in that field. Note that if this
        number is sufficiently close to some floating-point number
        (and, in particular, if this number is exactly representable in
        floating-point), then this will trigger exact computation, which
        may be very slow.

        EXAMPLES::

            sage: x = AA(2).sqrt()
            sage: y = x*x
            sage: x.interval(RIF)
            1.414213562373095?
            sage: x.interval_exact(RIF)
            1.414213562373095?
            sage: y.interval(RIF)
            2.000000000000000?
            sage: y.interval_exact(RIF)
            2
            sage: z = 1 + AA(2).sqrt() / 2^200
            sage: z.interval(RIF)
            1.000000000000001?
            sage: z.interval_exact(RIF)
            1.000000000000001?
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
        Given a ``RealField``, compute a good approximation to self in
        that field. The approximation will be off by at most two
        ulp's, except for numbers which are very close to 0, which
        will have an absolute error at most
        ``2**(-(field.prec()-1))``. Also, the rounding mode of the
        field is respected.

        EXAMPLES::

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
        return field(v)

    _mpfr_ = real_number

    def __float__(self):
        r"""
        Compute a good float approximation to self.

        EXAMPLES::

            sage: AA(golden_ratio).__float__()
            1.618033988749895
            sage: float(AA(sqrt(11)))
            3.3166247903554
        """
        return float(RR(self))

    def _complex_mpfr_field_(self, field):
        r"""
        Compute an approximation to this ``AlgebraicReal`` in the given field,
        which may be an interval field (in which case ``self.interval()`` is
        called) or any other real number field (in which case
        ``self.real_number()`` is called.

        Note that the field ``field`` should be a *complex* field (whose
        ``_real_field()`` method will be called to obtain a real subfield.)

        EXAMPLE::

            sage: AA(golden_ratio)._complex_mpfr_field_(ComplexIntervalField(100))
            1.618033988749894848204586834365?
            sage: AA(golden_ratio)._complex_mpfr_field_(ComplexField(100))
            1.6180339887498948482045868344
        """
        if is_ComplexIntervalField(field):
            return field(self.interval(field._real_field()))
        else:
            return field(self.real_number(field._real_field()))

    def real_exact(self, field):
        r"""
        Given a ``RealField``, compute the best possible approximation of
        this number in that field. Note that if this number is
        sufficiently close to the halfway point between two
        floating-point numbers in the field (for the default
        round-to-nearest mode) or if the number is sufficiently close
        to a floating-point number in the field (for directed rounding
        modes), then this will trigger exact computation, which may be
        very slow.

        The rounding mode of the field is respected.

        EXAMPLES::

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
            sage: y = (x-2).real_exact(RR).abs()
            sage: y == 0.0 or y == -0.0 # the sign of 0.0 is not significant in MPFI
            True
            sage: y = (x-2).real_exact(RealField(53, rnd='RNDD'))
            sage: y == 0.0 or y == -0.0 # same as above
            True
            sage: y = (x-2).real_exact(RealField(53, rnd='RNDU'))
            sage: y == 0.0 or y == -0.0 # idem
            True
            sage: y = (x-2).real_exact(RealField(53, rnd='RNDZ'))
            sage: y == 0.0 or y == -0.0 # ibidem
            True
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

        # Call the largest floating-point number <= self 'x'. Then
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
    r"""
    The subclass of ``ANDescr`` that represents an arbitrary
    rational. This class is private, and should not be used directly.
    """

    def __init__(self, x):
        """
        TESTS::

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
            raise TypeError("Illegal initializer for algebraic number rational")

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = AA(5/2); type(t._descr)
            <class 'sage.rings.qqbar.ANRational'>
            sage: loads(dumps(t)) == t
            True
        """
        return (ANRational, (self._value, ))

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLES::

            sage: QQbar(2/3)._repr_()
            '2/3'
        """
        return repr(self._value)

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always
        False, for rationals).

        EXAMPLES::

            sage: sage_input(QQbar(22/7), verify=True)
            # Verified
            QQbar(22/7)
            sage: sage_input(-AA(3)/5, verify=True)
            # Verified
            AA(-3/5)
            sage: sage_input(vector(AA, (0, 1/2, 1/3)), verify=True)
            # Verified
            vector(AA, [0, 1/2, 1/3])
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: rat = ANRational(9/10)
            sage: rat.handle_sage_input(sib, False, True)
            ({call: {atomic:QQbar}({binop:/ {atomic:9} {atomic:10}})}, False)
        """
        v = sib(self._value, True)
        if not coerce:
            v = sib.name('QQbar' if is_qqbar else 'AA')(v)
        return (v, False)

    def _interval_fast(self, prec):
        r"""
        Return an approximation to self in a real interval field of precision prec.

        EXAMPLE::

            sage: QQbar(355/113)._descr._interval_fast(30)
            3.14159292?
        """
        return RealIntervalField(prec)(self._value)

    def generator(self):
        r"""
        Return an :class:`AlgebraicGenerator` object associated to this
        element. Returns the trivial generator, since self is rational.

        EXAMPLE::

            sage: QQbar(0)._descr.generator()
            Trivial generator
        """
        return qq_generator

    def is_complex(self):
        r"""
        Return False, since rational numbers are real

        EXAMPLE::

            sage: QQbar(1/7)._descr.is_complex()
            False
        """
        return False

    def exactify(self):
        r"""
        Calculate self exactly. Since self is a rational number, return self.

        EXAMPLE::

            sage: a = QQbar(1/3)._descr
            sage: a.exactify() is a
            True
        """
        return self

    def is_simple(self):
        """
        Checks whether this descriptor represents a value with the same
        algebraic degree as the number field associated with the descriptor.

        This is always true for rational numbers.

        EXAMPLES::

            sage: AA(1/2)._descr.is_simple()
            True
        """
        return True

    def minpoly(self):
        r"""
        Return the min poly of self over `\QQ`.

        EXAMPLE::

            sage: QQbar(7)._descr.minpoly()
            x - 7
        """
        return QQx_x - self._value

    def neg(self, n):
        r"""
        Negation of self.

        EXAMPLE::

            sage: a = QQbar(3)
            sage: b = a._descr
            sage: type(b)
            <class 'sage.rings.qqbar.ANRational'>
            sage: b.neg(a)
            -3
        """
        return ANRational(-self._value)

    def invert(self, n):
        r"""
        1/self.

        EXAMPLE::

            sage: a = QQbar(3)
            sage: b = a._descr
            sage: b.invert(a)
            1/3
        """
        return ANRational(~self._value)

    def abs(self, n):
        r"""
        Absolute value of self.

        EXAMPLE::

            sage: a = QQbar(3)
            sage: b = a._descr
            sage: b.abs(a)
            3
        """
        return ANRational(abs(self._value))

    def rational_argument(self, n):
        r"""
        Return the argument of self divided by `2 \pi`, or ``None`` if this
        element is 0.

        EXAMPLE::

            sage: QQbar(3)._descr.rational_argument(None)
            0
            sage: QQbar(-3)._descr.rational_argument(None)
            1/2
            sage: QQbar(0)._descr.rational_argument(None) is None
            True
        """
        if self._value > 0:
            return QQ(0)
        if self._value < 0:
            return QQ(1)/2
        return None

    def angle(self):
        r"""
        Return a rational number `q \in (-1/2, 1/2]` such that ``self`` is a rational multiple of
        `e^{2\pi i q}`. Always returns 0, since this element is rational.

        EXAMPLE::

            sage: QQbar(3)._descr.angle()
            0
            sage: QQbar(-3)._descr.angle()
            0
            sage: QQbar(0)._descr.angle()
            0
        """
        return QQ_0

    def scale(self):
        r"""
        Return a rational number `r` such that ``self`` is equal to `r e^{2 \pi
        i q}` for some `q \in (-1/2, 1/2]`.  In other words, just return self
        as a rational number.

        EXAMPLE::

            sage: QQbar(-3)._descr.scale()
            -3
        """
        return self._value

def is_AlgebraicReal(x):
    r"""
    Test if ``x`` is an instance of :class:`~AlgebraicReal`. For internal use.

    EXAMPLE::

        sage: from sage.rings.qqbar import is_AlgebraicReal
        sage: is_AlgebraicReal(AA(sqrt(2)))
        True
        sage: is_AlgebraicReal(QQbar(sqrt(2)))
        False
        sage: is_AlgebraicReal("spam")
        False
    """
    return isinstance(x, AlgebraicReal)

def is_AlgebraicNumber(x):
    r"""
    Test if ``x`` is an instance of :class:`~AlgebraicNumber`. For internal use.

    EXAMPLE::

        sage: from sage.rings.qqbar import is_AlgebraicNumber
        sage: is_AlgebraicNumber(AA(sqrt(2)))
        False
        sage: is_AlgebraicNumber(QQbar(sqrt(2)))
        True
        sage: is_AlgebraicNumber("spam")
        False
    """
    return isinstance(x, AlgebraicNumber)

QQbarPoly = PolynomialRing(QQbar, 'x')
AAPoly = PolynomialRing(AA, 'x')

class AlgebraicPolynomialTracker(SageObject):
    r"""
    Keeps track of a polynomial used for algebraic numbers.

    If multiple algebraic numbers are created as roots of a single
    polynomial, this allows the polynomial and information about
    the polynomial to be shared. This reduces work if the polynomial
    must be recomputed at higher precision, or if it must be factored.

    This class is private, and should only be constructed by
    ``AA.common_polynomial()`` or ``QQbar.common_polynomial()``, and should
    only be used as an argument to ``AA.polynomial_root()`` or
    ``QQbar.polynomial_root()``. (It doesn't matter whether you create
    the common polynomial with ``AA.common_polynomial()`` or
    ``QQbar.common_polynomial()``.)

    EXAMPLES::

        sage: x = polygen(QQbar)
        sage: P = QQbar.common_polynomial(x^2 - x - 1)
        sage: P
        x^2 - x - 1
        sage: QQbar.polynomial_root(P, RIF(1, 2))
        1.618033988749895?
    """

    def __init__(self, poly):
        r"""
        Initialize this AlgebraicPolynomialTracker object.

        EXAMPLE::

            sage: x = polygen(QQbar)
            sage: P = QQbar.common_polynomial(x^2 - x - 1)
            sage: type(P) # indirect doctest
            <class 'sage.rings.qqbar.AlgebraicPolynomialTracker'>
        """
        if not is_Polynomial(poly):
            raise ValueError("Trying to create AlgebraicPolynomialTracker on non-Polynomial")
        if isinstance(poly.base_ring(), AlgebraicField_common):
            complex = is_AlgebraicField(poly.base_ring())
        else:
            try:
                poly = poly.change_ring(AA)
                complex = False
            except (TypeError, ValueError):
                poly = poly.change_ring(QQbar)
                complex = True
        self._poly = poly
        self._complex = complex
        self._exact = False
        self._roots_cache = {}

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]
            sage: type(v._descr._poly)
            <class 'sage.rings.qqbar.AlgebraicPolynomialTracker'>
            sage: loads(dumps(v)) == v
            True
        """
        return (AlgebraicPolynomialTracker, (self._poly, ))

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: sage_input(AA.common_polynomial(x^3 - 7))
            R.<x> = AA[]
            AA.common_polynomial(x^3 - 7)
            sage: x = polygen(AA)
            sage: p = sqrt(AA(2)) * x^2 - sqrt(AA(3))
            sage: cp = AA.common_polynomial(p)
            sage: sage_input((cp, cp))
            R.<x> = AA[]
            cp = AA.common_polynomial(AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))*x^2 - AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774))))
            (cp, cp)
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: cp._sage_input_(sib, False)
            {call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:* {call: {getattr: {atomic:AA}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:2}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.4142135623730949})}, {call: {atomic:RR}({atomic:1.4142135623730951})})})} {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}}} {call: {getattr: {atomic:AA}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:3}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.7320508075688772})}, {call: {atomic:RR}({atomic:1.7320508075688774})})})}})}
        """
        # XXX It would be nicer to skip the "AA.common_polynomial()"
        # wrapper if the polynomial is not actually shared. But
        # sage_input.py isn't quite that generic.
        v = sib.name('AA').common_polynomial(self._poly)
        sib.id_cache(self, v, 'cp')
        return v

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: x = polygen(QQ)
            sage: AA.common_polynomial(x^3 - 7)._repr_()
            'x^3 - 7'
        """
        return repr(self._poly)

    def poly(self):
        r"""
        Return the underlying polynomial of self.

        EXAMPLE::

            sage: x = polygen(QQ); f = x^3 - 7
            sage: g = AA.common_polynomial(f)
            sage: g.poly() == f
            True
        """
        return self._poly

    def is_complex(self):
        r"""
        Return True if the coefficients of this polynomial are non-real.

        EXAMPLE::

            sage: x = polygen(QQ); f = x^3 - 7
            sage: g = AA.common_polynomial(f)
            sage: g.is_complex()
            False
            sage: QQbar.common_polynomial(x^3 - QQbar(I)).is_complex()
            True
        """
        return self._complex

    def complex_roots(self, prec, multiplicity):
        """
        Find the roots of self in the complex field to precision prec.

        EXAMPLES::

            sage: x = polygen(ZZ)
            sage: cp = AA.common_polynomial(x^4 - 2)

        Note that the precision is not guaranteed to find the tightest
        possible interval since complex_roots() depends on the
        underlying BLAS implementation. ::

            sage: cp.complex_roots(30, 1)
            [-1.18920711500272...?,
             1.189207115002721?,
             -1.189207115002721?*I,
             1.189207115002721?*I]
        """
        if multiplicity in self._roots_cache:
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

        EXAMPLE::

            sage: x = polygen(AA)
            sage: p = sqrt(AA(2)) * x^2 - sqrt(AA(3))
            sage: cp = AA.common_polynomial(p)
            sage: cp._exact
            False
            sage: cp.exactify()
            sage: cp._exact
            True
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
        r"""
        EXAMPLE::

            sage: x=polygen(QQ); f=QQbar.common_polynomial(x^4 + 4)
            sage: f.factors()
            [y^2 - 2*y + 2, y^2 + 2*y + 2]
        """
        self.exactify()
        return self._factors

    def generator(self):
        r"""
        Return an :class:`AlgebraicGenerator` for a number field containing all
        the coefficients of self.

        EXAMPLE::

            sage: x = polygen(AA)
            sage: p = sqrt(AA(2)) * x^2 - sqrt(AA(3))
            sage: cp = AA.common_polynomial(p)
            sage: cp.generator()
            Number Field in a with defining polynomial y^4 - 4*y^2 + 1 with a in 1.931851652578137?
        """
        self.exactify()
        return self._gen

class ANRoot(ANDescr):
    """
    The subclass of ``ANDescr`` that represents a particular
    root of a polynomial with algebraic coefficients.
    This class is private, and should not be used directly.
    """
    def __init__(self, poly, interval, multiplicity=1):
        r"""
        Initialize this ``ANRoot`` object.

        EXAMPLE::

            sage: x = polygen(QQ); f = (x^3 + x + 1).roots(AA,multiplicities=False)[0]._descr
            sage: type(f) # indirect doctest
            <class 'sage.rings.qqbar.ANRoot'>
        """
        if not isinstance(poly, AlgebraicPolynomialTracker):
            poly = AlgebraicPolynomialTracker(poly)
        self._poly = poly
        self._multiplicity = multiplicity
        self._complex = is_ComplexIntervalFieldElement(interval)
        self._complex_poly = poly.is_complex()
        self._interval = self.refine_interval(interval, 64)

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]
            sage: type(v._descr)
            <class 'sage.rings.qqbar.ANRoot'>
            sage: loads(dumps(v)) == v
            True
        """
        return (ANRoot, (self._poly, self._interval, self._multiplicity))

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: x=polygen(QQ); v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]
            sage: v._descr._repr_()
            'Root 1.618033988749894849? of x^2 - x - 1'
        """
        return 'Root %s of %s'%(self._interval, self._poly)

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always True,
        for ``ANRoot``).

        EXAMPLES::

            sage: sage_input((AA(3)^(1/2))^(1/3), verify=True)
            # Verified
            R.<x> = AA[]
            AA.polynomial_root(AA.common_polynomial(x^3 - AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))), RIF(RR(1.2009369551760025), RR(1.2009369551760027)))

        These two examples are too big to verify quickly. (Verification
        would create a field of degree 28.)::

            sage: sage_input((sqrt(AA(3))^(5/7))^(9/4))
            R.<x> = AA[]
            v1 = AA.polynomial_root(AA.common_polynomial(x^2 - 3), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))
            v2 = v1*v1
            v3 = AA.polynomial_root(AA.common_polynomial(x^7 - v2*v2*v1), RIF(RR(1.4804728524798112), RR(1.4804728524798114)))
            v4 = v3*v3
            v5 = v4*v4
            AA.polynomial_root(AA.common_polynomial(x^4 - v5*v5*v3), RIF(RR(2.4176921938267877), RR(2.4176921938267881)))
            sage: sage_input((sqrt(QQbar(-7))^(5/7))^(9/4))
            R.<x> = QQbar[]
            v1 = QQbar.polynomial_root(AA.common_polynomial(x^2 + 7), CIF(RIF(RR(0)), RIF(RR(2.6457513110645903), RR(2.6457513110645907))))
            v2 = v1*v1
            v3 = QQbar.polynomial_root(AA.common_polynomial(x^7 - v2*v2*v1), CIF(RIF(RR(0.8693488875796217), RR(0.86934888757962181)), RIF(RR(1.8052215661454434), RR(1.8052215661454436))))
            v4 = v3*v3
            v5 = v4*v4
            QQbar.polynomial_root(AA.common_polynomial(x^4 - v5*v5*v3), CIF(RIF(-RR(3.8954086044650791), -RR(3.8954086044650786)), RIF(RR(2.7639398015408925), RR(2.7639398015408929))))
            sage: x = polygen(QQ)
            sage: sage_input(AA.polynomial_root(x^2-x-1, RIF(1, 2)), verify=True)
            # Verified
            R.<x> = AA[]
            AA.polynomial_root(AA.common_polynomial(x^2 - x - 1), RIF(RR(1.6180339887498947), RR(1.6180339887498949)))
            sage: sage_input(QQbar.polynomial_root(x^3-5, CIF(RIF(-3, 0), RIF(0, 3))), verify=True)
            # Verified
            R.<x> = AA[]
            QQbar.polynomial_root(AA.common_polynomial(x^3 - 5), CIF(RIF(-RR(0.85498797333834853), -RR(0.85498797333834842)), RIF(RR(1.4808826096823642), RR(1.4808826096823644))))
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: rt = ANRoot(x^3 - 2, RIF(0, 4))
            sage: rt.handle_sage_input(sib, False, True)
            ({call: {getattr: {atomic:QQbar}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:3}} {atomic:2}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.259921049894873})}, {call: {atomic:RR}({atomic:1.2599210498948732})})})}, True)
        """
        parent = sib.name('QQbar' if is_qqbar else 'AA')
        poly = sib(self._poly)
        intv = self._interval
        # Check whether a 53-bit interval actually isolates the root.
        # If so, use it, because 53-bit intervals print prettier.
        if is_ComplexIntervalFieldElement(intv):
            loose_intv = CIF(intv)
        else:
            loose_intv = RIF(intv)
        # If the derivative of the polynomial is bounded away from 0
        # over this interval, then it definitely isolates a root.
        if self._poly._poly.derivative()(loose_intv) != 0:
            good_intv = loose_intv
        else:
            good_intv = intv
        return (parent.polynomial_root(poly, sib(good_intv)), True)

    def is_complex(self):
        r"""
        Whether this is a root in `\overline{\QQ}` (rather than `\mathbf{A}`).
        Note that this may return True even if the root is actually real, as
        the second example shows; it does *not* trigger exact computation to
        see if the root is real.

        EXAMPLE::

            sage: x = polygen(QQ)
            sage: (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.is_complex()
            False
            sage: (x^2 - x - 1).roots(ring=QQbar, multiplicities=False)[1]._descr.is_complex()
            True
        """
        return self._complex

    def conjugate(self, n):
        r"""
        Complex conjugate of this ANRoot object.

        EXAMPLE::

            sage: a = (x^2 + 23).roots(ring=QQbar, multiplicities=False)[0]
            sage: b = a._descr
            sage: type(b)
            <class 'sage.rings.qqbar.ANRoot'>
            sage: c = b.conjugate(a); c
            <class 'sage.rings.qqbar.ANUnaryExpr'>
            sage: c.exactify()
            -2*a + 1 where a^2 - a + 6 = 0 and a in 0.50000000000000000? - 2.397915761656360?*I
        """
        if not self._complex:
            return self
        if not self._complex_poly:
            return ANRoot(self._poly, self._interval.conjugate(), self._multiplicity)

        return ANUnaryExpr(n, 'conjugate')

    def refine_interval(self, interval, prec):
        r"""
        Takes an interval which is assumed to enclose exactly one root
        of the polynomial (or, with multiplicity=`k`, exactly one root
        of the `k-1`-st derivative); and a precision, in bits.

        Tries to find a narrow interval enclosing the root using
        interval arithmetic of the given precision. (No particular
        number of resulting bits of precision is guaranteed.)

        Uses a combination of Newton's method (adapted for interval
        arithmetic) and bisection. The algorithm will converge very
        quickly if started with a sufficiently narrow interval.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(AA)
            sage: rt2 = ANRoot(x^2 - 2, RIF(0, 2))
            sage: rt2.refine_interval(RIF(0, 2), 75)
            1.4142135623730950488017?
        """
        if self._complex or self._complex_poly:
            v = self._complex_refine_interval(interval, prec)
            if self._complex:
                return v
            else:
                return v.real()
        else:
            return self._real_refine_interval(interval, prec)

    def _real_refine_interval(self, interval, prec):
        r"""
        Does the calculation for ``refine_interval``.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(AA)
            sage: rt2 = ANRoot(x^2 - 2, RIF(0, 2))
            sage: rt2.refine_interval(RIF(0, 2), 75) # indirect doctest
            1.4142135623730950488017?
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
        # this case here means we do not have to worry about iterating too
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
                raise ValueError("Refining interval that does not bound unique root!")

            # Use a simple algorithm:
            # Try an interval Newton-Raphson step. If this does not add at
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
        r"""
        Takes an interval which is assumed to enclose exactly one root
        of the polynomial (or, with multiplicity=`k`, exactly one root
        of the `k-1`-st derivative); and a precision, in bits.

        Tries to find a narrow interval enclosing the root using
        interval arithmetic of the given precision. (No particular
        number of resulting bits of precision is guaranteed.)

        Uses Newton's method (adapted for interval arithmetic). The
        algorithm will converge very quickly if started with a
        sufficiently narrow interval. If Newton's method fails, then
        we falls back on computing all the roots of the polynomial
        numerically, and select the appropriate root.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(QQbar)
            sage: intv = CIF(RIF(0, 1), RIF(0.1, 1))
            sage: rt = ANRoot(x^5 - 1, intv)
            sage: new_intv = rt.refine_interval(intv, 53); new_intv # indirect doctest
            0.3090169943749474241? + 0.951056516295153573?*I
            sage: rt.refine_interval(new_intv, 70)
            0.30901699437494742411? + 0.95105651629515357212?*I
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
        # this case here means we do not have to worry about iterating too
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
                    # Wow; we nailed it exactly. (This may happen
                    # whenever the root is exactly equal to some
                    # floating-point number, and cannot happen
                    # if the root is not equal to a floating-point
                    # number.)  We just return the perfect answer.
                    return interval

                if new_diam == diam:
                    # We're not getting any better. There are two
                    # possible reasons for this. Either we have
                    # refined as much as possible given the imprecision
                    # of our interval polynomial, and we have the best
                    # answer we're going to get at this precision;
                    # or we started with a poor approximation to the
                    # root, resulting in a broad range of possible
                    # slopes in this interval, and Newton-Raphson refining
                    # is not going to help.

                    # I do not have a formal proof, but I believe the
                    # following test differentiates between these two
                    # behaviors. (If I'm wrong, we might get bad behavior
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
                    # not to be a root. This should let us divide
                    # the interval in half, and improve on our previous
                    # estimates. I can only think of two reasons why
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
                    # be bounded away from zero. But we compare
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
        polynomial, and checking which one is in interval. Slow but sure.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(QQbar)
            sage: intv = CIF(RIF(0, 1), RIF(0.1, 1))
            sage: rt = ANRoot(x^5 - 1, intv)
            sage: rt._complex_isolate_interval(intv, 53)
            0.3090169943749474241? + 0.951056516295153573?*I
        """
        rts = self._poly.complex_roots(prec, self._multiplicity)

        # Find all the roots that overlap interval.
        our_root = [rt for rt in rts if rt.overlaps(interval)]

        if len(our_root) == 1:
            return our_root[0]

        if len(our_root) == 0:
            raise ValueError("Complex root interval does not include any roots")

        # We have more than one root that overlap the current interval.
        # Technically, this might not be an error; perhaps the actual
        # root is just outside our interval, even though the (presumably
        # tight) interval containing that root touches our interval.

        # But it seems far more likely that the provided interval is
        # just too big.

        raise ValueError("Complex root interval probably includes multiple roots")

    def exactify(self):
        """
        Returns either an ``ANRational`` or an
        ``ANExtensionElement`` with the same value as this number.

        EXAMPLES::

            sage: from sage.rings.qqbar import ANRoot
            sage: x = polygen(QQbar)
            sage: two = ANRoot((x-2)*(x-sqrt(QQbar(2))), RIF(1.9, 2.1))
            sage: two.exactify()
            2
            sage: strange = ANRoot(x^2 + sqrt(QQbar(3))*x - sqrt(QQbar(2)), RIF(-0, 1))
            sage: strange.exactify()
            a where a^8 - 6*a^6 + 5*a^4 - 12*a^2 + 4 = 0 and a in 0.6051012265139511?

        TESTS:

        Verify that :trac:`12727` is fixed::

            sage: m = sqrt(sin(pi/5)); a = QQbar(m); b = AA(m)
            sage: a.minpoly()
            x^8 - 5/4*x^4 + 5/16
            sage: b.minpoly()
            x^8 - 5/4*x^4 + 5/16
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
                v = [c.polynomial()(gen_val) for c in factor]
                # This try/except can be triggered if ifield is Real
                # but the entries in v have some imaginary part that
                # is only known to be 0 to very low precision, e.g.,
                # as in Trac #12727.  In such cases, we instead create
                # the polynomial over the appropriate complex interval
                # field, which is mathematically safe, unlike taking
                # real parts would be.
                try:
                    ip = if_poly(v)
                except TypeError:
                    if_poly = ComplexIntervalField(prec)['x']
                    ip = if_poly(v)
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
        Recompute the interval enclosing this ``ANRoot`` object at higher
        precision.

        EXAMPLE::

            sage: x = polygen(QQ); y = (x^3 + x + 1).roots(AA,multiplicities=False)[0]
            sage: z = y._descr
            sage: z._interval.prec()
            64
            sage: z._more_precision()
            sage: z._interval.prec()
            128
        """
        prec = self._interval.prec()
        self._interval = self.refine_interval(self._interval, prec*2)

    def _interval_fast(self, prec):
        """
        Given a RealIntervalField, compute the value of this number
        using interval arithmetic of at least the precision of the field,
        and return the value in that field. (More precision may be used
        in the computation.)

        EXAMPLE::

            sage: x = polygen(QQ); y = (x^3 + x + 1).roots(AA,multiplicities=False)[0]._descr
            sage: y._interval_fast(128)
            -0.68232780382801932736948373971104825689?

        Check that :trac:`15493` is fixed::

            sage: y._interval_fast(20).parent() is RealIntervalField(20)
            True
        """
        if prec == self._interval.prec():
            return self._interval
        if prec < self._interval.prec():
            return self._interval.parent().to_prec(prec)(self._interval)
        self._more_precision()
        return self._interval_fast(prec)

class ANExtensionElement(ANDescr):
    r"""
    The subclass of ``ANDescr`` that represents a number field
    element in terms of a specific generator. Consists of a polynomial
    with rational coefficients in terms of the generator, and the
    generator itself, an ``AlgebraicGenerator``.
    """

    def __new__(self, generator, value):
        if value.is_rational():
            return ANRational(value._rational_())
        else:
            return ANDescr.__new__(self)

    def __init__(self, generator, value):
        self._generator = generator
        self._value = value
        self._exactly_real = not generator.is_complex()

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]
            sage: v.exactify()
            sage: type(v._descr)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: loads(dumps(v)) == v
            True
        """
        return (ANExtensionElement, (self._generator, self._value))

    def _repr_(self):
        fgen = self._generator._field.gen()
        sgen = str(fgen)
        return '%s where %s = 0 and %s in %s'%(self._value,
                                               self._generator.field().polynomial()._repr(name=sgen),
                                               sgen,
                                               self._generator._interval_fast(53))

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always True,
        for ``ANExtensionElement``).

        EXAMPLES::

            sage: I = QQbar(I)
            sage: sage_input(3+4*I, verify=True)
            # Verified
            QQbar(3 + 4*I)
            sage: v = QQbar.zeta(3) + QQbar.zeta(5)
            sage: v - v == 0
            True
            sage: sage_input(vector(QQbar, (4-3*I, QQbar.zeta(7))), verify=True)
            # Verified
            R.<x> = AA[]
            vector(QQbar, [4 - 3*I, QQbar.polynomial_root(AA.common_polynomial(x^6 + x^5 + x^4 + x^3 + x^2 + x + 1), CIF(RIF(RR(0.62348980185873348), RR(0.62348980185873359)), RIF(RR(0.7818314824680298), RR(0.78183148246802991))))])
            sage: sage_input(v, verify=True)
            # Verified
            R.<x> = AA[]
            v = QQbar.polynomial_root(AA.common_polynomial(x^8 - x^7 + x^5 - x^4 + x^3 - x + 1), CIF(RIF(RR(0.91354545764260087), RR(0.91354545764260098)), RIF(RR(0.40673664307580015), RR(0.40673664307580021))))
            v^5 + v^3
            sage: v = QQbar(sqrt(AA(2)))
            sage: v.exactify()
            sage: sage_input(v, verify=True)
            # Verified
            R.<x> = AA[]
            QQbar(AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951))))
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: extel = ANExtensionElement(QQbar_I_generator, QQbar_I_generator.field().gen() + 1)
            sage: extel.handle_sage_input(sib, False, True)
            ({call: {atomic:QQbar}({binop:+ {atomic:1} {atomic:I}})}, True)
        """
        if self._generator is QQbar_I_generator:
            assert(is_qqbar)
            re, im = self._value.list()
            im_part = sib.prod([sib(im, True), sib.name('I')], simplify=True)
            v = sib.sum([sib(re, True), im_part], simplify=True)
            if coerce != 2:
                v = sib.name('QQbar')(v)
                return (v, True)
            return (v, False)

        result_is_qqbar = self._generator.is_complex()

        rt = sib(self._generator.root_as_algebraic())
        # For the best fidelity, we really ought to somehow ensure
        # that rt is exactified, but sage_input doesn't support that
        # nicely. Skip it for now.
        # The following is copied with slight mods from polynomial_element.pyx
        coeffs = [sib(c, True) for c in self._value.list()]
        terms = []
        for i in range(len(coeffs)-1, -1, -1):
            if i > 0:
                if i > 1:
                    rt_pow = rt**sib.int(i)
                else:
                    rt_pow = rt
                terms.append(sib.prod((coeffs[i], rt_pow), simplify=True))
            else:
                terms.append(coeffs[i])
        v = sib.sum(terms, simplify=True)
        if result_is_qqbar != is_qqbar:
            v = sib.name('QQbar' if is_qqbar else 'AA')(v)
        return (v, True)

    def is_complex(self):
        r"""
        Return True if the number field that defines this element is not real.
        This does not imply that the element itself is definitely non-real, as
        in the example below.

        EXAMPLE::

            sage: rt2 = QQbar(sqrt(2))
            sage: rtm3 = QQbar(sqrt(-3))
            sage: x = rtm3 + rt2 - rtm3
            sage: x.exactify()
            sage: y = x._descr
            sage: type(y)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: y.is_complex()
            True
            sage: x.imag() == 0
            True
        """
        return not self._exactly_real

    def is_simple(self):
        r"""
        Checks whether this descriptor represents a value with the same
        algebraic degree as the number field associated with the descriptor.

        For ``ANExtensionElement`` elements, we check this by
        comparing the degree of the minimal polynomial to the degree
        of the field.

        EXAMPLES::

            sage: rt2 = AA(sqrt(2))
            sage: rt3 = AA(sqrt(3))
            sage: rt2b = rt3 + rt2 - rt3
            sage: rt2.exactify()
            sage: rt2._descr
            a where a^2 - 2 = 0 and a in 1.414213562373095?
            sage: rt2._descr.is_simple()
            True

            sage: rt2b.exactify()
            sage: rt2b._descr
            a^3 - 3*a where a^4 - 4*a^2 + 1 = 0 and a in 1.931851652578137?
            sage: rt2b._descr.is_simple()
            False
        """
        try:
            return self._is_simple
        except AttributeError:
            self._is_simple = (self.minpoly().degree() == self.generator().field().degree())
            return self._is_simple

    def generator(self):
        r"""
        Return the :class:`~AlgebraicGenerator` object corresponding to self.

        EXAMPLE::

            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.exactify()
            sage: v.generator()
            Number Field in a with defining polynomial y^2 - y - 1 with a in 1.618033988749895?
        """
        return self._generator

    def exactify(self):
        r"""
        Return an exact representation of self. Since self is already exact,
        just return self.

        EXAMPLE::

            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.exactify()
            sage: type(v)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: v.exactify() is v
            True
        """
        return self

    def field_element_value(self):
        r"""
        Return the underlying number field element.

        EXAMPLE::

            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.exactify()
            sage: v.field_element_value()
            a
        """
        return self._value

    def minpoly(self):
        """
        Compute the minimal polynomial of this algebraic number.

        EXAMPLES::

            sage: v = (x^2 - x - 1).roots(ring=AA, multiplicities=False)[1]._descr.exactify()
            sage: type(v)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: v.minpoly()
            x^2 - x - 1
        """
        try:
            return self._minpoly
        except AttributeError:
            self._minpoly = self._value.minpoly()
            return self._minpoly

    def simplify(self, n):
        """
        Compute an exact representation for this descriptor, in the
        smallest possible number field.

        INPUT:

        - ``n`` -- The element of ``AA`` or ``QQbar`` corresponding
          to this descriptor.

        EXAMPLES::

            sage: rt2 = AA(sqrt(2))
            sage: rt3 = AA(sqrt(3))
            sage: rt2b = rt3 + rt2 - rt3
            sage: rt2b.exactify()
            sage: rt2b._descr
            a^3 - 3*a where a^4 - 4*a^2 + 1 = 0 and a in 1.931851652578137?
            sage: rt2b._descr.simplify(rt2b)
            a where a^2 - 2 = 0 and a in 1.414213562373095?
        """

        if self.is_simple():
            return self

        # This is very inefficient...
        # for instance, the .exactify() call will try to factor poly,
        # even though we know that poly is irreducible
        poly = self.minpoly()
        intv = isolating_interval(lambda prec: n._interval_fast(prec), poly)
        new_v = QQbar.polynomial_root(poly, intv)
        new_v.exactify()
        return new_v._descr

    def _interval_fast(self, prec):
        gen_val = self._generator._interval_fast(prec)
        v = self._value.polynomial()(gen_val)
        if self._exactly_real and is_ComplexIntervalFieldElement(v):
            return v.real()
        return v

    # for these three functions the argument n is not used (but it is there
    # anyway for compatibility)

    def neg(self, n):
        r"""
        Negation of self.

        EXAMPLE::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))
            sage: a.exactify()
            sage: b = a._descr
            sage: type(b)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: b.neg(a)
            1/3*a^3 - 2/3*a^2 + 4/3*a - 2 where a^4 - 2*a^3 + a^2 - 6*a + 9 = 0 and a in -0.7247448713915890? - 1.573132184970987?*I
            sage: b.neg("ham spam and eggs")
            1/3*a^3 - 2/3*a^2 + 4/3*a - 2 where a^4 - 2*a^3 + a^2 - 6*a + 9 = 0 and a in -0.7247448713915890? - 1.573132184970987?*I
        """
        return ANExtensionElement(self._generator, -self._value)

    def invert(self, n):
        r"""
        1/self.

        EXAMPLE::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))
            sage: a.exactify()
            sage: b = a._descr
            sage: type(b)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: b.invert(a)
            7/3*a^3 - 2/3*a^2 + 4/3*a - 12 where a^4 - 2*a^3 + a^2 - 6*a + 9 = 0 and a in -0.7247448713915890? - 1.573132184970987?*I
            sage: b.invert("ham spam and eggs")
            7/3*a^3 - 2/3*a^2 + 4/3*a - 12 where a^4 - 2*a^3 + a^2 - 6*a + 9 = 0 and a in -0.7247448713915890? - 1.573132184970987?*I
        """
        return ANExtensionElement(self._generator, ~self._value)

    def conjugate(self, n):
        r"""
        Negation of self.

        EXAMPLE::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))
            sage: a.exactify()
            sage: b = a._descr
            sage: type(b)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: b.conjugate(a)
            -1/3*a^3 + 2/3*a^2 - 4/3*a + 2 where a^4 - 2*a^3 + a^2 - 6*a + 9 = 0 and a in -0.7247448713915890? + 1.573132184970987?*I
            sage: b.conjugate("ham spam and eggs")
            -1/3*a^3 + 2/3*a^2 - 4/3*a + 2 where a^4 - 2*a^3 + a^2 - 6*a + 9 = 0 and a in -0.7247448713915890? + 1.573132184970987?*I
        """
        if self._exactly_real:
            return self
        else:
            return ANExtensionElement(self._generator.conjugate(), self._value)

    # The rest of these unary operations do actually use n, which is an
    # AlgebraicNumber pointing to self.

    def norm(self, n):
        r"""
        Norm of self (square of complex absolute value)

        EXAMPLE::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))
            sage: a.exactify()
            sage: b = a._descr
            sage: type(b)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: b.norm(a)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        if self._exactly_real:
            return (n*n)._descr
        elif self._generator is QQbar_I_generator:
            return ANRational(self._value.norm())
        else:
            return ANUnaryExpr(n, 'norm')

    def abs(self, n):
        r"""
        Return the absolute value of self (square root of the norm).

        EXAMPLE::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(-3))
            sage: a.exactify()
            sage: b = a._descr
            sage: type(b)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: b.abs(a)
            Root 3.146264369941972342? of x^2 - 9.89897948556636?
        """
        return AlgebraicReal(self.norm(n)).sqrt()._descr

    def rational_argument(self, n):
        r"""
        If the argument of self is `2\pi` times some rational number in `[1/2,
        -1/2)`, return that rational; otherwise, return ``None``.

        EXAMPLE::

            sage: a = QQbar(sqrt(-2)) + QQbar(sqrt(3))
            sage: a.exactify()
            sage: b = a._descr
            sage: type(b)
            <class 'sage.rings.qqbar.ANExtensionElement'>
            sage: b.rational_argument(a) is None
            True
            sage: x = polygen(QQ)
            sage: a = (x^4 + 1).roots(QQbar, multiplicities=False)[0]
            sage: a.exactify()
            sage: b = a._descr
            sage: b.rational_argument(a)
            -3/8
        """
        # If the argument of self is 2*pi times some rational number a/b,
        # then self/abs(self) is a root of the b'th cyclotomic polynomial.
        # This implies that the algebraic degree of self is at least
        # phi(b). Working backward, we know that the algebraic degree
        # of self is at most the degree of the generator, so that gives
        # an upper bound on phi(b). According to
        # http://mathworld.wolfram.com/TotientFunction.html,
        # phi(b) >= sqrt(b) for b > 6; this gives us an upper bound on b.
        # We then check to see if this is possible; if so, we test
        # to see if it actually holds.

        if self._exactly_real:
            if n > 0:
                return QQ.zero()
            else:
                return QQ((1,2))

        gen_degree = self._generator._field.degree()
        if gen_degree <= 2:
            max_b = 6
        else:
            max_b = gen_degree*gen_degree
        rat_arg_fl = n._interval_fast(128).argument() / RealIntervalField(128).pi() / 2
        rat_arg = rat_arg_fl.simplest_rational()
        if rat_arg.denominator() > max_b:
            return None
        n_exp = n ** rat_arg.denominator()
        if n_exp.real() > AA.zero() and n_exp.imag().is_zero():
            return rat_arg
        # Strictly speaking, we need to look for the second-simplest
        # rational in rat_arg_fl and make sure its denominator is > max_b.
        # For now, we just punt.
        raise NotImplementedError

class ANUnaryExpr(ANDescr):
    def __init__(self, arg, op):
        r"""
        Initialize this ANUnaryExpr.

        EXAMPLE::

            sage: t = ~QQbar(sqrt(2)); type(t._descr) # indirect doctest
            <class 'sage.rings.qqbar.ANUnaryExpr'>
        """
        self._arg = arg
        self._op = op
        self._complex = True

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = ~QQbar(sqrt(2)); type(t._descr)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
            sage: loads(dumps(t)) == 1/QQbar(sqrt(2))
            True
        """
        return (ANUnaryExpr, (self._arg, self._op))

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always
        True for ``ANUnaryExpr``).

        EXAMPLES::

            sage: sage_input(-sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            -AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(~sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            ~AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(sqrt(QQbar(-3)).conjugate(), verify=True)
            # Verified
            R.<x> = QQbar[]
            QQbar.polynomial_root(AA.common_polynomial(x^2 + 3), CIF(RIF(RR(0)), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))).conjugate()
            sage: sage_input(QQbar.zeta(3).real(), verify=True)
            # Verified
            R.<x> = AA[]
            QQbar.polynomial_root(AA.common_polynomial(x^2 + x + 1), CIF(RIF(-RR(0.50000000000000011), -RR(0.49999999999999994)), RIF(RR(0.8660254037844386), RR(0.86602540378443871)))).real()
            sage: sage_input(QQbar.zeta(3).imag(), verify=True)
            # Verified
            R.<x> = AA[]
            QQbar.polynomial_root(AA.common_polynomial(x^2 + x + 1), CIF(RIF(-RR(0.50000000000000011), -RR(0.49999999999999994)), RIF(RR(0.8660254037844386), RR(0.86602540378443871)))).imag()
            sage: sage_input(abs(sqrt(QQbar(-3))), verify=True)
            # Verified
            R.<x> = QQbar[]
            abs(QQbar.polynomial_root(AA.common_polynomial(x^2 + 3), CIF(RIF(RR(0)), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))))
            sage: sage_input(sqrt(QQbar(-3)).norm(), verify=True)
            # Verified
            R.<x> = QQbar[]
            QQbar.polynomial_root(AA.common_polynomial(x^2 + 3), CIF(RIF(RR(0)), RIF(RR(1.7320508075688772), RR(1.7320508075688774)))).norm()
            sage: sage_input(QQbar(QQbar.zeta(3).real()), verify=True)
            # Verified
            R.<x> = AA[]
            QQbar(QQbar.polynomial_root(AA.common_polynomial(x^2 + x + 1), CIF(RIF(-RR(0.50000000000000011), -RR(0.49999999999999994)), RIF(RR(0.8660254037844386), RR(0.86602540378443871)))).real())
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: unexp = ANUnaryExpr(sqrt(AA(2)), '~')
            sage: unexp.handle_sage_input(sib, False, False)
            ({unop:~ {call: {getattr: {atomic:AA}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:2}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.4142135623730949})}, {call: {atomic:RR}({atomic:1.4142135623730951})})})}},
             True)
            sage: unexp.handle_sage_input(sib, False, True)
            ({call: {atomic:QQbar}({unop:~ {call: {getattr: {atomic:AA}.polynomial_root}({call: {getattr: {atomic:AA}.common_polynomial}({binop:- {binop:** {gen:x {constr_parent: {subscr: {atomic:AA}[{atomic:'x'}]} with gens: ('x',)}} {atomic:2}} {atomic:2}})}, {call: {atomic:RIF}({call: {atomic:RR}({atomic:1.4142135623730949})}, {call: {atomic:RR}({atomic:1.4142135623730951})})})}})},
             True)
        """
        arg_is_qqbar = self._arg.parent() is QQbar
        v = sib(self._arg)
        op = self._op
        if op == '-':
            v = -v
        elif op == '~':
            v = ~v
        elif op == 'conjugate':
            v = v.conjugate()
        elif op == 'real':
            v = v.real()
        elif op == 'imag':
            v = v.imag()
        elif op == 'abs':
            v = abs(v)
        elif op == 'norm':
            v = v.norm()
        else:
            raise NotImplementedError

        result_is_qqbar = arg_is_qqbar
        if op in ('real', 'imag', 'abs', 'norm'):
            result_is_qqbar = False
        if result_is_qqbar != is_qqbar:
            # The following version is not safe with respect to caching;
            # with the current sage_input.py, anything that gets entered
            # into the cache must be safe at all coercion levels.
#             if is_qqbar and not coerce:
#                 v = sib.name('QQbar')(v)
#             if not is_qqbar and coerce != 2:
#                 v = sib.name('AA')(v)
            v = sib.name('QQbar' if is_qqbar else 'AA')(v)

        return (v, True)

    def is_complex(self):
        r"""
        Return whether or not this element is complex. Note that this is a data
        type check, and triggers no computations -- if it returns False, the
        element might still be real, it just doesn't know it yet.

        EXAMPLE::

            sage: t = AA(sqrt(2))
            sage: s = (-t)._descr
            sage: s
            <class 'sage.rings.qqbar.ANUnaryExpr'>
            sage: s.is_complex()
            False
            sage: QQbar(-sqrt(2))._descr.is_complex()
            True
        """
        return self._complex

    def _interval_fast(self, prec):
        r"""
        Calculate an approximation to this ``ANUnaryExpr`` object in an interval field of precision ``prec``.

        EXAMPLE::

            sage: t = AA(sqrt(2))
            sage: s = (-t)._descr
            sage: s
            <class 'sage.rings.qqbar.ANUnaryExpr'>
            sage: s._interval_fast(150)
            -1.414213562373095048801688724209698078569671876?
        """
        op = self._op

        v = self._arg._interval_fast(prec)

        if not is_ComplexIntervalFieldElement(v):
            self._complex = False

        if op == '-':
            return -v

        if op == '~':
            return ~v

        if op == 'conjugate':
            if is_ComplexIntervalFieldElement(v):
                return v.conjugate()
            else:
                return v

        self._complex = False

        if op == 'real':
            if is_ComplexIntervalFieldElement(v):
                return v.real()
            else:
                return v

        if op == 'imag':
            if is_ComplexIntervalFieldElement(v):
                return v.imag()
            else:
                return RealIntervalField(prec)(0)

        if op == 'abs':
            return abs(v)

        if op == 'norm':
            if is_ComplexIntervalFieldElement(v):
                return v.norm()
            else:
                return v.square()

        raise NotImplementedError

    def exactify(self):
        r"""
        Trigger exact computation of self.

        EXAMPLE::

            sage: v = (-QQbar(sqrt(2)))._descr
            sage: type(v)
            <class 'sage.rings.qqbar.ANUnaryExpr'>
            sage: v.exactify()
            -a where a^2 - 2 = 0 and a in 1.414213562373095?
        """
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
        r"""
        Initialize this ANBinaryExpr.

        EXAMPLE::

            sage: t = QQbar(sqrt(2)) + QQbar(sqrt(3)); type(t._descr) # indirect doctest
            <class 'sage.rings.qqbar.ANBinaryExpr'>
        """
        self._left = left
        self._right = right
        self._op = op
        self._complex = True

    def __reduce__(self):
        """
        Add customized pickling support.

        EXAMPLES::

            sage: t = QQbar(sqrt(2)) + QQbar(sqrt(3)); type(t._descr)
            <class 'sage.rings.qqbar.ANBinaryExpr'>
            sage: loads(dumps(t)) == QQbar(sqrt(2)) + QQbar(sqrt(3))
            True
        """
        return (ANBinaryExpr, (self._left, self._right, self._op))

    def handle_sage_input(self, sib, coerce, is_qqbar):
        r"""
        Produce an expression which will reproduce this value when evaluated,
        and an indication of whether this value is worth sharing (always
        True for ``ANBinaryExpr``).

        EXAMPLES::

            sage: sage_input(2 + sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            2 + AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(sqrt(AA(2)) + 2, verify=True)
            # Verified
            R.<x> = AA[]
            AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951))) + 2
            sage: sage_input(2 - sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            2 - AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(2 / sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            2/AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(2 + (-1*sqrt(AA(2))), verify=True)
            # Verified
            R.<x> = AA[]
            2 - AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: sage_input(2*sqrt(AA(2)), verify=True)
            # Verified
            R.<x> = AA[]
            2*AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            sage: rt2 = sqrt(AA(2))
            sage: one = rt2/rt2
            sage: n = one+3
            sage: sage_input(n)
            R.<x> = AA[]
            v = AA.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951)))
            v/v + 3
            sage: one == 1
            True
            sage: sage_input(n)
            1 + AA(3)
            sage: rt3 = QQbar(sqrt(3))
            sage: one = rt3/rt3
            sage: n = sqrt(AA(2))+one
            sage: one == 1
            True
            sage: sage_input(n)
            R.<x> = AA[]
            QQbar.polynomial_root(AA.common_polynomial(x^2 - 2), RIF(RR(1.4142135623730949), RR(1.4142135623730951))) + 1
            sage: from sage.rings.qqbar import *
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: binexp = ANBinaryExpr(AA(3), AA(5), operator.mul)
            sage: binexp.handle_sage_input(sib, False, False)
            ({binop:* {atomic:3} {call: {atomic:AA}({atomic:5})}}, True)
            sage: binexp.handle_sage_input(sib, False, True)
            ({call: {atomic:QQbar}({binop:* {atomic:3} {call: {atomic:AA}({atomic:5})}})}, True)
        """
        arg1 = self._left
        arg2 = self._right
        op = self._op

        # We want 2+QQbar.zeta(3) and QQbar.zeta(3)+2, not
        # QQbar(2)+QQbar.zeta(3). So we want to pass coerced=True to
        # an argument if it is rational (but if both arguments are
        # rational, we only want to set it for one of them).

        arg1_coerced = False
        arg2_coerced = False

        if isinstance(arg1._descr, ANRational):
            arg1_coerced = True
        elif isinstance(arg2._descr, ANRational):
            arg2_coerced = True

        arg1_is_qqbar = arg1.parent() is QQbar
        arg2_is_qqbar = arg2.parent() is QQbar

        result_is_qqbar = \
            (arg1_is_qqbar and not arg1_coerced) or \
            (arg2_is_qqbar and not arg2_coerced)

        v1 = sib(arg1, arg1_coerced)
        v2 = sib(arg2, arg2_coerced)

        if op is operator.add:
            v = sib.sum([v1, v2], simplify=True)
        elif op is operator.sub:
            v = sib.sum([v1, -v2], simplify=True)
        elif op is operator.mul:
            v = sib.prod([v1, v2], simplify=True)
        elif op is operator.div:
            v = v1 / v2
        else:
            raise RuntimeError("op is {}".format(op))

        if result_is_qqbar != is_qqbar:
            # The following version is not safe with respect to caching;
            # with the current sage_input.py, anything that gets entered
            # into the cache must be safe at all coercion levels.
#             if is_qqbar and not coerce:
#                 v = sib.name('QQbar')(v)
#             if not is_qqbar and coerce != 2:
#                 v = sib.name('AA')(v)
            v = sib.name('QQbar' if is_qqbar else 'AA')(v)

        return (v, True)

    def is_complex(self):
        r"""
        Whether this element is complex. Does not trigger exact computation, so
        may return True even if the element is real.

        EXAMPLE::

            sage: x = (QQbar(sqrt(-2)) / QQbar(sqrt(-5)))._descr
            sage: x.is_complex()
            True
        """
        return self._complex

    def _interval_fast(self, prec):
        r"""
        Calculate an approximation to self in an interval field of precision prec.

        EXAMPLE::

            sage: x = (QQbar(sqrt(-2)) / QQbar(sqrt(-5)))._descr
            sage: y= x._interval_fast(64); y
            0.632455532033675867?
            sage: y.parent()
            Complex Interval Field with 64 bits of precision
        """
        op = self._op

        lv = self._left._interval_fast(prec)
        rv = self._right._interval_fast(prec)

        if not (is_ComplexIntervalFieldElement(lv) or is_ComplexIntervalFieldElement(rv)):
            self._complex = False

        return op(lv, rv)

    def exactify(self):
        """
        TESTS::

            sage: rt2c = QQbar.zeta(3) + AA(sqrt(2)) - QQbar.zeta(3)
            sage: rt2c.exactify()

        We check to make sure that this method still works even. We
        do this by increasing the recursion level at each step and
        decrease it before we return::

            sage: import sys; sys.getrecursionlimit()
            1000
            sage: s = SymmetricFunctions(QQ).schur()
            sage: a=s([3,2]).expand(8)(flatten([[QQbar.zeta(3)^d for d in range(3)], [QQbar.zeta(5)^d for d in range(5)]]))
            sage: a.exactify(); a # long time
            0
            sage: sys.getrecursionlimit()
            1000

        """
        import sys
        old_recursion_limit = sys.getrecursionlimit()
        sys.setrecursionlimit(old_recursion_limit + 10)
        try:
            left = self._left
            right = self._right
            left.exactify()
            right.exactify()
            gen = left._exact_field().union(right._exact_field())
            left_value = gen(left._exact_value())
            right_value = gen(right._exact_value())

            value = self._op(left_value, right_value)

            if gen.is_trivial():
                return ANRational(value)
            else:
                return ANExtensionElement(gen, value)
        finally:
            sys.setrecursionlimit(old_recursion_limit)

# These are the functions used to add, subtract, multiply, and divide
# algebraic numbers. Basically, we try to compute exactly if both
# arguments are already known to be in the same number field. Otherwise
# we fall back to floating-point computation to be backed up by exact
# symbolic computation only as required.

def an_binop_rational(a, b, op):
    r"""
    Used to add, subtract, multiply or divide algebraic numbers.

    Used when both are actually rational.

    EXAMPLES::

        sage: from sage.rings.qqbar import an_binop_rational
        sage: f = an_binop_rational(QQbar(2), QQbar(3/7), operator.add)
        sage: f
        17/7
        sage: type(f)
        <class 'sage.rings.qqbar.ANRational'>

        sage: f = an_binop_rational(QQbar(2), QQbar(3/7), operator.mul)
        sage: f
        6/7
        sage: type(f)
        <class 'sage.rings.qqbar.ANRational'>
    """
    return ANRational(op(a._descr._value, b._descr._value))

def an_binop_expr(a, b, op):
    r"""
    Add, subtract, multiply or divide algebraic numbers represented as
    binary expressions.

    INPUT:

    - ``a``, ``b`` -- two elements

    - ``op`` -- an operator

    EXAMPLE::

        sage: a = QQbar(sqrt(2)) + QQbar(sqrt(3))
        sage: b = QQbar(sqrt(3)) + QQbar(sqrt(5))
        sage: type(a._descr); type(b._descr)
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        sage: from sage.rings.qqbar import an_binop_expr
        sage: x = an_binop_expr(a, b, operator.add); x
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        sage: x.exactify()
        -6/7*a^7 + 2/7*a^6 + 71/7*a^5 - 26/7*a^4 - 125/7*a^3 + 72/7*a^2 + 43/7*a - 47/7 where a^8 - 12*a^6 + 23*a^4 - 12*a^2 + 1 = 0 and a in 3.12580...?

        sage: a = QQbar(sqrt(2)) + QQbar(sqrt(3))
        sage: b = QQbar(sqrt(3)) + QQbar(sqrt(5))
        sage: type(a._descr)
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        sage: x = an_binop_expr(a, b, operator.mul); x
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        sage: x.exactify()
        2*a^7 - a^6 - 24*a^5 + 12*a^4 + 46*a^3 - 22*a^2 - 22*a + 9 where a^8 - 12*a^6 + 23*a^4 - 12*a^2 + 1 = 0 and a in 3.1258...?
    """
    return ANBinaryExpr(a, b, op)

def an_binop_element(a, b, op):
    r"""
    Add, subtract, multiply or divide two elements represented as elements of
    number fields.

    EXAMPLES::

        sage: sqrt2 = QQbar(2).sqrt()
        sage: sqrt3 = QQbar(3).sqrt()
        sage: sqrt5 = QQbar(5).sqrt()

        sage: a = sqrt2 + sqrt3; a.exactify()
        sage: b = sqrt3 + sqrt5; b.exactify()
        sage: type(a._descr)
        <class 'sage.rings.qqbar.ANExtensionElement'>
        sage: from sage.rings.qqbar import an_binop_element
        sage: an_binop_element(a, b, operator.add)
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        sage: an_binop_element(a, b, operator.sub)
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        sage: an_binop_element(a, b, operator.mul)
        <class 'sage.rings.qqbar.ANBinaryExpr'>
        sage: an_binop_element(a, b, operator.div)
        <class 'sage.rings.qqbar.ANBinaryExpr'>

    The code tries to use existing unions of number fields::

        sage: sqrt17 = QQbar(17).sqrt()
        sage: sqrt19 = QQbar(19).sqrt()
        sage: a = sqrt17 + sqrt19
        sage: b = sqrt17 * sqrt19 - sqrt17 + sqrt19 * (sqrt17 + 2)
        sage: b, type(b._descr)
        (40.53909377268655?, <class 'sage.rings.qqbar.ANBinaryExpr'>)
        sage: a.exactify()
        sage: b = sqrt17 * sqrt19 - sqrt17 + sqrt19 * (sqrt17 + 2)
        sage: b, type(b._descr)
        (40.53909377268655?, <class 'sage.rings.qqbar.ANExtensionElement'>)
    """
    ad = a._descr
    bd = b._descr
    adg = ad.generator()
    bdg = bd.generator()
    if adg == qq_generator or adg == bdg:
        return ANExtensionElement(bdg, op(ad._value, bd._value))

    if bdg == qq_generator:
        return ANExtensionElement(adg, op(ad._value, bd._value))

    if adg in bdg._unions or bdg in adg._unions:
        p = bdg._unions[adg] if adg in bdg._unions else adg._unions[bdg]
        p = p.parent
        adg2 = adg.super_poly(p)
        bdg2 = bdg.super_poly(p)
        av = ad._value.polynomial()(adg2)
        bv = bd._value.polynomial()(bdg2)
        v = op(av, bv)
        return ANExtensionElement(p, op(av, bv))

    adg2 = adg.super_poly(bdg)
    if adg2 is not None:
        av = ad._value.polynomial()(adg2)
        return ANExtensionElement(bdg, op(av, bd._value))

    bdg2 = bdg.super_poly(adg)
    if bdg2 is not None:
        bv = bd._value.polynomial()(bdg2)
        return ANExtensionElement(adg, op(ad._value, bv))

    return ANBinaryExpr(a, b, op)

# instanciation of the multimethod dispatch
_binop_algo[ANRational, ANRational] = an_binop_rational
_binop_algo[ANRational, ANExtensionElement] = \
_binop_algo[ANExtensionElement, ANRational] = \
_binop_algo[ANExtensionElement, ANExtensionElement ] = an_binop_element
for t1 in (ANRational, ANRoot, ANExtensionElement, ANUnaryExpr, ANBinaryExpr):
    for t2 in (ANUnaryExpr, ANBinaryExpr, ANRoot):
        _binop_algo[t1, t2] = _binop_algo[t2, t1] = an_binop_expr

qq_generator = AlgebraicGenerator(QQ, ANRoot(AAPoly.gen() - 1, RIF(1)))


def _init_qqbar():
    """
    This code indirectly uses a huge amount of sage, despite the fact
    that qqbar is imported rather early on in the sage loading. This function
    is called at the end of sage.all.

    EXAMPLE::

        sage: sage.rings.qqbar.QQbar_I_generator # indirect doctest
        Number Field in I with defining polynomial x^2 + 1 with a in 1*I
    """
    global ZZX_x, AA_0, QQbar_I, AA_hash_offset, QQbar_hash_offset, QQbar_I_generator, QQbar_I_nf
    global QQ_0, QQ_1, QQ_1_2, QQ_1_4, RR_1_10

    RR_1_10 = RR(1)/10
    QQ_0 = QQ.zero()
    QQ_1 = QQ.one()
    QQ_1_2 = QQ((1,2))
    QQ_1_4 = QQ((1,4))

    AA_0 = AA.zero()

    QQbar_I_nf = QuadraticField(-1, 'I', embedding=CC.gen())
    QQbar_I_generator = AlgebraicGenerator(QQbar_I_nf, ANRoot(AAPoly.gen()**2 + 1, CIF(0, 1)))
    QQbar_I = AlgebraicNumber(ANExtensionElement(QQbar_I_generator, QQbar_I_nf.gen()))

    AA_hash_offset = AA(~ZZ(123456789))

    QQbar_hash_offset = AlgebraicNumber(ANExtensionElement(QQbar_I_generator, ~ZZ(123456789) + QQbar_I_nf.gen()/ZZ(987654321)))

    ZZX_x = ZZ['x'].gen()

# This is used in the _algebraic_ method of the golden_ratio constant,
# in sage/symbolic/constants.py
AA_golden_ratio = None

def get_AA_golden_ratio():
    r"""
    Return the golden ratio as an element of the algebraic real field. Used by
    :meth:`sage.symbolic.constants.golden_ratio._algebraic_`.

    EXAMPLE::

        sage: AA(golden_ratio) # indirect doctest
        1.618033988749895?
    """
    global AA_golden_ratio
    if AA_golden_ratio is None:
        AA_golden_ratio_nf = NumberField(ZZX_x**2 - ZZX_x - 1, 'phi')
        AA_golden_ratio_generator = AlgebraicGenerator(AA_golden_ratio_nf, ANRoot(AAPoly.gen()**2 - AAPoly.gen() - 1, RIF(1.618, 1.6181)))
        AA_golden_ratio = AlgebraicReal(ANExtensionElement(AA_golden_ratio_generator, AA_golden_ratio_nf.gen()))
    return AA_golden_ratio

class ANRootOfUnity(ANExtensionElement):
    r"""
    Deprecated class to support unpickling

    TESTS::

        sage: from sage.rings.qqbar import ANRootOfUnity
        sage: ANRootOfUnity(1/5, 3/2)
        doctest:...: DeprecationWarning: ANRootOfUnity is deprecated
        See http://trac.sagemath.org/19954 for details.
        3/2*zeta5 where zeta5^4 + zeta5^3 + zeta5^2 + zeta5 + 1 = 0
        and zeta5 in 0.3090169943749474? + 0.9510565162951536?*I
    """
    def __new__(self, a, b):
        from sage.misc.superseded import deprecation
        deprecation(19954, "ANRootOfUnity is deprecated")
        descr = QQbar.zeta(a.denominator())._descr
        generator = descr._generator
        value = b * descr._value ** (a.numerator())
        return ANExtensionElement(generator, value)
