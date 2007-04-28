r"""
Computation of Frobenius matrix on Monsky-Washnitzer cohomology.

The most interesting functions to be exported here are matrix_of_frobenius()
and adjusted_prec(), and matrix_of_frobenius_alternate().

Currently this code is limited to the case $p \geq 5$ (no $GF(p^n)$
for $n > 1$), and only handles the elliptic curve case (not more
general hyperelliptic curves).

REFERENCES:
  -- Kedlaya, K., ``Counting points on hyperelliptic curves using
     Monsky-Washnitzer cohomology'', J. Ramanujan Math. Soc. 16 (2001)
     no 4, 323--338
  -- Edixhoven, B., ``Point counting after Kedlaya'', EIDMA-Stieltjes graduate
     course, Lieden (lecture notes?).
  -- Harvey, D. ``Kedlaya's algorithm in larger characteristic'',
     http://arxiv.org/abs/math.NT/0610973

AUTHORS:
    -- David Harvey and Robert Bradshaw (initial code developed at the 2006
       MSRI graduate workshop, working with Jennifer Balakrishnan and Liang
       Xiao)
    -- David Harvey (Aug/Sep 2006): cleaned up, rewrote some chunks, lots
       more documentation, added Newton iteration method, added more complete
       "trace trick", integrated better into SAGE.
    -- David Harvey (Feb 2007): added algorithm with sqrt(p) complexity
    -- Robert Bradshaw (Mar 2007): keep track of exact form in reduction algorithms
    -- Robert Bradshaw (Apr 2007): generalization to hyperelliptic curves

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 Robert Bradshaw <robertwb@math.washington.edu>
#                     2006 David Harvey <dmharvey@math.harvard.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.all import Integers, Integer, PolynomialRing, is_Polynomial, PowerSeriesRing, Rationals, Rational, LaurentSeriesRing
from sage.algebras.all import FreeAlgebra, FreeAlgebraQuotient, Algebra, AlgebraElement
from sage.modules.module import Module
from sage.structure.element import ModuleElement
from sage.algebras.free_algebra_quotient_element import FreeAlgebraQuotientElement
from sage.matrix.all import matrix
from sage.modules.all import vector
from sage.rings.ring import CommutativeAlgebra
from sage.structure.element import CommutativeAlgebraElement
from sage.matrix.matrix_space import MatrixSpace

from sage.rings.arith import binomial, floor
from sage.misc.functional import log, ceil, sqrt

from ell_generic import is_EllipticCurve


class SpecialCubicQuotientRing(CommutativeAlgebra):
  r""" Specialised class for representing the quotient ring
  $R[x,T]/(T - x^3 - ax - b)$, where $R$ is an arbitrary commutative base ring
  (in which 2 and 3 are invertible), $a$ and $b$ are elements of that ring.

  Polynomials are represented internally in the form $p_0 + p_1 x + p_2 x^2$
  where the $p_i$ are polynomials in $T$. Multiplication of polynomials
  always reduces high powers of $x$ (i.e. beyond $x^2$) to powers of $T$.

  Hopefully this ring is faster than a general quotient ring because it uses
  the special structure of this ring to speed multiplication (which is the
  dominant operation in the frobenius matrix calculation). I haven't actually
  tested this theory though...

  TODO:
    -- Eventually we will want to run this in characteristic 3, so we need to:
       (a) Allow Q(x) to contain an $x^2$ term, and
       (b) Remove the requirement that 3 be invertible. Currently this is used
           in the Toom-Cook algorithm to speed multiplication.

  EXAMPLES:
      sage: B.<t> = PolynomialRing(Integers(125))
      sage: R = monsky_washnitzer.SpecialCubicQuotientRing(t^3 - t + B(1/4))
      sage: R
      SpecialCubicQuotientRing over Ring of integers modulo 125 with polynomial T = x^3 + 124*x + 94

    Get generators:
      sage: x, T = R.gens()
      sage: x
      (0) + (1)*x + (0)*x^2
      sage: T
      (T) + (0)*x + (0)*x^2

    Coercions:
      sage: R(7)
      (7) + (0)*x + (0)*x^2

    Create elements directly from polynomials:
      sage: A, z = R.poly_ring().objgen()
      sage: A
      Univariate Polynomial Ring in T over Ring of integers modulo 125
      sage: R.create_element(z^2, z+1, 3)
      (T^2) + (T + 1)*x + (3)*x^2

    Some arithmetic:
      sage: x^3
      (T + 31) + (1)*x + (0)*x^2
      sage: 3 * x**15 * T**2 + x - T
      (3*T^7 + 90*T^6 + 110*T^5 + 20*T^4 + 58*T^3 + 26*T^2 + 124*T) + (15*T^6 + 110*T^5 + 35*T^4 + 63*T^2 + 1)*x + (30*T^5 + 40*T^4 + 8*T^3 + 38*T^2)*x^2

    Retrieve coefficients (output is zero-padded):
      sage: x^10
      (3*T^2 + 61*T + 8) + (T^3 + 93*T^2 + 12*T + 40)*x + (3*T^2 + 61*T + 9)*x^2
      sage: (x^10).coeffs()
      [[8, 61, 3, 0], [40, 12, 93, 1], [9, 61, 3, 0]]

    TODO: write an example checking multiplication of these polynomials
    against SAGE's ordinary quotient ring arithmetic. I can't seem to get
    the quotient ring stuff happening right now...
  """

  def __init__(self, Q, laurent_series = False):
    """
    Constructor.

    INPUT:
        Q -- a polynomial of the form Q(x) = x^3 + ax + b, where a, b
             belong to a ring in which 2, 3 are invertible.
        laurent_series -- whether or not to allow negative powers of T (default=False)
    """
    if not is_Polynomial(Q):
      raise TypeError, "Q (=%s) must be a polynomial" % Q

    if Q.degree() != 3 or not Q[2].is_zero():
      raise ValueError, "Q (=%s) must be of the form x^3 + ax + b" % Q

    base_ring = Q.parent().base_ring()

    if not base_ring(6).is_unit():
      raise ArithmeticError, \
            "2 and 3 must be invertible in the coefficient ring (=%s) of Q" % \
            base_ring

    CommutativeAlgebra.__init__(self, base_ring)
    self._a = Q[1]
    self._b = Q[0]
    if laurent_series:
      self._poly_ring = LaurentSeriesRing(base_ring, 'T')    # R[T]
    else:
      self._poly_ring = PolynomialRing(base_ring, 'T')    # R[T]
    self._poly_generator = self._poly_ring.gen(0)    # the generator T

    # Precompute a matrix that is used in the Toom-Cook multiplication.
    # This is where we need 2 and 3 invertible.

    # (a good description of Toom-Cook is online at:
    # http://www.gnu.org/software/gmp/manual/html_node/Toom-Cook-3-Way-Multiplication.html )

    self._speedup_matrix = \
        (matrix(Integers(), 3, 3, [2, 4, 8,
                                   1, 1, 1,
                                   8, 4, 2])**(-1)
         ).change_ring(base_ring).list()

    # todo: get rid of the next line.
    # It's currently necessary in __mul__ for the elements of this list to
    # be in poly_ring. But in an ideal world they should be elements of the
    # base_ring. Unfortunately SAGE's binop stuff is screwed up somewhere,
    # and __mul__ breaks for certain base rings.
    self._speedup_matrix = [self._poly_ring(x) for x in self._speedup_matrix]


  def __repr__(self):
    return "SpecialCubicQuotientRing over %s with polynomial T = %s" % \
           (self.base_ring(), PolynomialRing(self.base_ring(), 'x')(
                                                [self._b, self._a, 0, 1]))

  def poly_ring(self):
    """ Return the underlying polynomial ring in T. """
    return self._poly_ring

  def gens(self):
    """ Return a list [x, T] where x and T are the generators of the ring
    (as element *of this ring*).

    NOTE:
        I have no idea if this is compatible with the usual SAGE
        "gens" interface.
    """
    return [SpecialCubicQuotientRingElement(self,
                 self._poly_ring(0), self._poly_ring(1), self._poly_ring(0),
                 check=False),
            SpecialCubicQuotientRingElement(self,
                 self._poly_generator, self._poly_ring(0), self._poly_ring(0),
                 check=False)]

  def create_element(self, p0, p1, p2, check=True):
    """ Creates the element $p_0 + p_1*x + p_2*x^2$, where pi's are
    polynomials in T.

    INPUT:
        p0, p1, p2 -- coefficients; must be coercable into poly_ring()
        check -- bool (default True): whether to carry out coercion

    """
    return SpecialCubicQuotientRingElement(self, p0, p1, p2, check)


  # todo: work out exactly where the following coercions are being
  # called, and exactly how they should operate

  def __call__(self, value):
    return self._coerce_(value)


  def _coerce_impl(self, value):
    # todo: I don't understand why the direct _poly_ring.__call__()
    # doesn't work....

    # try coercing to base ring
    try:
      value = self.base_ring()._coerce_(value)
      value = self._poly_ring._coerce_(value)

    except TypeError:
      # try coercing to underlying polynomial ring
      value = self._poly_ring._coerce_(value)

    return SpecialCubicQuotientRingElement(self,
                 value, self._poly_ring(0), self._poly_ring(0), check=False)



class SpecialCubicQuotientRingElement(CommutativeAlgebraElement):
  """ An element of a SpecialCubicQuotientRing. """

  def __init__(self, parent, p0, p1, p2, check=True):
    """ Constructs the element $p_0 + p_1*x + p_2*x^2$, where pi's are
    polynomials in T.

    INPUT:
        parent -- a SpecialCubicQuotientRing
        p0, p1, p2 -- coefficients; must be coercable into parent.poly_ring()
        check -- bool (default True): whether to carry out coercion

    """
    if not isinstance(parent, SpecialCubicQuotientRing):
      raise TypeError, \
              "parent (=%s) must be a SpecialCubicQuotientRing" % parent

    CommutativeAlgebraElement.__init__(self, parent)

    if check:
      poly_ring = parent.poly_ring()
      p0 = poly_ring(p0)
      p1 = poly_ring(p1)
      p2 = poly_ring(p2)

    self._triple = (p0, p1, p2)


  def coeffs(self):
    """ Returns list of three lists of coefficients, corresponding to the
    $x^0$, $x^1$, $x^2$ coefficients. The lists are zero padded to the same
    length. The list entries belong to the base ring.

    """
    coeffs = [column.coeffs() for column in self._triple]
    degree = max([len(x) for x in coeffs])
    base_ring = self.parent().base_ring()
    for column in coeffs:
      column.extend([base_ring(0)] * (degree - len(column)))
    return coeffs


  def _repr_(self):
    return "(%s) + (%s)*x + (%s)*x^2" % self._triple


  def _latex_(self):
    return "(%s) + (%s)x + (%s)x^2" % \
           tuple([column._latex_() for column in self._triple])


  def _add_(self, other):
    return SpecialCubicQuotientRingElement(self.parent(),
                    self._triple[0] + other._triple[0],
                    self._triple[1] + other._triple[1],
                    self._triple[2] + other._triple[2],
                    check=False)


  def _sub_(self, other):
    return SpecialCubicQuotientRingElement(self.parent(),
                    self._triple[0] - other._triple[0],
                    self._triple[1] - other._triple[1],
                    self._triple[2] - other._triple[2],
                    check=False)


  def shift(self, n):
    """ Returns this element multiplied by $T^n$. """
    return SpecialCubicQuotientRingElement(self.parent(),
                    self._triple[0].shift(n),
                    self._triple[1].shift(n),
                    self._triple[2].shift(n),
                    check=False)


  def scalar_multiply(self, scalar):
    """ Multiplies this element by a "scalar".

    i.e. just multiply each coefficient of $x^j$ by the scalar.

    INPUT:
        scalar -- either an element of base_ring,
                  or an element of poly_ring.
    """
    # try to coerce scalar into underlying polynomial ring
    # todo: why bother with this coercion here?
    scalar = self.parent()._poly_ring(scalar)
    return SpecialCubicQuotientRingElement(self.parent(),
                                           scalar * self._triple[0],
                                           scalar * self._triple[1],
                                           scalar * self._triple[2],
                                           check=False)


  def square(self):
    # todo: we can maybe do this faster with the toom-cook squaring algorithm,
    # should be particularly effective for very large input
    return self * self


  def _mul_(self, other):
    # todo: cache results of toom-cook splitting, i.e. often we multiply
    # the same polynomial by a bunch of other things, and we can save
    # on part of the repetitve work

    # todo: if the degree is small, perhaps just use the naive
    # algorithm instead?

    # todo: I did a bit of simple profiling on this code and it looks to
    # be spending WAY too much time in the additions/subtraction/scalar
    # operations. It should be spending most of its time in the polynomial
    # multiplications. So this all needs to be revisited when the underlying
    # polynomial code has been pyrexified and all that.

    if not isinstance(other, SpecialCubicQuotientRingElement):
      return self.scalar_multiply(other)

    # Here we do Toom-Cook three-way multiplication, which reduces the
    # naive 9 polynomial multiplications to only 5 polynomial multiplications.

    a0, a1, a2 = self._triple
    b0, b1, b2 = other._triple
    M = self.parent()._speedup_matrix

    p0 = a0 * b0
    p1 = (a0 + 2*a1 + 4*a2) * (b0 + 2*b1 + 4*b2)
    p2 = (a0 + a1 + a2) * (b0 + b1 + b2)
    p3 = (4*a0 + 2*a1 + a2) * (4*b0 + 2*b1 + b2)
    p4 = a2 * b2

    q1 = p1 - p0 - 16*p4
    q2 = p2 - p0 - p4
    q3 = p3 - 16*p0 - p4

    c0 = p0
    c1 = M[0]*q1 + M[1]*q2 + M[2]*q3
    c2 = M[3]*q1 + M[4]*q2 + M[5]*q3
    c3 = M[6]*q1 + M[7]*q2 + M[8]*q3
    c4 = p4

    # Now the product is c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4.
    # We need to reduce mod y = x^3 + ax + b and return result.

    parent = self.parent()
    T = parent._poly_generator
    b = parent._b
    a = parent._a

    # todo: These lines are necessary to get binop stuff working
    # for certain base rings, e.g. when we compute b*c3 in the
    # final line. They shouldn't be necessary. Need to fix this
    # somewhere else in SAGE.
    a = parent._poly_ring(a)
    b = parent._poly_ring(b)

    return SpecialCubicQuotientRingElement(parent,
                                           -b*c3 + c0 + c3*T,
                                           -b*c4 - a*c3 + c1 + c4*T,
                                           -a*c4 + c2,
                                           check=False)

def transpose_list(input):
    """
    INPUT:
        input -- a list of lists, each list of the same length

    OUTPUT:
        output -- a list of lists such that output[i][j] = input[j][i]

    """

    h = len(input)
    w = len(input[0])

    output = []
    for i in range(w):
        row = []
        for j in range(h):
            row.append(input[j][i])
        output.append(row)
    return output



def helper_matrix(Q):
    """
    Computes the (constant) matrix used to calculate the linear combinations
    of the $d(x^i y^j)$ needed to eliminate the negative powers of $y$
    in the cohomology (i.e. in reduce_negative()).

    INPUT:
        Q -- cubic polynomial
    """

    a = Q[1]
    b = Q[0]

    # Discriminant (should be invertible for a curve of good reduction)
    D = 4*a**3 + 27*b**2

    # This is the inverse of the matrix
    #   [  a,  -3b,    0 ]
    #   [  0,  -2a,  -3b ]
    #   [  3,    0,  -2a ]

    return (1/D) * matrix(Q.base_ring(), 3, 3,
                          [  4*a**2 , -6*b*a  , 9*b**2,
                             -9*b   , -2*a**2 , 3*b*a,
                              6*a   , -9*b    , -2*a**2 ])


def lift(x):
    r"""
    Tries to call x.lift(), presumably from the p-adics to ZZ.

    If this fails, it assumes the input is a power series, and tries to
    lift it to a power series over QQ.

    This function is just a very kludgy solution to the problem of trying
    to make the reduction code (below) work over both Zp and Zp[[t]].
    """
    try:
        return x.lift()
    except AttributeError:
        return PowerSeriesRing(Rationals(), "t")(x.list(), x.prec())



def reduce_negative(Q, p, coeffs, offset, exact_form=None):
    """
    Applies cohomology relations to incorporate negative powers of $y$
    into the $y^0$ term.

    INPUT:
        p -- prime
        Q -- cubic polynomial
        coeffs -- list of length 3 lists. The i^th list [a, b, c]
                  represents $y^{2(i - offset)} (a + bx + cx^2) dx/y$.
        offset -- nonnegative integer

    OUTPUT:
        The reduction is performed in-place. The output is placed in
        coeffs[offset]. Note that coeffs[i] will be meaningless for
        i < offset after this function is finished.

    EXAMPLE:
        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: coeffs = [[10, 15, 20], [1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_negative(Q, 5, coeffs, 3)
        sage: coeffs[3]
         [28, 52, 9]

        sage: R.<x> = Integers(7^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: coeffs = [[7, 14, 21], [1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_negative(Q, 7, coeffs, 3)
        sage: coeffs[3]
         [245, 332, 9]

    """

    m = helper_matrix(Q).list()
    base_ring = Q.base_ring()
    next_a = coeffs[0]

    if exact_form is not None:
        x = exact_form.parent().gen(0)
        y = exact_form.parent().base_ring().gen(0)

    try:
        three_j_plus_5 = 5 - base_ring(6*offset)
        three_j_plus_7 = 7 - base_ring(6*offset)
        six = base_ring(6)

        for i in range(0, offset):

            j = 2*(i-offset)
            a = next_a
            next_a = coeffs[i+1]

            # todo: the following divisions will sometimes involve
            # a division by (a power of) p. In all cases, we know (from
            # Kedlaya's estimates) that the answer should be p-integral.
            # However, since we're working over $Z/p^k Z$, we're not allowed
            # to "divide by p". So currently we lift to Q, divide, and coerce
            # back. Eventually, when pAdicInteger is implemented, and plays
            # nicely with pAdicField, we should reimplement this stuff
            # using pAdicInteger.

            if (p.divides(j+1)):
               # need to lift here to perform the division
               a[0] = base_ring(lift(a[0]) / (j+1))
               a[1] = base_ring(lift(a[1]) / (j+1))
               a[2] = base_ring(lift(a[2]) / (j+1))
            else:
               j_plus_1_inv = ~base_ring(j+1)
               a[0] = a[0] * j_plus_1_inv
               a[1] = a[1] * j_plus_1_inv
               a[2] = a[2] * j_plus_1_inv

            c1 = m[3]*a[0] + m[4]*a[1] + m[5]*a[2]
            c2 = m[6]*a[0] + m[7]*a[1] + m[8]*a[2]
            next_a[0] = next_a[0] - three_j_plus_5 * c1
            next_a[1] = next_a[1] - three_j_plus_7 * c2

            three_j_plus_7 = three_j_plus_7 + six
            three_j_plus_5 = three_j_plus_5 + six

            if exact_form is not None:
                c0 = m[0]*a[0] + m[1]*a[1] + m[2]*a[2]
                exact_form += (c0 + c1*x + c2 * x**2) * y**(j+1)


    except NotImplementedError:
        raise NotImplementedError, \
            "It looks like you've found a non-integral matrix of Frobenius! " \
            "(Q=%s, p=%s)\nTime to write a paper." % (Q, p)

    coeffs[int(offset)] = next_a

    return exact_form



def reduce_positive(Q, p, coeffs, offset, exact_form=None):
    """
    Applies cohomology relations to incorporate positive powers of $y$
    into the $y^0$ term.

    INPUT:
        Q -- cubic polynomial
        coeffs -- list of length 3 lists. The i^th list [a, b, c]
                  represents $y^{2(i - offset)} (a + bx + cx^2) dx/y$.
        offset -- nonnegative integer

    OUTPUT:
        The reduction is performed in-place. The output is placed in
        coeffs[offset]. Note that coeffs[i] will be meaningless for
        i > offset after this function is finished.

    EXAMPLE:
        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)

        sage: coeffs = [[1, 2, 3], [10, 15, 20]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_positive(Q, 5, coeffs, 0)
        sage: coeffs[0]
         [16, 102, 88]

        sage: coeffs = [[9, 8, 7], [10, 15, 20]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_positive(Q, 5, coeffs, 0)
        sage: coeffs[0]
         [24, 108, 92]

    """

    base_ring = Q.base_ring()
    next_a = coeffs[len(coeffs) - 1]

    Qa = Q[1]
    Qb = Q[0]

    A = 2*Qa
    B = 3*Qb

    offset = Integer(offset)


    if exact_form is not None:
        x = exact_form.parent().gen(0)
        y = exact_form.parent().base_ring().gen(0)

    for i in range(len(coeffs)-1, offset, -1):
        j = 2*(i-offset) - 2
        a = next_a
        next_a = coeffs[i-1]

        a[0] = a[0] - Qa*a[2]/3   # subtract d(y^j + 3)
        if exact_form is not None:
            exact_form += Q.base_ring()(a[2].lift() / (3*j+9)) * y**(j+3)

        # todo: see comments about pAdicInteger in reduceNegative()

        # subtract off c1 of d(x y^j + 1), and
        if p.divides(3*j + 5):
            c1 = base_ring(lift(a[0]) / (3*j + 5))
        else:
            c1 = a[0] / (3*j + 5)

        # subtract off c2 of d(x^2 y^j + 1)
        if p.divides(3*j + 7):
            c2 = base_ring(lift(a[1]) / (3*j + 7))
        else:
            c2 = a[1] / (3*j + 7)

        next_a[0] = next_a[0] + B*c1*(j+1)
        next_a[1] = next_a[1] + A*c1*(j+1) + B*c2*(j+1)
        next_a[2] = next_a[2]              + A*c2*(j+1)

        if exact_form is not None:
            exact_form += (c1*x + c2 * x**2) * y**(j+1)

    coeffs[int(offset)] = next_a

    return exact_form



def reduce_zero(Q, coeffs, offset, exact_form=None):
    """
    Applies cohomology relation to incorporate $x^2 y^0$ term into $x^0 y^0$
    and $x^1 y^0$ terms.

    INPUT:
        Q -- cubic polynomial
        coeffs -- list of length 3 lists. The i^th list [a, b, c]
                  represents $y^{2(i - offset)} (a + bx + cx^2) dx/y$.
        offset -- nonnegative integer

    OUTPUT:
        The reduction is performed in-place. The output is placed in
        coeffs[offset]. This method completely ignores coeffs[i] for
        i != offset.

    EXAMPLE:
        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: coeffs = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_zero(Q, coeffs, 1)
        sage: coeffs[1]
         [6, 5, 0]

    """

    a = coeffs[int(offset)]
    if a[2] == 0:
      return exact_form

    Qa = Q[1]

    a[0] = a[0] - a[2]*Qa/3    # $3x^2 dx/y = -a dx/y$

    coeffs[int(offset)] = a

    if exact_form is not None:
        x = exact_form.parent().gen(0)
        y = exact_form.parent().base_ring().gen(0)
        exact_form += Q.base_ring()(a[2] / 3) * y

    a[2] = 0

    coeffs[int(offset)] = a
    return exact_form



def reduce_all(Q, p, coeffs, offset, compute_exact_form=False):
    """
    Applies cohomology relations to reduce all terms to a linear combination
    of $dx/y$ and $x dx/y$.

    INPUT:
        Q -- cubic polynomial
        coeffs -- list of length 3 lists. The i^th list [a, b, c]
                  represents $y^{2(i - offset)} (a + bx + cx^2) dx/y$.
        offset -- nonnegative integer

    OUTPUT:
        A, B -- pair such that the input differential is cohomologous to
                (A + Bx) dx/y.

    NOTE:
        The algorithm operates in-place, so the data in coeffs is destroyed.

    EXAMPLE:
        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: coeffs = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_all(Q, 5, coeffs, 1)
         (21, 106)

    """

    R = Q.base_ring()

    if compute_exact_form:
#        exact_form = SpecialCubicQuotientRing(Q, laurent_series=True)(0)
        exact_form = PolynomialRing(LaurentSeriesRing(Q.base_ring(), 'y'), 'x')(0)
#        t = (Q.base_ring().order().factor())[0]
#        from sage.rings.padics.qp import pAdicField
#        exact_form = PolynomialRing(LaurentSeriesRing(pAdicField(p, t[1]), 'y'), 'x')(0)
    else:
        exact_form = None

    while len(coeffs) <= offset:
        coeffs.append([R(0), R(0), R(0)])

    exact_form = reduce_negative(Q, p, coeffs, offset, exact_form)
    exact_form = reduce_positive(Q, p, coeffs, offset, exact_form)
    exact_form = reduce_zero(Q, coeffs, offset, exact_form)

    if exact_form is None:
        return coeffs[int(offset)][0], coeffs[int(offset)][1]
    else:
        return (coeffs[int(offset)][0], coeffs[int(offset)][1]), exact_form



def frobenius_expansion_by_newton(Q, p, M):
  r"""
  Computes the action of Frobenius on $dx/y$ and on $x dx/y$, using
  Newton's method (as suggested in Kedlaya's paper).

  (This function does *not* yet use the cohomology relations
  -- that happens afterwards in the "reduction" step.)

  More specifically, it finds $F_0$ and $F_1$ in the quotient ring
  $R[x, T]/(T - Q(x))$, such that
  $$
     F(  dx/y) = T^{-r} F0 dx/y,   \text{ and }
     F(x dx/y) = T^{-r} F1 dx/y
  $$
  where
  $$
     r = ( (2M-3)p - 1 )/2.
  $$
  (Here $T$ is $y^2 = z^{-2}$, and $R$ is the coefficient ring of $Q$.)

  $F_0$ and $F_1$ are computed in the SpecialCubicQuotientRing
  associated to $Q$, so all powers of $x^j$ for $j \geq 3$ are reduced
  to powers of $T$.

  INPUT:
      Q -- cubic polynomial of the form Q(x) = x^3 + ax + b,
           whose coefficient ring is a Z/(p^M)Z-algebra
      p -- residue characteristic of the p-adic field
      M -- p-adic precision of the coefficient ring (this will be used
           to determine the number of Newton iterations)

  OUTPUT:
      F0, F1 -- elements of SpecialCubicQuotientRing(Q), as described above
      r -- non-negative integer, as described above

  """

  S = SpecialCubicQuotientRing(Q)
  x, _ = S.gens()     # T = y^2
  base_ring = S.base_ring()

  # When we compute Frob(1/y) we actually only need precision M-1, since
  # we're going to multiply by p at the end anyway.
  M = float(M - 1)

  # Kedlaya sets s = Q(x^p)/T^p = 1 + p T^{-p} E, where
  # E = (Q(x^p) - Q(x)^p) / p (has integral coefficients).
  # Then he computes s^{-1/2} in S, using Newton's method to find
  # successive approximations. We follow this plan, but we normalise our
  # approximations so that we only ever need positive powers of T.

  # Start by setting r = Q(x^p)/2 = 1/2 T^p s.
  # (The 1/2 is for convenience later on.)
  x_to_p_less_one = x**(p-1)
  x_to_p = x_to_p_less_one * x
  x_to_p_cubed = x_to_p.square() * x_to_p
  r = (base_ring(1) / base_ring(2)) * (x_to_p_cubed + Q[1]*x_to_p + S(Q[0]))

  # todo: this next loop would be clearer if it used the newton_method_sizes()
  # function

  # We will start with a hard-coded initial approximation, which we provide
  # up to precision 3. First work out what precision is best to start with.
  if M <= 3:
      initial_precision = M
  elif ceil(log(M/2, 2)) == ceil(log(M/3, 2)):
      # In this case there's no advantage to starting with precision three,
      # because we'll overshoot at the end. E.g. suppose the final precision
      # is 8. If we start with precision 2, we need two iterations to get us
      # to 8. If we start at precision 3, we will still need two iterations,
      # but we do more work along the way. So may as well start with only 2.
      initial_precision = 2
  else:
      initial_precision = 3

  # Now compute the first approximation. In the main loop below, X is the
  # normalised approximation, and k is the precision. More specifically,
  # X = T^{p(k-1)} x_i, where x_i is an approximation to s^{-1/2}, and the
  # approximation is correct mod p^k.
  if initial_precision == 1:
    k = 1
    X = S(1)
  elif initial_precision == 2:
    # approximation is 3/2 - 1/2 s
    k = 2
    X = S(base_ring(3) / base_ring(2)).shift(p) - r
  elif initial_precision == 3:
    # approximation is (15 - 10 s + 3 s^2) / 8
    k = 3
    X = (base_ring(1) / base_ring(8)) * (
           S(15).shift(2*p) - \
           (base_ring(20) * r).shift(p) + \
           (base_ring(12) * r.square()) \
        )

  # The key to the following calculation is that the T^{-m} coefficient
  # of every x_i is divisible by p^(ceil(m/p)) (for m >= 0). Therefore if
  # we are only expecting an answer correct mod p^k, we can truncate
  # beyond the T^{-(k-1)p} term without any problems.

  # todo: what would be really nice is to be able to work in a lower
  # precision *coefficient ring* when we start the iteration, and move up to
  # higher precision rings as the iteration proceeds. This would be feasible
  # over Integers(p**n), but quite complicated (maybe impossible) over a more
  # general base ring. This might give a decent constant factor speedup;
  # or it might not, depending on how much the last iteration dominates the
  # whole runtime. My guess is that it isn't worth the effort.

  three_halves = base_ring(3) / base_ring(2)

  # Newton iteration loop
  while k < M:
    # target_k = k' = precision we want our answer to be after this iteration
    target_k = 2*k

    # This prevents us overshooting. For example if the current precision
    # is 3 and we want to get to 10, we're better off going up to 5
    # instead of 6, because it's less work to get from 5 to 10 than it
    # is to get from 6 to 10.
    if ceil(log(M/target_k, 2)) == ceil(log(M/(target_k-1), 2)):
      target_k -= 1

    # temp = T^{p(3k-2)} 1/2 s x_i^3
    temp = X.square() * (X * r)

    # We know that the final result is only going to be correct mod
    # p^(target_k), so we might as well truncate the extraneous terms now.
    # temp = T^{p(k'-1)} 1/2 s x_i^3
    temp = temp.shift(-p*(3*k - target_k - 1))

    # X = T^{p(k'-1)} (3/2 x_i - 1/2 s x_i^3)
    #   = T^{p(k'-1)} x_{i+1}
    X = (three_halves * X).shift(p*(target_k - k)) - temp

    k = target_k

  # Now k should equal M, since we're up to the correct precision
  assert k == M, "Oops, something went wrong in the iteration"

  # We should have s^{-1/2} correct to precision M.
  # The following line can be uncommented to verify this.
  # (It's a slow verification though, can double the whole computation time.)

  #assert (p * X.square() * r * base_ring(2)).coeffs() == \
  #       R(p).shift(p*(2*M - 1)).coeffs()

  # Finally incorporate frobenius of dx and x dx, and choose offset that
  # compensates for our normalisations by powers of T.
  F0 = base_ring(p) * x_to_p_less_one * X
  F1 = F0 * x_to_p
  offset = ((2*k-1)*p - 1)/2

  return F0, F1, offset




def frobenius_expansion_by_series(Q, p, M):
  r"""
  Computes the action of Frobenius on dx/y and on x dx/y, using
  a series expansion.

  (This function computes the same thing as frobenius_expansion_by_newton(),
  using a different method. Theoretically the Newton method should be
  asymptotically faster, when the precision gets large. However, in practice,
  this functions seems to be marginally faster for moderate precision, so I'm
  keeping it here until I figure out exactly why it's faster.)

  (This function does *not* yet use the cohomology relations
  -- that happens afterwards in the "reduction" step.)

  More specifically, it finds F0 and F1 in the quotient ring
  $R[x, T]/(T - Q(x))$, such that
     $F(  dx/y) = T^{-r} F0 dx/y$,   and
     $F(x dx/y) = T^{-r} F1 dx/y$
  where
     $r = ( (2M-3)p - 1 )/2$.
  (Here T is $y^2 = z^{-2}$, and R is the coefficient ring of Q.)

  $F_0$ and $F_1$ are computed in the SpecialCubicQuotientRing associated
  to $Q$, so all powers of $x^j$ for $j \geq 3$ are reduced to powers of $T$.

  It uses the sum
     $$ F0 = \sum_{k=0}^{M-2} {-1/2 \choose k} p x^{p-1} E^k T^{(M-2-k)p}$$
  and
     $$ F1 = x^p F0,$$
  where $E = Q(x^p) - Q(x)^p$.

  INPUT:
      Q -- cubic polynomial of the form $Q(x) = x^3 + ax + b$,
           whose coefficient ring is a $\Z/(p^M)\Z$-algebra
      p -- residue characteristic of the p-adic field
      M -- p-adic precision of the coefficient ring (this will be used
           to determine the number of terms in the series)

  OUTPUT:
      F0, F1 -- elements of SpecialCubicQuotientRing(Q), as described above
      r -- non-negative integer, as described above

  """

  S = SpecialCubicQuotientRing(Q)
  x, _ = S.gens()
  base_ring = S.base_ring()

  x_to_p_less_1 = x**(p-1)
  x_to_p = x_to_p_less_1 * x

  # compute frobQ = Q(x^p)
  x_to_p_squared = x_to_p * x_to_p
  x_to_p_cubed = x_to_p_squared * x_to_p
  frobQ = x_to_p_cubed + Q[1]*x_to_p + Q[0]*S(1)
  # anticipating the day when p = 3 is supported:
  # frobQ = x_to_p_cubed + Q[2]*x_to_p_squared + Q[1]*x_to_p + Q[0]*S(1)

  E = frobQ - S(1).shift(p)    # E =  Q(x^p) - Q(x)^p

  offset = int( ((2*M-3)*p-1)/2 )
  term = p * x_to_p_less_1
  F0 = term.shift((M-2)*p)

  # todo: Possible speedup idea, perhaps by a factor of 2, but
  # it requires a lot of work:
  # Note that p divides E, so p^k divides E^k. So when we are
  # working with high powers of E, we're doing a lot more work
  # in the multiplications than we need to. To take advantage of
  # this we would need some protocol for "lowering the precision"
  # of a SpecialCubicQuotientRing. This would be quite messy to
  # do properly over an arbitrary base ring. Perhaps it is
  # feasible to do for the most common case (i.e. Z/p^nZ).
  # (but it probably won't save much time unless p^n is very
  # large, because the machine word size is probably pretty
  # big anyway.)

  for k in range(int(1), int(M-1)):
    term = term * E
    c = base_ring(binomial(-Integer(1)/2, k))
    F0 += (term * c).shift((M-k-2)*p)

  return F0, F0 * x_to_p, offset



def adjusted_prec(p, prec):
    r"""
    Computes how much precision is required in matrix_of_frobenius to get
    an answer correct to prec $p$-adic digits.

    The issue is that the algorithm used in matrix_of_frobenius sometimes
    performs divisions by $p$, so precision is lost during the algorithm.

    The estimate returned by this function is based on Kedlaya's result
    (Lemmas 2 and 3 of ``Counting Points on Hyperelliptic Curves...''), which
    implies that if we start with $M$ $p$-adic digits, the total precision
    loss is at most
       $1 + \lfloor \log_p(2M-3) \rfloor$
    $p$-adic digits. (This estimate is somewhat less than the amount you
    would expect by naively counting the number of divisions by $p$.)

    INPUT:
        p -- a prime >= 5
        prec -- integer, desired output precision, >= 1

    OUTPUT:
        adjusted precision (usually slightly more than prec)

    """

    # initial estimate:
    if prec <= 2:
      adjusted = 2
    else:
      adjusted = prec + int(log(2*prec - 3, p)) - 1

    # increase it until we have enough
    while adjusted - int(log(2*adjusted - 3, p)) - 1 < prec:
        adjusted += 1

    return adjusted



def matrix_of_frobenius(Q, p, M, trace=None, compute_exact_forms=False):
  """
  Computes the matrix of Frobenius on Monsky-Washnitzer cohomology,
  with respect to the basis $(dx/y, x dx/y)$.

  INPUT:
      Q -- cubic polynomial $Q(x) = x^3 + ax + b$ defining an elliptic
           curve E by $y^2 = Q(x)$. The coefficient ring of Q should be
           a $\Z/(p^M)\Z$-algebra in which the matrix of frobenius will
           be constructed.
      p -- prime >= 5 for which E has good reduction
      M -- integer >= 2; $p$-adic precision of the coefficient ring
      trace -- (optional) the trace of the matrix, if known in advance.
           This is easy to compute because it's just the $a_p$ of the
           curve. If the trace is supplied, matrix_of_frobenius will
           use it to speed the computation (i.e. we know the determinant
           is $p$, so we have two conditions, so really only column of
           the matrix needs to be computed. It's actually a little more
           complicated than that, but that's the basic idea.)
           If trace=None, then both columns will be computed
           independently, and you can get a strong indication of
           correctness by verifying the trace afterwards.

  WARNING:
      -- THE RESULT WILL NOT NECESSARILY BE CORRECT TO M p-ADIC DIGITS.
         If you want prec digits of precision, you need to use the function
         adjusted_prec(), and then you need to reduce the answer mod p^prec
         at the end.

  OUTPUT:
      2x2 matrix of frobenius on Monsky-Washnitzer cohomology, with entries
      in the coefficient ring of Q.

  EXAMPLES:
    A simple example:
      sage: p = 5
      sage: prec = 3
      sage: M = monsky_washnitzer.adjusted_prec(p, prec)
      sage: M
      5
      sage: R.<x> = PolynomialRing(Integers(p**M))
      sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4), p, M)
      sage: A
      [3090  187]
      [2945  408]

    But the result is only accurate to prec digits:
      sage: B = A.change_ring(Integers(p**prec))
      sage: B
       [90 62]
       [70 33]

    Check trace (123 = -2 mod 125) and determinant:
      sage: B.det()
       5
      sage: B.trace()
       123
      sage: EllipticCurve([-1, 1/4]).ap(5)
       -2

    Try using the trace to speed up the calculation:
      sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4),
      ...                                             p, M, -2)
      sage: A
       [2715  187]
       [1445  408]

    Hmmm... it looks different, but that's because the trace of our
    first answer was only -2 modulo $5^3$, not -2 modulo $5^5$. So the
    right answer is:
    sage: A.change_ring(Integers(p**prec))
       [90 62]
       [70 33]

    Check it works with only one digit of precision:
      sage: p = 5
      sage: prec = 1
      sage: M = monsky_washnitzer.adjusted_prec(p, prec)
      sage: R.<x> = PolynomialRing(Integers(p**M))
      sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4), p, M)
      sage: A.change_ring(Integers(p))
       [0 2]
       [0 3]

    Here's an example that's particularly badly conditioned for using the
    trace trick:
      sage: p = 11
      sage: prec = 3
      sage: M = monsky_washnitzer.adjusted_prec(p, prec)
      sage: R.<x> = PolynomialRing(Integers(p**M))
      sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 + 7*x + 8, p, M)
      sage: A.change_ring(Integers(p**prec))
       [1144  176]
       [ 847  185]

    The problem here is that the top-right entry is divisible by 11,
    and the bottom-left entry is divisible by $11^2$. So when you
    apply the trace trick, neither $F(dx/y)$ nor $F(x dx/y)$ is enough
    to compute the whole matrix to the desired precision, even if you
    try increasing the target precision by one. Nevertheless,
    \code{matrix_of_frobenius} knows how to get the right answer by
    evaluating $F((x+1) dx/y)$ instead:

      sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 + 7*x + 8, p, M, -2)
      sage: A.change_ring(Integers(p**prec))
       [1144  176]
       [ 847  185]

    The running time is about \code{O(p * prec**2)} (times some
    logarithmic factors), so it's feasible to run on fairly large
    primes, or precision (or both?!?!):

      sage: p = 10007
      sage: prec = 2
      sage: M = monsky_washnitzer.adjusted_prec(p, prec)
      sage: R.<x> = PolynomialRing(Integers(p**M))
      sage: A = monsky_washnitzer.matrix_of_frobenius(            # long time
      ...                             x^3 - x + R(1/4), p, M)     # long time
      sage: B = A.change_ring(Integers(p**prec)); B               # long time
       [74311982 57996908]
       [95877067 25828133]
      sage: B.det()                                               # long time
       10007
      sage: B.trace()                                             # long time
       66
      sage: EllipticCurve([-1, 1/4]).ap(10007)                    # long time
       66

      sage: p = 5
      sage: prec = 300
      sage: M = monsky_washnitzer.adjusted_prec(p, prec)
      sage: R.<x> = PolynomialRing(Integers(p**M))
      sage: A = monsky_washnitzer.matrix_of_frobenius(            # long time
      ...                             x^3 - x + R(1/4), p, M)     # long time
      sage: B = A.change_ring(Integers(p**prec))                  # long time
      sage: B.det()                                               # long time
       5
      sage: -B.trace()                                            # long time
       2
      sage: EllipticCurve([-1, 1/4]).ap(5)                        # long time
       -2

    Let's check consistency of the results for a range of precisions:
      sage: p = 5
      sage: max_prec = 60
      sage: M = monsky_washnitzer.adjusted_prec(p, max_prec)
      sage: R.<x> = PolynomialRing(Integers(p**M))
      sage: A = monsky_washnitzer.matrix_of_frobenius(            # long time
      ...                         x^3 - x + R(1/4), p, M)         # long time
      sage: A = A.change_ring(Integers(p**max_prec))              # long time
      sage: result = []                                           # long time
      sage: for prec in range(1, max_prec):                       # long time
      ...       M = monsky_washnitzer.adjusted_prec(p, prec)      # long time
      ...       R.<x> = PolynomialRing(Integers(p^M),'x')         # long time
      ...       B = monsky_washnitzer.matrix_of_frobenius(        # long time
      ...                         x^3 - x + R(1/4), p, M)         # long time
      ...       B = B.change_ring(Integers(p**prec))              # long time
      ...       result.append(B == A.change_ring(                 # long time
      ...                                Integers(p**prec)))      # long time
      sage: result == [True] * (max_prec - 1)                     # long time
       True


    The remaining examples discuss what happens when you take the coefficient
    ring to be a power series ring; i.e. in effect you're looking at a family
    of curves.

    The code does in fact work...
      sage: p = 11
      sage: prec = 3
      sage: M = monsky_washnitzer.adjusted_prec(p, prec)
      sage: S.<t> = PowerSeriesRing(Integers(p**M), default_prec=4)
      sage: a = 7 + t + 3*t^2
      sage: b = 8 - 6*t + 17*t^2
      sage: R.<x> = PolynomialRing(S)
      sage: Q = x**3 + a*x + b
      sage: A = monsky_washnitzer.matrix_of_frobenius(Q, p, M)    # long time
      sage: B = A.change_ring(PowerSeriesRing(                    # long time
      ...         Integers(p**prec), 't', default_prec=4))        # long time
      sage: B                                                     # long time
       [1144 + 264*t + 841*t^2 + 1025*t^3 + O(t^4)  176 + 1052*t + 216*t^2 + 523*t^3 + O(t^4)]
       [   847 + 668*t + 81*t^2 + 424*t^3 + O(t^4)   185 + 341*t + 171*t^2 + 642*t^3 + O(t^4)]

    The trace trick should work for power series rings too, even in the badly-
    conditioned case. Unfortunately I don't know how to compute the trace in
    advance, so I'm not sure exactly how this would help. Also, I suspect
    the running time will be dominated by the expansion, so the trace trick
    won't really speed things up anyway. Another problem is that the
    determinant is not always p:
      sage: B.det()                                               # long time
       11 + 484*t^2 + 451*t^3 + O(t^4)

    However, it appears that the determinant always has the property that if
    you substitute t -> 11t, you do get the constant series p (mod p**prec).
    Similarly for the trace. And since the parameter only really makes sense
    when it's divisible by p anyway, perhaps this isn't a problem after all.

  """

  M = int(M)
  if M < 2:
    raise ValueError, "M (=%s) must be at least 2" % M

  base_ring = Q.base_ring()

  # Expand out frobenius of dx/y and x dx/y.
  # (You can substitute frobenius_expansion_by_series here, that will work
  # as well. See its docstring for some performance notes.)
  F0, F1, offset = frobenius_expansion_by_newton(Q, p, M)
  #F0, F1, offset = frobenius_expansion_by_series(Q, p, M)

  if compute_exact_forms:
    # we need to do all the work to get the exact expressions f such that F(x^i dx/y) = df + \sum a_i x^i dx/y

    F0_coeffs = transpose_list(F0.coeffs())
    F0_reduced, f_0 = reduce_all(Q, p, F0_coeffs, offset, True)

    F1_coeffs = transpose_list(F1.coeffs())
    F1_reduced, f_1 = reduce_all(Q, p, F1_coeffs, offset, True)

  elif M == 2:
    # This implies that only one digit of precision is valid, so we only need
    # to reduce the second column. Also, the trace doesn't help at all.

    F0_reduced = [ base_ring(0), base_ring(0) ]

    F1_coeffs = transpose_list(F1.coeffs())
    F1_reduced = reduce_all(Q, p, F1_coeffs, offset)

  elif trace is None:
    # No trace provided, just reduce F(dx/y) and F(x dx/y) separately.

    F0_coeffs = transpose_list(F0.coeffs())
    F0_reduced = reduce_all(Q, p, F0_coeffs, offset)

    F1_coeffs = transpose_list(F1.coeffs())
    F1_reduced = reduce_all(Q, p, F1_coeffs, offset)

  else:
    # Trace has been provided.

    # In most cases this can be used to quickly compute F(dx/y) from
    # F(x dx/y). However, if we're unlucky, the (dx/y)-component of
    # F(x dx/y) (i.e. the top-right corner of the matrix) may be divisible
    # by p, in which case there isn't enough information to get the
    # (x dx/y)-component of F(dx/y) to the desired precision. When this
    # happens, it turns out that F((x+1) dx/y) always *does* give enough
    # information (together with the trace) to get both columns to the
    # desired precision.

    # First however we need a quick way of telling whether the top-right
    # corner is divisible by p, i.e. we want to compute the second column
    # of the matrix mod p. We could do this by just running the entire
    # algorithm with M = 2 (which assures precision 1). Luckily, we've
    # already done most of the work by computing F1 to high precision; so
    # all we need to do is extract the coefficients that would correspond
    # to the first term of the series, and run the reduction on them.

    # todo: actually we only need to do this reduction step mod p^2, not
    # mod p^M, which is what the code currently does. If the base ring
    # is Integers(p^M), then it's easy. Otherwise it's tricky to construct
    # the right ring, I don't know how to do it.

    F1_coeffs = transpose_list(F1.coeffs())
    F1_modp_coeffs = F1_coeffs[int((M-2)*p):]
    # make a copy, because reduce_all will destroy the coefficients:
    F1_modp_coeffs = [[cell for cell in row] for row in F1_modp_coeffs]
    F1_modp_offset = offset - (M-2)*p
    F1_modp_reduced = reduce_all(Q, p, F1_modp_coeffs, F1_modp_offset)

    if F1_modp_reduced[0].is_unit():
      # If the first entry is invertible mod p, then F(x dx/y) is sufficient
      # to get the whole matrix.

      F1_reduced = reduce_all(Q, p, F1_coeffs, offset)

      F0_reduced = [ base_ring(trace) - F1_reduced[1], None ]
      # using that the determinant is p:
      F0_reduced[1] = (F0_reduced[0] * F1_reduced[1] - base_ring(p)) \
                      / F1_reduced[0]

    else:
      # If the first entry is zero mod p, then F((x+1) dx/y) will be sufficient
      # to get the whole matrix. (Here we are using the fact that the second
      # entry *cannot* be zero mod p. This is guaranteed by some results in
      # section 3.2 of ``Computation of p-adic Heights and Log Convergence''
      # by Mazur, Stein, Tate. But let's quickly check it anyway :-))
      assert F1_modp_reduced[1].is_unit(), \
         "Hey that's impossible! The second entry in the second column " \
         "should be invertible mod p!"

      G0_coeffs = transpose_list( (F0 + F1).coeffs())
      G0_reduced = reduce_all(Q, p, G0_coeffs, offset)

      # Now G0_reduced expresses F((x+1) dx/y) in terms of dx/y and x dx/y.
      # Re-express this in terms of (x+1) dx/y and x dx/y.
      H0_reduced = [ G0_reduced[0], G0_reduced[1] - G0_reduced[0] ]

      # The thing we're about to divide by better be a unit.
      assert H0_reduced[1].is_unit(), \
         "Hey that's impossible! The second entry in this column " \
         "should be invertible mod p!"

      # Figure out the second column using the trace...
      H1_reduced = [ None, base_ring(trace) - H0_reduced[0] ]
      # ... and using that the determinant is p:
      H1_reduced[0] = (H0_reduced[0] * H1_reduced[1] - base_ring(p)) \
                      / H0_reduced[1]

      # Finally, change back to the usual basis (dx/y, x dx/y)
      F1_reduced = [ H1_reduced[0], \
                     H1_reduced[0] + H1_reduced[1] ]
      F0_reduced = [ H0_reduced[0] - F1_reduced[0],
                     H0_reduced[0] + H0_reduced[1] - F1_reduced[1] ]

    # One more sanity check: our final result should be congruent mod p
    # to the approximation we used earlier.
    assert not (
      (F1_reduced[0] - F1_modp_reduced[0]).is_unit() or \
      (F1_reduced[1] - F1_modp_reduced[1]).is_unit() or \
      F0_reduced[0].is_unit() or F0_reduced[1].is_unit()), \
      "Hey that's impossible! The output matrix is not congruent mod p " \
      "to the approximation found earlier!"

  if compute_exact_forms:
      return matrix(base_ring, 2, 2, [F0_reduced[0], F1_reduced[0],
                                      F0_reduced[1], F1_reduced[1]]), f_0, f_1
  else:
      return matrix(base_ring, 2, 2, [F0_reduced[0], F1_reduced[0],
                                      F0_reduced[1], F1_reduced[1]])



#*****************************************************************************
# Here is an implementation of D. Harvey's modification of Kedlaya's
# algorithm, ``Kedlaya's algorithm in larger characteristic'' (see arXiv).
# Actually the algorithm implemented here is better than the one in the paper;
# the $N^3$ contribution is reduced to $N^{5/2}$, and incorporates some ideas
# from ``Linear recurrences with polynomial coefficients and application to
# integer factorization and Cartier-Manin operator'' (Bostan, Gaudry, Schost)
# to reduce the running time by another factor of $\log p$.
#
# todo: I haven't yet done the precision analysis (i.e. the analogue of
# lemmas 2 and 3 from Kedlaya's original paper), so for now I'm playing it
# safe, using more precision than is presumably necessary.
#
# AUTHOR:
#    -- David Harvey (2007-02)
#
#*****************************************************************************


def shift_evaluation_points(R, values, b, r, p, cache):
    r"""
    Given the values of a polynomial $P(x)$ on an arithmetic progression,
    this function computes the values of $P(x)$ on a shift of that
    arithmetic progression.

    INPUT:
      p -- a prime
      R -- coefficient ring, of the form $\ZZ/p^m \ZZ$
      values -- a list of $d+1$ values of a polynomial $P(x)$ of degree $d$,
          at the $d+1$ points $x = 0, r, 2r, \ldots, dr$.
      b -- an integer
      r -- a nonzero integer
      cache -- a dictionary

    OUTPUT:
      A list of $d+1$ values $P(b)$, $P(b+r)$, ..., $P(b+dr)$.
      (Also the cache might get updated.)

    NOTE: The polynomial $P$ is *not* part of the input.

    PRECONDITIONS: $d$ must be less than $p$

    EXAMPLES:
    Run it on some random input, and compare to a naive implementation:
        sage: from sage.schemes.elliptic_curves.monsky_washnitzer import shift_evaluation_points
        sage: R = Integers(7**3)
        sage: S.<x> = PolynomialRing(R)
        sage: r = 3
        sage: for d in range(7):
        ...      for b in range(50):
        ...         p = S([R.random_element() for i in range(d+1)])
        ...         L = [p(i*r) for i in range(d+1)]
        ...         shifted = shift_evaluation_points(R, L, b, r, 7, {})
        ...         if shifted != [p(b + i*r) for i in range(d+1)]:
        ...            print "failed"

    ALGORITHM:
    The basic idea (as described in the Bostan/Gaudry/Schost paper) is as
    follows. The Lagrange interpolation formula says that
      $$ P(x) = \sum_{i=0}^d P(ir) \prod_{j \neq i} (x-jr)/(ir-jr). $$
    We substitute $x = b+kr$, and rearrange to eventually obtain
      $$ P(b+kr) = r^{-d} \prod_{j=0}^d (b+(k-j)r)
        \sum_{i=0}^d [P(ir) \prod_{j \neq i} (i-j)^{-1}] \frac{1}{b+(k-i)r}. $$
    The point is that the latter sum is a convolution and so can be
    evaluated for all $k$ simultaneously by a single polynomial multiplication.

    There are some complications related to non-invertible elements of $R$,
    but we just deal with it :-)

    """

    d = len(values) - 1
    assert p > d

    try:
        # try to retrieve some precomputed stuff from the cache
        bad, input_scale, S, kernel, output_scale = cache[d, b, r]

    except KeyError:
        # that didn't work; need to compute all of it

        # figure out for which i in [-d, d] is b+ir not invertible;
        # put these i into a "bad" list
        i = (-R(b)/r).lift() % p
        bad = []
        if i <= d:
            bad.append(i)
        if i >= p-d:
            bad.append(i-p)

        # compute inverse factorials [1/0!, 1/1!, 1/2!, ..., 1/d!] in R
        total = R(1)
        for k in range(2, d+1):
            total = total * R(k)
        inv = total**(-1)
        inv_fac = [inv]
        for k in range(d, 1, -1):
            inv_fac.append(inv_fac[-1] * R(k))
        inv_fac.append(R(1))
        inv_fac.reverse()

        # compute some coefficients to be used in the preprocessing step
        scale = R(r)**(-d)
        input_scale = [scale * inv_fac[i] * inv_fac[d-i] * \
                       (-1 if (d-i) & 1 else 1) for i in range(d+1)]

        # Let B_i = b+ir for each i in [-d, d], except we declare B_i = 1
        # if i is "bad".
        B = [(R(1) if i in bad else R(b + i*r)) for i in range(-d, d+1)]

        # compute partial products of B_i
        products = [R(1)]
        for z in B:
            products.append(products[-1] * z)

        # compute partial products of inverses of B_i
        inverse = ~(products[-1])
        inv_products = [inverse]
        for z in reversed(B):
            inv_products.append(inv_products[-1] * z)
        inv_products.reverse()

        # compute inverses of each B_i
        B_inv = [inv_products[i] * products[i-1] \
                 for i in range(1, len(products))]

        # replace the "bad" B_i with 0, so they don't contribute to the
        # convolution at all
        for i in bad:
            B_inv[i+d] = R(0)

        # kernel is the polynomial we need to multiply by in the convolution
        S = PolynomialRing(R, "x")
        kernel = S(B_inv, check=False)

        # coefficients for postprocessing step
        output_scale = [products[k+d+1] * inv_products[k] for k in range(d+1)]

        # save everything useful to the cache
        cache[d, b, r] = (bad, input_scale, S, kernel, output_scale)


    # preprocessing step:
    # multiply each P(ir) by \prod_{j \neq i} r (i-j)^{-1}
    values = [values[i] * input_scale[i] for i in range(d+1)]

    # do the convolution
    H = S(values, check=False) * kernel
    output = H.list()[d:(2*d+1)]
    # (need to zero-pad output up to length d+1)
    output.extend([R(0)] * (d+1 - len(output)))

    # postprocessing:
    # multiply each output by \prod_{j=0}^d B_{k-j}
    output = [output[k] * output_scale[k] for k in range(d+1)]


    # Now need to deal with "bad" indices. This shouldn't happen too often,
    # especially when $d$ is much smaller than $p$.

    # The first problem is that in the postprocessing step we left out
    # multiplying by the bad factors. Now we put them back in.
    for i in bad:
        factor = b + i*r
        for k in range(max(i, 0), min(d+i, d) + 1):
            output[k] = output[k] * factor

    # The second problem is that in the polynomial S(B_inv) the monomials
    # corresponding to the bad i were completely wrong.
    for i in bad:
        if i >= 0:
            correction = ([R(0)] * i) + [values[k-i] * output_scale[k] \
                          for k in range(i, d+1)]
        else:
            correction = [values[k-i] * output_scale[k] \
                          for k in range(0, i+d+1)] + ([R(0)] * (-i))

        # again, we left out all the bad factors; put them back in
        for ii in bad:
            if ii != i:
                factor = b + ii*r
                for k in range(max(ii, 0), min(d+ii, d) + 1):
                    correction[k] = correction[k] * factor

        # add correction term to output
        for k in range(d+1):
            output[k] = output[k] + correction[k]

    return output



def naive_shift_evaluation_points(R, values, b, r, p):
    r"""
    Same as shift_evaluation_points, but uses naive algorithm (naive
    Lagrange interpolation, followed by naive evaluation).

    For testing/debugging purposes only.
    """
    d = len(values) - 1
    S, x = PolynomialRing(R, "x").objgen()
    total = S(0)
    for j in range(d+1):
        prod = S(1)
        for i in range(d+1):
            if i != j:
                prod = prod * (x - i*r) / (r*(j-i))
        total = total + values[j] * prod

    return [total(b + i*r) for i in range(d+1)]



def matrix_polys(A, R, n, b, r, q, p, cache):
    r"""
    INPUT:
      p -- a prime
      R -- coefficient ring, of the form $\ZZ/p^m \ZZ$
      n -- a positive integer
      A -- an $n \times n$ matrix of linear polys over $R$ (stored as a list
           of lists)
      b -- an integer
      r -- a positive integer
      q -- a positive integer
      cache -- dictionary (passed through to shift_evaluation_points)

    OUTPUT:
      Let $B(x) = A(x+1) A(x+2) ... A(x+r-1) A(x+r)$, a matrix of degree
      $r$ polynomials.

      Output an $n \times n$ matrix, the $(i, j)$-entry is a list of the
      $(i, j)$-entries of the evaluations $B(b), B(b+q), \ldots, B(b+rq)$;
      i.e. it evaluates $B(x)$ at $r+1$ points spaced out by $q$.

    ALGORITHM:
      Umm, it's probably the one described in the Bostan/Gaudry/Schost paper,
      although I haven't looked at that part of the paper too closely yet.

    EXAMPLES:
    Run it on some random input, and compare to a naive implementation:
        sage: from sage.schemes.elliptic_curves.monsky_washnitzer import matrix_polys, naive_matrix_polys
        sage: R = Integers(7**3)
        sage: S.<x> = PolynomialRing(R, "x")
        sage: b = 5
        sage: for r in range(1, 14):
        ...       A = [[R.random_element() * x + R.random_element() \
        ...             for i in range(3)] for j in range(3)]
        ...       X = matrix_polys(A, R, 3, b, r, 3, 7, {})
        ...       Y = naive_matrix_polys(A, R, 3, b, r, 3, 7)
        ...       if X != Y: print "failed"

    """
    assert r >= 1
    assert r < 2*p

    if r == 1:
        # base case; i.e. evaluate original polys at x = b+1, b+q+1
        return [[[A[j][i](b+1), A[j][i](b+q+1)] for i in range(n)]
                for j in range(n)]

    # let r' = floor(r/2)
    half_r = r // 2

    # evaluate A(x+1) A(x+2) ... A(x+r')
    # at x = b, b + q, ..., b + r'q.
    X = matrix_polys(A, R, n, b, half_r, q, p, cache)

    # shift to obtain evaluations of A(x+r'+1) A(x+r'+2) ... A(x+2r')
    # at x = b, b + q, ..., b + r'q.
    shift1 = [[shift_evaluation_points(R, values, half_r, q, p, cache) \
               for values in row] for row in X]

    # shift to obtain evaluations of A(x+1) A(x+2) ... A(x+r')
    # at x = b + (r'+1)q, b + (r'+2)q, ..., b + (2r'+1)q
    shift2 = [[shift_evaluation_points(R, values, (half_r+1)*q, q, p, cache) \
               for values in row] for row in X]

    # shift to obtain evaluations of A(x+r'+1) A(x+r'+2) ... A(x+2r')
    # at x = b + (r'+1)q, b + (r'+2)q, ..., b + (2r'+1)q
    shift3 = [[shift_evaluation_points(R, values, half_r, q, p, cache) \
               for values in row] for row in shift2]

    # multiply matrices to obtain evaluations of A(x+1) A(x+2) ... A(x+2r')
    # at x = b, b + q, ..., b + r'q.
    prod1 = [[[sum([X[i][m][k] * shift1[m][j][k] for m in range(n)]) \
               for k in range(half_r+1)] for j in range(n)] for i in range(n)]

    # multiply matrices to obtain evaluations of A(x+1) A(x+2) ... A(x+2r')
    # at x = b + (r'+1)q, ..., b + (2r'+1)q.
    prod2 = [[[sum([shift2[i][m][k] * shift3[m][j][k] for m in range(n)]) \
               for k in range(half_r+1)] for j in range(n)] for i in range(n)]

    if r & 1 == 0:
        # if r is even, then the polys are the right degree, but we have one
        # extra evaluation point (namely x = b + (r+1)q). Throw it away.
        for i in range(n):
            for j in range(n):
                prod2[i][j].pop()
                prod1[i][j].extend(prod2[i][j])
        output = prod1

    else:
        # if r is odd, then we have precisely the right number of evaluation
        # points, but our polynomials are short by one term, namely A(x+r).
        # Add in the effect of that term.
        for i in range(n):
            for j in range(n):
                prod1[i][j].extend(prod2[i][j])

        output = \
               [[[sum([prod1[i][m][k] * A[m][j](b+k*q+r) for m in range(n)]) \
                  for k in range(r+1)] for j in range(n)] for i in range(n)]

    return output



def naive_matrix_polys(A, R, n, b, r, q, p):
    r"""
    Same as matrix_polys, but uses a naive algorithm.

    For testing/debugging purposes only.
    """
    result = [[[] for _ in range(n)] for _ in range(n)]
    for i in range(r+1):
        prod = MatrixSpace(R, n, n).identity_matrix()
        for j in range(1, r+1):
            A_eval = matrix(R, [[f(b+q*i+j) for f in row] for row in A])
            prod = prod * A_eval
        for k1 in range(n):
            for k2 in range(n):
                result[k1][k2].append(prod[k1][k2])
    return result




def reduce_one(T, s, n, M, D, R):
    r"""
    A single step reduction.

    INPUT:
      R -- coefficient ring, of the form $\ZZ/p^m \ZZ$
      n -- a positive integer
      T -- an $n$-tuple of elements of $R$.
      M -- an $n \times n$ matrix of linear polys over $R$ (stored as a list
           of lists)
      D -- a linear polynomial over $R$.

    OUTPUT:
      Let $A(x) = M(x)/D(x)$.
      Return value is the matrix-vector product $A(s) T$, as a list.

    """
    d = D(s).lift()
    return [R(sum([M[i][j](s) * T[j] for j in range(n)]).lift() / d) \
            for i in range(n)]



def reduce_sequence(T, d, s, n, M, D, p, R):
    r"""
    Reduces a sequence of terms.

    INPUT:
      p -- a prime
      R -- coefficient ring, of the form $\ZZ/p^m \ZZ$
      T -- a nonempty list of $n$-tuples of elements of $R$. Denote by $L$
           the length of $T$.
      d -- a list of positive integers, of length $L$; they must be in
           strictly ascending order
      s -- an integer >= 0
      n -- an integer >= 1
      M -- an $n \times n$ matrix of linear polynomials over $R$,
           stored as a list of lists
      D -- a linear polynomial over $R$

    PRECONDITIONS:
       The largest $d_i$ must be less than $\sqrt p$.
       (or something like that... todo: check this...)

    OUTPUT:
       Let $A(x) = M(x)/D(x)$, and let $B(t)$ denote the matrix product
         $$ A(s+1) A(s+2) \cdots A(s+t). $$

       This function returns the sum
         $$ \sum_{i=1}^L  B(d_i) T_i. $$

       In other words, treating each term $T_i$ as appearing in position
       $s + d_i$, it applies the reduction matrices to reduce them all down to
       a single $n$-tuple in position $s$.

       Running time is something like
          $O(M(e^{1/2}) + L M(e^{1/4}))$
       where $e = \max(d_i)$, and where $M(k)$ is the time for degree $k$
       polynomial multiplication over $R$.

    """
    # L = number of terms to reduce
    L = len(T)
    # e = total distance to reduce
    e = d[-1]

    if e <= 30:
        # naive algorithm when the reduction length is small enough
        output = [R(0)] * n
        for i in range(L):
            term = T[i]
            for t in range(d[i], 0, -1):
                term = reduce_one(term, s + t, n, M, D, R)
            for j in range(n):
                output[j] = output[j] + term[j]
        return output

    # We do the reduction by splitting into segments. There are $w+1$
    # segments of length $w$, where $w$ is the smallest integer such
    # that $w(w+1) \geq e$ (so in particular $w \approx \sqrt e$).
    w = ceil(sqrt(e))
    while w*(w+1) >= e:
        w = w - 1
    if w*(w+1) < e:
        w = w + 1

    # The BGS method only works if the terms are on segment boundaries.
    # So first we (recursively) reduce each term separately to make it
    # lie on a segment boundary.
    # This is the step that takes time $O(L M(e^{1/4}))$.
    new_T = []
    new_d = []
    for i in range(L):
        this_d = d[i]
        reduced_d = w * (this_d // w)
        if reduced_d == this_d:
            reduced_T = T[i]
        else:
            reduced_T = reduce_sequence([T[i]], [this_d - reduced_d],
                                        s + reduced_d, n, M, D, p, R)
        new_T.append(reduced_T)
        new_d.append(reduced_d)
    T = new_T
    d = new_d

    # Now all terms are on segment boundaries. Next we compute a reduction
    # matrix for each segment, using BGS method.
    cache = {}
    M_segments = matrix_polys([[M[i][j] for j in range(n)] for i in range(n)],
                              R, n, s, w, w, p, cache)
    D_segments = matrix_polys([[D]], R, 1, s, w, w, p, cache)[0][0]

    # Now we use the reduction matrices to sweep up all the terms.
    total = [R(0), R(0), R(0)]
    for segment in range(w+1, -1, -1):
        # add in terms at the current segment boundary
        while len(d) and (d[-1] == segment*w):
            for j in range(n):
                total[j] = total[j] + T[-1][j]
            T.pop()
            d.pop()

        # apply reduction matrix
        if segment > 0:
            M_here = [[M_segments[i][j][segment-1] \
                       for j in range(n)] for i in range(n)]
            D_here = D_segments[segment-1].lift()
            total = [R(sum([M_here[i][j] * total[j]
                            for j in range(n)]).lift() / D_here) \
                     for i in range(n)]

    return total



def starting_differentials(R, a, b, p, M):
    r"""
    Computes the initial terms in the expansion of $\sigma(x^i dx/y)$
    (where $\sigma$ is frobenius) to which the cohomology reduction
    algorithm will be applied.

    INPUT:
      p -- a prime $\geq 5$
      R -- coefficient ring, of the form $\ZZ/p^m \ZZ$
      a, b -- elements of $R$, coefficients of an elliptic curve
              $y^2 = x^3 + ax + b$
      M -- number of terms being computed (see ``Kedlaya's algorithm in
           larger characteristic'' for the exact meaning)

    OUTPUT:
      A list of length $M$; the $j$th entry is a list of length $3j+1$;
      its $m$th entry is the coefficient of
         $ x^{p(m+i+1)-1} y^{-p(2j+1)+1} dx/y $
      in the expansion of $\sigma(x^i dx/y)$.

      Complexity is $O(M^2)$ operations in $R$.

    """
    # compute B[k][j] = binomial(k, j) in R, for 0 <= j <= k < 2M
    B = [[R(1)]]
    for k in range(1, 2*M):
        B.append([R(1)] + [B[-1][j] + B[-1][j-1] for j in range(1,k)] + [R(1)])

    # compute temp[k] = binomial(2k, k) / 4^k in R, for 0 <= k < M
    four_to_minus_1 = R(4)**(-1)
    four_to_minus_k = R(1)
    temp = [R(1)]
    for k in range(1, M):
        four_to_minus_k = four_to_minus_k * four_to_minus_1
        temp.append(four_to_minus_k * B[2*k][k])

    # compute constants
    # C[j] = (-1)^j p \sum_{k=j}^{M-1} 4^{-k} (2k \choose k) (k \choose j)
    # in R, for 0 <= j < M
    C = [p * sum([temp[k] * B[k][j] for k in range(j, M)]) for j in range(M)]
    for j in range(1, M, 2):
        C[j] = -C[j]

    # compute A[j][m] = u^m coefficient of (u^3 + a*u + b)^j
    # for 0 <= j < M
    Ru, u = PolynomialRing(R, "u").objgen()
    Q = u**3 + a*u + b
    Q_to_j = Ru(1)
    A = [[R(1)]]
    for j in range(1, M):
        Q_to_j = Q_to_j * Q
        A.append(Q_to_j.list())

    # compute products C[j] * A[j][m]
    return [[C[j] * x for x in A[j]] for j in range(M)]




def matrix_of_frobenius_alternate(a, b, p, N):
    r"""
    Computes the matrix of Frobenius on Monsky-Washnitzer cohomology,
    with respect to the basis $(dx/y, x dx/y)$.

    INPUT:
        a, b -- rational numbers, coefficients of elliptic curve
                $y^2 = x^3 + ax + b$
        p -- a prime for which the curve has good reduction
        N -- desired output precision for matrix (i.e. modulo $p^N$).

    OUTPUT:
        2x2 matrix of frobenius, with coefficients in Integers(p^N)

    TODO:
    Haven't yet checked yet exactly which combinations of $p$ and $N$ are
    valid. If $p$ is much larger than $N$, there shouldn't be a problem.

    EXAMPLES:
        sage: from monsky_washnitzer import matrix_of_frobenius_alternate
        sage: M = matrix_of_frobenius_alternate(-1, 1/4, 101, 3); M
        [824362 606695]
        [641451 205942]
        sage: M.det()
        101
        sage: M.trace()
        3

    """
    assert p > N    # todo: check this is the right condition....

    # choose minimal M such that M >= floor(log_p(2M+1)) + N
    M = N
    while M < floor(log(2*M+1, p)) + N:
      M = M + 1

    # I haven't done the precision analysis properly yet, so for now
    # we work to precision 3M, with M digits at the bottom and M spare
    # digits at the top, which should be pretty safe.
    R = Integers(p**(3*M))
    R_output = Integers(p**N)
    S, x = PolynomialRing(R, "x").objgen()
    a = R(a)
    b = R(b)

    differentials = starting_differentials(R, a, b, p, M)

    # multiply through by p^M for safety, need to come back and fix this...
    differentials = [[item * p**M for item in row] for row in differentials]

    output = []

    for i in range(2):
        column = []

        for j in range(M):
            t = p*(2*j+1)-1

            terms = [R(0)] if i == 1 else []
            terms.extend(differentials[j])
            terms = [[y, R(0), R(0)] for y in terms]

            # reduce this row to the left
            triple = reduce_sequence(
              terms, [p*(m+1)-1 for m in range(len(terms))],
              0, 3,
              [[S(0),            S(0),            -2*b*x          ],
               [2*x - 3*t + 3,   S(0),            -a*(2*x - t + 1)],
               [S(0),            2*x - 3*t + 3,   S(0)            ]],
              2*x - 3*t + 3,
              p, R)

            # do one reduction to make it only two terms instead of three
            triple = reduce_one(triple, t//2, 3,
              [[ 9*b*(6*x-5), 2*a*a*(6*x-5), -3*a*b*(6*x-5)],
               [-6*a*(6*x-7),   9*b*(6*x-7),  2*a*a*(6*x-7)],
               [        S(0),          S(0),           S(0)]],
              (27*b**2 + 4*a**3)*(2*x-1), R)

            column.append([triple[0], triple[1]])

        # reduce the column to the bottom
        result = reduce_sequence(
              column, [(p*(2*m+1)-3)//2 for m in range(len(column))],
              0, 2,
              [[ 9*b*(6*x-5), 2*a*a*(6*x-5)],
               [-6*a*(6*x-7),   9*b*(6*x-7)]],
              (27*b**2 + 4*a**3)*(2*x-1),
              p, R)

        # undo stupid p^M adjustment, and reduce down to Z/p^N
        for k in range(2):
            result[k] = R_output(result[k].lift() / p**M)

        output.append(result)

    return matrix(R_output, 2, 2, [output[0][0], output[1][0],
                                   output[0][1], output[1][1]])





#*****************************************************************************
# This is a generalization of the above functionality for hyperelliptic curves.
#
# THIS IS A WORK IN PROGRESS.
#
# I tried to embed must stuff into the rings themselves rather than
# just extract and manipulate lists of coefficents. Hence the implementations
# below are much less optimized, so are much slower, but should hopefully be
# easier to follow. (E.g. one can print/make sense of intermediate results.)
#
# AUTHOR:
#    -- Robert Bradshaw (2007-04)
#
#*****************************************************************************


import weakref

from sage.schemes.hyperelliptic_curves.all import is_HyperellipticCurve
from sage.rings.padics.all import pAdicField
from sage.rings.all import QQ, is_LaurentSeries, is_LaurentSeriesRing, is_IntegralDomain
from sage.modules.all import FreeModule, is_FreeModuleElement

from sage.misc.profiler import Profiler
from sage.misc.misc import repr_lincomb

def matrix_of_frobenius_hyperelliptic(Q, p=None, prec=None, M=None):
    prof = Profiler()
    prof("setup")
    if p is None:
        try:
            K = Q.base_ring()
            p = K.prime()
            prec = K.precision_cap()
            #prec -= (adjusted_prec(p, prec) - prec)
        except AttributeError:
            raise ValueError, "p and prec must be specified if Q is not defined over a p-adic ring"
    if M is None:
        M = adjusted_prec(p, prec)
    extra_prec_ring = Integers(p**M) # pAdicField(p, M) # SLOW!
    real_prec_ring = pAdicField(p, prec) # pAdicField(p, prec) # To capped absolute?
    S = SpecialHyperellipticQuotientRing(Q, extra_prec_ring, True)
    MW = S.monsky_washnitzer()
    prof("frob basis elements")
    F = MW.frob_basis_elements(M, p)

    prof("rationalize")
    # do reduction over Q in case we have non-integral entries (and it's so much faster than padics)
    rational_S = S.change_ring(QQ)
    # this is a hack until pAdics are fast
    # (They are in the latest development bundle, but its not standard and I'd need to merge.
    # (it will periodically cast into this ring to reduce coefficent size)
    rational_S._prec_cap = p**M
    rational_S._p = p
#    S._p = p
#    rational_S(F[0]).reduce_fast()
#    prof("reduce others")
    F = [rational_S(F_i) for F_i in F]

    prof("reduce")
    reduced = [F_i.reduce_fast() for F_i in F]
#    reduced = [F_i.reduce() for F_i in F]

    #print reduced[0][0].diff() - F[0]

    # but the coeffs are WAY more precision than they need to be
    # print reduced[0][1].H1_vector()

    prof("make matrix")
    # now take care of precision capping
    M = matrix(real_prec_ring, [a.H1_vector() for f, a in reduced])
    M += real_prec_ring(0).add_bigoh(prec)
    print prof
#    print len(S._monomials)
    return M.transpose(), [f for f, a in reduced]




# For uniqueness (as many of the non-trivial calculations are cached along the way).

_special_ring_cache = {}
_mw_cache = {}

def SpecialHyperellipticQuotientRing(*args):
    if _special_ring_cache.has_key(args):
        R = _special_ring_cache[args]()
        if R is not None:
            return R
    R = SpecialHyperellipticQuotientRing_class(*args)
    _special_ring_cache[args] = weakref.ref(R)
    return R

def MonskyWashnitzerDifferentialRing(base_ring):
    if _mw_cache.has_key(base_ring):
        R = _mw_cache[base_ring]()
        if R is not None:
            return R

    R = MonskyWashnitzerDifferentialRing_class(base_ring)
    _mw_cache[base_ring] = weakref.ref(R)
    return R


class SpecialHyperellipticQuotientRing_class(CommutativeAlgebra):
    def __init__(self, Q, R=None, invert_y=True):
        if R is None:
            R = Q.base_ring()

        CommutativeAlgebra.__init__(self, R)

        x = PolynomialRing(R, 'x').gen(0)
        if is_EllipticCurve(Q):
            E = Q
            if E.a1() != 0 or E.a2() != 0:
                raise NotImplementedError, "Curve must be in Weierstrass normal form."
            Q = (-E.defining_polynomial()).change_ring(R)(x,0,1)

        elif is_HyperellipticCurve(Q):
            C = Q
            if C.hyperelliptic_polynomials()[1] != 0:
                raise NotImplementedError, "Curve must be of form y^2 = Q(x)."
            Q = C.hyperelliptic_polynomials()[0].change_ring(R)

        if is_Polynomial(Q):
            self._Q = Q.change_ring(R)
            self._coeffs = self._Q.coeffs()
            if self._coeffs.pop() != 1:
                raise NotImplementedError, "Polynomial must be monic."

        else:
            raise NotImplementedError, "Must be an elliptic curve or polynomial Q for y^2 = Q(x)"

        self._n = degree = int(Q.degree())

        self._series_ring = (LaurentSeriesRing if invert_y else PolynomialRing)(R, 'y')
        self._series_ring_y = self._series_ring.gen(0)
        self._series_ring_0 = self._series_ring(0)

        self._poly_ring = PolynomialRing(self._series_ring, 'x')

        self._x = self(self._poly_ring.gen(0))
        self._y = self(self._series_ring.gen(0))

        self._Q_coeffs = Q.change_ring(self._series_ring).list()
        self._dQ = Q.derivative().change_ring(self)(self._x)
        self._monsky_washnitzer = MonskyWashnitzerDifferentialRing(self)

        self._monomial_diffs = {}
        self._monomial_diff_coeffs = {}

    def _repr_(self):
        y_inverse = ",y^-1" if is_LaurentSeriesRing(self._series_ring) else ""
        return "SpecialHyperellipticQuotientRing K[x,y%s] / (y^2 = %s) over %s"%(y_inverse, self._Q, self.base_ring())

    def base_extend(self, R):
        if R.has_coerce_map_from(self.base_ring()):
            self.change_ring(R)
        else:
            raise TypeError, "no such base extension"

    def change_ring(self, R):
        return SpecialHyperellipticQuotientRing(self._Q, R, is_LaurentSeriesRing(self._series_ring))

    def __call__(self, val, offset=0):
        if isinstance(val, SpecialHyperellipticQuotientElement) and val.parent() is self:
            if offset == 0:
                return val
            else:
                return val << offset
        elif isinstance(val, MonskyWashnitzerDifferential):
            return self._monsky_washnitzer(val)
        return SpecialHyperellipticQuotientElement(self, val, offset)

    def gens(self):
        return self._x, self._y

    def x(self):
        return self._x

    def y(self):
        return self._y

    def monomial(self, i, j, b=None):
        """
        Returns $b y^j x^i$, computed quickly.
        """
        i = int(i)
        j = int(j)

        if 0 < i and i < self._n:
            if b is None:
                by_to_j = self._series_ring_y << (j-1)
            else:
                by_to_j = self._series_ring(b) << j
            v = [self._series_ring_0] * self._n
            v[i] = by_to_j
            return self(v)
        else:
            return (self._x ** i) << j if b is None else self.base_ring()(b) * (self._x ** i) << j

    def monomial_diff_coeffs(self, i, j):
        """
        The key here is that the formula for $d(x^iy^j)$ is messy
        in terms of i, but varies nicely with j.

        $d(x^iy^j) = y^{j-1} (2ix^{i-1}y^2 + j (A_i(x) + B_i(x)y^2)) \frac{dx}{2y}$

        Where $A,B$ have degree at most $n-1$ for each $i$.
        Pre-compute $A_i, B_i$ for each $i$ the "hard" way, and
        the rest are easy.
        """
        try:
            return self._monomial_diff_coeffs[i,j]
        except KeyError:
            pass
        if i < self._n:
            try:
                A, B, two_i_x_to_i = self._precomputed_diff_coeffs[i]
            except AttributeError:
                self._precomputed_diff_coeffs = self._precompute_monomial_diffs()
                A, B, two_i_x_to_i = self._precomputed_diff_coeffs[i]
            if i == 0:
                return j*A, j*B
            else:
                return j*A, j*B + two_i_x_to_i
        else:
            dg = self.monomial(i, j).diff()
            coeffs = [dg.extract_pow_y(j-1), dg.extract_pow_y(j+1)]
            self._monomial_diff_coeffs[i,j] = coeffs
            return coeffs

    def monomial_diff_coeffs_matrices(self):
        self.monomial_diff_coeffs(0, 0) # precompute stuff
        R = self.base_ring()
        mat_1 = matrix(R, self._n, self._n)
        mat_2 = matrix(R, self._n, self._n)
        for i in range(self._n):
            mat_1[i] = self._precomputed_diff_coeffs[i][1]
            mat_2[i] = self._precomputed_diff_coeffs[i][2]
        return mat_1.transpose(), mat_2.transpose()

    def _precompute_monomial_diffs(self):
        x, y = self.gens()
        R = self.base_ring()
        V = FreeModule(R, self.degree())
        As = []
        for i in range(self.degree()):
            dg = self.monomial(i, 1).diff()
            two_i_x_to_i = R(2*i) * x**(i-1) * y*y if i > 0 else self(0)
            A = dg - self._monsky_washnitzer(two_i_x_to_i)
            As.append( (V(A.extract_pow_y(0)), V(A.extract_pow_y(2)), V(two_i_x_to_i.extract_pow_y(2))) )
        return As


    def Q(self):
        return self._Q

    def degree(self):
        return self._n

    def prime(self):
        return self._p

    def monsky_washnitzer(self):
        return self._monsky_washnitzer

    def is_field(self):
        return False



class SpecialHyperellipticQuotientElement(CommutativeAlgebraElement):

    def __init__(self, parent, val=0, offset=0):
        CommutativeAlgebraElement.__init__(self, parent)
        if isinstance(val, tuple):
            val, offset = val
        if isinstance(val, SpecialHyperellipticQuotientElement):
            val, offset = val.coeffs()
        if isinstance(val, list) and len(val) > 0 and is_FreeModuleElement(val[0]):
            val = transpose_list(val)
        self._f = parent._poly_ring(val)
        if offset != 0:
            self._f = self._f.parent()([a << offset for a in self._f])

    def change_ring(self, R):
        return self.parent().change_ring(R)(self.coeffs())

    def __call__(self, *x):
        return self._f(*x)

    def __invert__(self):
        """
        The general element in our ring is not invertible, but y may be.
        We do not want to pass to the fraction field.
        """
        if self._f.degree() == 0 and self._f[0].is_unit():
            return SpecialHyperellipticQuotientElement(self.parent(), ~self._f[0])
        else:
            raise ZeroDivisionError, "Element not invertible"

    def is_zero(self):
        return self._f.is_zero()

    def __eq__(self, other):
        if not isinstance(other, SpecialHyperellipticQuotientElement):
            other = self.parent()(other)
        return self._f == other._f

    def _add_(self, other):
        return SpecialHyperellipticQuotientElement(self.parent(), self._f + other._f)

    def _sub_(self, other):
        return SpecialHyperellipticQuotientElement(self.parent(), self._f - other._f)

    def _mul_(self, other):
        # over laurent series, addition and subtraction can be expensive,
        # and the degree of this poly is small enough that Karatsuba actually hurts
        prod = self._f._mul_generic(other._f)
#        prod = self._f * other._f
        v = prod.list()
        parent = self.parent()
        Q_coeffs = parent._Q_coeffs
        n = len(Q_coeffs) - 1
        y2 = self.parent()._series_ring_y << 1
        for i in range(len(v)-1, n-1, -1):
            for j in range(n):
                v[i-n+j] -= Q_coeffs[j] * v[i]
            v[i-n] += y2 * v[i]
        return SpecialHyperellipticQuotientElement(parent, v[0:n])

    def _rmul_(self, c):
        coeffs = self._f.list()
        return self.parent()([c*a for a in coeffs])

    def _lmul_(self, c):
        coeffs = self._f.list()
        return self.parent()([a*c for a in coeffs])

    def __lshift__(self, k):
        coeffs = self._f.list()
        return self.parent()([a << k for a in coeffs])

    def __rshift__(self, k):
        coeffs = self._f.list()
        return self.parent()([a >> k for a in coeffs])

    def _repr_(self):
        x = PolynomialRing(QQ, 'x').gen(0)
        coeffs = self._f.list()
        return repr_lincomb([x**i for i in range(len(coeffs))], coeffs)

    def _latex_(self):
        x = PolynomialRing(QQ, 'x').gen(0)
        coeffs = self._f.list()
        return repr_lincomb([x**i for i in range(len(coeffs))], coeffs, is_latex=True)


    def diff(self):

#        try:
#            return self._diff_x
#        except AttributeError:
#            pass

        # d(self) = A dx + B dy
        #         = (2y A + BQ') dx/2y
        parent = self.parent()
        R = parent.base_ring()
        x, y = parent.gens()
        v = self._f.list()
        n = len(v)
        A = parent([R(i) * v[i] for i in range(1,n)])
        B = parent([a.derivative() for a in v])
        dQ = parent._dQ
        return parent._monsky_washnitzer( (R(2) * A << 1) + dQ * B )
#        self._diff = self.parent()._monsky_washnitzer( two_y * A + dQ * B )
#        return self._diff

    def extract_pow_y(self, k):
        v = [a[k] for a in self._f.list()]
        while len(v) < self.parent()._n:
            v.append(0)
        return v

    def min_pow_y(self):
        if self._f.degree() == -1:
            return 0
        return min([a.valuation() for a in self._f.list()])

    def max_pow_y(self):
        if self._f.degree() == -1:
            return 0
        return max([a.degree() for a in self._f.list()])

    def coeffs(self, R=None):
        zero = self.base_ring()(0) if R is None else R(0)
        y_offset = min(self.min_pow_y(), 0)
        y_degree = max(self.max_pow_y(), 0)
        coeffs = []
        for a in self._f.list():
            coeffs.append( [zero] * (a.valuation()-y_offset) + a.list() + [zero]*(y_degree - a.degree()) )
        while len(coeffs) < self.parent().degree():
            coeffs.append( [zero] * (y_degree - y_offset + 1) )
        V = FreeModule(self.base_ring() if R is None else R, self.parent().degree())
        coeffs = transpose_list(coeffs)
        return [V(a) for a in coeffs], y_offset



class MonskyWashnitzerDifferentialRing_class(Module):

    def __init__(self, base_ring):
        Module.__init__(self, base_ring)
        self._cache = {}

    def invariant_differential(self):
        return self(1)

    def __call__(self, val, offset=0):
        return MonskyWashnitzerDifferential(self, val, offset)

    def base_extend(self, R):
        return MonskyWashnitzerDifferentialRing(self.base_ring().base_extend(R))

    def change_ring(self, R):
        return MonskyWashnitzerDifferentialRing(self.base_ring().change_ring(R))

    def degree(self):
        return self.base_ring().degree()

    def Q(self):
        return self.base_ring().Q()

    def x_to_p(self, p):
        try:
            return self._cache["x_to_p", p]
        except KeyError:
            x_to_p = self.base_ring().x() ** p
            self._cache["x_to_p", p] = x_to_p
            return x_to_p

    def frob_Q(self, p):
        try:
            return self._cache["frobQ", p]
        except KeyError:
            x_to_p = self.x_to_p(p)
            frobQ = self.base_ring()._Q.change_ring(self.base_ring())(x_to_p)
            self._cache["frobQ", p] = frobQ
            return frobQ

    def frob_invariant_differential(self, prec, p):
        """
        $F_p(dx/y) = px^{p-1} y(F_py)^{-1} dx/y
                   = px^{p-1} y^{1-p} (1+pEy^{-2p})^{-1/2} dx/y
                   = px^{p-1} y^{1-p} (F_pQ y^{-p})^{-1/2} dx/y$

        Use Newton's method to calculate the square root.
        """
        prof = Profiler()
        prof("setup")
        # TODO, would it be useful to be able to take Frobenius of any element? Less efficient?
        x, y = self.base_ring().gens()
        prof("x_to_p")
        x_to_p_less_1 = x**(p-1)
        x_to_p = x*x_to_p_less_1

        # cache for future use
        self._cache["x_to_p", p] = x_to_p

        prof("frob_Q")
        a = self.frob_Q(p) >> 2*p  # frobQ * y^{-2p}

        prof("sqrt")
        Q = self.base_ring()._Q

#        three_halves = Q.parent().base_ring()(Rational((3,2)))
#        one_half = Q.parent().base_ring()(Rational((1,2)))
        three_halves = self.base_ring()._series_ring(Rational((3,2)))
        one_half     = self.base_ring()._series_ring(Rational((1,2)))

        # We are solving for t = a^{-1/2} = (F_pQ y^{-p})^{-1/2}
        # This converges because we know the root is in the same residue class as 1.

        t = self.base_ring()(1)
#        t = self.base_ring()(three_halves) - a._rmul_(one_half)

        print prec, "->", ceil(log(prec, 2))
        for _ in range(ceil(log(prec, 2))):
            t = t._rmul_(three_halves) - (t*t*t * a)._rmul_(one_half)  # t = (3/2) t - (1/2) a t^3
#            t = three_halves * t - one_half * t*t*t * a

#        print "a =", a
#        print "t =", t

#        prof("verify")
#        print "a*t^2 =", a * t**2

        prof("compose")
        F_dx_y = (p * x_to_p_less_1 * t) >> (p-1)  # px^{p-1} sqrt(a) * y^{-p+1}

#        print "-----", F_dx_y
#        print "-----", x_to_p * F_dx_y
        prof("done")
        print prof
        return MonskyWashnitzerDifferential(self, F_dx_y)

    def frob_basis_elements(self, prec, p):
        F_i = self.frob_invariant_differential(prec, p)
        x_to_p = self.x_to_p(p)
        F = [F_i]
        for i in range(1, self.degree()-1):
            F_i *= x_to_p
            F.append(F_i)
        return F

    def helper_matrix(self):
        """
        We use this to solve for the linear combination of $x^i y^j$ needed
        to clear all terms with $y^{j-1}$.
        """
        try:
            return self._helper_matrix
        except:
            AttributeError
        # The smallest y term of (1/j) d(x^i y^j) is constant for all j.
        L = []
        x, y = self.base_ring().gens()
        n = self.degree()
        for i in range(n):
            L.append( (y*x**i).diff().extract_pow_y(0) )
        A = matrix(L).transpose()
        if not is_IntegralDomain(A.base_ring()):
            # must be using integer_mod or something to approximate
            self._helper_matrix = (~A.change_ring(QQ)).change_ring(A.base_ring())
        else:
            self._helper_matrix = ~A
        return self._helper_matrix


class MonskyWashnitzerDifferential(ModuleElement):
    """
    Represents an element of the form F dx/2y
    """
    def __init__(self, parent, val=0, offset=0):
        ModuleElement.__init__(self, parent)
        if isinstance(val, MonskyWashnitzerDifferential):
            val = val._coeff.coeffs()
        self._coeff = self.parent().base_ring()(val, offset)

    def _add_(left, right):
        return MonskyWashnitzerDifferential(left.parent(),
                                            left._coeff + right._coeff)

    def _sub_(left, right):
        return MonskyWashnitzerDifferential(left.parent(),
                                            left._coeff - right._coeff)

    def __neg__(self):
        return MonskyWashnitzerDifferential(self.parent(), -self._coeff)

    def _lmul_(self, a):
        return MonskyWashnitzerDifferential(self.parent(), self._coeff * a)

    def _rmul_(self, a):
        return MonskyWashnitzerDifferential(self.parent(), a * self._coeff)

    def coeff(self):
        """
        This is a one-dimensional module over the base ring, generated by dx/2y.
        Return $A$ where $A dx/2y = self$.
        """
        return self._coeff

    def _repr_(self):
        s = self._coeff._repr_()
        if s.find("+") != -1 or s.find("-") != -1:
            s = "(%s)"%s
        return s + " dx/2y"

    def _latex_(self):
        s = self._coeff._latex_()
        if s.find("+") != -1 or s.find("-") != -1:
            s = "\\left(%s\\right)"%s
        return s + " \\frac{dx}{2y}"

    def extract_pow_y(self, k):
        """
        Really the power of y in A where self = A dx/2y.
        """
        return self._coeff.extract_pow_y(k)

    def min_pow_y(self):
        """
        Really the minimum power of y in A where self = A dx/2y.
        """
        return self._coeff.min_pow_y()

    def max_pow_y(self):
        """
        Really the maximum power of y in A where self = A dx/2y.
        """
        return self._coeff.max_pow_y()

    def reduce_neg_y(self):
        """
        Use homology relations to eliminate negative powers of y.
        """
        S = self.parent().base_ring()
        R = S.base_ring()
        M = self.parent().helper_matrix()
        p = S._p
        n = S.degree()
        x, y = S.gens()
        f = S(0)
        reduced = self
        for j in range(self.min_pow_y()+1, 0):
            if p.divides(j):
                cs = [R(QQ(a)/j) for a in reduced.extract_pow_y(j-1)]
            else:
                j_inverse = ~R(j)
                cs = [a*j_inverse for a in reduced.extract_pow_y(j-1)]
            lin_comb = M * vector(M.base_ring(), cs)
#            print "j =", j, "b =", cs, "lin_comb =", lin_comb
            g = self.parent().base_ring()(0)
            if not lin_comb.is_zero():
                for i in range(n):
                    if lin_comb[i] != 0:
                        g += S.monomial(i, j, lin_comb[i])
                if not g.is_zero():
                    f += g
                    reduced -= g.diff()
#                    print g, g.diff()
#                    print "reduced", reduced

        return f, reduced

    def reduce_neg_y_fast(self):
        """
        Use homology relations to eliminate negative powers of y.
        """
#        prof = Profiler()
#        prof("reduce setup")
        S = self.parent().base_ring()
        R = S.base_ring()
        M = self.parent().helper_matrix()

#        prof("extract coeffs")
        coeffs, offset = self.coeffs(R)
        V = coeffs[0].parent()

        if offset == 0:
            return S(0), self

#        prof("loop %s"%self.min_pow_y())
        forms = []
        for j in range(self.min_pow_y()+1, 0):
            if coeffs[j-offset-1].is_zero():
                forms.append(V(0))
            else:
                # this is a total hack to deal with the fact that we're using
                # rational numbers to approximate fixed precision p-adics
                if j % 3 == 1:
                    try:
                        v = coeffs[j-offset-1]
                        for kk in range(len(v)):
                            a = v[kk]
                            ppow = S._p**max(-a.valuation(S._p), 0)
                            v[kk] = ((a * ppow) % S._prec_cap) / ppow
                    except AttributeError:
                        pass
                lin_comb = ~R(j) * (M * coeffs[j-offset-1])
                forms.append(lin_comb)
                for i in lin_comb.nonzero_positions():
                    # g = lin_comb[i] x^i y^j
                    # self -= dg
                    coeffs[j-offset+1] -= lin_comb[i] * S.monomial_diff_coeffs(i, j)[1]

#        prof("recreate forms")
        f = S(forms, offset+1)
        reduced = S._monsky_washnitzer(coeffs[-1-offset:], -1)
#        print self - f.diff() - reduced
#        prof("done")
#        print prof
        return f, reduced

    def reduce_neg_y_faster(self):
        """
        Use homology relations to eliminate negative powers of y.
        """
        # Timings indicate this isn't any faster after all...

        S = self.parent().base_ring()
        R = S.base_ring()
        M = self.parent().helper_matrix()

        coeffs, offset = self.coeffs(R)
        V = coeffs[0].parent()

        if offset == 0:
            return S(0), self

        # See monomial_diff_coeffs
        # this is the B_i and x_to_i contributions respectively for all i
        d_mat_1, d_mat_2 = S.monomial_diff_coeffs_matrices()

        forms = []
        for j in range(self.min_pow_y()+1, 0):
            if coeffs[j-offset-1].is_zero():
                forms.append(V(0))
            else:
                # this is a total hack to deal with the fact that we're using
                # rational numbers to approximate fixed precision p-adics
                if j % 3 == 0:
                    try:
                        v = coeffs[j-offset-1]
                        for kk in range(len(v)):
                            a = v[kk]
                            ppow = S._p**max(-a.valuation(S._p), 0)
                            v[kk] = ((a * ppow) % S._prec_cap) / ppow
                    except AttributeError:
                        pass
                j_inverse =  ~R(j)
                lin_comb = (M * coeffs[j-offset-1])
                forms.append(j_inverse * lin_comb)
                coeffs[j-offset+1] -= (d_mat_1 + j_inverse * d_mat_2) * lin_comb

        f = S(forms, offset+1)
        reduced = S._monsky_washnitzer(coeffs[-1-offset:], -1)
#        reduced = self - f.diff()
        return f, reduced


    def reduce_pos_y(self):
        """
        Use homology relations to eliminate positive powers of y.
        """
        S = self.parent().base_ring()
        series = S.base_ring()
        n = S.Q().degree()
        p = S._p
        x, y = S.gens()
        f = S(0)
        reduced = self
        for j in range(self.max_pow_y(), 0, -1):
            for i in range(n-1, -1, -1):
                c = reduced.extract_pow_y(j)[i]
#                print "x^%s y^%s"%(i,j), c
                if c != 0:
                    g = S.monomial(0, j+1) if i == n-1 else S.monomial(i+1, j-1)
                    dg = g.diff()
#                    print reduced, " - ", dg
                    denom = dg.extract_pow_y(j)[i]
                    if p.divides(denom):
                        R = c.parent()
                        c = R(QQ(c)/QQ(denom))
                    else:
                        c /= denom
                    c = g.parent()(c)
                    f += c * g
                    reduced -= c * dg

        return f, reduced


    def reduce_pos_y_fast(self):
        """
        Use homology relations to eliminate positive powers of y.
        """
        S = self.parent().base_ring()
        R = S.base_ring()
        n = S.Q().degree()

        coeffs, offset = self.coeffs(R)
        V = coeffs[0].parent()
        forms = [V(0), V(0)]

        for j in range(self.max_pow_y(), -1, -1):

            form = V(0)
            i = n-1
            c = coeffs[j-offset][i]
            if c != 0:
                dg_coeffs = S.monomial_diff_coeffs(0, j+1)[0]
                c /= dg_coeffs[i]
                forms[len(forms)-2][0] = c
                # self -= c d(y^{j+1})
                coeffs[j-offset] -= c*dg_coeffs

            if j == 0:
                # the others are basis elements
                break

            for i in range(n-2, -1, -1):
                c = coeffs[j-offset][i]
                if c != 0:
                    dg_coeffs = S.monomial_diff_coeffs(i+1, j-1)
                    denom = dg_coeffs[1][i]
                    c /= denom
                    form[i+1] = c
                    # self -= c d(x^{i+1} y^{j-1})
                    coeffs[j-offset] -= c*dg_coeffs[1]
                    coeffs[j-offset-2] -= c*dg_coeffs[0]
            forms.append(form)

        forms.reverse()
        f = S(forms)
        reduced = self.parent()(coeffs[:1-offset], offset)
        return f, reduced

    def reduce(self):
        """
        Use homology relations to find $a$ and $f$ such that
        $self = a + df$ where $a$ is given in terms of the $x^i dx/2y$.
        """
#        print "max_pow_y = ", self.max_pow_y(), "min_pow_y = ", self.min_pow_y()
        n = self.parent().base_ring().Q().degree()
        f1, a = self.reduce_neg_y()
        f2, a = a.reduce_pos_y()
        f = f1 + f2

        c = a.extract_pow_y(0)[n-1]
        if c != 0:
            x, y = self.parent().base_ring().gens()
            g = y
            dg = g.diff()
            c = g.parent()(QQ(c)/QQ(dg.extract_pow_y(0)[n-1])) # TODO: fix when we have fast p-adics
            f += c * g
            a -= c * dg
#            print g, dg

        return f, a

    def reduce_fast(self):
        """
        Use homology relations to find $a$ and $f$ such that
        $self = a + df$ where $a$ is given in terms of the $x^i dx/2y$.
        """
#        print "max_pow_y = ", self.max_pow_y(), "min_pow_y = ", self.min_pow_y()
        f1, reduced = self.reduce_neg_y_fast()
        f2, reduced = reduced.reduce_pos_y_fast()

        return f1+f2, reduced


    def H1_vector(self):
        """
        Returns self as a element of $H^1$ under the basis $dx/2y, x dx/2y, ..., x^{n-1} dx/2y$.
        """
        f, reduced = self.reduce_fast()
        v = reduced.extract_pow_y(0)
        v.pop()
        return v

    def coeffs(self, R=None):
        return self._coeff.coeffs(R)

### end of file
