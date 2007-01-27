r"""
Computation of Frobenius matrix on Monsky-Washnitzer cohomology.

The most interesting functions to be exported here are matrix_of_frobenius()
and adjusted_prec().

Currently this code is limited to the case $p \geq 5$ (no $GF(p^n)$
for $n > 1$), and only handles the elliptic curve case (not more
general hyperelliptic curves).

REFERENCES:
  -- Kedlaya, K., ``Counting points on hyperelliptic curves using
     Monsky-Washnitzer cohomology'', J. Ramanujan Math. Soc. 16 (2001)
     no 4, 323--338
  -- Edixhoven, B., ``Point counting after Kedlaya'', EIDMA-Stieltjes graduate
     course, Lieden (lecture notes?).

AUTHORS:
    -- David Harvey and Robert Bradshaw (initial code developed at the 2006
       MSRI graduate workshop, working with Jennifer Balakrishnan and Liang
       Xiao)
    -- David Harvey (Aug/Sep 2006): cleaned up, rewrote some chunks, lots
       more documentation, added Newton iteration method, added more complete
       "trace trick", integrated better into SAGE.

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 Robert Bradshaw <robertwb@math.washington.edu>
#                     2006 David Harvey <dmharvey@math.harvard.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.all import Integers, Integer, PolynomialRing, is_Polynomial
from sage.matrix.all import matrix
from sage.rings.ring import CommutativeAlgebra
from sage.structure.element import CommutativeAlgebraElement

from sage.rings.arith import binomial
from sage.misc.functional import log, ceil


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

  def __init__(self, Q):
    """
    Constructor.

    INPUT:
        Q -- a polynomial of the form Q(x) = x^3 + ax + b, where a, b
             belong to a ring in which 2, 3 are invertible.
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
    self._poly_ring = PolynomialRing(base_ring, 'T')    # R[T]
    self._poly_generator = self._poly_ring.gen(0)    # the generator T

    # Precompute a matrix that is used in the Toom-Cook multiplication.
    # This is where we need 2 and 3 invertible.

    # (a good description of Toom-Cook is online at:
    # http://www.gnu.org/software/gmp/manual/html_node/Toom-Cook-3-Way-Multiplication.html)

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


  def __mul__(self, other):
    # todo: cache results of toom-cook splitting, i.e. often we multiply
    # the same polynomial by a bunch of other things, and we can save
    # on part of the repetitve work

    # todo: this should really be _mul_, not __mul__. But currently of the
    # _mul_ dispatch code is a bit messed up and it doesn't seem to work
    # if I make it _mul_. This is very frustrating. Maybe this will work
    # properly when PolynomialRings become algebras rather than just
    # commutative rings.

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
    y = parent._poly_generator
    b = parent._b
    a = parent._a

    # todo: These lines are necessary to get binop stuff working
    # for certain base rings, e.g. when we compute b*c3 in the
    # final line. They shouldn't be necessary. Need to fix this
    # somewhere else in SAGE.
    a = parent._poly_ring(a)
    b = parent._poly_ring(b)

    return SpecialCubicQuotientRingElement(parent,
                                           -b*c3 + c0 + c3*y,
                                           -b*c4 - a*c3 + c1 + c4*y,
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



def reduce_negative(Q, p, coeffs, offset):
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

            a[0] = base_ring(a[0].lift() / (j+1))
            a[1] = base_ring(a[1].lift() / (j+1))
            a[2] = base_ring(a[2].lift() / (j+1))

            c1 = m[3]*a[0] + m[4]*a[1] + m[5]*a[2]
            c2 = m[6]*a[0] + m[7]*a[1] + m[8]*a[2]
            next_a[0] = next_a[0] - three_j_plus_5 * c1
            next_a[1] = next_a[1] - three_j_plus_7 * c2

            three_j_plus_7 = three_j_plus_7 + six
            three_j_plus_5 = three_j_plus_5 + six

    except NotImplementedError:
        raise NotImplementedError, \
            "It looks like you've found a non-integral matrix of Frobenius! " \
            "(Q=%s, p=%s)\nTime to write a paper." % (Q, p)

    coeffs[int(offset)] = next_a



def reduce_positive(Q, coeffs, offset):
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
        sage: monsky_washnitzer.reduce_positive(Q, coeffs, 0)
        sage: coeffs[0]
         [16, 102, 88]

        sage: coeffs = [[9, 8, 7], [10, 15, 20]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_positive(Q, coeffs, 0)
        sage: coeffs[0]
         [24, 108, 92]

    """

    base_ring = Q.base_ring()
    next_a = coeffs[len(coeffs) - 1]

    Qa = Q[1]
    Qb = Q[0]

    A = 2*Qa
    B = 3*Qb

    for i in range(len(coeffs)-1, offset, -1):
        j = 2*(i-offset) - 2
        a = next_a
        next_a = coeffs[i-1]

        a[0] = a[0] - Qa*a[2]/3   # subtract d(y^j + 1)

        # todo: see comments about pAdicInteger in reduceNegative()

        # subtract off c1 of d(x y^j + 1)
        c1 = base_ring(a[0].lift() * (j+1) / (3*j + 5))
        # subtract off c2 of d(x^2 y^j + 1)
        c2 = base_ring(a[1].lift() * (j+1) / (3*j + 7))

        next_a[0] = next_a[0] + B*c1
        next_a[1] = next_a[1] + A*c1 + B*c2
        next_a[2] = next_a[2]        + A*c2

    coeffs[int(offset)] = next_a



def reduce_zero(Q, coeffs, offset):
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
      return

    Qa = Q[1]

    a[0] = a[0] - a[2]*Qa/3    # $3x^2 dx/y = -a dx/y$
    a[2] = 0

    coeffs[int(offset)] = a



def reduce_all(Q, p, coeffs, offset):
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

    while len(coeffs) <= offset:
        coeffs.append([R(0), R(0), R(0)])

    reduce_negative(Q, p, coeffs, offset)
    reduce_positive(Q, coeffs, offset)
    reduce_zero(Q, coeffs, offset)

    return coeffs[int(offset)][0], coeffs[int(offset)][1]



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



def matrix_of_frobenius(Q, p, M, trace=None):
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

  if M == 2:
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

  return matrix(base_ring, 2, 2, [F0_reduced[0], F1_reduced[0],
                                  F0_reduced[1], F1_reduced[1]])


### end of file
