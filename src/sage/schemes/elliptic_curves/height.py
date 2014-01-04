##############################################################################
#       Copyright (C) 2010 Robert Bradshaw <robertwb@math.washington.edu>
#                     2014 John Cremona <john.cremona@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

import numpy
import math, bisect

from sage.rings.all import (ZZ, QQ, RR, RDF, RIF, CC, CDF, CIF,
    infinity, RealField, ComplexField)

from sage.misc.all import cached_method, cartesian_product_iterator
from sage.rings.arith import primes, lcm, factor, factorial
from sage.ext.fast_callable import fast_callable
from sage.functions.log import log, exp
from sage.symbolic.all import SR

from period_lattice_region import PeriodicRegion

class UnionOfIntervals:
    r"""
    A class which represents a union of closed intervals in
    `\RR` which can be scaled, shifted, intersected, etc.

    The intervals are represented as an ordered list of their
    endpoints, which may include $-\infty$ and $+\infty$.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
        sage: R = UnionOfIntervals([1,2,3,infinity]); R
        ([1, 2] U [3, +Infinity])
        sage: R + 5
        ([6, 7] U [8, +Infinity])
        sage: ~R
        ([-Infinity, 1] U [2, 3])
        sage: ~R | (10*R + 100)
        ([-Infinity, 1] U [2, 3] U [110, 120] U [130, +Infinity])
    """
    def __init__(self, endpoints):
        r"""
        An union of intervals is initialized by giving an increasing list
        of endpoints, the first of which may be $-\infty$ and the last of
        which may be $+\infty$.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: UnionOfIntervals([0,1])
            ([0, 1])
            sage: UnionOfIntervals([-infinity, pi, 17, infinity])
            ([-Infinity, pi] U [17, +Infinity])
            sage: UnionOfIntervals([])
            ()

            sage: UnionOfIntervals([1])
            Traceback (most recent call last):
            ...
            ValueError: an even number of endpoints must be given (got 1)
            sage: UnionOfIntervals([3,2,1,0])
            Traceback (most recent call last):
            ...
            ValueError: endpoints must be given in order
        """
        if len(endpoints) % 2:
            raise ValueError, "an even number of endpoints must be given (got %s)" % len(endpoints)
        if endpoints != sorted(endpoints):
            raise ValueError, "endpoints must be given in order"
        self._endpoints = endpoints

    def finite_endpoints(self):
        r"""
        Returns the finite endpoints of this union of intervals.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: UnionOfIntervals([0,1]).finite_endpoints()
            [0, 1]
            sage: UnionOfIntervals([-infinity, 0, 1, infinity]).finite_endpoints()
            [0, 1]
        """
        return [e for e in self._endpoints if -infinity < e < infinity]

    def intervals(self):
        r"""
        Returns the intervals in self, as a list of 2-tuples.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: UnionOfIntervals(range(10)).intervals()
            [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]
            sage: UnionOfIntervals([-infinity, pi, 17, infinity]).intervals()
            [(-Infinity, pi), (17, +Infinity)]
        """
        return zip(self._endpoints[::2], self._endpoints[1::2])

    def is_empty(self):
        r"""
        Returns whether self is empty.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: UnionOfIntervals([3,4]).is_empty()
            False
            sage: all = UnionOfIntervals([-infinity, infinity])
            sage: all.is_empty()
            False
            sage: (~all).is_empty()
            True
            sage: A = UnionOfIntervals([0,1]) & UnionOfIntervals([2,3])
            sage: A.is_empty()
            True
        """
        return not self._endpoints

    def __add__(left, right):
        r"""
        If both left an right are unions of intervals, take their union,
        otherwise treat the non-union of intervals as a scalar and shift.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([0, 1/2, 2, infinity]); A
            ([0, 1/2] U [2, +Infinity])
            sage: A + 1
            ([1, 3/2] U [3, +Infinity])
            sage: pi + A
            ([pi, pi + 1/2] U [pi + 2, +Infinity])
            sage: A + UnionOfIntervals([-infinity, -1])
            ([-Infinity, -1] U [0, 1/2] U [2, +Infinity])
        """
        if not isinstance(left, UnionOfIntervals):
            left, right = right, left
        elif not isinstance(right, UnionOfIntervals):
            return UnionOfIntervals([right + e for e in left._endpoints])
        else:
            return left.union([left, right])

    def __mul__(left, right):
        r"""
        Scale a union of inervals on the left or right.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([0, 1/2, 2, infinity]); A
            ([0, 1/2] U [2, +Infinity])
            sage: 2 * A
            ([0, 1] U [4, +Infinity])
            sage: A * 100
            ([0, 50] U [200, +Infinity])
            sage: 1.5 * A
            ([0.000000000000000, 0.750000000000000] U [3.00000000000000, +Infinity])
        """
        if not isinstance(right, UnionOfIntervals):
            return UnionOfIntervals([e*right for e in left._endpoints])
        elif not isinstance(left, UnionOfIntervals):
            return UnionOfIntervals([left*e for e in right._endpoints])
        else:
            return NotImplemented

    def __rmul__(self, other):
        r"""
        Scale by an operand on the left.

        TESTS::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([0, 1/2, 2, infinity]); A
            ([0, 1/2] U [2, +Infinity])
            sage: pi * A
            ([0, 1/2*pi] U [2*pi, +Infinity])
        """
        return self * other

    def __radd__(self, other):
        r"""
        Add a scalar operand on the left.

        TESTS::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([0, 1/2, 2, infinity]); A
            ([0, 1/2] U [2, +Infinity])
            sage: 100 + A
            ([100, 201/2] U [102, +Infinity])
        """
        return self + other

    def __invert__(self):
        r"""
        Return the closure of the compliment of self.

        .. NOTE::

        We take the closure because open intervals are not supported.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([0, 1/2, 2, infinity]); A
            ([0, 1/2] U [2, +Infinity])
            sage: ~A
            ([-Infinity, 0] U [1/2, 2])
            sage: A | ~A
            ([-Infinity, +Infinity])
            sage: A & ~A
            ([0, 0] U [1/2, 1/2] U [2, 2])
        """
        endpoints = list(self._endpoints)
        if endpoints[0] == -infinity:
            del endpoints[0]
        else:
            endpoints.insert(0, -infinity)
        if endpoints[-1] == infinity:
            endpoints.pop()
        else:
            endpoints.append(infinity)
        return UnionOfIntervals(endpoints)

    @staticmethod
    def join(L, condition):
        r"""
        Utility function to form the union or intersection of a list of UnionOfIntervals.

        INPUT:

        - ``L`` (list) -- a list of UnionOfIntervals instances

        - ``condition`` (function) -- either ``any`` or ``all``, or
          some other boolean function of a list of boolean values.

        OUTPUT:

        A new UnionOfIntervals instance representing the subset of
        '\R' equal to those reals in any/all/condition of the
        UnionOfIntervals in the list.

        .. NOTE::

        This is a static method for the class.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([1,3,5,7]); A
            ([1, 3] U [5, 7])
            sage: B = A+1; B
            ([2, 4] U [6, 8])
            sage: A.join([A,B],any) # union
            ([1, 4] U [5, 8])
            sage: A.join([A,B],all) # intersection
            ([2, 3] U [6, 7])
            sage: A.join([A,B],sum) # symmetric difference
            ([1, 2] U [3, 4] U [5, 6] U [7, 8])
        """
        all = []
        for ix, region in enumerate(L):
            for i, e in enumerate(region._endpoints):
                all.append((e, -(not (i % 2)), ix))
        all.sort()
        join = []
        in_join = False
        in_L = [False] * len(L)
        for e, start, ix in all:
            in_L[ix] = start
            if condition(in_L) != in_join:
                join.append(e)
                in_join = not in_join
        return UnionOfIntervals(join)

    @classmethod
    def union(cls, L):
        r"""
        Return the union of a list of UnionOfIntervals.

        INPUT:

        - ``L`` (list) -- a list of UnionOfIntervals instances

        OUTPUT:

        A new UnionOfIntervals instance representing the union of the
        UnionOfIntervals in the list.

        .. NOTE::

        This is a class method.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([1,3,5,7]); A
            ([1, 3] U [5, 7])
            sage: B = A+1; B
            ([2, 4] U [6, 8])
            sage: A.union([A,B])
            ([1, 4] U [5, 8])
        """
        return cls.join(L, any)

    @classmethod
    def intersection(cls, L):
        r"""
        Return the intersection of a list of UnionOfIntervals.

        INPUT:

        - ``L`` (list) -- a list of UnionOfIntervals instances

        OUTPUT:

        A new UnionOfIntervals instance representing the intersection
        of the UnionOfIntervals in the list.

        .. NOTE::

        This is a class method.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([1,3,5,7]); A
            ([1, 3] U [5, 7])
            sage: B = A+1; B
            ([2, 4] U [6, 8])
            sage: A.intersection([A,B])
            ([2, 3] U [6, 7])
        """
        for R in L:
            if R.is_empty():
                return R
        return cls.join(L, all)

    def __or__(left, right):
        r"""
        Return the union of a two UnionOfIntervals instances.

        INPUT:

        - ``left``, ``right`` (UnionOfIntervals) -- two UnionOfIntervals instances

        OUTPUT:

        A new UnionOfIntervals instance representing the union of ``left`` and ``right``.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([1,3,5,7]); A
            ([1, 3] U [5, 7])
            sage: B = A+1; B
            ([2, 4] U [6, 8])
            sage: A | B
            ([1, 4] U [5, 8])
        """
        return left.union([left, right])

    def __and__(left, right):
        r"""
        Return the intersection of a two UnionOfIntervals instances.

        INPUT:

        - ``left``, ``right`` (UnionOfIntervals) -- two UnionOfIntervals instances

        OUTPUT:

        A new UnionOfIntervals instance representing the intersection of ``left`` and ``right``.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([1,3,5,7]); A
            ([1, 3] U [5, 7])
            sage: B = A+1; B
            ([2, 4] U [6, 8])
            sage: A & B
            ([2, 3] U [6, 7])
        """
        return left.intersection([left, right])

    def __contains__(self, x):
        r"""
        Return True if ``x`` is in the UnionOfIntervals.

        INPUT:

        - ``x`` (real) -- a real number

        OUTPUT:

        Boolean: True if and only if ``x`` is in the union of intervals.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([1,3,5,7]); A
            ([1, 3] U [5, 7])
            sage: 1 in A
            True
            sage: 4 in A
            False
            sage: -infinity in A
            False
            sage: 'a' in A
            False
        """
        return x in self._endpoints or bisect.bisect_left(self._endpoints, x) % 2 == 1

    def __str__(self):
        r"""
        Return the string representation of this UnionOfIntervals.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([1,3,5,7])
            sage: str(A)
            '([1, 3] U [5, 7])'
        """
        return repr(self)

    def __repr__(self):
        r"""
        Return the string representation of this UnionOfIntervals.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import UnionOfIntervals
            sage: A = UnionOfIntervals([1,3,5,7]); A
            ([1, 3] U [5, 7])
        """
        return "(%s)" % " U ".join(str(list(I)) for I in self.intervals())

def nonneg_region(f):
    r"""
    Returns the UnionOfIntervals representing the region where ``f`` is non-negative.

    INPUT:

    - ``f`` (polynomial) -- a univariate polynomial over `\RR`.

    OUTPUT:

    A UnionOfIntervals representing the set `\{x \in\RR mid f(x) \ge 0\}`.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.height import nonneg_region
        sage: x = polygen(RR)
        sage: nonneg_region(x^2-1)
        ([-Infinity, -1.00000000000000] U [1.00000000000000, +Infinity])
        sage: nonneg_region(1-x^2)
        ([-1.00000000000000, 1.00000000000000])
        sage: nonneg_region(1-x^3)
        ([-Infinity, 1.00000000000000])
        sage: nonneg_region(x^3-1)
        ([1.00000000000000, +Infinity])
        sage: nonneg_region((x-1)*(x-2))
        ([-Infinity, 1.00000000000000] U [2.00000000000000, +Infinity])
        sage: nonneg_region(-(x-1)*(x-2))
        ([1.00000000000000, 2.00000000000000])
        sage: nonneg_region((x-1)*(x-2)*(x-3))
        ([1.00000000000000, 2.00000000000000] U [3.00000000000000, +Infinity])
        sage: nonneg_region(-(x-1)*(x-2)*(x-3))
        ([-Infinity, 1.00000000000000] U [2.00000000000000, 3.00000000000000])
        sage: nonneg_region(x^4+1)
        ([-Infinity, +Infinity])
        sage: nonneg_region(-x^4-1)
        ()
    """
    roots = f.roots()
    roots.sort()
    sign_changes = [r for r,e in roots if e%2 == 1]
    if (f.leading_coefficient() * (-1)**f.degree()) > 0:
        sign_changes = [-infinity] + sign_changes
    if f.leading_coefficient() > 0:
        sign_changes += [infinity]
    return UnionOfIntervals(sign_changes)

def inf_max_abs(f, g, D):
    r"""
    Returns `\inf_D(\max(|f|, |g|))`.

    INPUT:

    - ``f``, ``g`` (polynomials) -- real univariate polynomaials

    - ``D`` (UnionOfIntervals) -- a subset of `\RR`

    OUTPUT:

    A real number approximating the value of `\inf_D(\max(|f|, |g|))`.

    ALGORITHM:

    The extreme values must occur at an endpoint of a subinterval of
    `D` or at a point where one of `f`, `f'`, `g`, `g'`, `f\pm g` is
    zero.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.height import inf_max_abs, UnionOfIntervals
        sage: x = polygen(RR)
        sage: f = (x-10)^4+1
        sage: g = 2*x^3+100
        sage: inf_max_abs(f,g,UnionOfIntervals([1,2,3,4,5,6]))
        425.638201706391
        sage: r0 = (f-g).roots()[0][0]
        sage: r0
        5.46053402234697
        sage: max(abs(f(r0)),abs(g(r0)))
        425.638201706391

    """
    xs =  f.roots() + f.derivative().roots()
    xs += g.roots() + g.derivative().roots()
    xs += (f-g).roots() + (f+g).roots()
    xs = [r for r,e in xs if r in D]  # ignore multiplicities and points outside D
    xs += D.finite_endpoints()        # include endpoints of intervals
    if xs:
        return min([max(abs(f(r)), abs(g(r))) for r in xs])
    return infinity

def min_on_disk(f, tol, max_iter=10000):
    r"""
    Returns the minimum of a real-valued complex function on a square.

    INPUT:

    - ``f`` -- a function from CIF to RIF

    - ``tol`` (real) -- a positive real number

    - ``max_iter`` (integer, default 10000) -- a positive integer
      bounding the number of iterations to be used

    OUTPUT:

    A 2-tuple `(s,t)`, where `t=f(s)` and `s` is a CIF element
    contained in the disk `|z|\le1`, at which `f` takes its minumum
    value.
    """
    # L holds a list of 4-tuples (mfs, ds, s, in_disk) where s is a
    # subregion of the initial square, ds its relative diameter,
    # mfs=-f(s) (actually minus the lower bound on f(s)) and in_disk
    # is a flag indicating whether or not s is a subset of the unit
    # disk.

    # We store the negative of the lower bound on f(s) so that we can
    # use the bisect module to sort these 4-tuples.

    # Initially L contains one element, the whole unit box, which is
    # not contained in the unit square.

    s = CIF(RIF(-1,1), RIF(-1,1))
    fs = f(s)
    L = [(-fs.lower(), fs.relative_diameter(), s, False)]

    # min_max holds the minumum over L of fs.upper().

    min_max = fs.upper()

    # We iterate at most max_iter times.  At each step we look at the
    # best so far and return it if is good enough, meaning that its
    # relative diameter is less than the given tolerance; otherwise we
    # bisect this best region (into 4 pieces) and replace the entry in
    # L with at most 4 smaller entries.

    for k in range(max_iter):
        value, err, region, in_disk = L.pop()
        if err < tol:       # reached desired tolerance, so return
            return region, -value
        for s in region.bisection(): # 4 sub-regions
            if in_disk:
                s_in_disk = True     # if the original region si in the disk so are all its children
            else:
                 r = abs(s)          # otherwise we test each one
                 if r > 1:
                     continue        # skip this subregion if it is entirely outside the disk
                 s_in_disk = r < 1   # meaning it is entirely inside the disk

            fs = f(s)

            if fs.upper() < min_max: # we definitely beat the record
                min_max = fs.upper()
                unneeded = bisect.bisect(L, (-min_max,))
                if unneeded > 100:   # discard the worse entries (if there are many)
                    L = L[unneeded:]

            if fs.lower() < min_max: # we may beat the record, cannot yet tell: insert this region
                                     # into the list at the appropriate palce to maintain sorting
                bisect.insort(L, (-fs.lower(), fs.relative_diameter(), s, s_in_disk))

    # If we get here, then even after max_iter iterations the tolerance has not been reached.
    raise ValueError, "too many iterations"

two_pi_i_CDF = CDF(0, 2*RDF.pi())
two_pi_i_CIF = CIF(0, 2*RIF.pi())
i_CIF = CIF.gen()

# Ideas: We know tau, so we know the direction of the diagonal.
#        We can solve for x in p1, will this allow us to find the maxima exactly?

def rat_term_CIF(z, try_strict=True):
    r"""
    $rat_term(z) = u/(1-u)^2$ where $u=exp(2*pi*i*z)$.
    """
    two_pi_i_z = two_pi_i_CIF * z
    r = (two_pi_i_z.real()).exp() # = |u|
    x, y = two_pi_i_z.imag().cos(), two_pi_i_z.imag().sin()

    real_part = imag_part = None

    # If there are no local minima the intervals are
    # strictly determined by their values at the endpoints.
    if try_strict:

        corner_reals = []
        corner_imags = []
        for a, b in cartesian_product_iterator([z.real().endpoints(), z.imag().endpoints()]):
            zz = CDF(a,b)
            u = (two_pi_i_CDF*zz).exp()
            f = u/(1-u)**2
            corner_reals.append(f.real())
            corner_imags.append(f.imag())

        p1 = (((((r+2*x)*r - 6)*r + 2*x) * r) + 1)  # r^4 + 2*r^3*x - 6*r^2 + 2*r*x + 1
        p2 = (r*(x*(r+2*x)-4)+x)                    # r^2*x + 2*r*x^2 - 4*r + x

        df_dr = (r**2-1) * p2
        df_dx = p1 * r

        dg_dr = p1 * y
        dg_dx = r * df_dr / y

        if not dg_dr.contains_zero() or not dg_dx.contains_zero():
            real_part = RIF(min(corner_reals), max(corner_reals))

        if not dg_dr.contains_zero() or not dg_dx.contains_zero():
            imag_part = RIF(min(corner_imags), max(corner_imags))

    if real_part is None or imag_part is None:
        denom = (1-r*(2*x-r))**2
    if real_part is None:
        real_part = r*(x*(1+r**2)-2*r)/denom
    if imag_part is None:
        imag_part = -(r**2-1)*y*r/denom

    return CIF(real_part, imag_part)

def eps(err, is_real):
    r"""
    Return a Real or Complex interval centered on 0 with radius err.

    INPUT:

    - ``err`` (real) -- a positive real number, the radius of the interval

    - ``is_real`` (boolean) -- if True, returns a real interval in
      RIF, else a complex interval in CIF

    OUTPUT:

    An element of RIF or CIF (as specified), centered on 0, with given radius.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.height import eps
        sage: eps(0.01, True)
        0.0?
        sage: eps(0.01, False)
        0.0? + 0.0?*I
    """
    e = RIF(-err, err)
    if is_real:
        return e
    else:
        return CIF(e, e)


def prime_ideals_of_bounded_norm(K, B):
    if K is QQ:
        return primes(B+1)
    all = []
    for p in primes(B+1):
        for pp, e in K.factor(p):
            if pp.norm() <= B:
                all.append(pp)
    return all

class EllipticCurveCanonicalHeight:
    r""" Class for computing canonical heights of points on elliptic
    curves defined over number fields, including rigorous lower bounds for
    the canonical height of non-torsion points.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.height import EllipticCurveCanonicalHeight
        sage: E = EllipticCurve([0,0,0,0,1])
        sage: EllipticCurveCanonicalHeight(E)
        EllipticCurveCanonicalHeight object associated to Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field
    """

    def __init__(self, E):
        r"""
        Initialize the class with an elliptic curve.

        INPUT:

        - `E` -- an elliptic curve defined over a number field

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import EllipticCurveCanonicalHeight
            sage: E = EllipticCurve([0,0,0,0,1])
            sage: EllipticCurveCanonicalHeight(E)
            EllipticCurveCanonicalHeight object associated to Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field

        TESTS:

        The base field must be a number field::

            sage: from sage.schemes.elliptic_curves.height import EllipticCurveCanonicalHeight
            sage: E = EllipticCurve(GF(7),[0,0,0,0,1])
            sage: EllipticCurveCanonicalHeight(E)
            Traceback (most recent call last):
            ...
            ValueError: EllipticCurveCanonicalHeight class can only be created from an elliptic curve defined over a number field
        """
        from sage.schemes.elliptic_curves.ell_generic import is_EllipticCurve
        if is_EllipticCurve(E):
            self.E = E
            from sage.rings.number_field.number_field_base import is_NumberField
            K = E.base_ring()
            if is_NumberField(K):
                self.K = K
            else:
                raise ValueError, "EllipticCurveCanonicalHeight class can only be created from an elliptic curve defined over a number field"
        else:
            raise ValueError, "EllipticCurveCanonicalHeight class can only be created from an elliptic curve"

    def __repr__(self):
        r"""
        Return the string representation.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import EllipticCurveCanonicalHeight
            sage: E = EllipticCurve([0,0,0,0,1])
            sage: EllipticCurveCanonicalHeight(E)
            EllipticCurveCanonicalHeight object associated to Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field
        """
        return "EllipticCurveCanonicalHeight object associated to %s" % self.E

    def curve(self):
        r"""
        Return the elliptic curve.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import EllipticCurveCanonicalHeight
            sage: E = EllipticCurve([0,0,0,0,1])
            sage: H = EllipticCurveCanonicalHeight(E)
            sage: H.curve()
            Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field
        """
        return self.E

    def base_field(self):
        r"""
        Return the base field.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import EllipticCurveCanonicalHeight
            sage: E = EllipticCurve([0,0,0,0,1])
            sage: H = EllipticCurveCanonicalHeight(E)
            sage: H.base_field()
            Rational Field
        """
        return self.K

    def __call__(self, P):
        r"""
        Return the canonical height of the point ``P``.

        INPUT:

        - ``P`` -- a point on the elliptic curve

        OUTPUT:

        The canonical height of ``P``.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.height import EllipticCurveCanonicalHeight
            sage: E = EllipticCurve([0,0,1,-1,0])
            sage: P = E(0,0)
            sage: P.height()
            0.0511114082399688
            sage: H = EllipticCurveCanonicalHeight(E)
            sage: H(P)
            0.0511114082399688
            sage: H([0,0])
            0.0511114082399688
            sage: H((0,0))
            0.0511114082399688

        Over a number field other than `\QQ`::

            sage: from sage.schemes.elliptic_curves.height import EllipticCurveCanonicalHeight
            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve(K, [0,0,0,1,-27])
            sage: H = EllipticCurveCanonicalHeight(E)
            sage: H.base_field()
            Number Field in i with defining polynomial x^2 + 1
            sage: H((1,5*i))
            1.22257115164148
        """
        return self.E(P).height()

    @cached_method
    def alpha(self, v, tol=0.01):
        r"""
        Return the constant `\alpha_v` associated to the embedding ``v``.

        INPUT:

        - ``v`` -- an embedding of the base field into `\RR` or `\CC`

        OUTPUT:

        The constant `\alpha_v`.  In the notation of Cremona,
        Prickett, Siksek (2006), `\alpha_v^3=\epsilon_v`.
        """
        b2, b4, b6, b8 = [v(b) for b in self.E.b_invariants()]
        x = v.codomain()['x'].gen()
        f = 4*x**3 + b2*x**2 + 2*b4*x + b6
        g = x**4 - b4*x**2 - 2*b6*x - b8
        F = f.reverse() << (4-f.degree())
        G = g.reverse() << (4-g.degree())

        if v(self.K.gen()) in RR:
            I = UnionOfIntervals([-1,1])
            min_fg = inf_max_abs(f, g, nonneg_region(f) & I)
            min_FG = inf_max_abs(F, G, nonneg_region(F) & I)
            return min(min_fg, min_FG) ** (-1/QQ(3))

        else:
            def pair_max(f, g):
                f = f.change_ring(CIF)
                g = g.change_ring(CIF)
                max = type(RIF(0)).max
                def max_f_g(z):
                    return max(abs(f(z)), abs(g(z)))
                return max_f_g
            pair_max_old = pair_max
            def pair_max(f, g):
                f = f.change_ring(CDF)
                g = g.change_ring(CDF)
                dfn = [fast_callable(f.derivative(n)/factorial(n), CDF) for n in range(f.degree()+1)]
                dgn = [fast_callable(g.derivative(n)/factorial(n), CDF) for n in range(g.degree()+1)]
                def max_f_g(s):
                    (a,b),(c,d) = s.real().endpoints(), s.imag().endpoints()
                    dx = a-b; dy = c-d
                    eta = RDF(dx*dx + dy*dy).sqrt()
                    z = CDF(s.center())
                    err_f = sum(eta ** n * abs(df(z)) for n, df in enumerate(dfn) if n)
                    err_g = sum(eta ** n * abs(dg(z)) for n, dg in enumerate(dgn) if n)
                    return RIF(max(abs(f(z)), abs(g(z)))) + eps(max(err_f, err_g), True)
                return max_f_g
            _, min_fg = min_on_disk(pair_max(f, g), tol)
            _, min_FG = min_on_disk(pair_max(F, G), tol)
            return min(min_fg, min_FG) ** (-1/QQ(3))

    @cached_method
    def e_p(self, p):
        kp = self.K.residue_field(p)
        if self.E.has_bad_reduction(p):
            if self.E.has_additive_reduction(p):
                ep = kp.characteristic()
            elif self.E.has_split_multiplicative_reduction(p):
                ep = len(kp) - 1
            else:
                ep = len(kp) + 1
        else:
            ep = self.E.reduction(p).abelian_group().exponent()
        return ZZ(ep)

    @cached_method
    def DE(self, n):
        s = 0
        N = self.E.conductor()
        for p in prime_ideals_of_bounded_norm(self.K, (n+1) ** max(2, self.K.degree())):
            ep = self.e_p(p)
            if ep.divides(n):
                kp = self.K.residue_field(p)
                s += 2*(1+(n/ep).valuation(kp.characteristic())) * log(len(kp))
        return s

    @cached_method
    def ME(self):
        from sage.misc.misc import prod
        return prod([p ** (e - self.E.local_data(p).discriminant_valuation()) for p, e in self.E.discriminant().factor()], self.K.one_element())

    def B(self, n, mu, verbose=False):
        K = self.K
        B = exp(K.degree() * n**2 * mu - self.DE(n)) / self.ME().norm() ** 6
        for v in K.places():
            if v(K.gen()) in RR:
                B *= self.alpha(v)
            else:
                B *= self.alpha(v) ** 2
        return B

    # Empty real intersection detection.

    def psi(self, xi, v):
        if xi > 1e9:
            return 1
        L = self.E.period_lattice(v)
        w1, w2 = L.basis()
        t = L(self.E.lift_x(xi)) / w1
        if t < 0.5:
            t = 1 - t
        return t

    def S(self, xi1, xi2, v):
        L = self.E.period_lattice(v)
        w1, w2 = L.basis()
        beta = L.elliptic_exponential(w1/2)[0]
        if xi2 < beta:
            return UnionOfIntervals([])
        elif xi1 < beta <= xi2:
            a = self.psi(xi2, v)
            return UnionOfIntervals([1-a, a])
        else:
            a, b = self.psi(xi1, v), self.psi(xi2, v)
            return UnionOfIntervals([1-b, 1-a, a, b])

    def Sn(self, xi1, xi2, n, v):
        SS = 1/ZZ(n) * self.S(xi1, xi2, v)
        return UnionOfIntervals.union([t/ZZ(n) + SS for t in range(n)])

    def real_intersection_is_empty(self, Bk, v):
        return UnionOfIntervals.intersection([self.Sn(-B, B, k+1, v) for k,B in enumerate(Bk)]).is_empty()

    # Empty complex intersection detection.

#    @cached_method
    def tau(self, v):
        w1, w2 = self.E.period_lattice(v).basis()
        return w2/w1

    def wp_c(self, v):
        r"""
        Given the recurance relations for the laurent series expansion of wp,
        it is easy to see that there is a constant c such that

            $$|\wp(z) - z^-2| \le \frac{c^2|z|^2}{1-c|z|^2}$$

        whenever $c|z|^2<1$. We return this constant.
        """
        K = self.E.base_ring()
        w1, w2 = self.E.period_lattice(v).basis()
        series = list(self.E._pari_().ellwp(n=5))
        c2, c3 = v(K(series[2])) * abs(w1) ** 4, v(K(series[4])) * abs(w1) ** 6
        return max(RR(abs(c2)) ** 0.5, RR(abs(c3)) ** (1.0/3))

    def fk_intervals(self, v=None, N=20, domain=CIF):
        r"""
        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: L = E.period_lattice()
            sage: w1, w2 = L.basis()
            sage: z = CDF(0.3, 0.4)
            sage: L.elliptic_exponential(z*w1, to_curve=False)[0] * w1 ** 2
            -0.284513926792076 - 0.269698913503063*I
            sage: fk, err = E.height_function().fk_intervals(N=10)
            sage: fk(CIF(z))
            -0.284513926792076? - 0.269698913503064?*I
            sage: err(CIF(z))
            2.58152201135278e-20

            sage: fk, err = E.height_function().fk_intervals(N=10, domain=CDF)
            sage: fk(z)
            -0.284513926792 - 0.269698913503*I
        """
        if v is None:
            if self.K is QQ:
                v = QQ.hom(RR)
            else:
                raise ValueError, "must specify embedding"
        # pre-compute some constants
        tau = self.tau(v)
        const_term = 1/CC(12)
        qn = q = (2 * CC(-1).sqrt() * CC.pi() * tau).exp()
        for n in range(1, N):
            const_term -= 2 * qn/(1-qn) ** 2
            qn *= q

        two_pi_i = 2 * domain(-1).sqrt() * domain.pi()
        neg_four_pi2 = -4 * domain.pi() ** 2
        const_term = domain(const_term)
        tau = domain(tau)

        abs_q = abs(domain(q))
        abs_qN = abs(domain(qn))
        err_factor = abs(neg_four_pi2) / (1-abs_q)
        err_term = 2*abs_qN/(1-abs_qN) ** 2

        # choose u/(1-u)^2 evaluation method
        if domain is CIF:
            rat_term = rat_term_CIF
        else:
            def rat_term(z):
                u = (two_pi_i*z).exp()
                return u/(1-u)**2

        # the actual series
        def fk(z):
            s = const_term + rat_term(z)
            for n in range(1, N):
                s += rat_term(n*tau+z)
                s += rat_term(n*tau-z)
            return s * neg_four_pi2

        # the error function
        def err(z):
            alpha = z.imag() / tau.imag()
            qNa = abs_q**(N+alpha)
            qNai = abs_q**(N-alpha)
            return (err_factor * (qNa/(1-qNa) ** 2 + qNai/(1-qNai) ** 2 + err_term)).upper()

        return fk, err

    @cached_method
    def wp_intervals(self, v=None, N=20, abs_only=False):
        r"""
        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: wp = E.height_function().wp_intervals()
            sage: z = CDF(0.3, 0.4)
            sage: wp(CIF(z))
            -0.28451392679208? - 0.26969891350307?*I

            sage: L = E.period_lattice()
            sage: w1, w2 = L.basis()
            sage: L.elliptic_exponential(z*w1, to_curve=False)[0] * w1^2
            -0.284513926792076 - 0.269698913503063*I

            sage: z = CDF(0.3, 0.1)
            sage: wp(CIF(z))
            9.1896706142338? - 5.2087270125357?*I
            sage: L.elliptic_exponential(z*w1, to_curve=False)[0] * w1^2
            9.18967061423377 - 5.20872701253569*I

        """
        if v is None:
            if self.K is QQ:
                v = QQ.hom(RR)
            else:
                raise ValueError, "must specify embedding"
        tau = self.tau(v)
        fk, fk_err = self.fk_intervals(v, N)
        c = self.wp_c(v)

        def wp(z):

            # center around origin
            offset = (z.imag().lower() / tau.imag()).round()
            if offset:
                z -= CIF(offset * tau)
            offset = z.real().lower().round()
            if offset:
                z -= offset

            # estimate using the series
            approx = fk(z)
            err = fk_err(z)
            if abs_only:
                approx = abs(approx)
            approx += eps(err, abs_only)
    #        print "fk_approx", approx

            # refine using an esitmate that's better near the pole
            z_bound = abs(z).upper()
            cz2 = c * z_bound ** 2
            if cz2 < 1:
                err = (c * cz2) / (1 - cz2)
                if abs_only:
                    pole_approx = abs(z) ** -2
                else:
                    pole_approx = z ** -2
    #            print "pole approx", pole_approx + eps(err, abs_only)
    #            print approx in approx.intersection(pole_approx + eps(err, abs_only))
                approx = approx.intersection(pole_approx + eps(err, abs_only))

            return approx

        return wp

    @cached_method
    def wp_on_grid(self, v, N, half=False):
        tau = self.tau(v)
        fk, err = self.fk_intervals(v, 15, CDF)
        var_z = SR.var('z')
        ff = fast_callable(fk(var_z), CDF, [var_z])
        N_or_half = N // (1+half)
        vals = numpy.empty((N,N_or_half))
        for i in range(N):
            for j in range(N_or_half):
                vals[i,j] = abs(ff((i+.5)/N + (j+.5)*tau/N))
        return vals

    def complex_intersection_is_empty(self, Bk, v, verbose=False, use_half=True):

        b2 = v(self.E.b2())
        w1, w2 = self.E.period_lattice(v).basis()
        tau = w2/w1
        bounds = [RR((B.sqrt() + abs(b2)/12) * abs(w1) ** 2) for B in Bk]
        vals = self.wp_on_grid(v, 30, half=use_half)
        wp = self.wp_intervals(v, abs_only=True)

        k = len(bounds)

        # First try and prove a negative result (cheap).
        if verbose:
            print "trying to prove negative result..."
        intersection = None
        for B, n in sorted(zip(bounds, ZZ.range(1, k+1))):
            T = PeriodicRegion(CDF(1), CDF(tau), vals < B, full=not use_half)
            if intersection is None:
                intersection = PeriodicRegion(CDF(1), CDF(tau), vals < B, full=not use_half)
            else:
                intersection &= T/n
                if intersection.is_empty():
                    break
        else:
            z = CIF(intersection.innermost_point())
            if all(wp((k+1)*z) < B for k, B in enumerate(bounds)):
                return False

        # Now try to prove a positive result.
        if verbose:
            print "trying to prove positive result..."
        intersection = None
        for B, n in sorted(zip(bounds, ZZ.range(1, k+1))):

            T = PeriodicRegion(CDF(1), CDF(tau), vals < B, full=not use_half).expand().refine()
            leaning_right = tau.real() / tau.imag() >= 0
            def check_line(z):
                wpz = wp(z)
                if wpz > B:
                    return True
                # Try refining once before we declare failure.
                z00, z01, z10, z11 = z.bisection()
                if leaning_right:
                    start, end = z00, z11
                else:
                    start, end = z01, z10
                if wp(start) > B and wp(end) > B:
                    return True
                return False

            # This step here is the bottleneck.
            while not T.verify(check_line):
                print "bad"
                T = T.expand()
            if intersection is None:
                intersection = T
            else:
                intersection &= T/n
                if intersection.is_empty():
                    return True

        return False

    def test_mu(self, mu, N, verbose=True):
        Bk = []
        for n in ZZ.range(1, N+1):
            b = self.B(n, mu)
            b = RR(b)
            if verbose:
                print "B_%s(%s) = %s" % (n, mu, b)
            if b < 1:
               return True
            Bk.append(b)

        for v in self.K.places():
            if v(self.K.gen()) in RR:
                if self.real_intersection_is_empty(Bk, v):
                    return True
            else:
                if self.complex_intersection_is_empty(Bk, v):
                    return True
        return False # Couldn't prove it...

    def min_gr(self, tol, n_max, verbose=False):
        r"""
        EXAMPLES::

        Example from Cremona and Siksek "Computing a Lower Bound for the
        Canonical Height on Elliptic Curves over Q"::

            sage: E = EllipticCurve([1, 0, 1, 421152067, 105484554028056]) # 60490d1
            sage: E.height_function().min_gr(.0001, 5)
            1.98684388147
        """
        test = self.test_mu
        if test(1, n_max, verbose):
            mu = 2
            while test(mu, n_max, False):
                mu *= 2
            mu /= 2
        else:
            mu = .5
            while not test(mu, n_max, False):
                mu /= 2
        # The true value lies between mu and eps * mu.
        eps = 2.0
        while eps > tol+1:
            if verbose:
                print "height bound in [%s, %s]" % (mu, mu*eps)
            eps = math.sqrt(eps)
            if test(mu*eps, n_max, False):
                mu = mu*eps
        return RDF(mu)

    def min(self, tol, n_max, verbose=False):
        r"""
        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: h = E.height_function()
            sage: h.min(.01, 5)
            0.0398731805749
            sage: E.gen(0).height()
            0.0511114082399688
            sage: K.<a> = QuadraticField(-5)
            sage: E.change_ring(K).height_function().min(0.1, 5)
            0.0170391833229

            sage: E = EllipticCurve('389a')
            sage: h = E.height_function()
            sage: h.min(0.1, 5)
            0.0573127527003
            sage: [P.height() for P in E.gens()]
            [0.686667083305587, 0.327000773651605]

            sage: E = EllipticCurve([1, 0, 1, 421152067, 105484554028056]) # 60490d1
            sage: E.height_function().min(.01, 5)
            0.0011215737111

            sage: K.<i> = QuadraticField(-1)
            sage: E = EllipticCurve([91 - 26*i, -144 - 323*i])
            sage: E.height_function().min(.1, 4)
            0.162104944331
        """
        tp = [L.tamagawa_number() for L in self.E.local_data()]
        if self.K == QQ:
            tp.append(self.E.real_components())
        elif self.K.real_places():
            tp.append([2]) # only sometimes?
        return self.min_gr(tol, n_max, verbose) / lcm(tp) ** 2
