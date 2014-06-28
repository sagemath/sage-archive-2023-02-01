r"""
Hyperelliptic curves over a finite field

AUTHORS:

- David Kohel (2006)

- Robert Bradshaw (2007)

- Alyson Deines, Marina Gresham, Gagan Sekhon, (2010)

- Jean-Pierre Flori, Jan Tuitman (2013)

EXAMPLES::

    sage: K.<a> = GF(9, 'a')
    sage: x = polygen(K)
    sage: C = HyperellipticCurve(x^7 - x^5 - 2, x^2 + a)
    sage: C._points_fast_sqrt()
    [(0 : 1 : 0), (a + 1 : a : 1), (a + 1 : a + 1 : 1), (2 : a + 1 : 1), (2*a : 2*a + 2 : 1), (2*a : 2*a : 1), (1 : a + 1 : 1)]
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#  Copyright (C) 2010 Alyson Deines <aly.deines@gmail.com>, Marina Gresham
#  <marina.gresham@coloradocollege.edu>, Gagan Sekhon <gagan.d.sekhon@gmail.com>
#  Copyright (C) 2013 Jean-Pierre Flori <jean-pierre.flori@ssi.gouv.fr>,
#  Jan Tuitman <jan.tuitman@wis.kuleuven.be>
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
#*****************************************************************************

from sage.rings.all import ZZ, RR, QQ, GF
from sage.rings.arith import binomial
from sage.rings.power_series_ring import PowerSeriesRing
import hyperelliptic_generic
from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import identity_matrix, matrix
from sage.misc.functional import rank


class HyperellipticCurve_finite_field(hyperelliptic_generic.HyperellipticCurve_generic):
    def _frobenius_coefficient_bound_charpoly(self):
        """
        Computes bound on number of ``p``-adic digits needed to recover
        frobenius polynomial computing the characteristic polynomial
        of the frobenius matrix, i.e. returns ``B`` so that knowledge of
        ``a_1``, ..., ``a_g`` modulo ``p^B`` determine frobenius polynomial
        uniquely.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound_charpoly()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound_charpoly()
            2
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound_charpoly()
            3

            sage: R.<t> = PolynomialRing(GF(next_prime(10^9)))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound_charpoly()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound_charpoly()
            2
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound_charpoly()
            2
            sage: HyperellipticCurve(t^9 + t + 1)._frobenius_coefficient_bound_charpoly()
            3
            sage: HyperellipticCurve(t^11 + t + 1)._frobenius_coefficient_bound_charpoly()
            3
            sage: HyperellipticCurve(t^13 + t + 1)._frobenius_coefficient_bound_charpoly()
            4
        """
        assert self.base_ring().is_finite()
        p = self.base_ring().characteristic()
        q = self.base_ring().order()
        sqrtq = RR(q).sqrt()
        g = self.genus()

        # note: this bound is from Kedlaya's paper, but he tells me it's not
        # the best possible
        M = 2 * binomial(2*g, g) * sqrtq**g
        B = ZZ(M.ceil()).exact_log(p)
        if p**B < M:
            B += 1
        return B

    def _frobenius_coefficient_bound_traces(self, n=1):
        """
        Computes bound on number of ``p``-adic digits needed to recover
        the number of rational points on ``n`` extensions computing
        traces of the frobenius matrix powers, i.e. returns ``B`` so that
        knowledge of ``\tr(M^1)``, ..., ``\tr(M^n)`` modulo ``p^B`` determine
        ``N_1``, ..., ``N_n`` uniquely.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound_traces()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound_traces()
            1
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound_traces()
            1

            sage: R.<t> = PolynomialRing(GF(next_prime(10^9)))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound_traces()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound_traces()
            1
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound_traces()
            1
            sage: HyperellipticCurve(t^9 + t + 1)._frobenius_coefficient_bound_traces(n=3)
            2
            sage: HyperellipticCurve(t^11 + t + 1)._frobenius_coefficient_bound_traces(n=3)
            2
            sage: HyperellipticCurve(t^13 + t + 1)._frobenius_coefficient_bound_traces(n=5)
            3
        """
        g = self.genus()
        p = self.base_ring().characteristic()
        return (ZZ(4*g).exact_log(p) + n//2).floor() + 1

    def frobenius_matrix_hypellfrob(self, N=None):
        """
        Compute ``p``-adic frobenius matrix to precision ``p^N``.
        If ``N`` not supplied, a default value is selected, which is the
        minimum needed to recover the charpoly unambiguously.

        .. note::

            Currently only implemented using hypellfrob, which means only works
            over the prime field ``GF(p)``, and must have ``p > (2g+1)(2N-1)``.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_matrix_hypellfrob()
            [1258 + O(37^2)  925 + O(37^2)  132 + O(37^2)  587 + O(37^2)]
            [1147 + O(37^2)  814 + O(37^2)  241 + O(37^2) 1011 + O(37^2)]
            [1258 + O(37^2) 1184 + O(37^2) 1105 + O(37^2)  482 + O(37^2)]
            [1073 + O(37^2)  999 + O(37^2)  772 + O(37^2)  929 + O(37^2)]

        The `hypellfrob` program doesn't support non-prime fields::

            sage: K.<z> = GF(37**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + z*t^3 + 1)
            sage: H.frobenius_matrix_hypellfrob()
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of Frobenius matrix only implemented for hyperelliptic curves defined over prime fields.

        nor too small characteristic::

            sage: K = GF(7)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^3 + 1)
            sage: H.frobenius_matrix_hypellfrob()
            Traceback (most recent call last):
            ...
            ValueError: In the current implementation, p must be greater than (2g+1)(2N-1) = 81
        """
        p = self.base_ring().characteristic()
        e = self.base_ring().degree()
        if e != 1:
            raise NotImplementedError("Computation of Frobenius matrix only implemented for hyperelliptic curves defined over prime fields.")

        f, h = self.hyperelliptic_polynomials()
        if h != 0:
            # need y^2 = f(x)
            raise NotImplementedError("only implemented for curves y^2 = f(x)")

        sign = 1
        if not f.is_monic():
            # at this time we need a monic f
            c = f.leading_coefficient()
            f = f / c
            if c.is_square():
                # solutions of $y^2 = c * f(x)$ correspond naturally to
                # solutions of $(sqrt(c) y)^2 = f(x)$
                pass
            else:
                # we'll count points on a twist and then correct by untwisting...
                sign = -1
        assert f.is_monic()

        # By default, use precision enough to be able to compute the
        # frobenius minimal polynomial
        if N is None:
            N = self._frobenius_coefficient_bound_charpoly()

        matrix_of_frobenius = hypellfrob(p, N, f)
        matrix_of_frobenius = sign * matrix_of_frobenius
        return matrix_of_frobenius

    def frobenius_matrix(self, N=None, algorithm='hypellfrob'):
        """
        Compute ``p``-adic frobenius matrix to precision ``p^N``.
        If ``N`` not supplied, a default value is selected, which is the
        minimum needed to recover the charpoly unambiguously.

        .. note::

            Currently only implemented using hypellfrob, which means only works
            over the prime field ``GF(p)``, and must have ``p > (2g+1)(2N-1)``.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_matrix()
            [1258 + O(37^2)  925 + O(37^2)  132 + O(37^2)  587 + O(37^2)]
            [1147 + O(37^2)  814 + O(37^2)  241 + O(37^2) 1011 + O(37^2)]
            [1258 + O(37^2) 1184 + O(37^2) 1105 + O(37^2)  482 + O(37^2)]
            [1073 + O(37^2)  999 + O(37^2)  772 + O(37^2)  929 + O(37^2)]

        The `hypellfrob` program doesn't support non-prime fields::

            sage: K.<z> = GF(37**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + z*t^3 + 1)
            sage: H.frobenius_matrix(algorithm='hypellfrob')
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of Frobenius matrix only implemented for hyperelliptic curves defined over prime fields.

        nor too small characteristic::

            sage: K = GF(7)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^3 + 1)
            sage: H.frobenius_matrix(algorithm='hypellfrob')
            Traceback (most recent call last):
            ...
            ValueError: In the current implementation, p must be greater than (2g+1)(2N-1) = 81
        """
        if algorithm != 'hypellfrob':
            raise ValueError("Unknown algorithm")

        # By default, use precision enough to be able to compute the
        # frobenius minimal polynomial
        if N is None:
            N = self._frobenius_coefficient_bound_charpoly()

        return self.frobenius_matrix_hypellfrob(N=N)

    def frobenius_polynomial_cardinalities(self, a=None):
        """
        Compute the charpoly of frobenius, as an element of ``\ZZ[x]``,
        by computing the number of points on the curve over ``g`` extensions
        of the base field where ``g`` is the genus of the curve.

        .. WARNING::

            This is highly inefficient when the base field or the genus of the
            curve are large.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_polynomial_cardinalities()
            x^4 + x^3 - 52*x^2 + 37*x + 1369

        A quadratic twist::

            sage: H = HyperellipticCurve(2*t^5 + 2*t + 4)
            sage: H.frobenius_polynomial_cardinalities()
            x^4 - x^3 - 52*x^2 - 37*x + 1369

        Over a non-prime field, this currently does not work for ``g \geq 2``::

            sage: K.<z> = GF(17**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + z*t + z^2)
            sage: H.frobenius_polynomial_cardinalities()
            Traceback (most recent call last):
            ...
            NotImplementedError: Naive method only implemented over base field or extensions of prime fields in odd characteristic.

        This method may actually be useful when `hypellfrob` does not work::

            sage: K = GF(7)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^3 + 1)
            sage: H.frobenius_polynomial_matrix(algorithm='hypellfrob')
            Traceback (most recent call last):
            ...
            ValueError: In the current implementation, p must be greater than (2g+1)(2N-1) = 81
            sage: H.frobenius_polynomial_cardinalities()
            x^8 - 5*x^7 + 7*x^6 + 36*x^5 - 180*x^4 + 252*x^3 + 343*x^2 - 1715*x + 2401
        """
        g = self.genus()
        q = self.base_ring().cardinality()

        if a is None:
            # this may actually call frobenius_polynomial()
            a = self.count_points(g)
            # maybe calling count_points_exhaustive() would make more sense
            # but the method is currently only called with a precomputed list
            # of number of points so it does not really matter

        # computation of the reciprocal polynomial
        s = [ai - q**(i+1) - 1 for i, ai in enumerate(a)]
        coeffs = [1]
        for i in xrange(1,g+1):
            c = 0
            for j in xrange(i):
                c += s[i-1-j]*coeffs[j]
            coeffs.append(c/i)
        coeffs = coeffs + [coeffs[g-i] * q**(i) for i in xrange(1,g+1)]

        return ZZ['x'](coeffs).reverse()

    def frobenius_polynomial_matrix(self, M=None, algorithm='hypellfrob'):
        """
        Compute the charpoly of frobenius, as an element of ``\ZZ[x]``,
        by computing the charpoly of the frobenius matrix.

        This is currently only supported when the base field is prime
        and large enough using the `hypellfrob` library.

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_polynomial_matrix()
            x^4 + x^3 - 52*x^2 + 37*x + 1369

        A quadratic twist::

            sage: H = HyperellipticCurve(2*t^5 + 2*t + 4)
            sage: H.frobenius_polynomial_matrix()
            x^4 - x^3 - 52*x^2 - 37*x + 1369

        Curves defined over larger prime fields::

            sage: K = GF(49999)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^5 + 1)
            sage: H.frobenius_polynomial_matrix()
            x^8 + 281*x^7 + 55939*x^6 + 14144175*x^5 + 3156455369*x^4 + 707194605825*x^3 + 139841906155939*x^2 + 35122892542149719*x + 6249500014999800001
            sage: H = HyperellipticCurve(t^15 + t^5 + 1)
            sage: H.frobenius_polynomial_matrix() # long time
            x^14 - 76*x^13 + 220846*x^12 - 12984372*x^11 + 24374326657*x^10 - 1203243210304*x^9 + 1770558798515792*x^8 - 74401511415210496*x^7 + 88526169366991084208*x^6 - 3007987702642212810304*x^5 + 3046608028331197124223343*x^4 - 81145833008762983138584372*x^3 + 69007473838551978905211279154*x^2 - 1187357507124810002849977200076*x + 781140631562281254374947500349999

        This `hypellfrob` program doesn't support non-prime fields::

            sage: K.<z> = GF(37**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + z*t^3 + 1)
            sage: H.frobenius_polynomial_matrix(algorithm='hypellfrob')
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of Frobenius matrix only implemented for hyperelliptic curves defined over prime fields.
        """
        K = self.base_ring()
        p = K.characteristic()
        q = K.cardinality()
        g = self.genus()
        N = self._frobenius_coefficient_bound_charpoly()
        # compute chapoly over ZZ and then reduce back
        # (because charpoly of p-adic matrices sometimes loses precision)
        M = self.frobenius_matrix(N=N, algorithm=algorithm).change_ring(ZZ)

        # get a_g, ..., a_0 in ZZ (i.e. with correct signs)
        f = M.charpoly().list()[g:2*g+1]
        ppow = p**N
        f = [x % ppow for x in f]
        f = [x if 2*x < ppow else x - ppow for x in f]

        # get a_{2g}, ..., a_{g+1}
        f = [f[g-i] * q**(g-i) for i in range(g)] + f

        return ZZ['x'](f)

    @cached_method
    def frobenius_polynomial(self):
        """
        Compute the charpoly of frobenius, as an element of ZZ[x].

        EXAMPLES::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_polynomial()
            x^4 + x^3 - 52*x^2 + 37*x + 1369

        A quadratic twist::

            sage: H = HyperellipticCurve(2*t^5 + 2*t + 4)
            sage: H.frobenius_polynomial()
            x^4 - x^3 - 52*x^2 - 37*x + 1369

        Slightly larger example::

            sage: K = GF(2003)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^7 + 487*t^5 + 9*t + 1)
            sage: H.frobenius_polynomial()
            x^6 - 14*x^5 + 1512*x^4 - 66290*x^3 + 3028536*x^2 - 56168126*x + 8036054027

        Curves defined over a non-prime field are only supported for ``g < 2``
        in odd characteristic and a naive algorithm is used
        (and we should actually use fast point counting on elliptic curves)::

            sage: K.<z> = GF(23**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^3 + z*t + 4)
            sage: H.frobenius_polynomial() # long time
            x^2 - 15*x + 12167

        To support higher genera, we need a way of easily chacking for
        squareness in quotient rings of polynomial rings over finite fields
        of odd characteristic::

            sage: K.<z> = GF(23**3)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + z*t + z**3)
            sage: H.frobenius_polynomial()
            Traceback (most recent call last):
            ...
            NotImplementedError: Naive method only implemented over base field or extensions of prime fields in odd characteristic.

        Over prime fields of odd characteristic, when ``h`` is non-zero,
        this naive algorithm is currently used as well, whereas we should
        rather use another defining equation::

            sage: K = GF(101)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + 27*t + 3, t)
            sage: H.frobenius_polynomial()
            x^4 + 2*x^3 - 58*x^2 + 202*x + 10201

        In even characteristic, the naive algorithm could cover all cases
        because we can easily check for squareness in quotient rings of
        polynomial rings over finite fields but these rings unfortunately
        do not support iteration::

            sage: K.<z> = GF(2**5)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + z*t + z**3, t)
            sage: H.frobenius_polynomial()
            Traceback (most recent call last):
            ...
            NotImplementedError: object does not support iteration
        """
        K = self.base_ring()
        e = K.degree()
        q = K.cardinality()

        g = self.genus()
        f, h = self.hyperelliptic_polynomials()

        if e == 1 and q >= (2*g+1)*(2*self._frobenius_coefficient_bound_charpoly()-1) and h == 0:
            return self.frobenius_polynomial_matrix()
        else:
            return self.frobenius_polynomial_cardinalities()

    def _points_fast_sqrt(self):
        """
        List points by enumerating over x and solving the resulting
        quadratic for y.

        EXAMPLES::

            sage: K.<a> = GF(9, 'a')
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^7 - 1, x^2 + a)
            sage: C._points_fast_sqrt()
            [(0 : 1 : 0), (a : 2*a + 1 : 1), (2 : a + 1 : 1), (2*a + 2 : 2*a : 1), (2*a + 2 : 1 : 1), (1 : 2*a + 2 : 1), (1 : 0 : 1)]
            sage: K.<a> = GF(49, 'a')
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 - x^2 - 1, x^2 + a)
            sage: len(C._points_fast_sqrt())
            31

        TESTS::

            sage: x = polygen(GF(16, 'a'))
            sage: C = HyperellipticCurve(x^5 - x + 1, x^2 + x + 1)
            sage: set(C._points_fast_sqrt()) == set(C._points_cache_sqrt())
            True
            sage: x = polygen(GF(19))
            sage: C = HyperellipticCurve(x^5 + 5*x^2 + 1, x + 1)
            sage: set(C._points_fast_sqrt()) == set(C._points_cache_sqrt())
            True
            sage: x = polygen(GF(13))
            sage: C = HyperellipticCurve(x^3 + x^2 - 1)
            sage: C._points_fast_sqrt()
            [(0 : 1 : 0), (0 : 5 : 1), (0 : 8 : 1), (1 : 1 : 1), (1 : 12 : 1), (3 : 3 : 1), (3 : 10 : 1), (4 : 1 : 1), (4 : 12 : 1), (6 : 2 : 1), (6 : 11 : 1), (7 : 1 : 1), (7 : 12 : 1), (8 : 4 : 1), (8 : 9 : 1), (9 : 4 : 1), (9 : 9 : 1), (12 : 5 : 1), (12 : 8 : 1)]
            sage: set(C._points_fast_sqrt()) == set(C._points_cache_sqrt())
            True
        """
        # For givaro finite fields, taking square roots is very fast
        # so no need to cache as in prime case
        K = self.base_ring()
        f, h = self.hyperelliptic_polynomials()
        one = K(1)

        # start with the points at infinity
        P = self.defining_polynomial()
        if not P(K(0), K(1), K(0)):
            # (0:1:0) is a point on the curve
            points = [self.point([K(0), K(1), K(0)], check=True)]
        else:
            points=[]
        if P.degree() > 2:
            # P(1, y, 0) = r*y + s
            s = P(K(1), K(0), K(0))
            r = P(K(1), K(1), K(0)) - s
            if r: # r not zero
                points.append(self.point([K(1), -s/r, K(0)], check=True))
            # the case r = 0 need not be considered
        elif K.characteristic() == 2: # deg(P) = 2 and char(K) = 2
            # quadratic equation doesn't work in characteristic 2 so use brute
            # force
            points += [self.point([K(1), y, K(0)], check=True) for y in K \
                       if not P(K(1), y, K(0))]
        else: # deg(P) = 2 and char(K) not 2
            # P(1, y, 0) = y^2 + r*y + s
            s = -f[2]
            r = h[1]
            d = r**2/4 - s
            if not d: # d = 0
                points.append(self.point([K(1), -r/2, K(0)], check=True))
            elif d.is_square():
                sqrtd = d.sqrt()
                points.append(self.point([K(1), -r/2+sqrtd, K(0)], check=True))
                points.append(self.point([K(1), -r/2-sqrtd, K(0)], check=True))

        if K.characteristic() == 2:
            # quadratic equation doesn't work in characteristic 2
            if h.is_zero():
                for x in K:
                    points.append(self.point([x, f(x).sqrt(), one], check=True))
            else:
                R = f.base_ring()
                a_sqrts = { } # Artin-Schreier 2-roots
                for x in K:
                    a_sqrts[x**2 + x] = x  # char 2 => x^2 - x == x^2 + x
                for x in K:
                    b = h(x)
                    c = f(x)
                    if b:
                        try:
                            r = a_sqrts[c / b**2]
                            points.append(self.point([x, r*b, one], check=True))
                            points.append(self.point([x, r*b+b, one], check=True))
                        except KeyError:
                            # y^2 + by + c irreducible, so yields no points
                            pass
                    else:  # b == 0
                        points.append(self.point([x, c.sqrt(), one], check=True))
        elif h.is_zero():
            # special case to save work if we are of the form y^2 = f(x)
            for x in K:
                y2 = f(x)
                if not y2: # y = 0
                    points.append(self.point([x, y2, one], check=True))
                elif y2.is_square():
                    y = y2.sqrt()
                    points.append(self.point([x, y, one], check=True))
                    points.append(self.point([x, -y, one], check=True))
        else:
            b = -h/2
            D = b*b + f
            for x in K:
                Dval = D(x)
                if not Dval:  # D(x) = 0
                    points.append(self.point([x, b(x), one], check=True))
                elif Dval.is_square():
                    sqrtD = Dval.sqrt()
                    v = b(x)
                    points.append(self.point([x, v+sqrtD, one], check=True))
                    points.append(self.point([x, v-sqrtD, one], check=True))
        return points

    def _points_cache_sqrt(self, brute_force=False):
        """
        List points by enumerating over x and solving the resulting
        quadratic for y.

        Caches all square roots ahead of time by squaring every element of
        the field. Elements must have an __index__ method.

        EXAMPLES::

            sage: x = polygen(GF(7))
            sage: C = HyperellipticCurve(x^3 + x^2 - 1)
            sage: C._points_cache_sqrt()
            [(0 : 1 : 0), (1 : 6 : 1), (1 : 1 : 1), (2 : 5 : 1), (2 : 2 : 1), (3 : 0 : 1), (4 : 4 : 1), (4 : 3 : 1), (5 : 4 : 1), (5 : 3 : 1)]
            sage: set(C._points_cache_sqrt()) == set(C._points_cache_sqrt(brute_force=True))
            True
        """
        K = self.base_ring()
        if K.characteristic() != 2:
            # cache the squares (faster than O(p) sqrts)
            square_roots = [None] * len(K)
            for x in K:
                square_roots[x*x] = x
        f, h = self.hyperelliptic_polynomials()
        one = K(1)

        # start with the points at infinity
        P = self.defining_polynomial()
        if not P(K(0), K(1), K(0)):
            # (0:1:0) is a point on the curve
            points = [self.point([K(0), K(1), K(0)], check=True)]
        else:
            points=[]
        if P.degree() > 2:
            # P(1, y, 0) = r*y + s
            s = P(K(1), K(0), K(0))
            r = P(K(1), K(1), K(0)) - s
            if r: # r not zero
                points.append(self.point([K(1), -s/r, K(0)], check=True))
            # the case r = 0 need not be considered
        elif K.characteristic() == 2: # deg(P) = 2 and char(K) = 2
            # quadratic equation doesn't work in characteristic 2 so use brute
            # force
            points += [self.point([K(1), y, K(0)], check=True) for y in K \
                       if not P(K(1), y, K(0))]
        else: # deg(P) = 2 and char(K) not 2
            # P(1, y, 0) = y^2 + r*y + s
            s = -f[2]
            r = h[1]
            d = r**2/4 - s
            sqrtd = square_roots[d]
            if not d: # d = 0
                points.append(self.point([K(1), -r/2, K(0)], check=True))
            elif sqrtd is not None:
                points.append(self.point([K(1), -r/2+sqrtd, K(0)], check=True))
                points.append(self.point([K(1), -r/2-sqrtd, K(0)], check=True))

        if K.characteristic() == 2 or brute_force:
            # quadratic equation doesn't work in characteristic 2
            # but there are only 4 affine points, so just test them
            f = self.defining_polynomial()
            points += [self.point([x, y, one], check=True) for x in K for y in K if not f(x, y, one)]
        elif h.is_zero():
            # special case to save work if we are of the form y^2 = f(x)
            for x in K:
                y2 = f(x)
                y = square_roots[y2]
                if not y2: # y = 0
                    points.append(self.point([x, y2, one], check=True))
                elif y is not None:
                    points.append(self.point([x, y, one], check=True))
                    points.append(self.point([x, -y, one], check=True))
        else:
            b = -h/2
            D = b*b + f  # this is really disc/4
            for x in K:
                Dval = D(x)
                sqrtD = square_roots[Dval]
                if not Dval:  # D(x) = 0
                    points.append(self.point([x, b(x), one], check=True))
                elif sqrtD is not None:
                    v = b(x)
                    points.append(self.point([x, v+sqrtD, one], check=True))
                    points.append(self.point([x, v-sqrtD, one], check=True))
        return points

    def points(self):
        r"""
        All the points on this hyperelliptic curve.

        EXAMPLES::

            sage: x = polygen(GF(7))
            sage: C = HyperellipticCurve(x^7 - x^2 - 1)
            sage: C.points()
            [(0 : 1 : 0), (2 : 5 : 1), (2 : 2 : 1), (3 : 0 : 1), (4 : 6 : 1), (4 : 1 : 1), (5 : 0 : 1), (6 : 5 : 1), (6 : 2 : 1)]

        ::

            sage: x = polygen(GF(121, 'a'))
            sage: C = HyperellipticCurve(x^5 + x - 1, x^2 + 2)
            sage: len(C.points())
            122

        Conics are allowed (the issue reported at #11800 has been resolved)::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurve(3*x^2 + 5*x + 1)
            sage: H.points()
            [(0 : 6 : 1), (0 : 1 : 1), (1 : 4 : 1), (1 : 3 : 1), (2 : 4 : 1), (2 : 3 : 1), (3 : 6 : 1), (3 : 1 : 1)]

        The method currently lists points on the plane projective model, that
        is the closure in $\mathbb{P}^2$ of the curve defined by $y^2+hy=f$.
        This means that one point $(0:1:0)$ at infinity is returned if the
        degree of the curve is at least 4 and $\deg(f)>\deg(h)+1$. This point
        is a singular point of the plane model. Later implementations may
        consider a smooth model instead since that would be a more relevant
        object. Then, for a curve whose only singularity is at $(0:1:0)$, the
        point at infinity would be replaced by a number of rational points of
        the smooth model. We illustrate this with an example of a genus 2
        hyperelliptic curve::

            sage: R.<x>=GF(11)[]
            sage: H = HyperellipticCurve(x*(x+1)*(x+2)*(x+3)*(x+4)*(x+5))
            sage: H.points()
            [(0 : 1 : 0), (0 : 0 : 1), (1 : 7 : 1), (1 : 4 : 1), (5 : 7 : 1), (5 : 4 : 1), (6 : 0 : 1), (7 : 0 : 1), (8 : 0 : 1), (9 : 0 : 1), (10 : 0 : 1)]

        The plane model of the genus 2 hyperelliptic curve in the above example
        is the curve in $\mathbb{P}^2$ defined by $y^2z^4=g(x,z)$ where
        $g(x,z)=x(x+z)(x+2z)(x+3z)(x+4z)(x+5z).$ This model has one point at
        infinity $(0:1:0)$ which is also the only singular point of the plane
        model. In contrast, the hyperelliptic curve is smooth and imbeds via
        the equation $y^2=g(x,z)$ into weighted projected space
        $\mathbb{P}(1,3,1)$. The latter model has two points at infinity:
        $(1:1:0)$ and $(1:-1:0)$.
        """
        from sage.rings.finite_rings.constructor import zech_log_bound
        try:
            return self.__points
        except AttributeError: pass

        if self.base_ring().is_prime_field():
            self.__points = self._points_cache_sqrt()
        else:
            if self.base_ring().order() < zech_log_bound:
                self.__points = self._points_fast_sqrt() # this is fast using Zech logarithms
            else:
                self.__points = self._points_cache_sqrt()

        return self.__points

    def count_points_matrix_traces(self, n=1, M=None, N=None):
        r"""
        Count the number of points on the curve over the first ``n`` extensions
        of the base field by computing traces of powers of the frobenius
        matrix.
        This requires less ``p``-adic precision than computing the charpoly
        of the matrix when ``n < g`` where ``g`` is the genus of the curve.

        EXAMPLES::

            sage: K = GF(49999)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^19 + t + 1)
            sage: H.count_points_matrix_traces(3)
            [49491, 2500024375, 124992509154249]
        """
        if N is None:
            N = self._frobenius_coefficient_bound_traces(n=n)

        if M is None:
            M = self.frobenius_matrix(N=N)

        K = self.base_ring()
        p = K.characteristic()
        q = K.cardinality()
        ppow = p**N

        t = []
        Mpow = 1
        for i in xrange(n):
            Mpow *= M
            t.append(Mpow.trace())

        t = [x.lift() for x in t]
        t = [x if 2*x < ppow else x - ppow for x in t]

        return [q**(i+1) + 1 - t[i] for i in xrange(n)]

    def count_points_frobenius_polynomial(self, n=1, f=None):
        r"""
        Count the number of points on the curve over the first n extensions
        of the base field by computing the frobenius polynomial.

        EXAMPLES::

            sage: K = GF(49999)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^19 + t + 1)

        The following computation takes a long time as the complete
        characteristic polynomial of the frobenius is computed::

            sage: H.count_points_frobenius_polynomial(3) # long time, 20s on a Corei7
            [49491, 2500024375, 124992509154249]

        As the polynomial is cached, further computations of number of points
        are really fast::

            sage: H.count_points_frobenius_polynomial(19)
            [49491,
            2500024375,
            124992509154249,
            6249500007135192947,
            312468751250758776051811,
            15623125093747382662737313867,
            781140631562281338861289572576257,
            39056250437482500417107992413002794587,
            1952773465623687539373429411200893147181079,
            97636720507718753281169963459063147221761552935,
            4881738388665429945305281187129778704058864736771824,
            244082037694882831835318764490138139735446240036293092851,
            12203857802706446708934102903106811520015567632046432103159713,
            610180686277519628999996211052002771035439565767719719151141201339,
            30508424133189703930370810556389262704405225546438978173388673620145499,
            1525390698235352006814610157008906752699329454643826047826098161898351623931,
            76268009521069364988723693240288328729528917832735078791261015331201838856825193,
            3813324208043947180071195938321176148147244128062172555558715783649006587868272993991,
            190662397077989315056379725720120486231213267083935859751911720230901597698389839098903847]
        """
        if f is None:
            f = self.frobenius_polynomial()

        g = self.genus()
        q = self.base_ring().cardinality()
        S = PowerSeriesRing(QQ, default_prec=n+1, names='t')
        frev = f.reverse()
        # the coefficients() method of power series only returns non-zero coefficients
        # so let's use the list() method
        # but this does not work for zero which gives the empty list
        flog = S(frev).log()
        return [q**(i+1) + 1 + ZZ((i+1)*flog[i+1]) for i in xrange(n)]

    def count_points_exhaustive(self, n=1, naive=False):
        r"""
        Count the number of points on the curve over the first ``n`` extensions
        of the base field by exhaustive search if ``n`` if smaller than ``g``,
        the genus of the curve, and by computing the frobenius polynomial
        after performing exhaustive search on the first ``g`` extensions if
        ``n > g`` (unless `naive == True`).

        EXAMPLES::

            sage: K = GF(5)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^3 + 1)
            sage: H.count_points_exhaustive(n=5)
            [9, 27, 108, 675, 3069]

        When ``n > g``, the frobenius polynomial is computed from the numbers
        of points of the curve over the first ``g`` extension, so that computing
        the number of points on extensions of degree ``n > g`` is not much more
        expensive than for ``n == g``::

            sage: H.count_points_exhaustive(n=15)
            [9,
            27,
            108,
            675,
            3069,
            16302,
            78633,
            389475,
            1954044,
            9768627,
            48814533,
            244072650,
            1220693769,
            6103414827,
            30517927308]

        This behavior can be disabled by passing `naive=True`::

           sage: H.count_points_exhaustive(n=6, naive=True)
           [9, 27, 108, 675, 3069, 16302]
        """
        g = self.genus()
        a = []
        for i in xrange(1, min(n, g)+1):
            a.append(self.cardinality_exhaustive(extension_degree=i))

        if n <= g:
            return a

        if naive:
            for i in xrange(g+1,n+1):
                a.append(self.cardinality_exhaustive(extension_degree=i))

        # let's not be too naive and compute the frobenius polynomial
        f = self.frobenius_polynomial_cardinalities(a=a)
        return self.count_points_frobenius_polynomial(n=n, f=f)

    def count_points_hypellfrob(self, n=1, N=None, algorithm=None):
        r"""
        Count the number of points on the curve over the first ``n``` extensions
        of the base field using the `hypellfrob` program.

        This only supports prime fields of large enough characteristic.

        EXAMPLES::

            sage: K = GF(49999)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^21 + 3*t^5 + 5)
            sage: H.count_points_hypellfrob()
            [49804]
            sage: H.count_points_hypellfrob(2)
            [49804, 2499799038]

            sage: K = GF(2**7-1)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^11 + 3*t^5 + 5)
            sage: H.count_points_hypellfrob()
            [127]
            sage: H.count_points_hypellfrob(n=5)
            [127, 16335, 2045701, 260134299, 33038098487]

            sage: K = GF(2**7-1)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^13 + 3*t^5 + 5)
            sage: H.count_points(n=6)
            [112, 16360, 2045356, 260199160, 33038302802, 4195868633548]

        The base field should be prime::

            sage: K.<z> = GF(19**10)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + (z+1)*t^5 + 1)
            sage: H.count_points_hypellfrob()
            Traceback (most recent call last):
            ...
            ValueError: hypellfrob does not support non-prime fields

        and the characteristic should be large enough::

            sage: K = GF(7)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + t^3 + 1)
            sage: H.count_points_hypellfrob()
            Traceback (most recent call last):
            ...
            ValueError: p=7 should be greater than (2*g+1)(2*N-1)=27
        """
        K = self.base_ring()
        e = K.degree()

        if e != 1:
            raise ValueError("hypellfrob does not support non-prime fields")

        # K is a prime field
        p = K.cardinality()
        g = self.genus()

        if algorithm is None:
            if n < g:
                algorithm = 'traces'
            else:
                algorithm = 'charpoly'

        if N is None:
            if algorithm == 'traces':
                N = self._frobenius_coefficient_bound_traces(n)
            elif algorithm == 'charpoly':
                N = self._frobenius_coefficient_bound_charpoly()
            else:
                raise ValueError("Unknown algorithm")

        if p <= (2*g+1)*(2*N-1):
            raise ValueError("p=%d should be greater than (2*g+1)(2*N-1)=%d"%(p,(2*g+1)*(2*N-1)))

        if algorithm == 'traces':
            M = self.frobenius_matrix(N=N, algorithm='hypellfrob')
            return self.count_points_matrix_traces(n=n,M=M,N=N)
        elif algorithm == 'charpoly':
            f = self.frobenius_polynomial_matrix(algorithm='hypellfrob')
            return self.count_points_frobenius_polynomial(n=n,f=f)
        else:
            raise ValueError("Unknown algorithm")

    def count_points(self, n=1):
        r"""
        Count points over finite fields.

        INPUT:

        - ``n`` -- integer.

        OUTPUT:

        An integer. The number of points over `\GF{q}, \ldots,
        \GF{q^n}` on a hyperelliptic curve over a finite field `\GF{q}`.

        .. WARNING::

           This is currently using exhaustive search for hyperelliptic curves
           over non-prime fields, which can be awfully slow.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(GF(3))
            sage: C = HyperellipticCurve(x^3+x^2+1)
            sage: C.count_points(4)
            [6, 12, 18, 96]
            sage: C.base_extend(GF(9,'a')).count_points(2)
            [12, 96]

            sage: K = GF(2**31-1)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^5 + 3*t + 5)
            sage: H.count_points() # long time, 2.4 sec on a Corei7
            [2147464821]
            sage: H.count_points(n=2) # long time, 30s on a Corei7
            [2147464821, 4611686018988310237]

            sage: K = GF(2**7-1)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^13 + 3*t^5 + 5)
            sage: H.count_points(n=6)
            [112, 16360, 2045356, 260199160, 33038302802, 4195868633548]
        """
        K = self.base_ring()
        q = K.cardinality()
        e = K.degree()
        g = self.genus()
        f, h = self.hyperelliptic_polynomials()

        if e == 1 and h == 0:
            N1 = self._frobenius_coefficient_bound_traces(n)
            N2 = self._frobenius_coefficient_bound_charpoly()
            if n < g and q > (2*g+1)*(2*N1-1):
                return self.count_points_hypellfrob(n, N=N1, algorithm='traces')
            elif q > (2*g+1)*(2*N2-1):
                return self.count_points_hypellfrob(n, N=N2, algorithm='charpoly')

        # No smart method available
        return self.count_points_exhaustive(n)

    def cardinality_exhaustive(self, extension_degree=1, algorithm=None):
        """
        Count points on a single extension of the base field
        by enumerating over x and solving the resulting quadratic
        equation for y.

        EXAMPLES::

            sage: K.<a> = GF(9, 'a')
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^7 - 1, x^2 + a)
            sage: C.cardinality_exhaustive()
            7

            sage: K = GF(next_prime(1<<10))
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^7 + 3*t^5 + 5)
            sage: H.cardinality_exhaustive()
            1025
        """
        K = self.base_ring()
        g = self.genus()
        n = extension_degree

        if g == 0:
            # here is the projective line
            return K.cardinality()**n + 1

        f, h = self.hyperelliptic_polynomials()
        a = 0

        # begin with points at infinity (on the smooth model)
        if g == 1:
            # elliptic curves always have one smooth point at infinity
            a += 1
        else:
            # g > 1
            # solve y^2 + y*h[g+1] == f[2*g+2], i.e., y^2 + r*y - s == 0
            s = f[2*g+2]
            r = h[g+1]
            if K.characteristic() == 2:
                if r == 0:
                    # one unique solution
                    a += 1
                elif n % 2 == 0 or (s/r**2).trace() == 0:
                    # artin-schreier equation t^2 + t = s/r^2
                    # it always has a solution in extensions of even degree
                    a += 2
            else:
                # compute the usual discriminant
                d = r**2 + 4*s
                if d == 0:
                    a += 1
                elif n % 2 == 0 or d.is_square():
                    a += 2

        # affine points
        if K.characteristic() == 2:
            if n == 1:
                # the base field
                L = K
            elif K.degree() == 1:
                # extension of the prime field
                L = GF(2**n, name='t')
            else:
                # general extension
                # this might be costly, but everything here is anyway
                # WARNING:
                # 1. beware that if K is actually defined as a quotient of
                # a polynomial ring, then polynomial_ring() will return
                # will return the polynomial ring it is a quotient of,
                # and not a polynomial ring on top of the quotient ring...
                # IMHO this is a bug
                # 2. if we rather used the PolynomialRing constructor, we would
                # get the correct construction, but we then need to make sure
                # we call trace enough below, e.g. with a while loop where we
                # repeatedly apply trace() till the result lives in the prime
                # field
                # 3. for our application we don't really care as hyperelliptic
                # curves defined over finite field represented as quotient rings
                # won't end up in the HyperellipticCurve_finite_field class
                # 4. and this does not work as a quotient ring of a polynomial
                # ring over a finite field is not iterable, so the following
                # for loop will fail at the moment
                L = K.extension(K.polynomial_ring().irreducible_element(n))
            for x in L:
                s = f(x)
                r = h(x)
                if r == 0:
                    a += 1
                elif (s/r**2).trace().trace() == 0:
                    # currently extension of finite fields may be represented
                    # as quotient rings of polynomial rings,
                    # so for these we need:
                    # - one trace to go to the base field,
                    # - another trace to go to F_2.
                    # if L is actually represented directly as a finite field,
                    # the second trace does not hurt as tr(F_p) = id
                    a += 2
        else:
            if n == 1:
                # L is the base field
                L = K
            elif K.degree() == 1:
                # L is an extension of the prime field
                L = GF(K.cardinality()**n, name='t')
            else:
                # there is no easy way to construct extensions of finite fields
                # of odd characteristic in such a way that we can check for
                # squares easily
                raise NotImplementedError("Naive method only implemented over base field or extensions of prime fields in odd characteristic.")
            for x in L:
                s = f(x)
                r = h(x)
                d = r**2 + 4*s
                if d.is_zero():
                    a += 1
                elif d.is_square():
                    a += 2
        return a

    def cardinality_hypellfrob(self, extension_degree=1, algorithm=None):
        r"""
        Count points on a single extension of the base field
        using the hypellfrob prgoram.

        EXAMPLES::

            sage: K = GF(next_prime(1<<10))
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^7 + 3*t^5 + 5)
            sage: H.cardinality_hypellfrob()
            1025

            sage: K = GF(49999)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^7 + 3*t^5 + 5)
            sage: H.cardinality_hypellfrob()
            50162
            sage: H.cardinality_hypellfrob(3)
            124992471088310
        """
        # the following actually computes the cardinality for several extensions
        # but the overhead is negligible
        return self.count_points_hypellfrob(n=extension_degree, algorithm=algorithm)[-1]

    @cached_method
    def cardinality(self, extension_degree=1):
        r"""
        Count points on a single extension of the base field.

        EXAMPLES::

            sage: K = GF(101)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + 3*t^5 + 5)
            sage: H.cardinality()
            106
            sage: H.cardinality(15)
            1160968955369992567076405831000
            sage: H.cardinality(100)
            270481382942152609326719471080753083367793838278100277689020104911710151430673927943945601434674459120495370826289654897190781715493352266982697064575800553229661690000887425442240414673923744999504000

            sage: K = GF(37)
            sage: R.<t> = PolynomialRing(K)
            sage: H = HyperellipticCurve(t^9 + 3*t^5 + 5)
            sage: H.cardinality()
            40
            sage: H.cardinality(2)
            1408
            sage: H.cardinality(3)
            50116
        """
        K = self.base_ring()
        q = K.cardinality()
        e = K.degree()
        g = self.genus()
        f, h = self.hyperelliptic_polynomials()
        n = extension_degree

        # We may:
        # - check for actual field of definition of the curve (up to isomorphism)
        if e == 1 and h == 0:
            N1 = self._frobenius_coefficient_bound_traces(n)
            N2 = self._frobenius_coefficient_bound_charpoly()
            if n < g and q > (2*g+1)*(2*N1-1):
                return self.cardinality_hypellfrob(n, algorithm='traces')
            elif q > (2*g+1)*(2*N2-1):
                return self.cardinality_hypellfrob(n, algorithm='charpoly')

        # No smart method available
        return self.cardinality_exhaustive(n)

    #This where Cartier Matrix is actually computed. This is either called by E.Cartier_matrix, E.a_number, or E.Hasse_Witt
    @cached_method
    def _Cartier_matrix_cached(self):
        r"""
        INPUT:

        - 'E' - Hyperelliptic Curve of the form `y^2 = f(x)` over a
          finite field, `\mathbb{F}_q`

        OUTPUT:

        - matrix(Fq,M)' The matrix $M = (c_(pi-j)), f(x)^((p-1)/2) = \sum c_i x^i$
        - 'Coeff' List of Coeffs of F, this is needed for Hasse-Witt function.
        - 'g' genus of the curve self, this is needed by a-number.
        - 'Fq' is the base field of self, and it is needed for Hasse-Witt
        - 'p' is the char(Fq), this is needed for Hasse-Witt.
        - 'E' The initial elliptic curve to check some caching conditions.

        EXAMPLES::

            sage: K.<x>=GF(9,'x')[]
            sage: C=HyperellipticCurve(x^7-1,0)
            sage: C._Cartier_matrix_cached()
            (
            [0 0 2]
            [0 0 0]
            [0 1 0], [2, 0, 0, 0, 0, 0, 0, 1, 0], 3, Finite Field in x of size 3^2, 3, Hyperelliptic Curve over Finite Field in x of size 3^2 defined by y^2 = x^7 + 2
            )
            sage: K.<x>=GF(49,'x')[]
            sage: C=HyperellipticCurve(x^5+1,0)
            sage: C._Cartier_matrix_cached()
            (
            [0 3]
            [0 0], [1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 1], 2, Finite Field in x of size 7^2, 7, Hyperelliptic Curve over Finite Field in x of size 7^2 defined by y^2 = x^5 + 1
            )

            sage: P.<x>=GF(9,'a')[]
            sage: C=HyperellipticCurve(x^29+1,0)
            sage: C._Cartier_matrix_cached()
            (
            [0 0 1 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 1 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 1 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 14, Finite Field in a of size 3^2, 3, Hyperelliptic Curve over Finite Field in a of size 3^2 defined by y^2 = x^29 + 1
            )

        TESTS::

            sage: K.<x>=GF(2,'x')[]
            sage: C=HyperellipticCurve(x^7-1,x)
            sage: C._Cartier_matrix_cached()
            Traceback (most recent call last):
            ...
            ValueError: p must be odd

            sage: K.<x>=GF(5,'x')[]
            sage: C=HyperellipticCurve(x^7-1,4)
            sage: C._Cartier_matrix_cached()
            Traceback (most recent call last):
            ...
            ValueError: E must be of the form y^2 = f(x)

            sage: K.<x>=GF(5,'x')[]
            sage: C=HyperellipticCurve(x^8-1,0)
            sage: C._Cartier_matrix_cached()
            Traceback (most recent call last):
            ...
            ValueError: In this implementation the degree of f must be odd

            sage: K.<x>=GF(5,'x')[]
            sage: C=HyperellipticCurve(x^5+1,0,check_squarefree=False)
            sage: C._Cartier_matrix_cached()
            Traceback (most recent call last):
            ...
            ValueError: curve is not smooth
        """
        #Compute the finite field and prime p.
        Fq=self.base_ring();
        p=Fq.characteristic()
        #checks

        if p == 2:
            raise ValueError("p must be odd");


        g = self.genus()

        #retrieve the function f(x) ,where y^2=f(x)
        f,h = self.hyperelliptic_polynomials()
        #This implementation only deals with h=0
        if h!=0:
            raise ValueError("E must be of the form y^2 = f(x)")

        d = f.degree()
        #this implementation is for odd degree only, even degree will be handled later.
        if d%2 == 0:
            raise ValueError("In this implementation the degree of f must be odd")
        #Compute resultant to make sure no repeated roots
        df=f.derivative()
        R=df.resultant(f)
        if R == 0:
            raise ValueError("curve is not smooth")

        #computing F, since the entries of the matrix are c_i where F= \sum c_i x^i

        F = f**((p-1)/2)

        #coefficients returns a_0, ... , a_n where f(x) = a_n x^n + ... + a_0

        Coeff = F.list()

        #inserting zeros when necessary-- that is, when deg(F) < p*g-1, (simplified if p <2g-1)
        #which is the highest powered coefficient needed for our matrix
        #So we don't have enough coefficients we add extra zeros to have the same poly,
        #but enough coeff.

        zeros = [0 for i in range(p*g-len(Coeff))]
        Coeff = Coeff + zeros

        # compute each row of matrix as list and then M=list of lists(rows)

        M=[];
        for j in range(1,g+1):
            H=[Coeff[i] for i in range((p*j-1), (p*j-g-1),-1)]
            M.append(H);
        return matrix(Fq,M), Coeff, g, Fq,p, self

    #This is what is called from command line
    def Cartier_matrix(self):
        r"""
        INPUT:

        - ``E`` : Hyperelliptic Curve of the form `y^2 = f(x)` over a finite field, `\mathbb{F}_q`

        OUTPUT:


        - ``M``: The matrix `M = (c_{pi-j})`, where `c_i` are the coeffients of  `f(x)^{(p-1)/2} = \sum c_i x^i`

        Reference-N. Yui. On the Jacobian varieties of hyperelliptic curves over fields of characteristic `p > 2`.

        EXAMPLES::

            sage: K.<x>=GF(9,'x')[]
            sage: C=HyperellipticCurve(x^7-1,0)
            sage: C.Cartier_matrix()
            [0 0 2]
            [0 0 0]
            [0 1 0]

            sage: K.<x>=GF(49,'x')[]
            sage: C=HyperellipticCurve(x^5+1,0)
            sage: C.Cartier_matrix()
            [0 3]
            [0 0]

            sage: P.<x>=GF(9,'a')[]
            sage: E=HyperellipticCurve(x^29+1,0)
            sage: E.Cartier_matrix()
            [0 0 1 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 1 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 1 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 1 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 1 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 1 0]

        TESTS::

            sage: K.<x>=GF(2,'x')[]
            sage: C=HyperellipticCurve(x^7-1,x)
            sage: C.Cartier_matrix()
            Traceback (most recent call last):
            ...
            ValueError: p must be odd

            sage: K.<x>=GF(5,'x')[]
            sage: C=HyperellipticCurve(x^7-1,4)
            sage: C.Cartier_matrix()
            Traceback (most recent call last):
            ...
            ValueError: E must be of the form y^2 = f(x)

            sage: K.<x>=GF(5,'x')[]
            sage: C=HyperellipticCurve(x^8-1,0)
            sage: C.Cartier_matrix()
            Traceback (most recent call last):
            ...
            ValueError: In this implementation the degree of f must be odd

            sage: K.<x>=GF(5,'x')[]
            sage: C=HyperellipticCurve(x^5+1,0,check_squarefree=False)
            sage: C.Cartier_matrix()
            Traceback (most recent call last):
            ...
            ValueError: curve is not smooth
        """
        #checking first that Cartier matrix is not already cached. Since
        #it can be called by either Hasse_Witt or a_number.
        #This way it does not matter which function is called first
        #in the code.
        # Trac Ticket #11115: Why shall we waste time by studying
        # the cache manually? We only need to check whether the cached
        # data belong to self.
        M, Coeffs,g, Fq, p, E= self._Cartier_matrix_cached()
        if E!=self:
            self._Cartier_matrix_cached.clear_cache()
            M, Coeffs,g, Fq, p, E= self._Cartier_matrix_cached()
        return M

    #This where Hasse_Witt is actually computed. This is either called by E.Hasse_Witt or p_rank
    @cached_method
    def _Hasse_Witt_cached(self):
        r"""
        INPUT:

        - ``E`` : Hyperelliptic Curve of the form `y^2 = f(x)` over a finite field, `\mathbb{F}_q`

        OUTPUT:

        - ``N`` : The matrix `N = M M^p \dots M^{p^{g-1}}` where `M = c_{pi-j}, f(x)^{(p-1)/2} = \sum c_i x^i`
        - ``E`` : The initial curve to check some caching conditions.

        EXAMPLES::

            sage: K.<x>=GF(9,'x')[]
            sage: C=HyperellipticCurve(x^7-1,0)
            sage: C._Hasse_Witt_cached()
            (
            [0 0 0]
            [0 0 0]
            [0 0 0], Hyperelliptic Curve over Finite Field in x of size 3^2 defined by y^2 = x^7 + 2
            )

            sage: K.<x>=GF(49,'x')[]
            sage: C=HyperellipticCurve(x^5+1,0)
            sage: C._Hasse_Witt_cached()
            (
            [0 0]
            [0 0], Hyperelliptic Curve over Finite Field in x of size 7^2 defined by y^2 = x^5 + 1
            )

            sage: P.<x>=GF(9,'a')[]
            sage: C=HyperellipticCurve(x^29+1,0)
            sage: C._Hasse_Witt_cached()
            (
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0], Hyperelliptic Curve over Finite Field in a of size 3^2 defined by y^2 = x^29 + 1
            )
        """
        # If Cartier Matrix is already cached for this curve, use that or evaluate it to get M,
        #Coeffs, genus, Fq=base field of self, p=char(Fq). This is so we have one less matrix to
        #compute.

        #We use caching here since Cartier matrix is needed to compute Hasse Witt. So if the Cartier
        #is already computed it is stored in list A. If it was not cached (i.e. A is empty), we simply
        #compute it. If it is cached then we need to make sure that we have the correct one. So check
        #which curve does the matrix correspond to. Since caching stores a lot of stuff, we only check
        #the last entry in A. If it does not match, clear A and compute Cartier.
        #
        #Since Trac Ticket #11115, there is a different cache for methods
        #that don't  accept arguments. Anyway, the easiest is to call
        #the cached method and simply see whether the data belong to self.
        M, Coeffs,g, Fq, p, E= self._Cartier_matrix_cached()
        if E!=self:
            self._Cartier_matrix_cached.clear_cache()
            M, Coeffs,g, Fq, p, E= self._Cartier_matrix_cached()

        #This compute the action of p^kth Frobenius  on list of coefficients
        def frob_mat(Coeffs, k):
            a = p**k
            mat = []
            Coeffs_pow = [c**a for c in Coeffs]
            for i in range(1,g+1):
                H=[(Coeffs[j]) for j in range((p*i-1), (p*i - g-1), -1)]
                mat.append(H);
            return matrix(Fq,mat)

        #Computes all the different possible action of frobenius on matrix M and stores in list Mall
        Mall = [M] + [frob_mat(Coeffs,k) for k in range(1,g)]

        #initial N=I, so we can go through Mall and multiply all matrices with I and
        #get the Hasse-Witt matrix.
        N = identity_matrix(Fq,g)
        for l in Mall:
            N = N*l;
        return N, E

    #This is the function which is actually called by command line
    def Hasse_Witt(self):
        r"""
        INPUT:

        - ``E`` : Hyperelliptic Curve of the form `y^2 = f(x)` over a finite field, `\mathbb{F}_q`

        OUTPUT:

        - ``N`` : The matrix `N = M M^p \dots M^{p^{g-1}}` where `M = c_{pi-j}`, and `f(x)^{(p-1)/2} = \sum c_i x^i`



        Reference-N. Yui. On the Jacobian varieties of hyperelliptic curves over fields of characteristic `p > 2`.

        EXAMPLES::

            sage: K.<x>=GF(9,'x')[]
            sage: C=HyperellipticCurve(x^7-1,0)
            sage: C.Hasse_Witt()
            [0 0 0]
            [0 0 0]
            [0 0 0]

            sage: K.<x>=GF(49,'x')[]
            sage: C=HyperellipticCurve(x^5+1,0)
            sage: C.Hasse_Witt()
            [0 0]
            [0 0]

            sage: P.<x>=GF(9,'a')[]
            sage: E=HyperellipticCurve(x^29+1,0)
            sage: E.Hasse_Witt()
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        """
        # Since Trac Ticket #11115, there is a special
        # type of cached for those methods that don't
        # accept arguments. We want to get Hasse-Witt
        # from the cache - but apparently it could be
        # that the cached value does not belong to self.
        # So, the easiest is:
        N, E= self._Hasse_Witt_cached()
        if E!=self:
            self._Hasse_Witt_cached.clear_cache()
            N, E= self._Hasse_Witt_cached()
        return N

    def a_number(self):
        r"""
        INPUT:

        - ``E``: Hyperelliptic Curve of the form `y^2 = f(x)` over a finite field, `\mathbb{F}_q`

        OUTPUT:

        - ``a`` : a-number


         EXAMPLES::

            sage: K.<x>=GF(49,'x')[]
            sage: C=HyperellipticCurve(x^5+1,0)
            sage: C.a_number()
            1

            sage: K.<x>=GF(9,'x')[]
            sage: C=HyperellipticCurve(x^7-1,0)
            sage: C.a_number()
            1

            sage: P.<x>=GF(9,'a')[]
            sage: E=HyperellipticCurve(x^29+1,0)
            sage: E.a_number()
            5
        """
        #We use caching here since Cartier matrix is needed to compute a_number. So if the Cartier
        #is already computed it is stored in list A. If it was not cached (i.e. A is empty), we simply
        #compute it. If it is cached then we need to make sure that we have the correct one. So check
        #which curve does the matrix correspond to. Since caching stores a lot of stuff, we only check
        #the last entry in A. If it does not match, clear A and compute Cartier.
        # Since Trac Ticket #11115, there is a special cache for methods
        # that don't accept arguments. The easiest is: Call the cached
        # method, and test whether the last entry is self.
        M,Coeffs,g, Fq, p,E= self._Cartier_matrix_cached()
        if E != self:
            self._Cartier_matrix_cached.clear_cache()
            M,Coeffs,g, Fq, p,E= self._Cartier_matrix_cached()
        a=g-rank(M);
        return a;

    def p_rank(self):
        r"""
        INPUT:

        - ``E`` : Hyperelliptic Curve of the form `y^2 = f(x)` over a finite field, `\mathbb{F}_q`

        OUTPUT:

        - ``pr`` :p-rank


        EXAMPLES::

            sage: K.<x>=GF(49,'x')[]
            sage: C=HyperellipticCurve(x^5+1,0)
            sage: C.p_rank()
            0

            sage: K.<x>=GF(9,'x')[]
            sage: C=HyperellipticCurve(x^7-1,0)
            sage: C.p_rank()
            0

            sage: P.<x>=GF(9,'a')[]
            sage: E=HyperellipticCurve(x^29+1,0)
            sage: E.p_rank()
            0
        """
        #We use caching here since Hasse Witt is needed to compute p_rank. So if the Hasse Witt
        #is already computed it is stored in list A. If it was not cached (i.e. A is empty), we simply
        #compute it. If it is cached then we need to make sure that we have the correct one. So check
        #which curve does the matrix correspond to. Since caching stores a lot of stuff, we only check
        #the last entry in A. If it does not match, clear A and compute Hasse Witt.
        # However, it seems a waste of time to manually analyse the cache
        # -- See Trac Ticket #11115
        N,E= self._Hasse_Witt_cached()
        if E!=self:
            self._Hasse_Witt_cached.clear_cache()
            N,E= self._Hasse_Witt_cached()
        pr=rank(N);
        return pr
