"""
Hyperelliptic curves over a finite field

EXAMPLES::

    sage: K.<a> = GF(9, 'a')
    sage: x = polygen(K)
    sage: C = HyperellipticCurve(x^7 - x^5 - 2, x^2 + a)
    sage: C._points_fast_sqrt()
    [(0 : 1 : 0), (2*a : 2*a + 2 : 1), (2*a : 2*a : 1), (a + 1 : a : 1), (a + 1 : a + 1 : 1), (2 : a + 1 : 1), (1 : a + 1 : 1)]
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ, RR, binomial
import hyperelliptic_generic
from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob


class HyperellipticCurve_finite_field(hyperelliptic_generic.HyperellipticCurve_generic):

    def _frobenius_coefficient_bounds(self):
        """
        Computes bounds on coefficients of frobenius polynomial, from
        Weil conjectures, via Kedlaya's paper:
           | a_i | <= (2g choose i) * q**(i/2)

        Return value is a list of integers [B_0, ..., B_2g] so that knowledge
        of a_i mod B_i determines a_i uniquely.

        AUTHORS::
          -- Nick Alexander, massaged by David Harvey (2009-03)
        """
        q = self.base_ring().order()
        sqrtq = RR(q).sqrt()
        g = self.genus()
        Bs = []
        for i in range(2*g + 1):
            B = 2 * binomial(2*g, i) * sqrtq**i
            Bs.append(ZZ(B.ceil()))
        Bs.reverse()
        return Bs


    def _frobenius_coefficient_bound(self):
        """
        Computes bound on number of p-adic digits needed to recover
        frobenius polynomial, i.e. returns B so that knowledge of a_i
        modulo p^B determines a_i uniquely.

        AUTHORS::
          -- Nick Alexander, massaged by David Harvey (2009-03)
        """
        assert self.base_ring().is_finite()
        p = self.base_ring().characteristic()

        Bs = self._frobenius_coefficient_bounds()
        M = ZZ(max(Bs))
        B = M.exact_log(p)
        if p**B < M:
            B += 1
        assert p**B >= M
        return B


    def _frobenius_matrix(self, N=None):
        """
        Compute p-adic frobenius matrix to precision p^N. If N not supplied,
        a default value is selected, which is the minimum needed to recover
        the charpoly unambiguously.

        Currently only implemented using hypellfrob, which means only works
        over GF(p^1), and must have p > (2g+1)(2N-1).

        AUTHORS::
          -- Nick Alexander, massaged by David Harvey (2009-03)
        """
        p = self.base_ring().characteristic()
        f, h = self.hyperelliptic_polynomials()
        if h != 0:
            # need y^2 = f(x)
            raise NotImplementedError, "only implemented for curves y^2 = f(x)"

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

        if N is None:
            N = self._frobenius_coefficient_bound()

        matrix_of_frobenius = hypellfrob(p, N, f)
        matrix_of_frobenius = sign * matrix_of_frobenius
        return matrix_of_frobenius


    def frobenius_polynomial(self):
        """
        Charpoly of frobenius, as an element of ZZ[x].

        TODO::
            -- use naive point counting for small problems
            -- use BSGS on jacobian for some parameter ranges
            -- currently only works for p > (2g-1)(2N-1), where N is the
               working precision (requirement of hypellfrob)
            -- doesn't work over non-prime fields yet
            -- doesn't handle equations y^2 + yh = f for h != 0
            -- depending on genus, can be faster to compute frobenius matrix
               modulo a lower power of p, and then get rest of data via
               group operations on jacobian. For example when g = 3, instead
               of computing charpoly mod p^2, do it only mod p^1 and recover
               remaining 1/2 digit via BSGS. When g = 4, instead of doing it
               mod p^3, do it mod p^2 and then only finitely many candidates
               to test. See Andrew Sutherland's papers for more ideas along
               these lines.

        TESTS::
            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_polynomial()
            x^4 + x^3 - 52*x^2 + 37*x + 1369

        A quadratic twist:
            sage: H = HyperellipticCurve(2*t^5 + 2*t + 4)
            sage: H.frobenius_polynomial()
            x^4 - x^3 - 52*x^2 - 37*x + 1369

        AUTHORS::
            -- David Harvey (2009-03)
        """
        p = self.base_ring().characteristic()
        N = self._frobenius_coefficient_bound()
        # compute chapoly over ZZ and then reduce back
        # (because charpoly of p-adic matrices sometimes loses precision)
        M = self._frobenius_matrix(N=N).change_ring(ZZ)
        f = M.charpoly().list()
        ppow = p**N
        f = [x % ppow for x in f]
        f = [x if 2*x < ppow else x - ppow for x in f]
        return ZZ['x'](f)



    def _points_fast_sqrt(self):
        """
        Count points by enumerating over x and solving the resulting
        quadratic for y.

        EXAMPLES::

            sage: K.<a> = GF(9, 'a')
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^7 - 1, x^2 + a)
            sage: C._points_fast_sqrt()
            [(0 : 1 : 0), (2 : a + 1 : 1), (a : 2*a + 1 : 1), (2*a + 2 : 2*a : 1), (2*a + 2 : 1 : 1), (1 : 2*a + 2 : 1), (1 : 0 : 1)]
            sage: K.<a> = GF(49, 'a')
            sage: x = polygen(K)
            sage: C = HyperellipticCurve(x^5 - x^2 - 1, x^2 + a)
            sage: len(C._points_fast_sqrt())
            31

        TESTS::

            sage: x = polygen(GF(16, 'a'))
            sage: C = HyperellipticCurve(x^5 - x + 1, x^2 + x)
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
        points = [self.point([K(0), one, K(0)], check=True)]

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
                            # y^2 + by + c irreducable, so yields no points
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
        Count points by enumerating over x and solving the resulting
        quadratic for y.

        Caches all square roots ahead of time by sqaring every element of
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
        points = [self.point([K(0), one, K(0)], check=True)]

        if K.characteristic() == 2 or brute_force:
            # quadratic equation doesn't work in characteristic 2
            # but there are only 4 affinte points, so just test them
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
        """
        from sage.rings.finite_field import zech_log_bound
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


