r"""
Hyperelliptic curves over a finite field

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
#  Copyright (C) 2010  Alyson Deines <aly.deines@gmail.com>, Marina Gresham
#  <marina.gresham@coloradocollege.edu>, Gagan Sekhon <gagan.d.sekhon@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ, RR
from sage.rings.arith import binomial
import hyperelliptic_generic
from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import identity_matrix, matrix
from sage.misc.functional import rank


class HyperellipticCurve_finite_field(hyperelliptic_generic.HyperellipticCurve_generic):

    def _frobenius_coefficient_bound(self):
        """
        Computes bound on number of p-adic digits needed to recover
        frobenius polynomial, i.e. returns B so that knowledge of
        a_1, ..., a_g modulo p^B determine frobenius polynomial uniquely.

        TESTS::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound()
            2
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound()
            3

            sage: R.<t> = PolynomialRing(GF(next_prime(10^9)))
            sage: HyperellipticCurve(t^3 + t + 1)._frobenius_coefficient_bound()
            1
            sage: HyperellipticCurve(t^5 + t + 1)._frobenius_coefficient_bound()
            2
            sage: HyperellipticCurve(t^7 + t + 1)._frobenius_coefficient_bound()
            2
            sage: HyperellipticCurve(t^9 + t + 1)._frobenius_coefficient_bound()
            3
            sage: HyperellipticCurve(t^11 + t + 1)._frobenius_coefficient_bound()
            3
            sage: HyperellipticCurve(t^13 + t + 1)._frobenius_coefficient_bound()
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


    def _frobenius_matrix(self, N=None):
        """
        Compute p-adic frobenius matrix to precision p^N. If N not supplied,
        a default value is selected, which is the minimum needed to recover
        the charpoly unambiguously.

        Currently only implemented using hypellfrob, which means only works
        over GF(p^1), and must have p > (2g+1)(2N-1).

        TESTS::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H._frobenius_matrix()
            [1258 + O(37^2)  925 + O(37^2)  132 + O(37^2)  587 + O(37^2)]
            [1147 + O(37^2)  814 + O(37^2)  241 + O(37^2) 1011 + O(37^2)]
            [1258 + O(37^2) 1184 + O(37^2) 1105 + O(37^2)  482 + O(37^2)]
            [1073 + O(37^2)  999 + O(37^2)  772 + O(37^2)  929 + O(37^2)]

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

        TESTS::

            sage: R.<t> = PolynomialRing(GF(37))
            sage: H = HyperellipticCurve(t^5 + t + 2)
            sage: H.frobenius_polynomial()
            x^4 + x^3 - 52*x^2 + 37*x + 1369

        A quadratic twist:

        ::

            sage: H = HyperellipticCurve(2*t^5 + 2*t + 4)
            sage: H.frobenius_polynomial()
            x^4 - x^3 - 52*x^2 - 37*x + 1369

        """
        p = self.base_ring().characteristic()
        g = self.genus()
        N = self._frobenius_coefficient_bound()
        # compute chapoly over ZZ and then reduce back
        # (because charpoly of p-adic matrices sometimes loses precision)
        M = self._frobenius_matrix(N=N).change_ring(ZZ)

        # get a_g, ..., a_0 in ZZ (i.e. with correct signs)
        f = M.charpoly().list()[g:2*g+1]
        ppow = p**N
        f = [x % ppow for x in f]
        f = [x if 2*x < ppow else x - ppow for x in f]

        # get a_{2g}, ..., a_{g+1}
        f = [f[g-i] * p**(g-i) for i in range(g)] + f

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
        Count points by enumerating over x and solving the resulting
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
            raise ValueError, "p must be odd";


        g = self.genus()

        #retrieve the function f(x) ,where y^2=f(x)
        f,h = self.hyperelliptic_polynomials()
        #This implementation only deals with h=0
        if h!=0:
            raise ValueError, "E must be of the form y^2 = f(x)"

        d = f.degree()
        #this implementation is for odd degree only, even degree will be handled later.
        if d%2 == 0:
            raise ValueError, "In this implementation the degree of f must be odd"
        #Compute resultant to make sure no repeated roots
        df=f.derivative()
        R=df.resultant(f)
        if R == 0:
            raise ValueError, "curve is not smooth"

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

