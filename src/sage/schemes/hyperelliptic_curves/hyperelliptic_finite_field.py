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

import hyperelliptic_generic

class HyperellipticCurve_finite_field(hyperelliptic_generic.HyperellipticCurve_generic):

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


