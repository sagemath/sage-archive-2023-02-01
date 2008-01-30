r"""
Formal groups of elliptic curves.

AUTHORS:
    -- William Stein: original implementations
    -- David Harvey: improved asymptotics of some methods
    -- Nick Alexander: separation from ell_generic.py, bugfixes and docstrings
"""

from sage.structure.sage_object import SageObject

import sage.misc.misc as misc
import sage.rings.all as rings
from sage.rings.all import O

import ell_generic

class EllipticCurveFormalGroup(SageObject):
    r"""
    The formal group associated to an elliptic curve.
    """
    def __init__(self, E):
        self.__E = E

    def _repr_(self):
        return "Formal Group associated to the %s"%self.__E

    def curve(self):
        r"""
        The elliptic curve this formal group is associated to.

        EXAMPLES:
            sage: E = EllipticCurve("37a")
            sage: F = E.formal_group()
            sage: F.curve()
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        return self.__E

    def w(self, prec=20):
        r"""
        The formal group power series w.

        INPUT:
            prec -- integer

        OUTPUT:
            a power series with given precision

        DETAILS:
            Return the formal power series
            $$
                   w(t) = t^3 + a_1 t^4 + (a_2 + a_1^2) t^5 + \cdots
            $$
            to precision $O(t^prec)$ of Proposition IV.1.1 of [Silverman
            AEC1].  This is the formal expansion of $w = -1/y$ about the
            formal parameter $t = -x/y$ at $\\infty$.

            The result is cached, and a cached version is returned if
            possible.

        WARNING:
            The resulting power series will have precision prec, but its
            parent PowerSeriesRing will have default precision 20 (or whatever
            the default default is).

        ALGORITHM:
            Uses Newton's method to solve the elliptic curve equation
            at the origin. Complexity is roughly $O(M(n))$ where
            $n$ is the precision and $M(n)$ is the time required to multiply
            polynomials of length $n$ over the coefficient ring of $E$.

        AUTHOR:
            -- David Harvey (2006-09-09): modified to use Newton's method
            instead of a recurrence formula.

        EXAMPLES:
            sage: e = EllipticCurve([0, 0, 1, -1, 0])
            sage: e.formal_group().w(10)
             t^3 + t^6 - t^7 + 2*t^9 + O(t^10)

          Check that caching works:
            sage: e = EllipticCurve([3, 2, -4, -2, 5])
            sage: e.formal_group().w(20)
             t^3 + 3*t^4 + 11*t^5 + 35*t^6 + 101*t^7 + 237*t^8 + 312*t^9 - 949*t^10 - 10389*t^11 - 57087*t^12 - 244092*t^13 - 865333*t^14 - 2455206*t^15 - 4366196*t^16 + 6136610*t^17 + 109938783*t^18 + 688672497*t^19 + O(t^20)
            sage: e.formal_group().w(7)
             t^3 + 3*t^4 + 11*t^5 + 35*t^6 + O(t^7)
            sage: e.formal_group().w(35)
             t^3 + 3*t^4 + 11*t^5 + 35*t^6 + 101*t^7 + 237*t^8 + 312*t^9 - 949*t^10 - 10389*t^11 - 57087*t^12 - 244092*t^13 - 865333*t^14 - 2455206*t^15 - 4366196*t^16 + 6136610*t^17 + 109938783*t^18 + 688672497*t^19 + 3219525807*t^20 + 12337076504*t^21 + 38106669615*t^22 + 79452618700*t^23 - 33430470002*t^24 - 1522228110356*t^25 - 10561222329021*t^26 - 52449326572178*t^27 - 211701726058446*t^28 - 693522772940043*t^29 - 1613471639599050*t^30 - 421817906421378*t^31 + 23651687753515182*t^32 + 181817896829144595*t^33 + 950887648021211163*t^34 + O(t^35)
        """
        prec = max(prec, 0)
        k = self.curve().base_ring()

        try:
            # Try cached version
            w = self.__w
            cached_prec = w.prec()
            R = w.parent()
        except AttributeError:
            # No cached version available
            R = rings.PowerSeriesRing(k, "t")
            w = R([k(0), k(0), k(0), k(1)], 4)
            cached_prec = 4
            self.__w = w

        if prec < cached_prec:
            return R(w, prec)

        # We use the following iteration, which doubles the precision
        # at each step:
        #
        #              z^3 - a_3 w^2 - a_4 z w^2 - 2 a_6 w^3
        # w' = -----------------------------------------------------
        #      1 - a_1 z - a_2 z^2 - 2 a_3 w - 2 a_4 z w - 3 a_6 w^2

        a1, a2, a3, a4, a6 = self.curve().ainvs()
        current_prec = cached_prec
        w = w.truncate()   # work with polynomials instead of power series

        numerator_const = w.parent()([0, 0, 0, 1])      # z^3
        denominator_const = w.parent()([1, -a1, -a2])   # 1 - a_1 z - a_2 z^2

        last_prec = 0
        for next_prec in misc.newton_method_sizes(prec):
            if next_prec > current_prec:
                if w.degree() - 1 > last_prec:
                    # Here it's best to throw away some precision to get us
                    # in sync with the sizes recommended by
                    # newton_method_sizes(). This is especially counter-
                    # intuitive when we throw away almost half of our
                    # cached data!

                    # todo: this might not actually be true, depending on
                    # the overhead of truncate(), which is currently very
                    # high e.g. for NTL based polynomials (but this is
                    # slated to be fixed...)

                    w = w.truncate(last_prec)

                w_squared = w.square()
                w_cubed = (w_squared * w).truncate(next_prec)

                numerator = numerator_const                \
                            -  a3 * w_squared              \
                            -  a4 * w_squared.shift(1)     \
                            -  (2*a6) * w_cubed

                denominator = denominator_const           \
                              - (2*a3) * w                \
                              - (2*a4) * w.shift(1)       \
                              - (3*a6) * w_squared

                # todo: this is quite inefficient, because it gets
                # converted to a power series, then the power series
                # inversion works with polynomials again, and then
                # it gets converted *back* to a power series, and
                # then we convert it to a polynomial again! That's four
                # useless conversions!!!

                inverse = ~R(denominator, prec=next_prec)
                inverse = inverse.truncate(next_prec)

                w = (numerator * inverse).truncate(next_prec)

            last_prec = next_prec

        # convert back to power series
        w = R(w, prec)
        self.__w = (prec, w)
        return self.__w[1]

    def x(self, prec=20):
        r"""
        Return the formal series $x(t) = t/w(t)$ in terms of the local
        parameter $t = -x/y$ at infinity.

        INPUT:
            prec -- integer

        OUTPUT:
            a laurent series with given precision

        DETAILS:
            Return the formal series
            $$
                   x(t) = t^{-2} - a_1 t^{-1} - a_2 - a_3 t - \cdots
            $$
            to precision $O(t^prec)$ of page 113 of [Silverman
            AEC1].

        WARNING:
            The resulting series will have precision prec, but its
            parent PowerSeriesRing will have default precision 20 (or whatever
            the default default is).

        EXAMPLES:
            sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().x(10)
             t^-2 - t + t^2 - t^4 + 2*t^5 - t^6 - 2*t^7 + 6*t^8 - 6*t^9 + O(t^10)
        """
        prec = max(prec, 0)
        y = self.y(prec)
        t = y.parent().gen()
        return -t*y + O(t**prec)

    def y(self, prec=20):
        r"""
        Return the formal series $y(t) = -1/w(t)$ in terms of the local
        parameter $t = -x/y$ at infinity.

        INPUT:
            prec -- integer

        OUTPUT:
            a laurent series with given precision

        DETAILS:
            Return the formal series
            $$
                   y(t) = - t^{-3} + a_1 t^{-2} + a_2 t + a_3 + \cdots
            $$
            to precision $O(t^prec)$ of page 113 of [Silverman
            AEC1].

            The result is cached, and a cached version is returned if
            possible.

        WARNING:
            The resulting series will have precision prec, but its
            parent PowerSeriesRing will have default precision 20 (or whatever
            the default default is).

        EXAMPLES:
            sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().y(10)
             -t^-3 + 1 - t + t^3 - 2*t^4 + t^5 + 2*t^6 - 6*t^7 + 6*t^8 + 3*t^9 + O(t^10)
        """
        prec = max(prec,0)
        try:
            pr, y = self.__y
        except AttributeError:
            pr = -1
        if prec <= pr:
            t = y.parent().gen()
            return y + O(t**prec)
        w = self.w(prec+6) # XXX why 6?
        t = w.parent().gen()
        y = -1/w + O(t**prec)
        self.__y = (prec, y)
        return self.__y[1]

    def differential(self, prec=20):
        r"""
        Returns the power series $f(t) = 1 + \cdots$ such that $f(t) dt$ is
        the usual invariant differential $dx/(2y + a_1 x + a_3)$.

        INPUT:
           prec -- nonnegative integer, answer will be returned $O(t^{\var{prec}})$

        OUTPUT:
            a power series with given precision

        DETAILS:
            Return the formal series
            $$
                   f(t) = 1 + a_1 t + ({a_1}^2 + a_2) t^2 + \cdots
            $$
            to precision $O(t^prec)$ of page 113 of [Silverman
            AEC1].

            The result is cached, and a cached version is returned if
            possible.

        WARNING:
            The resulting series will have precision prec, but its
            parent PowerSeriesRing will have default precision 20 (or whatever
            the default default is).

        EXAMPLES:
           sage: EllipticCurve([-1, 1/4]).formal_group().differential(15)
            1 - 2*t^4 + 3/4*t^6 + 6*t^8 - 5*t^10 - 305/16*t^12 + 105/4*t^14 + O(t^15)
           sage: EllipticCurve(Integers(53), [-1, 1/4]).formal_group().differential(15)
            1 + 51*t^4 + 14*t^6 + 6*t^8 + 48*t^10 + 24*t^12 + 13*t^14 + O(t^15)

        AUTHOR:
           -- David Harvey (2006-09-10): factored out of log
        """
        prec = max(prec,0)
        try:
            cached_prec, omega = self.__omega
        except AttributeError:
            cached_prec = -1
        if prec <= cached_prec:
            return omega.add_bigoh(prec)

        a = self.curve().ainvs()
        x = self.x(prec+1)
        y = self.y(prec+1)
        xprime = x.derivative()
        g = xprime / (2*y + a[0]*x + a[2])
        self.__omega = (prec, g.power_series().add_bigoh(prec))
        return self.__omega[1]

    def log(self, prec=20):
        r"""
        Returns the power series $f(t) = t + \cdots$ which is an isomorphism
        to the additive formal group.

        Generally this only makes sense in characteristic zero, although the
        terms before $t^p$ may work in characteristic $p$.

        INPUT:
           prec -- nonnegative integer

        OUTPUT:
            a power series with given precision

        EXAMPLES:
           sage: EllipticCurve([-1, 1/4]).formal_group().log(15)
            t - 2/5*t^5 + 3/28*t^7 + 2/3*t^9 - 5/11*t^11 - 305/208*t^13 + O(t^15)

        AUTHOR:
           -- David Harvey (2006-09-10): rewrote to use differential
        """
        return self.differential(prec-1).integral().add_bigoh(prec)

    def inverse(self, prec=20):
        r"""
        The formal group inverse law i(t), which satisfies F(t, i(t)) = 0.

        INPUT:
            prec -- integer

        OUTPUT:
            a power series with given precision

        DETAILS:
            Return the formal power series
            $$
                   i(t) = - t + a_1 t^2 + \cdots
            $$
            to precision $O(t^prec)$ of page 114 of [Silverman AEC1].

            The result is cached, and a cached version is returned if
            possible.

        WARNING:
            The resulting power series will have precision prec, but its
            parent PowerSeriesRing will have default precision 20 (or whatever
            the default default is).

        EXAMPLES:
            sage: e = EllipticCurve([1, 2])
            sage: F = e.formal_group().group_law(5)
            sage: i = e.formal_group().inverse(5)
            sage: Fx = F.base_base_extend(i.parent())
            sage: Fx (i) (i.parent().gen())
            O(t^5)
        """
        prec = max(prec,0)
        try:
            pr, inv = self.__inverse
        except AttributeError:
            pr = -1
        if prec <= pr:
            t = inv.parent().gen()
            return inv + O(t**prec)
        x = self.x(prec)
        y = self.y(prec)
        a1, _, a3, _, _ = self.curve().ainvs()
        inv = x / ( y + a1*x + a3)          # page 114 of Silverman, AEC I
        inv = inv.power_series().add_bigoh(prec)
        self.__inverse = (prec, inv)
        return inv

    def group_law(self, prec=10):
        r"""
        The formal group law.

        INPUT:
            prec -- integer

        OUTPUT:
            a power series with given precision in ZZ[[ ZZ[['t1']], 't2']]

        DETAILS:
            Return the formal power series
            $$
                   F(t1, t2) = t1 + t2 - a1 t1 t2 - \cdots
            $$
            to precision $O(t^prec)$ of page 115 of [Silverman AEC1].

            The result is cached, and a cached version is returned if
            possible.

        WARNING:
            The resulting power series will have precision prec, but its
            parent PowerSeriesRing will have default precision 20 (or whatever
            the default default is).

        AUTHOR:
            -- Nick Alexander: minor fixes, docstring

        EXAMPLES:
            sage: e = EllipticCurve([1, 2])
            sage: F = e.formal_group().group_law(5); F
             t1 + O(t1^5) + (1 - 2*t1^4 + O(t1^5))*t2 + (-4*t1^3 + O(t1^5))*t2^2 + (-4*t1^2 - 30*t1^4 + O(t1^5))*t2^3 + (-2*t1 - 30*t1^3 + O(t1^5))*t2^4 + O(t2^5)
            sage: i = e.formal_group().inverse(5)
            sage: Fx = F.base_base_extend(i.parent())
            sage: Fx (i.parent().gen()) (i)
             O(t^5)

            Let's ensure caching with changed precision is working:
            sage: e.formal_group().group_law(4)
             t1 + O(t1^4) + (1 + O(t1^4))*t2 + (-4*t1^3 + O(t1^4))*t2^2 + (-4*t1^2 + O(t1^4))*t2^3 + O(t2^4)
        """
        prec = max(prec,0)
        R1 = rings.PowerSeriesRing(self.curve().base_ring(),"t1")
        R2 = rings.PowerSeriesRing(R1,"t2")
        t1 = R1.gen().add_bigoh(prec)
        t2 = R2.gen().add_bigoh(prec)

        def fix_prec(F, final_prec):
            return R2([ c + O(t1**final_prec) for c in F ]) + O(t2 ** final_prec)

        try:
            pr, F = self.__group_law
        except AttributeError:
            pr = -1
        if prec <= pr:
            # we have to 'fix up' coefficient precisions
            return fix_prec(F, prec)

        w = self.w(prec)
        def tsum(n):
            return sum([t2**m * t1**(n-m-1) for m in range(n)])
        lam = sum([tsum(n)*w[n] for n in range(3,prec)])
        w1 = R1(w, prec)
        nu = w1 - lam*t1 + O(t1**prec)
        a1, a2, a3, a4, a6 = self.curve().ainvs()
        lam2 = lam*lam
        lam3 = lam2*lam
        t3 = -t1 - t2 - \
             (a1*lam + a3*lam2 + a2*nu + 2*a4*lam*nu + 3*a6*lam2*nu)/  \
             (1 + a2*lam + a4*lam2 + a6*lam3)
        inv = self.inverse(prec)

        # we have to 'fix up' precision issues
        F = fix_prec(inv(t3), prec)
        self.__group_law = (prec, F)
        return F

    def mult_by_n(self, n, prec=10):
        """
        The formal 'multiplication by n' endomorphism $[n]$.

        INPUT:
            prec -- integer

        OUTPUT:
            a power series with given precision

        DETAILS:
            Return the formal power series
            $$
                   [n](t) = n t + \cdots
            $$
            to precision $O(t^prec)$ of Proposition 2.3 of [Silverman AEC1].

        WARNING:
            The resulting power series will have precision prec, but its
            parent PowerSeriesRing will have default precision 20 (or whatever
            the default default is).

        AUTHOR:
            -- Nick Alexander: minor fixes, docstring
            -- David Harvey (2007-03): faster algorithm for char 0 field case

        EXAMPLES:
            sage: e = EllipticCurve([1, 2, 3, 4, 6])
            sage: e.formal_group().mult_by_n(0, 5)
             O(t^5)
            sage: e.formal_group().mult_by_n(1, 5)
             t + O(t^5)

        We verify an identity of low degree:

            sage: none = e.formal_group().mult_by_n(-1, 5)
            sage: two = e.formal_group().mult_by_n(2, 5)
            sage: ntwo = e.formal_group().mult_by_n(-2, 5)
            sage: ntwo - none(two)
             O(t^5)
            sage: ntwo - two(none)
             O(t^5)

        It's quite fast:
            sage: E = EllipticCurve("37a"); F = E.formal_group()
            sage: F.mult_by_n(100, 20)
            100*t - 49999950*t^4 + 3999999960*t^5 + 14285614285800*t^7 - 2999989920000150*t^8 + 133333325333333400*t^9 - 3571378571674999800*t^10 + 1402585362624965454000*t^11 - 146666057066712847999500*t^12 + 5336978000014213190385000*t^13 - 519472790950932256570002000*t^14 + 93851927683683567270392002800*t^15 - 6673787211563812368630730325175*t^16 + 320129060335050875009191524993000*t^17 - 45670288869783478472872833214986000*t^18 + 5302464956134111125466184947310391600*t^19 + O(t^20)
        """
        if self.curve().base_ring().is_field() and self.curve().base_ring().characteristic() == 0 and n != 0:
            # The following algorithm only works over a field of
            # characteristic zero. I don't know whether something similar
            # can be done over a general ring. It would be nice if it did,
            # since it's much faster than using the formal group law.
            # -- dmharvey

            # Create a "formal point" on the original curve E.
            # Our anwer only needs prec-1 coefficients (since lowest term
            # is t^1), and x(t) = t^(-2) + ... and y(t) = t^(-3) + ...,
            # so we only need x(t) mod t^(prec-3) and y(t) mod t^(prec-4)
            x = self.x(prec-3)
            y = self.y(prec-4)
            R = x.parent()    # the laurent series ring over the base ring
            X = self.curve().change_ring(R)
            P = X(x, y)

            # and multiply it by n, using the group law on E
            Q = n*P

            # express it in terms of the formal parameter
            return -Q[0] / Q[1]


        # Now the general case, not necessarily over a field.

        # For unknown reasons, this seems to lose one place of precision:
        # the coefficient of t**(prec-1) seems off.  So we apply the easy fix.
        orig_prec = max(prec, 0)
        prec = orig_prec + 1

        R = rings.PowerSeriesRing(self.curve().base_ring(),"t")
        t = R.gen()
        if n == 1:
            return t + O(t**orig_prec)
        if n == 0:
            return R(0) + O(t**orig_prec)
        if n == -1:
            return R(self.inverse(orig_prec))
        if n < 0:
            F = self.inverse(prec)(self.mult_by_n(-n,orig_prec))
            return R(F.add_bigoh(orig_prec))
        F = self.group_law(prec)
        g = F.parent().base_ring().gen()
        for m in range(2,n+1):
            g = F(g)
        return R(g.add_bigoh(orig_prec))

    def sigma(self, prec=10):
        a1,a2,a3,a4,a6 = self.curve().ainvs()

        k = self.curve().base_ring()
        fl = self.log(prec)
        R = rings.PolynomialRing(k,'c'); c = R.gen()
        F = fl.reversion()

        S = rings.LaurentSeriesRing(R,'z')
        c = S(c)
        z = S.gen()
        F = F(z + O(z**prec))
        wp = self.x()(F)
        e2 = 12*c - a1**2 - 4*a2
        g = (1/z**2 - wp + e2/12).power_series()
        h = g.integral().integral()
        sigma_of_z = z.power_series() * h.exp()

        T = rings.PowerSeriesRing(R,'t')
        fl = fl(T.gen()+O(T.gen()**prec))
        sigma_of_t = sigma_of_z(fl)
        return sigma_of_t
