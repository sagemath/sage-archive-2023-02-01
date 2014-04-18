r"""
Series constructor for modular forms for Hecke triangle groups

AUTHORS:

- Based on the thesis of John Garrett Leo (2008)
- Jonas Jermann (2013): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ, QQ, infinity, rising_factorial, PolynomialRing, LaurentSeries, PowerSeriesRing, FractionField
from sage.functions.all import exp

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from hecke_triangle_groups import HeckeTriangleGroup


class MFSeriesConstructor(SageObject,UniqueRepresentation):
    r"""
    Constructor for the Fourier expansion of some
    (specific, basic) modular forms.

    The constructor is used by forms elements in case
    their Fourier expansion is needed or requested.
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, prec=ZZ(10), fix_d=False, set_d=None, d_num_prec=ZZ(53)):
        r"""
        Return a (cached) instance with canonical parameters.

        In particular in case ``fix_d = True`` or if ``set_d`` is
        set then the ``base_ring`` is replaced by the common parent
        of ``base_ring`` and the parent of ``set_d`` (resp. the
        numerical value of ``d`` in case ``fix_d=True``).

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor() == MFSeriesConstructor(3, ZZ, 10, False, None, 53)
            True
            sage: MFSeriesConstructor(base_ring = CC, set_d=CC(1)) == MFSeriesConstructor(set_d=CC(1))
            True
            sage: MFSeriesConstructor(group=4, fix_d=True).base_ring() == QQ
            True
            sage: MFSeriesConstructor(group=5, fix_d=True).base_ring() == RR
            True
        """

        if (group==infinity):
            group = HeckeTriangleGroup(infinity)
        else:
            try:
                group = HeckeTriangleGroup(ZZ(group))
            except TypeError:
                group = HeckeTriangleGroup(group.n())
        prec=ZZ(prec)
        #if (prec<1):
        #    raise Exception("prec must be an Integer >=1")

        fix_d = bool(fix_d)
        if (fix_d):
            n = group.n()
            d = group.dvalue()
            if (group.is_arithmetic()):
                d_num_prec = None
                set_d = 1/base_ring(1/d)
            else:
                d_num_prec = ZZ(d_num_prec)
                set_d = group.dvalue().n(d_num_prec)
        else:
            d_num_prec = None

        if (set_d is not None):
            base_ring=(base_ring(1)*set_d).parent()
        #elif (not base_ring.is_exact()):
        #    raise NotImplementedError

        return super(MFSeriesConstructor,cls).__classcall__(cls, group, base_ring, prec, fix_d, set_d, d_num_prec)

    def __init__(self, group, base_ring, prec, fix_d, set_d, d_num_prec):
        r"""
        Constructor for the Fourier expansion of some
        (specific, basic) modular forms.

        INPUT:

        - ``group``       - A Hecke triangle group (default: HeckeTriangleGroup(3)).

        - ``base_ring``   - The base ring (default: ZZ)

        - ``prec``        - An integer (default: 10), the default precision used
                            in calculations in the LaurentSeriesRing or PowerSeriesRing.

        - ``fix_d``       - ``True`` or ``False`` (default: ``False``).

                            If ``fix_d == False`` the base ring of the power series
                            is (the fraction field) of the polynomial ring over the base
                            ring in one formal parameter ``d``.

                            If ``fix_d == True`` the formal parameter ``d`` is replaced
                            by its numerical value with numerical precision at least ``d_num_prec``
                            (or exact in case n=3, 4, 6). The base ring of the PowerSeriesRing
                            or LaurentSeriesRing is changed to a common parent of
                            ``base_ring`` and the parent of the mentioned value ``d``.

        - ``set_d``       - A number which replaces the formal parameter ``d``.
                            The base ring of the PowerSeriesRing or LaurentSeriesRing is
                            changed to a common parent of ``base_ring``
                            and the parent of the specified value for ``d``.
                            Note that in particular ``set_d=1`` will produce
                            rational Fourier expansions.

        - ``d_num_prec``  - An integer, a lower bound for the precision of the
                            numerical value of ``d``.

        OUTPUT:

        The constructor for Fourier expansion with the specified settings.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFC = MFSeriesConstructor()
            sage: MFC
            Power series constructor for Hecke modular forms for n=3, base ring=Integer Ring
            with (basic series) precision 10 with formal parameter d
            sage: MFC.group()
            Hecke triangle group for n = 3
            sage: MFC.prec()
            10
            sage: MFC.d().parent()
            Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFC._ZZseries_ring
            Power Series Ring in q over Rational Field

            sage: MFSeriesConstructor(set_d=CC(1))
            Power series constructor for Hecke modular forms for n=3, base ring=Complex Field with 53 bits of precision
            with (basic series) precision 10 with parameter d=1.00000000000000
            
            sage: MFSeriesConstructor(group=4, fix_d=True)
            Power series constructor for Hecke modular forms for n=4, base ring=Rational Field
            with (basic series) precision 10 with parameter d=1/256

            sage: MFSeriesConstructor(group=5, fix_d=True)
            Power series constructor for Hecke modular forms for n=5, base ring=Real Field with 53 bits of precision
            with (basic series) precision 10 with parameter d=0.00705223418128563
        """

        self._group          = group
        self._base_ring      = base_ring
        self._prec           = prec
        self._fix_d          = fix_d
        self._set_d          = set_d
        self._d_num_prec     = d_num_prec

        if (set_d):
            self._coeff_ring = FractionField(base_ring)
            self._d          = set_d
        else:
            self._coeff_ring = FractionField(PolynomialRing(base_ring,"d"))
            self._d          = self._coeff_ring.gen()

        self._ZZseries_ring  = PowerSeriesRing(QQ,'q',default_prec=self._prec)
        self._qseries_ring   = PowerSeriesRing(self._coeff_ring,'q',default_prec=self._prec)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, fix_d=True)
            Power series constructor for Hecke modular forms for n=4, base ring=Rational Field
            with (basic series) precision 10 with parameter d=1/256

            sage: MFSeriesConstructor(group=5)
            Power series constructor for Hecke modular forms for n=5, base ring=Integer Ring
            with (basic series) precision 10 with formal parameter d
        """

        if (self._set_d):
            return "Power series constructor for Hecke modular forms for n={}, base ring={} with (basic series) precision {} with parameter d={}".\
                format(self._group.n(), self._base_ring, self._prec, self._d)
        else:
            return "Power series constructor for Hecke modular forms for n={}, base ring={} with (basic series) precision {} with formal parameter d".\
                format(self._group.n(), self._base_ring, self._prec)

    def group(self):
        r"""
        Return the (Hecke triangle) group of ``self``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, fix_d=True).group()
            Hecke triangle group for n = 4
        """

        return self._group

    def hecke_n(self):
        r"""
        Return the parameter ``n`` of the (Hecke triangle) group of ``self``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, fix_d=True).hecke_n()
            4
        """

        return self._group.n()

    def base_ring(self):
        r"""
        Return base ring of ``self``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=5, fix_d=True).base_ring()
            Real Field with 53 bits of precision
            sage: MFSeriesConstructor(group=5, fix_d=True, d_num_prec=100).base_ring()
            Real Field with 100 bits of precision
        """

        return self._base_ring

    def prec(self):
        r"""
        Return the used default precision for the PowerSeriesRing or LaurentSeriesRing.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=5, fix_d=True).prec()
            10
            sage: MFSeriesConstructor(group=5, prec=20).prec()
            20
        """

        return self._prec

    def fix_d(self):
        r"""
        Return whether the numerical value for the parameter
        ``d`` will be substituted or not.
        
        Note: Depending on whether ``set_d`` is ``None`` or
        not ``d`` might still be substituted despite ``fix_d``
        being ``False``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=5, fix_d=True, set_d=1).fix_d()
            True
            sage: MFSeriesConstructor(group=5, fix_d=True, set_d=1).set_d()
            0.00705223418128563
            sage: MFSeriesConstructor(group=5, set_d=1).fix_d()
            False
        """

        return self._fix_d

    def set_d(self):
        r"""
        Return the numerical value which is substituted for
        the parameter ``d``. Default: ``None``, meaning
        the formal parameter ``d`` is used.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=5, fix_d=True, set_d=1).set_d()
            0.00705223418128563
            sage: MFSeriesConstructor(group=5, set_d=1).set_d()
            1
            sage: MFSeriesConstructor(group=5, set_d=1).set_d().parent()
            Integer Ring
        """

        return self._set_d

    def is_exact(self):
        r"""
        Return whether used ``base_ring`` is exact.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, fix_d=True).is_exact()
            True
            sage: MFSeriesConstructor(group=5, fix_d=True).is_exact()
            False
            sage: MFSeriesConstructor(group=5, set_d=1).is_exact()
            True
        """

        return self._base_ring.is_exact()

    def d(self):
        r"""
        Return the formal parameter ``d`` respectively
        its (possibly numerical) value in case ``set_d``
        is not ``None``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, fix_d=True).d()
            1/256
            sage: MFSeriesConstructor(group=4).d()
            d
            sage: MFSeriesConstructor(group=4).d().parent()
            Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, fix_d=True).d()
            0.00705223418128563
            sage: MFSeriesConstructor(group=5, set_d=1).d()
            1
        """

        return self._d

    def q(self):
        r"""
        Return the generator of the used PowerSeriesRing.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, fix_d=True).q()
            q
            sage: MFSeriesConstructor(group=4, fix_d=True).q().parent()
            Power Series Ring in q over Rational Field
            sage: MFSeriesConstructor(group=5, fix_d=True).q().parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        return self._qseries_ring.gen()

    def coeff_ring(self):
        r"""
        Return coefficient ring of ``self``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, fix_d=True).coeff_ring()
            Rational Field
            sage: MFSeriesConstructor(group=4).coeff_ring()
            Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, fix_d=True).coeff_ring()
            Real Field with 53 bits of precision
            sage: MFSeriesConstructor(group=5).coeff_ring()
            Fraction Field of Univariate Polynomial Ring in d over Integer Ring
        """

        return self._coeff_ring

    def qseries_ring(self):
        r"""
        Return the used PowerSeriesRing.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, fix_d=True).qseries_ring()
            Power Series Ring in q over Rational Field
            sage: MFSeriesConstructor(group=4).qseries_ring()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, fix_d=True).qseries_ring()
            Power Series Ring in q over Real Field with 53 bits of precision
            sage: MFSeriesConstructor(group=5).qseries_ring()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
        """

        return self._qseries_ring

    @cached_method
    def J_inv_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``J_inv``,
        where ``d`` is replaced by ``1``.

        This is the main function used to determine all Fourier expansions!

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).J_inv_ZZ()
            q^-1 + 31/72 + 1823/27648*q + O(q^2)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).J_inv_ZZ()
            q^-1 + 79/200 + 42877/640000*q + O(q^2)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).J_inv_ZZ().parent()
            Laurent Series Ring in q over Rational Field
        """

        F1       = lambda a,b:   self._ZZseries_ring(\
                       [ ZZ(0) ] + [\
                           rising_factorial(a,k) * rising_factorial(b,k) / (ZZ(k).factorial())**2 * sum([\
                               ZZ(1)/(a+j)+ZZ(1)/(b+j)-ZZ(2)/ZZ(1+j) for j in range(ZZ(0),ZZ(k))\
                           ]) for k in range(ZZ(1),ZZ(self._prec+1))
                       ], ZZ(self._prec+1)\
                   )
        F        = lambda a,b,c: self._ZZseries_ring([\
                       rising_factorial(a,k) * rising_factorial(b,k) / rising_factorial(c,k) / (ZZ(k).factorial())\
                       for k in range(ZZ(0),ZZ(self._prec+1))\
                   ], ZZ(self._prec+1))
        a        = self._group.alpha()
        b        = self._group.beta()
        Phi      = F1(a,b) / F(a,b,ZZ(1))
        q        = self._ZZseries_ring.gen()
        J_inv_ZZ = ZZ(1) / ((q*Phi.exp()).reversion())
        return J_inv_ZZ

    @cached_method
    def J_inv(self):
        r"""
        Return the Fourier expansion of ``J_inv``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3, fix_d=True).J_inv()
            1/1728*q^-1 + 31/72 + 1823/16*q + O(q^2)
            sage: MFSeriesConstructor(prec=3).J_inv_ZZ() == MFSeriesConstructor(prec=3, set_d=1).J_inv()
            True

            sage: MFSeriesConstructor(group=5, prec=3).J_inv()
            d*q^-1 + 79/200 + 42877/(640000*d)*q + O(q^2)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).J_inv()
            0.00705223418128563*q^-1 + 0.395000000000000 + 9.49987064777062*q + O(q^2)

            sage: MFSeriesConstructor(group=5, prec=3).J_inv().parent()
            Laurent Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).J_inv().parent()
            Laurent Series Ring in q over Real Field with 53 bits of precision
        """

        return self.J_inv_ZZ()(self._qseries_ring.gen()/self._d)

    @cached_method
    def F_rho_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``F_rho``,
        where ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).F_rho_ZZ()
            1 + 5/36*q + 5/6912*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_rho_ZZ()
            1 + 7/100*q + 21/160000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_rho_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._ZZseries_ring.gen()
        n = self.hecke_n()
        temp_expr = ((-q*self.J_inv_ZZ().derivative())**2/(self.J_inv_ZZ()*(self.J_inv_ZZ()-1))).power_series()
        F_rho_ZZ = (temp_expr.log()/(n-2)).exp()
        return F_rho_ZZ

    @cached_method
    def F_rho(self):
        r"""
        Return the Fourier expansion of ``F_rho``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3, fix_d=True).F_rho()
            1 + 240*q + 2160*q^2 + O(q^3)

            sage: MFSeriesConstructor(prec=3).F_rho_ZZ() == MFSeriesConstructor(prec=3, set_d=1).F_rho()
            True

            sage: MFSeriesConstructor(group=5, prec=3).F_rho()
            1 + 7/(100*d)*q + 21/(160000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_rho()
            1.00000000000000 + 9.92593243510795*q + 2.63903932249093*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).F_rho().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_rho().parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        return self.F_rho_ZZ()(self._qseries_ring.gen()/self._d)

    @cached_method
    def F_i_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``F_i``,
        where ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).F_i_ZZ()
            1 - 7/24*q - 77/13824*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_i_ZZ()
            1 - 13/40*q - 351/64000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_i_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._ZZseries_ring.gen()
        n = self.hecke_n()
        temp_expr = ((-q*self.J_inv_ZZ().derivative())**n/(self.J_inv_ZZ()**(n-1)*(self.J_inv_ZZ()-1))).power_series()
        F_i_ZZ = (temp_expr.log()/(n-2)).exp()
        return F_i_ZZ

    @cached_method
    def F_i(self):
        r"""
        Return the Fourier expansion of ``F_i``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3, fix_d=True).F_i()
            1 - 504*q - 16632*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).F_i_ZZ() == MFSeriesConstructor(prec=3, set_d=1).F_i()
            True

            sage: MFSeriesConstructor(group=5, prec=3).F_i()
            1 - 13/(40*d)*q - 351/(64000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_i()
            1.00000000000000 - 46.0846863058583*q - 110.274143118371*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).F_i().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_i().parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        return self.F_i_ZZ()(self._qseries_ring.gen()/self._d)

    @cached_method
    def F_inf_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``F_inf``,
        where ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).F_inf_ZZ()
            q - 1/72*q^2 + 7/82944*q^3 + O(q^4)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_inf_ZZ()
            q - 9/200*q^2 + 279/640000*q^3 + O(q^4)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_inf_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._ZZseries_ring.gen()
        n = self.hecke_n()
        temp_expr  = ((-q*self.J_inv_ZZ().derivative())**(2*n)/(self.J_inv_ZZ()**(2*n-2)*(self.J_inv_ZZ()-1)**n)/q**(n-2)).power_series()
        F_inf_ZZ = (temp_expr.log()/(n-2)).exp()*q
        return F_inf_ZZ

    @cached_method
    def F_inf(self):
        r"""
        Return the Fourier expansion of ``F_inf``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3, fix_d=True).F_inf()
            q - 24*q^2 + 252*q^3 + O(q^4)
            sage: MFSeriesConstructor(prec=3).F_inf_ZZ() == MFSeriesConstructor(prec=3, set_d=1).F_inf()
            True

            sage: MFSeriesConstructor(group=5, prec=3).F_inf()
            q - 9/(200*d)*q^2 + 279/(640000*d^2)*q^3 + O(q^4)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_inf()
            0.000000000000000 + 1.00000000000000*q - 6.38095656542654*q^2 + 8.76538060684488*q^3 + O(q^4)

            sage: MFSeriesConstructor(group=5, prec=3).F_inf().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).F_inf().parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        return self._d*self.F_inf_ZZ()(self._qseries_ring.gen()/self._d)

    @cached_method
    def G_inv_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``G_inv``,
        where ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, prec=3).G_inv_ZZ()
            q^-1 - 3/32 - 955/16384*q + O(q^2)
            sage: MFSeriesConstructor(group=8, prec=3, fix_d=True).G_inv_ZZ()
            q^-1 - 15/128 - 15139/262144*q + O(q^2)
            sage: MFSeriesConstructor(group=8, prec=3, fix_d=True).G_inv_ZZ().parent()
            Laurent Series Ring in q over Rational Field
        """

        n = self.hecke_n()
        if (ZZ(2).divides(n)):
            return self.F_i_ZZ()*(self.F_rho_ZZ()**(ZZ(n/ZZ(2))))/self.F_inf_ZZ()
        else:
            #return self._qseries_ring([])
            raise Exception("G_inv doesn't exist for n={}.".format(self.hecke_n()))
    @cached_method
    def G_inv(self):
        r"""
        Return the Fourier expansion of ``G_inv``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, prec=3, fix_d=True).G_inv()
            1/16777216*q^-1 - 3/2097152 - 955/4194304*q + O(q^2)
            sage: MFSeriesConstructor(group=4, prec=3).G_inv_ZZ() == MFSeriesConstructor(group=4, prec=3, set_d=1).G_inv()
            True

            sage: MFSeriesConstructor(group=8, prec=3).G_inv()
            d^3*q^-1 - 15*d^2/128 - 15139*d/262144*q + O(q^2)
            sage: MFSeriesConstructor(group=8, prec=3, fix_d=True).G_inv()
            1.64838830030189e-6*q^-1 - 0.0000163526310530017 - 0.000682197999433738*q + O(q^2)

            sage: MFSeriesConstructor(group=8, prec=3).G_inv().parent()
            Laurent Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=8, prec=3, fix_d=True).G_inv().parent()
            Laurent Series Ring in q over Real Field with 53 bits of precision
        """

        return (self._d)**2*self.G_inv_ZZ()(self._qseries_ring.gen()/self._d)

    @cached_method
    def E4_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``E_4``,
        where ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).E4_ZZ()
            1 + 5/36*q + 5/6912*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E4_ZZ()
            1 + 21/100*q + 483/32000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E4_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._ZZseries_ring.gen()
        E4_ZZ = ((-q*self.J_inv_ZZ().derivative())**2/(self.J_inv_ZZ()*(self.J_inv_ZZ()-1))).power_series()
        return E4_ZZ

    @cached_method
    def E4(self):
        r"""
        Return the Fourier expansion of ``E_4``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3, fix_d=True).E4()
            1 + 240*q + 2160*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).E4_ZZ() == MFSeriesConstructor(prec=3, set_d=1).E4()
            True

            sage: MFSeriesConstructor(group=5, prec=3).E4()
            1 + 21/(100*d)*q + 483/(32000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E4()
            1.00000000000000 + 29.7777973053239*q + 303.489522086457*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).E4().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E4().parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        return self.E4_ZZ()(self._qseries_ring.gen()/self._d)

    @cached_method
    def E6_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``E_6``,
        where ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).E6_ZZ()
            1 - 7/24*q - 77/13824*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E6_ZZ()
            1 - 37/200*q - 14663/320000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E6_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._ZZseries_ring.gen()
        n = self.hecke_n()
        E6_ZZ = ((-q*self.J_inv_ZZ().derivative())**3/(self.J_inv_ZZ()**2*(self.J_inv_ZZ()-1))).power_series()
        return E6_ZZ

    @cached_method
    def E6(self):
        r"""
        Return the Fourier expansion of ``E_6``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3, fix_d=True).E6()
            1 - 504*q - 16632*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).E6_ZZ() == MFSeriesConstructor(prec=3, set_d=1).E6()
            True

            sage: MFSeriesConstructor(group=5, prec=3).E6()
            1 - 37/(200*d)*q - 14663/(320000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E6()
            1.00000000000000 - 26.2328214356424*q - 921.338894897250*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).E6().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E6().parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        return self.E6_ZZ()(self._qseries_ring.gen()/self._d)

    @cached_method
    def Delta_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``Delta``,
        where ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).Delta_ZZ()
            q - 1/72*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).Delta_ZZ()
            71/50*q + 28267/16000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).Delta_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        n = self.hecke_n()
        return self.E4_ZZ()**(2*n-6)*(self.E4_ZZ()**n-self.E6_ZZ()**2)

    @cached_method
    def Delta(self):
        r"""
        Return the Fourier expansion of ``Delta``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3, fix_d=True).Delta()
            q - 24*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).Delta_ZZ() == MFSeriesConstructor(prec=3, set_d=1).Delta()
            True

            sage: MFSeriesConstructor(group=5, prec=3).Delta()
            71/50*q + 28267/(16000*d)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).Delta()
            0.000000000000000 + 1.42000000000000*q + 250.514582270711*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).Delta().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).Delta().parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        return (self._d)*self.Delta_ZZ()(self._qseries_ring.gen()/self._d)

    @cached_method
    def E2_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``E2``,
        where ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).E2_ZZ()
            1 - 1/72*q - 1/41472*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E2_ZZ()
            1 - 9/200*q - 369/320000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E2_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._ZZseries_ring.gen()
        E2_ZZ = (q*self.F_inf_ZZ().derivative())/self.F_inf_ZZ()
        return E2_ZZ

    @cached_method
    def E2(self):
        r"""
        Return the Fourier expansion of ``E2``.

        EXAMPLES::

            sage: from sage.modular.hecke_mf.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3, fix_d=True).E2()
            1 - 24*q - 72*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).E2_ZZ() == MFSeriesConstructor(prec=3, set_d=1).E2()
            True

            sage: MFSeriesConstructor(group=5, prec=3).E2()
            1 - 9/(200*d)*q - 369/(320000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E2()
            1.00000000000000 - 6.38095656542654*q - 23.1858454761703*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).E2().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3, fix_d=True).E2().parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        return self.E2_ZZ()(self._qseries_ring.gen()/self._d)
