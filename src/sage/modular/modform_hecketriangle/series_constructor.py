r"""
Series constructor for modular forms for Hecke triangle groups

AUTHORS:

- Based on the thesis of John Garrett Leo (2008)
- Jonas Jermann (2013): initial version

.. NOTE:

   ``J_inv_ZZ`` is the main function used to determine all Fourier expansions.
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
from sage.rings.big_oh import O
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
    def __classcall__(cls, group = HeckeTriangleGroup(3), prec=ZZ(10)):
        r"""
        Return a (cached) instance with canonical parameters.

        .. NOTE:

            For each choice of group and precision the constructor is
            cached (only) once. Further calculations with different
            base rings and possibly numerical parameters are based on
            the same cached instance.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor() == MFSeriesConstructor(3, 10)
            True
            sage: MFSeriesConstructor(group=4).hecke_n()
            4
            sage: MFSeriesConstructor(group=5, prec=12).prec()
            12
        """

        if (group==infinity):
            group = HeckeTriangleGroup(infinity)
        else:
            try:
                group = HeckeTriangleGroup(ZZ(group))
            except TypeError:
                group = HeckeTriangleGroup(group.n())
        prec=ZZ(prec)
        # We don't need this assumption the precision may in principle also be negative.
        # if (prec<1):
        #     raise Exception("prec must be an Integer >=1")

        return super(MFSeriesConstructor,cls).__classcall__(cls, group, prec)

    def __init__(self, group, prec):
        r"""
        Constructor for the Fourier expansion of some
        (specific, basic) modular forms.

        INPUT:

        - ``group``      -- A Hecke triangle group (default: HeckeTriangleGroup(3)).

        - ``prec``       -- An integer (default: 10), the default precision used
                            in calculations in the LaurentSeriesRing or PowerSeriesRing.

        OUTPUT:

        The constructor for Fourier expansion with the specified settings.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFC = MFSeriesConstructor()
            sage: MFC
            Power series constructor for Hecke modular forms for n=3 with (basic series) precision 10
            sage: MFC.group()
            Hecke triangle group for n = 3
            sage: MFC.prec()
            10
            sage: MFC._series_ring
            Power Series Ring in q over Rational Field

            sage: MFSeriesConstructor(group=4)
            Power series constructor for Hecke modular forms for n=4 with (basic series) precision 10
            sage: MFSeriesConstructor(group=5, prec=12)
            Power series constructor for Hecke modular forms for n=5 with (basic series) precision 12
        """

        self._group          = group
        self._prec           = prec
        self._series_ring    = PowerSeriesRing(QQ,'q',default_prec=self._prec)

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4)
            Power series constructor for Hecke modular forms for n=4 with (basic series) precision 10

            sage: MFSeriesConstructor(group=5, prec=12)
            Power series constructor for Hecke modular forms for n=5 with (basic series) precision 12
        """

        return "Power series constructor for Hecke modular forms for n={} with (basic series) precision {}".\
                format(self._group.n(), self._prec)

    def group(self):
        r"""
        Return the (Hecke triangle) group of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4).group()
            Hecke triangle group for n = 4
        """

        return self._group

    def hecke_n(self):
        r"""
        Return the parameter ``n`` of the (Hecke triangle) group of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4).hecke_n()
            4
        """

        return self._group.n()

    def prec(self):
        r"""
        Return the used default precision for the PowerSeriesRing or LaurentSeriesRing.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=5).prec()
            10
            sage: MFSeriesConstructor(group=5, prec=20).prec()
            20
        """

        return self._prec

    @cached_method
    def J_inv_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``J_inv``,
        where the parameter ``d`` is replaced by ``1``.

        This is the main function used to determine all Fourier expansions!

        .. TODO:

          The functions that are used in this implementation are
          products of hypergeometric series with other, elementary,
          functions.  Implement them and clean up this representation.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).J_inv_ZZ()
            q^-1 + 31/72 + 1823/27648*q + O(q^2)
            sage: MFSeriesConstructor(group=5, prec=3).J_inv_ZZ()
            q^-1 + 79/200 + 42877/640000*q + O(q^2)
            sage: MFSeriesConstructor(group=5, prec=3).J_inv_ZZ().parent()
            Laurent Series Ring in q over Rational Field
        """

        F1       = lambda a,b:   self._series_ring(
                       [ ZZ(0) ]
                       + [
                           rising_factorial(a,k) * rising_factorial(b,k) / (ZZ(k).factorial())**2
                           * sum(ZZ(1)/(a+j) + ZZ(1)/(b+j) - ZZ(2)/ZZ(1+j)
                                  for j in range(ZZ(0),ZZ(k))
                             )
                           for k in range(ZZ(1), ZZ(self._prec+1))
                       ],
                       ZZ(self._prec+1)
                   )

        F        = lambda a,b,c: self._series_ring(
                       [
                         rising_factorial(a,k) * rising_factorial(b,k) / rising_factorial(c,k) / ZZ(k).factorial()
                         for k in range(ZZ(0), ZZ(self._prec+1))
                       ],
                       ZZ(self._prec+1)
                   )
        a        = self._group.alpha()
        b        = self._group.beta()
        Phi      = F1(a,b) / F(a,b,ZZ(1))
        q        = self._series_ring.gen()

        # the current implementation of power series reversion is slow
        # J_inv_ZZ = ZZ(1) / ((q*Phi.exp()).reversion())

        temp_f   = (q*Phi.exp()).polynomial()
        new_f    = temp_f.revert_series(temp_f.degree()+1)
        J_inv_ZZ = ZZ(1) / (new_f + O(q**(temp_f.degree()+1)))

        return J_inv_ZZ

    @cached_method
    def f_rho_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``f_rho``,
        where the parameter ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).f_rho_ZZ()
            1 + 5/36*q + 5/6912*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).f_rho_ZZ()
            1 + 7/100*q + 21/160000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).f_rho_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._series_ring.gen()
        n = self.hecke_n()
        temp_expr = ((-q*self.J_inv_ZZ().derivative())**2/(self.J_inv_ZZ()*(self.J_inv_ZZ()-1))).power_series()
        f_rho_ZZ = (temp_expr.log()/(n-2)).exp()
        return f_rho_ZZ

    @cached_method
    def f_i_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``f_i``,
        where the parameter ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).f_i_ZZ()
            1 - 7/24*q - 77/13824*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).f_i_ZZ()
            1 - 13/40*q - 351/64000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).f_i_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._series_ring.gen()
        n = self.hecke_n()
        temp_expr = ((-q*self.J_inv_ZZ().derivative())**n/(self.J_inv_ZZ()**(n-1)*(self.J_inv_ZZ()-1))).power_series()
        f_i_ZZ = (temp_expr.log()/(n-2)).exp()
        return f_i_ZZ

    @cached_method
    def f_inf_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``f_inf``,
        where the parameter ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).f_inf_ZZ()
            q - 1/72*q^2 + 7/82944*q^3 + O(q^4)
            sage: MFSeriesConstructor(group=5, prec=3).f_inf_ZZ()
            q - 9/200*q^2 + 279/640000*q^3 + O(q^4)
            sage: MFSeriesConstructor(group=5, prec=3).f_inf_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._series_ring.gen()
        n = self.hecke_n()
        temp_expr  = ((-q*self.J_inv_ZZ().derivative())**(2*n)/(self.J_inv_ZZ()**(2*n-2)*(self.J_inv_ZZ()-1)**n)/q**(n-2)).power_series()
        f_inf_ZZ = (temp_expr.log()/(n-2)).exp()*q
        return f_inf_ZZ

    @cached_method
    def G_inv_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``G_inv``,
        where the parameter ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, prec=3).G_inv_ZZ()
            q^-1 - 3/32 - 955/16384*q + O(q^2)
            sage: MFSeriesConstructor(group=8, prec=3).G_inv_ZZ()
            q^-1 - 15/128 - 15139/262144*q + O(q^2)
            sage: MFSeriesConstructor(group=8, prec=3).G_inv_ZZ().parent()
            Laurent Series Ring in q over Rational Field
        """

        n = self.hecke_n()
        if (ZZ(2).divides(n)):
            return self.f_i_ZZ()*(self.f_rho_ZZ()**(ZZ(n/ZZ(2))))/self.f_inf_ZZ()
        else:
            #return self._qseries_ring([])
            raise ValueError("G_inv doesn't exist for n={}.".format(self.hecke_n()))

    @cached_method
    def E4_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``E_4``,
        where the parameter ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).E4_ZZ()
            1 + 5/36*q + 5/6912*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).E4_ZZ()
            1 + 21/100*q + 483/32000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).E4_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._series_ring.gen()
        E4_ZZ = ((-q*self.J_inv_ZZ().derivative())**2/(self.J_inv_ZZ()*(self.J_inv_ZZ()-1))).power_series()
        return E4_ZZ

    @cached_method
    def E6_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``E_6``,
        where the parameter ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).E6_ZZ()
            1 - 7/24*q - 77/13824*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).E6_ZZ()
            1 - 37/200*q - 14663/320000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).E6_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._series_ring.gen()
        n = self.hecke_n()
        E6_ZZ = ((-q*self.J_inv_ZZ().derivative())**3/(self.J_inv_ZZ()**2*(self.J_inv_ZZ()-1))).power_series()
        return E6_ZZ

    @cached_method
    def Delta_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``Delta``,
        where the parameter ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).Delta_ZZ()
            q - 1/72*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).Delta_ZZ()
            71/50*q + 28267/16000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).Delta_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        n = self.hecke_n()
        return self.E4_ZZ()**(2*n-6)*(self.E4_ZZ()**n-self.E6_ZZ()**2)

    @cached_method
    def E2_ZZ(self):
        r"""
        Return the rational Fourier expansion of ``E2``,
        where the parameter ``d`` is replaced by ``1``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).E2_ZZ()
            1 - 1/72*q - 1/41472*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).E2_ZZ()
            1 - 9/200*q - 369/320000*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).E2_ZZ().parent()
            Power Series Ring in q over Rational Field
        """

        q = self._series_ring.gen()
        E2_ZZ = (q*self.f_inf_ZZ().derivative())/self.f_inf_ZZ()
        return E2_ZZ

    @cached_method
    def series_data(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Determine a set of useful series data associated to the specified parameters.

        INPUT:

        - ``base_ring``    -- The base ring (default: ZZ)

        - ``fix_d``        -- ``True`` or ``False`` (default: ``False``).

                              If ``fix_d == False`` the base ring of the power series
                              is (the fraction field) of the polynomial ring over the base
                              ring in one formal parameter ``d``.

                              If ``fix_d == True`` the formal parameter ``d`` is replaced
                              by its numerical value with numerical precision at least ``d_num_prec``
                              (or exact in case n=3, 4, 6). The base ring of the PowerSeriesRing
                              or LaurentSeriesRing is changed to a common parent of
                              ``base_ring`` and the parent of the mentioned value ``d``.

        - ``d``            -- A number which replaces the formal parameter ``d``.
                              The base ring of the PowerSeriesRing or LaurentSeriesRing is
                              changed to a common parent of ``base_ring``
                              and the parent of the specified value for ``d``.
                              Note that in particular ``d=1`` will produce
                              rational Fourier expansions.

        - ``d_num_prec``   -- An integer, a lower bound for the precision of the
                              numerical value of ``d``.

        OUTPUT:

        - ``base_ring``    -- The base ring used for the series construction.
                              Note that this can be used to check whether the
                              series expansion is exact.

        - ``coeff_ring``   -- The coefficient ring of the Fourier series.

        - ``qseries_ring`` -- The basic power series ring for series.
                              Note that the resulting series might instead
                              lie in a Laurent series ring.

        - ``d``            -- The either formal or explicit parameter ``d``
                              for the series.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFC = MFSeriesConstructor()
            sage: (base_ring, coeff_ring, qseries_ring, d) = MFC.series_data()
            sage: base_ring
            Integer Ring
            sage: coeff_ring
            Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: qseries_ring
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: d
            d
            sage: d.parent()
            Fraction Field of Univariate Polynomial Ring in d over Integer Ring

            sage: (base_ring, coeff_ring, qseries_ring, d) = MFSeriesConstructor().series_data(d=CC(1))
            sage: base_ring
            Complex Field with 53 bits of precision
            sage: coeff_ring
            Complex Field with 53 bits of precision
            sage: qseries_ring
            Power Series Ring in q over Complex Field with 53 bits of precision
            sage: d
            1.00000000000...

            sage: (base_ring, coeff_ring, qseries_ring, d) = MFSeriesConstructor(group=4).series_data(fix_d=True)
            sage: base_ring
            Rational Field
            sage: coeff_ring
            Rational Field
            sage: qseries_ring
            Power Series Ring in q over Rational Field
            sage: d
            1/256

            sage: (base_ring, coeff_ring, qseries_ring, d) = MFSeriesConstructor(group=5).series_data(fix_d=True)
            sage: base_ring
            Real Field with 53 bits of precision
            sage: coeff_ring
            Real Field with 53 bits of precision
            sage: qseries_ring
            Power Series Ring in q over Real Field with 53 bits of precision
            sage: d
            0.00705223418128...


        .. NOTE:

            This function may return different return data for different arguments
            but when calculations are done with the data the underlying cached instance
            of the series constructor remains the same in all cases.
        """

        fix_d = bool(fix_d)

        if (fix_d):
            d = self._group.dvalue()
            if (self._group.is_arithmetic()):
                d_num_prec = None
                d = 1 / base_ring(1/d)
            else:
                d_num_prec = ZZ(d_num_prec)
                d = self._group.dvalue().n(d_num_prec)
        else:
            d_num_prec = None

        if (d is not None):
            base_ring = (base_ring(1) * d).parent()

        if (d):
            coeff_ring = FractionField(base_ring)
        else:
            coeff_ring = FractionField(PolynomialRing(base_ring, "d"))
            d          = coeff_ring.gen()

        qseries_ring = PowerSeriesRing(coeff_ring, 'q', default_prec=self._prec)

        return (base_ring, coeff_ring, qseries_ring, d)


    @cached_method
    def J_inv(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Return the Fourier expansion of ``J_inv``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).J_inv(fix_d=True)
            1/1728*q^-1 + 31/72 + 1823/16*q + O(q^2)
            sage: MFSeriesConstructor(prec=3).J_inv_ZZ() == MFSeriesConstructor(prec=3).J_inv(d=1)
            True

            sage: MFSeriesConstructor(group=5, prec=3).J_inv()
            d*q^-1 + 79/200 + 42877/(640000*d)*q + O(q^2)
            sage: MFSeriesConstructor(group=5, prec=3).J_inv(fix_d=True)
            0.007052234181285...*q^-1 + 0.3950000000000... + 9.499870647770...*q + O(q^2)

            sage: MFSeriesConstructor(group=5, prec=3).J_inv().parent()
            Laurent Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3).J_inv(fix_d=True).parent()
            Laurent Series Ring in q over Real Field with 53 bits of precision
        """

        (base_ring, coeff_ring, qseries_ring, d) = self.series_data(base_ring, fix_d, d, d_num_prec)
        return self.J_inv_ZZ()(qseries_ring.gen()/d)

    @cached_method
    def f_rho(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Return the Fourier expansion of ``f_rho``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).f_rho(fix_d=True)
            1 + 240*q + 2160*q^2 + O(q^3)

            sage: MFSeriesConstructor(prec=3).f_rho_ZZ() == MFSeriesConstructor(prec=3).f_rho(d=1)
            True

            sage: MFSeriesConstructor(group=5, prec=3).f_rho()
            1 + 7/(100*d)*q + 21/(160000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).f_rho(fix_d=True)
            1.000000000000... + 9.925932435107...*q + 2.639039322490...*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).f_rho().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3).f_rho(fix_d=True).parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        (base_ring, coeff_ring, qseries_ring, d) = self.series_data(base_ring, fix_d, d, d_num_prec)
        return self.f_rho_ZZ()(qseries_ring.gen()/d)

    @cached_method
    def f_i(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Return the Fourier expansion of ``f_i``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).f_i(fix_d=True)
            1 - 504*q - 16632*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).f_i_ZZ() == MFSeriesConstructor(prec=3).f_i(d=1)
            True

            sage: MFSeriesConstructor(group=5, prec=3).f_i()
            1 - 13/(40*d)*q - 351/(64000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).f_i(fix_d=True)
            1.000000000000... - 46.08468630585...*q - 110.2741431183...*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).f_i().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3).f_i(fix_d=True).parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        (base_ring, coeff_ring, qseries_ring, d) = self.series_data(base_ring, fix_d, d, d_num_prec)
        return self.f_i_ZZ()(qseries_ring.gen()/d)

    @cached_method
    def f_inf(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Return the Fourier expansion of ``f_inf``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).f_inf(fix_d=True)
            q - 24*q^2 + 252*q^3 + O(q^4)
            sage: MFSeriesConstructor(prec=3).f_inf_ZZ() == MFSeriesConstructor(prec=3).f_inf(d=1)
            True

            sage: MFSeriesConstructor(group=5, prec=3).f_inf()
            q - 9/(200*d)*q^2 + 279/(640000*d^2)*q^3 + O(q^4)
            sage: MFSeriesConstructor(group=5, prec=3).f_inf(fix_d=True)
            0.0000000000000... + 1.000000000000...*q - 6.380956565426...*q^2 + 8.765380606844...*q^3 + O(q^4)

            sage: MFSeriesConstructor(group=5, prec=3).f_inf().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3).f_inf(fix_d=True).parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        (base_ring, coeff_ring, qseries_ring, d) = self.series_data(base_ring, fix_d, d, d_num_prec)
        return d*self.f_inf_ZZ()(qseries_ring.gen()/d)

    @cached_method
    def G_inv(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Return the Fourier expansion of ``G_inv``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(group=4, prec=3).G_inv(fix_d=True)
            1/16777216*q^-1 - 3/2097152 - 955/4194304*q + O(q^2)
            sage: MFSeriesConstructor(group=4, prec=3).G_inv_ZZ() == MFSeriesConstructor(group=4, prec=3).G_inv(d=1)
            True

            sage: MFSeriesConstructor(group=8, prec=3).G_inv()
            d^3*q^-1 - 15*d^2/128 - 15139*d/262144*q + O(q^2)
            sage: MFSeriesConstructor(group=8, prec=3).G_inv(fix_d=True)
            1.648388300301...e-6*q^-1 - 0.00001635263105300... - 0.0006821979994337...*q + O(q^2)

            sage: MFSeriesConstructor(group=8, prec=3).G_inv().parent()
            Laurent Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=8, prec=3).G_inv(fix_d=True).parent()
            Laurent Series Ring in q over Real Field with 53 bits of precision
        """

        (base_ring, coeff_ring, qseries_ring, d) = self.series_data(base_ring, fix_d, d, d_num_prec)
        return d**2*self.G_inv_ZZ()(qseries_ring.gen()/d)

    @cached_method
    def E4(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Return the Fourier expansion of ``E_4``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).E4(fix_d=True)
            1 + 240*q + 2160*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).E4_ZZ() == MFSeriesConstructor(prec=3).E4(d=1)
            True

            sage: MFSeriesConstructor(group=5, prec=3).E4()
            1 + 21/(100*d)*q + 483/(32000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).E4(fix_d=True)
            1.000000000000... + 29.77779730532...*q + 303.4895220864...*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).E4().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3).E4(fix_d=True).parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        (base_ring, coeff_ring, qseries_ring, d) = self.series_data(base_ring, fix_d, d, d_num_prec)
        return self.E4_ZZ()(qseries_ring.gen()/d)

    @cached_method
    def E6(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Return the Fourier expansion of ``E_6``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).E6(fix_d=True)
            1 - 504*q - 16632*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).E6_ZZ() == MFSeriesConstructor(prec=3).E6(d=1)
            True

            sage: MFSeriesConstructor(group=5, prec=3).E6()
            1 - 37/(200*d)*q - 14663/(320000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).E6(fix_d=True)
            1.000000000000... - 26.23282143564...*q - 921.3388948972...*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).E6().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3).E6(fix_d=True).parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        (base_ring, coeff_ring, qseries_ring, d) = self.series_data(base_ring, fix_d, d, d_num_prec)
        return self.E6_ZZ()(qseries_ring.gen()/d)

    @cached_method
    def Delta(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Return the Fourier expansion of ``Delta``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).Delta(fix_d=True)
            q - 24*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).Delta_ZZ() == MFSeriesConstructor(prec=3).Delta(d=1)
            True

            sage: MFSeriesConstructor(group=5, prec=3).Delta()
            71/50*q + 28267/(16000*d)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).Delta(fix_d=True)
            0.0000000000000... + 1.420000000000...*q + 250.5145822707...*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).Delta().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3).Delta(fix_d=True).parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        (base_ring, coeff_ring, qseries_ring, d) = self.series_data(base_ring, fix_d, d, d_num_prec)
        return d*self.Delta_ZZ()(qseries_ring.gen()/d)

    @cached_method
    def E2(self, base_ring = ZZ, fix_d=False, d=None, d_num_prec=ZZ(53)):
        r"""
        Return the Fourier expansion of ``E2``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.series_constructor import MFSeriesConstructor
            sage: MFSeriesConstructor(prec=3).E2(fix_d=True)
            1 - 24*q - 72*q^2 + O(q^3)
            sage: MFSeriesConstructor(prec=3).E2_ZZ() == MFSeriesConstructor(prec=3).E2(d=1)
            True

            sage: MFSeriesConstructor(group=5, prec=3).E2()
            1 - 9/(200*d)*q - 369/(320000*d^2)*q^2 + O(q^3)
            sage: MFSeriesConstructor(group=5, prec=3).E2(fix_d=True)
            1.000000000000... - 6.380956565426...*q - 23.18584547617...*q^2 + O(q^3)

            sage: MFSeriesConstructor(group=5, prec=3).E2().parent()
            Power Series Ring in q over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MFSeriesConstructor(group=5, prec=3).E2(fix_d=True).parent()
            Power Series Ring in q over Real Field with 53 bits of precision
        """

        (base_ring, coeff_ring, qseries_ring, d) = self.series_data(base_ring, fix_d, d, d_num_prec)
        return self.E2_ZZ()(qseries_ring.gen()/d)
