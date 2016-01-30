r"""
Common Asymptotic Expansions

.. WARNING::

    As this code is experimental, a warning is thrown when an
    asymptotic ring (or an associated structure) is created for the
    first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: AsymptoticRing(growth_group='z^ZZ * log(z)^QQ', coefficient_ring=ZZ)
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/17601 for details.
        Asymptotic Ring <z^ZZ * log(z)^QQ> over Integer Ring


Asymptotic expansions in SageMath can be built through the
``asymptotic_expansions`` object. It contains generators for common
asymptotic expressions. For example,
::

    sage: asymptotic_expansions.Stirling('n', precision=5)
    sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(1/2) +
    1/12*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-1/2) +
    1/288*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-3/2) +
    O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))

generates the first 5 summands of Stirling's approximation formula for
factorials.

To construct an asymptotic expression manually, you can use the class
:class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticRing`. See
:doc:`asymptotic ring <asymptotic_ring>` for more details and a lot of
examples.


**Asymptotic Expansions**

.. list-table::
   :class: contentstable
   :widths: 4 12
   :header-rows: 0

   * - :meth:`~AsymptoticExpansionGenerators.Stirling`
     - Stirling's approximation formula for factorials

   * - :meth:`~AsymptoticExpansionGenerators.log_Stirling`
     - the logarithm of Stirling's approximation formula for factorials

   * - :meth:`~AsymptoticExpansionGenerators.Binomial_kn_over_n`
     - an asymptotic expansion of the binomial coefficient


AUTHORS:

- Daniel Krenn (2015)


ACKNOWLEDGEMENT:

- Benjamin Hackl, Clemens Heuberger and Daniel Krenn are supported by the
  Austrian Science Fund (FWF): P 24644-N26.


Classes and Methods
===================
"""

#*****************************************************************************
# Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject

class AsymptoticExpansionGenerators(SageObject):
    r"""
    A collection of constructors for several common asymptotic expansions.

    A list of all asymptotic expansions in this database is available via tab
    completion. Type "``asymptotic_expansions.``" and then hit tab to see which
    expansions are available.

    The asymptotic expansions currently in this class include:

    - :meth:`~Stirling`
    - :meth:`~log_Stirling`
    - :meth:`~Binomial_kn_over_n`
    """

    @staticmethod
    def Stirling(var, precision=None, skip_constant_factor=False):
        r"""
        Return Stirling's approximation formula for factorials.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``precision`` -- (default: ``None``) an integer `\ge 3`. If ``None``, then
          the default precision of the asymptotic ring is used.

        - ``skip_constant_factor`` -- (default: ``False``) a
          boolean. If set, then the constant factor `\sqrt{2\pi}` is left out.
          As a consequence, the coefficient ring of the output changes
          from ``Symbolic Constants Subring`` (if ``False``) to
          ``Rational Field`` (if ``True``).

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: asymptotic_expansions.Stirling('n', precision=5)
            sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(1/2) +
            1/12*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-1/2) +
            1/288*sqrt(2)*sqrt(pi)*e^(n*log(n))*(e^n)^(-1)*n^(-3/2) +
            O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))
            sage: _.parent()
            Asymptotic Ring <(e^(n*log(n)))^QQ * (e^n)^QQ * n^QQ * log(n)^QQ>
            over Symbolic Constants Subring

        TESTS::

            sage: expansion = asymptotic_expansions.Stirling('n', precision=5)
            sage: n = expansion.parent().gen()
            sage: expansion.compare_with_values(n, lambda x: x.factorial(), [5, 10, 20]) # rel tol 1e-6
            [(5, 0.00675841118?), (10, 0.0067589306?), (20, 0.006744925?)]
            sage: asymptotic_expansions.Stirling('n', precision=5,
            ....:                                skip_constant_factor=True)
            e^(n*log(n))*(e^n)^(-1)*n^(1/2) +
            1/12*e^(n*log(n))*(e^n)^(-1)*n^(-1/2) +
            1/288*e^(n*log(n))*(e^n)^(-1)*n^(-3/2) +
            O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))
            sage: _.parent()
            Asymptotic Ring <(e^(n*log(n)))^QQ * (e^n)^QQ * n^QQ * log(n)^QQ>
            over Rational Field
            sage: asymptotic_expansions.Stirling('m', precision=4)
            sqrt(2)*sqrt(pi)*e^(m*log(m))*(e^m)^(-1)*m^(1/2) +
            O(e^(m*log(m))*(e^m)^(-1)*m^(-1/2))
            sage: asymptotic_expansions.Stirling('m', precision=3)
            O(e^(m*log(m))*(e^m)^(-1)*m^(1/2))
            sage: asymptotic_expansions.Stirling('m', precision=2)
            Traceback (most recent call last):
            ...
            ValueError: precision must be at least 3
        """
        if precision < 3:
            raise ValueError("precision must be at least 3")
        log_Stirling = AsymptoticExpansionGenerators.log_Stirling(
            var, precision=precision, skip_constant_summand=True)

        P = log_Stirling.parent().change_parameter(
            growth_group='(e^({n}*log({n})))^QQ * (e^{n})^QQ * {n}^QQ * log({n})^QQ'.format(n=var))
        from sage.functions.log import exp
        result = exp(P(log_Stirling))

        if not skip_constant_factor:
            from sage.symbolic.ring import SR
            SCR = SR.subring(no_variables=True)
            result *= (2*SCR('pi')).sqrt()

        return result


    @staticmethod
    def log_Stirling(var, precision=None, skip_constant_summand=False):
        r"""
        Return the logarithm of Stirling's approximation formula
        for factorials.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.

        - ``skip_constant_summand`` -- (default: ``False``) a
          boolean. If set, then the constant summand `\log(2\pi)/2` is left out.
          As a consequence, the coefficient ring of the output changes
          from ``Symbolic Constants Subring`` (if ``False``) to
          ``Rational Field`` (if ``True``).

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: asymptotic_expansions.log_Stirling('n', precision=7)
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + 1/12*n^(-1)
            - 1/360*n^(-3) + 1/1260*n^(-5) + O(n^(-7))

        TESTS::

            sage: expansion = asymptotic_expansions.log_Stirling('n', precision=7)
            sage: n = expansion.parent().gen()
            sage: expansion.compare_with_values(n, lambda x: x.factorial().log(), [5, 10, 20])  # rel tol 1e-6
            [(5, 0.000564287?), (10, 0.0005870?), (20, 0.0006?)]
            sage: asymptotic_expansions.log_Stirling('n')
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + 1/12*n^(-1)
            - 1/360*n^(-3) + 1/1260*n^(-5) - 1/1680*n^(-7) + 1/1188*n^(-9)
            - 691/360360*n^(-11) + 1/156*n^(-13) - 3617/122400*n^(-15)
            + 43867/244188*n^(-17) - 174611/125400*n^(-19) + 77683/5796*n^(-21)
            - 236364091/1506960*n^(-23) + 657931/300*n^(-25)
            - 3392780147/93960*n^(-27) + 1723168255201/2492028*n^(-29)
            - 7709321041217/505920*n^(-31) + O(n^(-33))
            sage: _.parent()
            Asymptotic Ring <n^ZZ * log(n)^ZZ> over Symbolic Constants Subring

        ::

            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=7, skip_constant_summand=True)
            n*log(n) - n + 1/2*log(n) + 1/12*n^(-1) - 1/360*n^(-3) +
            1/1260*n^(-5) + O(n^(-7))
            sage: _.parent()
            Asymptotic Ring <n^ZZ * log(n)^ZZ> over Rational Field
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=0)
            O(n*log(n))
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=1)
            n*log(n) + O(n)
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=2)
            n*log(n) - n + O(log(n))
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=3)
            n*log(n) - n + 1/2*log(n) + O(1)
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=4)
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + O(n^(-1))
            sage: asymptotic_expansions.log_Stirling(
            ....:     'n', precision=5)
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + 1/12*n^(-1)
            + O(n^(-3))
            sage: asymptotic_expansions.log_Stirling(
            ....:     'm', precision=7, skip_constant_summand=True)
            m*log(m) - m + 1/2*log(m) + 1/12*m^(-1) - 1/360*m^(-3) +
            1/1260*m^(-5) + O(m^(-7))
        """
        if not skip_constant_summand:
            from sage.symbolic.ring import SR
            coefficient_ring = SR.subring(no_variables=True)
        else:
            from sage.rings.rational_field import QQ
            coefficient_ring = QQ

        from asymptotic_ring import AsymptoticRing
        A = AsymptoticRing(growth_group='{n}^ZZ * log({n})^ZZ'.format(n=var),
                           coefficient_ring=coefficient_ring)
        n = A.gen()

        if precision is None:
            precision = A.default_prec

        from sage.functions.log import log
        result = A.zero()
        if precision >= 1:
            result += n * log(n)
        if precision >= 2:
            result += -n
        if precision >= 3:
            result += log(n) / 2
        if precision >= 4 and not skip_constant_summand:
            result += log(2*coefficient_ring('pi')) / 2

        result += AsymptoticExpansionGenerators._log_StirlingNegativePowers_(
            var, precision - 4)

        if precision < 1:
            result += (n * log(n)).O()
        elif precision == 1:
            result += n.O()
        elif precision == 2:
            result += log(n).O()
        elif precision == 3:
            result += A(1).O()

        return result


    @staticmethod
    def _log_StirlingNegativePowers_(var, precision):
        r"""
        Helper function to calculate the logarithm of Stirling's approximation
        formula from the negative powers of ``var`` on, i.e., it skips the
        summands `n \log n - n + (\log n)/2 + \log(2\pi)/2`.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``precision`` -- an integer specifying the number of exact summands.
          If this is negative, then the result is `0`.

        OUTPUT:

        An asymptotic expansion.

        TESTS::

            sage: asymptotic_expansions._log_StirlingNegativePowers_(
            ....:     'm', precision=-1)
            0
            sage: asymptotic_expansions._log_StirlingNegativePowers_(
            ....:     'm', precision=0)
            O(m^(-1))
            sage: asymptotic_expansions._log_StirlingNegativePowers_(
            ....:     'm', precision=3)
            1/12*m^(-1) - 1/360*m^(-3) + 1/1260*m^(-5) + O(m^(-7))
            sage: _.parent()
            Asymptotic Ring <m^ZZ> over Rational Field
        """
        from asymptotic_ring import AsymptoticRing
        from sage.rings.rational_field import QQ

        A = AsymptoticRing(growth_group='{n}^ZZ'.format(n=var),
                           coefficient_ring=QQ)
        if precision < 0:
            return A.zero()
        n = A.gen()

        from sage.arith.all import bernoulli
        from sage.misc.misc import srange

        result = sum((bernoulli(k) / k / (k-1) / n**(k-1)
                      for k in srange(2, 2*precision + 2, 2)),
                     A.zero())
        return result + (1 / n**(2*precision + 1)).O()


    @staticmethod
    def Binomial_kn_over_n(var, k, precision=None, skip_constant_factor=False):
        r"""
        Return the asymptotic expansion of the binomial coefficient
        `kn` choose `n`.

        INPUT:

        - ``var`` -- a string for the variable name.

        - ``k`` -- a number or symbolic constant.

        - ``precision`` -- (default: ``None``) an integer. If ``None``, then
          the default precision of the asymptotic ring is used.

        - ``skip_constant_factor`` -- (default: ``False``) a
          boolean. If set, then the constant factor `\sqrt{k/(2\pi(k-1))}`
          is left out.
          As a consequence, the coefficient ring of the output changes
          from ``Symbolic Constants Subring`` (if ``False``) to
          ``Rational Field`` (if ``True``).

        OUTPUT:

        An asymptotic expansion.

        EXAMPLES::

            sage: asymptotic_expansions.Binomial_kn_over_n('n', k=2, precision=3)
            1/sqrt(pi)*4^n*n^(-1/2)
            - 1/8/sqrt(pi)*4^n*n^(-3/2)
            + 1/128/sqrt(pi)*4^n*n^(-5/2)
            + O(4^n*n^(-7/2))
            sage: _.parent()
            Asymptotic Ring <QQ^n * n^QQ> over Symbolic Constants Subring

        ::

            sage: asymptotic_expansions.Binomial_kn_over_n('n', k=3, precision=3)
            1/2*sqrt(3)/sqrt(pi)*(27/4)^n*n^(-1/2)
            - 7/144*sqrt(3)/sqrt(pi)*(27/4)^n*n^(-3/2)
            + 49/20736*sqrt(3)/sqrt(pi)*(27/4)^n*n^(-5/2)
            + O((27/4)^n*n^(-7/2))

        ::

            sage: asymptotic_expansions.Binomial_kn_over_n('n', k=7/5, precision=3)
            1/2*sqrt(7)/sqrt(pi)*(7/10*7^(2/5)*2^(3/5))^n*n^(-1/2)
            - 13/112*sqrt(7)/sqrt(pi)*(7/10*7^(2/5)*2^(3/5))^n*n^(-3/2)
            + 169/12544*sqrt(7)/sqrt(pi)*(7/10*7^(2/5)*2^(3/5))^n*n^(-5/2)
            + O((7/10*7^(2/5)*2^(3/5))^n*n^(-7/2))
            sage: _.parent()
            Asymptotic Ring <(Symbolic Constants Subring)^n * n^QQ>
            over Symbolic Constants Subring

        TESTS::

            sage: expansion = asymptotic_expansions.Binomial_kn_over_n('n', k=7/5, precision=3)
            sage: n = expansion.parent().gen()
            sage: expansion.compare_with_values(n, lambda x: binomial(7/5*x, x), [5, 10, 20]) # rel tol 1e-6
            [(5, -0.0287383845047?), (10, -0.030845971026?), (20, -0.03162833549?)]
            sage: asymptotic_expansions.Binomial_kn_over_n(
            ....:     'n', k=5, precision=3, skip_constant_factor=True)
            (3125/256)^n*n^(-1/2)
            - 7/80*(3125/256)^n*n^(-3/2)
            + 49/12800*(3125/256)^n*n^(-5/2)
            + O((3125/256)^n*n^(-7/2))
            sage: _.parent()
            Asymptotic Ring <QQ^n * n^QQ> over Rational Field
            sage: asymptotic_expansions.Binomial_kn_over_n(
            ....:     'n', k=4, precision=1, skip_constant_factor=True)
            (256/27)^n*n^(-1/2) + O((256/27)^n*n^(-3/2))

        ::

            sage: S = asymptotic_expansions.Stirling('n', precision=5)
            sage: n = S.parent().gen()
            sage: all(  # long time
            ....:     SR(asymptotic_expansions.Binomial_kn_over_n(
            ....:         'n', k=k, precision=3)).canonicalize_radical() ==
            ....:     SR(S.subs(n=k*n) / (S.subs(n=(k-1)*n) * S)).canonicalize_radical()
            ....:     for k in [2, 3, 4])
            True
        """
        from sage.symbolic.ring import SR
        SCR = SR.subring(no_variables=True)
        try:
            SCR.coerce(k)
        except TypeError as e:
            from misc import combine_exceptions
            raise combine_exceptions(
                TypeError('Cannot use k={}.'.format(k)), e)

        S = AsymptoticExpansionGenerators._log_StirlingNegativePowers_(
                var, precision=max(precision - 2,0))
        n = S.parent().gen()
        result = (S.subs(n=k*n) - S.subs(n=(k-1)*n) - S).exp()

        from sage.rings.rational_field import QQ

        P = S.parent().change_parameter(
                growth_group='QQ^{n} * {n}^QQ'.format(n=var),
                coefficient_ring=QQ)
        n = P.gen()

        b = k**k / (k-1)**(k-1)
        if b.parent() is SR:
            b = SCR(b).canonicalize_radical()
        result *= n.rpow(b)
        result *= n**(-QQ(1)/QQ(2))
        if not skip_constant_factor:
            result *= (k/((k-1)*2*SCR('pi'))).sqrt()

        return result


# Easy access to the asymptotic expansions generators from the command line:
asymptotic_expansions = AsymptoticExpansionGenerators()
r"""
A collection of several common asymptotic expansions.

This is an instance of :class:`AsymptoticExpansionGenerators` whose documentation
provides more details.
"""
