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
    e^(1/2*log(2*pi))*e^(n*log(n))*(e^n)^(-1)*n^(1/2) +
    1/12*e^(1/2*log(2*pi))*e^(n*log(n))*(e^n)^(-1)*n^(-1/2) +
    1/288*e^(1/2*log(2*pi))*e^(n*log(n))*(e^n)^(-1)*n^(-3/2) +
    O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))

generates the first 5 summands of Stirling's approximation formula for
factorials.

To construct an asymptotic expression manually, you can use the class
:class:`~sage.rings.asymptotic.asymptotic_ring.AsymptoticRing`. See
:doc:`asymptotic ring <asymptotic_ring>` for more details and a lot of
examples.


**Asymptotic Expansions**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~TransducerGenerators.Identity` | Returns a transducer realizing the i

AUTHORS:

- Daniel Krenn (2015-09-03): initial version and Stirling

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


class AsymptoticExpansionGenerators(object):
    r"""
    A class consisting of constructors for several common asymptotic
    expansions.

    A list of all asymptotic expansions in this database is available via tab
    completion. Type "``asymptotic_expansions.``" and then hit tab to see which
    expansions are available.

    The asymptotic expansions currently in this class include:

    - :meth:`~Stirling`
    """

    @staticmethod
    def Stirling(var, precision=None):
        r"""

        EXAMPLES::

            sage: asymptotic_expansions.Stirling('n', precision=5)
            e^(1/2*log(2*pi))*e^(n*log(n))*(e^n)^(-1)*n^(1/2) +
            1/12*e^(1/2*log(2*pi))*e^(n*log(n))*(e^n)^(-1)*n^(-1/2) +
            1/288*e^(1/2*log(2*pi))*e^(n*log(n))*(e^n)^(-1)*n^(-3/2) +
            O(e^(n*log(n))*(e^n)^(-1)*n^(-5/2))
            sage: _.parent()
        """
        from sage.functions.log import exp
        log_Stirling = AsymptoticExpansionGenerators.log_Stirling(var, precision)
        P = log_Stirling.parent().change_parameter(
            growth_group='(e^(n*log(n)))^ZZ * (e^n)^ZZ * n^ZZ * log(n)^ZZ')
        return exp(P(log_Stirling))


    @staticmethod
    def log_Stirling(var, precision=None):
        r"""
        EXAMPLES::

            sage: asymptotic_expansions.log_Stirling('n', precision=7)
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + 1/12*n^(-1)
            - 1/360*n^(-3) + 1/1260*n^(-5) + O(n^(-7))

        TESTS::

            sage: asymptotic_expansions.log_Stirling('n')
            n*log(n) - n + 1/2*log(n) + 1/2*log(2*pi) + 1/12*n^(-1)
            - 1/360*n^(-3) + 1/1260*n^(-5) - 1/1680*n^(-7) + 1/1188*n^(-9)
            - 691/360360*n^(-11) + 1/156*n^(-13) - 3617/122400*n^(-15)
            + 43867/244188*n^(-17) - 174611/125400*n^(-19) + 77683/5796*n^(-21)
            - 236364091/1506960*n^(-23) + 657931/300*n^(-25)
            - 3392780147/93960*n^(-27) + 1723168255201/2492028*n^(-29)
            - 7709321041217/505920*n^(-31) + O(n^(-33))
            sage: _.parent()
            Asymptotic Ring <n^ZZ * log(n)^ZZ> over Symbolic Ring
        """
        from asymptotic_ring import AsymptoticRing
        from sage.symbolic.ring import SR
        A = AsymptoticRing(growth_group='%s^ZZ * log(%s)^ZZ' % ((var,)*2),
                           coefficient_ring=SR)
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
        if precision >= 4:
            result += log(2*SR('pi')) / 2

        from sage.misc.misc import srange
        from sage.rings.arith import bernoulli
        for k in srange(2, 2*precision - 6, 2):
            result += bernoulli(k) / k / (k-1) / n**(k-1)

        result += (1 / n**(2*precision - 7)).O()
        return result


# Easy access to the asymptotic expansions generators from the command line:
asymptotic_expansions = AsymptoticExpansionGenerators()
