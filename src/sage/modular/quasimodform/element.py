"""
Elements of quasimodular forms rings

AUTHORS:

- DAVID AYOTTE (2021-03-18): initial version
"""
# ****************************************************************************
#       Copyright (C) 2021 David Ayotte
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.modular.modform.eis_series import eisenstein_series_qexp
from sage.modular.modform.element import GradedModularFormElement, ModularFormElement

from sage.rings.all import Integer, QQ, ZZ

from sage.structure.element import Element
from sage.structure.richcmp import richcmp, op_NE, op_EQ

from sage.rings.polynomial.polynomial_element import Polynomial

class QuasiModularFormsElement(Element):
    r"""
    A quasimodular forms ring element

    .. TODO::

    Move this class somewhere else?
    """
    def __init__(self, parent, polynomial):
        r"""
        An element of a graded ring of quasimodular form.

        INPUTS:

        - ``parent`` - QuasiModularForms
        - ``polynomial`` - a polynomial `f_0 + f_1 E_2 + ... + f_n E_2^n` where each `f_i`
          are modular forms or base ring elements and `E_2` correspond to the weight 2 Eisenstein series.

        OUTPUT:

        - ``QuasiModularFormsElement``

        EXAMPLES::

            sage: QM = QuasiModularForms()
            sage: QM.gens()
            [1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6),
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6),
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)]
            sage: QM.0 + QM.1
            2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
        """
        if not isinstance(polynomial, Polynomial):
            raise TypeError("'polynomial' argument should be of type 'Polynomial'")
        for f in polynomial.coefficients():
            if not isinstance(f, GradedModularFormElement):
                raise ValueError("at least one coefficient is not a 'GradedModularFormElement'")
        self._polynomial = polynomial
        self._coefficients = polynomial.coefficients(sparse=False)
        self.__base_ring = parent.base_ring()
        Element.__init__(self, parent)

    def q_expansion(self, prec=6):
        r"""
        Computes the `q`-expansion of self to precision `prec`.

        EXAMPLES:::

            sage: QM = QuasiModularForms()
            sage: E2 = QM.0
            sage: E2.q_expansion()
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: E2.q_expansion(prec=10)
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 - 288*q^6 - 192*q^7 - 360*q^8 - 312*q^9 + O(q^10)
        """
        E2 = eisenstein_series_qexp(2, prec=prec, K=self.__base_ring, normalization='constant') #normalization -> to force integer coefficients
        return sum(f.q_expansion(prec=prec)*E2**idx for idx, f in enumerate(self._coefficients))

    def _repr_(self):
        r"""
        String representation of self.

        TESTS::

            sage: QM = QuasiModularForms()
            sage: QM.0
            1 - 24*q - 72*q^2 - 96*q^3 - 168*q^4 - 144*q^5 + O(q^6)
            sage: QM.1
            1 + 240*q + 2160*q^2 + 6720*q^3 + 17520*q^4 + 30240*q^5 + O(q^6)
            sage: QM.2
            1 - 504*q - 16632*q^2 - 122976*q^3 - 532728*q^4 - 1575504*q^5 + O(q^6)
        """
        return "%s" % (self.q_expansion())

    def _richcmp_(self, other, op):
        r"""
        Compare self with other.

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: QM.0 == QM.1
            False
            sage: QM.0 == QM.0
            True
            sage: QM.0 != QM.1
            True
            sage: QM.0 != QM.0
            False
            sage: QM.0 < QM.1
            Traceback (most recent call last):
            ...
            TypeError: invalid comparison between quasimodular forms ring elements
        """
        if op != op_EQ and op != op_NE:
            raise TypeError('invalid comparison between quasimodular forms ring elements')
        return richcmp(self._polynomial, other._polynomial, op)

    def _add_(self, other):
        r"""
        Addition of two ``QuasiModularFormElement``.

        INPUT:

        - ``other`` - ``QuasiModularFormElement``

        OUTPUT: a ``QuasiModularFormElement``

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: QM.0 + QM.1
            2 + 216*q + 2088*q^2 + 6624*q^3 + 17352*q^4 + 30096*q^5 + O(q^6)
            sage: QM.0 + (QM.1 + QM.2) == (QM.0 + QM.1) + QM.2
            True
        """
        return self.__class__(self.parent(), self._polynomial + other._polynomial)

    def __neg__(self):
        r"""
        The negation of ``self```

        TESTS::

            sage: -QuasiModularForms(1).0
            -1 + 24*q + 72*q^2 + 96*q^3 + 168*q^4 + 144*q^5 + O(q^6)
            sage: QuasiModularForms(1).0 - QuasiModularForms(1).0
            0
        """
        return self.__class__(self.parent(), -self._polynomial)

    def _mul_(self, other):
        r"""
        The multiplication of two ``QuasiModularFormElement``

        INPUT:

        - ``other`` - ``QuasiModularFormElement``

        OUTPUT: a ``QuasiModularFormElement``

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: QM.0 * QM.1
            1 + 216*q - 3672*q^2 - 62496*q^3 - 322488*q^4 - 1121904*q^5 + O(q^6)
            sage: (QM.0 * QM.1) * QM.2 == QM.0 * (QM.1 * QM.2)
            True
        """
        return self.__class__(self.parent(), self._polynomial * other._polynomial)

    def _lmul_(self, c):
        r"""
        The left action of the base ring on self.

        INPUT:

        - ``other`` - ``QuasiModularFormElement``

        OUTPUT: a ``QuasiModularFormElement``

        TESTS::

            sage: QM = QuasiModularForms(1)
            sage: (1/2) * QM.0
            1/2 - 12*q - 36*q^2 - 48*q^3 - 84*q^4 - 72*q^5 + O(q^6)
            sage: QM.0 * (3/2)
            3/2 - 36*q - 108*q^2 - 144*q^3 - 252*q^4 - 216*q^5 + O(q^6)
        """
        return self.__class__(self.parent(), c * self._polynomial)
