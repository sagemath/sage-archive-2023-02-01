r"""
Complete Discrete Valuation Rings (CDVR) and Fields (CDVF)
"""
from __future__ import absolute_import
#**************************************************************************
#  Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#**************************************************************************


from sage.misc.abstract_method import abstract_method

from sage.categories.category_singleton import Category_singleton
from .discrete_valuation import DiscreteValuationRings, DiscreteValuationFields
#from sage.misc.cachefunc import cached_method


class CompleteDiscreteValuationRings(Category_singleton):
    """
    The category of complete discrete valuation rings

    EXAMPLES::

        sage: Zp(7) in CompleteDiscreteValuationRings()
        True
        sage: QQ in CompleteDiscreteValuationRings()
        False
        sage: QQ[['u']] in CompleteDiscreteValuationRings()
        True
        sage: Qp(7) in CompleteDiscreteValuationRings()
        False
        sage: TestSuite(CompleteDiscreteValuationRings()).run()
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationRings().super_categories()
            [Category of discrete valuation rings]
        """
        return [DiscreteValuationRings()]

    class ElementMethods:
        @abstract_method
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: R = Zp(7)
                sage: x = R(7); x
                7 + O(7^21)
                sage: x.valuation()
                1
            """

        def denominator(self):
            """
            Return the denominator of this element normalized
            as a power of the uniformizer

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(1/21)
                sage: x.denominator()
                7 + O(7^21)

                sage: x = K(7)
                sage: x.denominator()
                1 + O(7^20)

            Note that the denominator lives in the ring of integers::

                sage: x.denominator().parent()
                7-adic Ring with capped relative precision 20

            An error is raised when the input is indistinguishable from 0::

                sage: x = K(0,5); x
                O(7^5)
                sage: x.denominator()
                Traceback (most recent call last):
                ...
                ValueError: Cannot determine the denominator of an element indistinguishable from 0
            """
            return self.parent()(1)

        def _matrix_determinant(self,M):
            """
            Return the determinant of this matrix

            ALGORITHM:

            We row-echenolize the matrix by always choosing the
            pivot of smallest valuation and allowing permutations
            of columns.

            We then compute separatedly the value of the determinant
            (as the product of the diagonal entries of the row-echelon
            form) and a bound on the precision on it.

            EXAMPLES::

                sage: R = Zp(5, 10)
                sage: M = matrix(R, 2, 2, [1, 6, 2, 7])
                sage: M.determinant()  # indirect doctest
                4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)

                sage: (5*M).determinant()  # indirect doctest
                4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12)

            TESTS:

            We check the stability of our algorithm::

                sage: R = Zp(5,10)
                sage: M = random_matrix(R,3) * diagonal_matrix([1,25,125]) * random_matrix(R,3)
                sage: d = M.determinant()
                sage: d.precision_absolute() >= 12
                True

                sage: for dim in range(3,10):
                ....:     M = matrix(dim, dim, [ R(1) for _ in range(dim^2) ])
                ....:     print M.determinant()
                O(5^20)
                O(5^30)
                O(5^40)
                O(5^50)
                O(5^60)
                O(5^70)
                O(5^80)
            """
            from sage.matrix.matrix_cdv_dense import determinant
            return determinant(M)


        @abstract_method
        def lift_to_precision(self, absprec=None):
            """
            Return another element of the same parent with absolute precision
            at least ``absprec``, congruent to this element modulo the
            precision of this element.

            INPUT:

            - ``absprec`` -- an integer or ``None`` (default: ``None``), the
              absolute precision of the result. If ``None``, lifts to the maximum
              precision allowed.

            .. NOTE::

                If setting ``absprec`` that high would violate the precision cap,
                raises a precision error.  Note that the new digits will not
                necessarily be zero.

            EXAMPLES::

                sage: R = ZpCA(17)
                sage: R(-1,2).lift_to_precision(10)
                16 + 16*17 + O(17^10)
                sage: R(1,15).lift_to_precision(10)
                1 + O(17^15)
                sage: R(1,15).lift_to_precision(30)
                Traceback (most recent call last):
                ...
                PrecisionError: Precision higher than allowed by the precision cap.
                sage: R(-1,2).lift_to_precision().precision_absolute() == R.precision_cap()
                True

                sage: R = Zp(5); c = R(17,3); c.lift_to_precision(8)
                2 + 3*5 + O(5^8)
                sage: c.lift_to_precision().precision_relative() == R.precision_cap()
                True

            """

class CompleteDiscreteValuationFields(Category_singleton):
    """
    The category of complete discrete valuation fields

    EXAMPLES::

        sage: Zp(7) in CompleteDiscreteValuationFields()
        False
        sage: QQ in CompleteDiscreteValuationFields()
        False
        sage: LaurentSeriesRing(QQ,'u') in CompleteDiscreteValuationFields()
        True
        sage: Qp(7) in CompleteDiscreteValuationFields()
        True
        sage: TestSuite(CompleteDiscreteValuationFields()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationFields().super_categories()
            [Category of discrete valuation fields]
        """
        return [DiscreteValuationFields()]

    class ParentMethods:
        def _matrix_determinant(self,M):
            """
            Return the determinant of this matrix

            ALGORITHM:

            We row-echenolize the matrix by always choosing the
            pivot of smallest valuation and allowing permutations
            of columns.

            We then compute separatedly the value of the determinant
            (as the product of the diagonal entries of the row-echelon
            form) and a bound on the precision on it.

            EXAMPLES::

                sage: R = Qp(5, 10)
                sage: M = matrix(R, 2, 2, [1, 6, 2, 7])
                sage: M.determinant()  # indirect doctest
                4*5 + 4*5^2 + 4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + O(5^10)

                sage: (5*M).determinant()  # indirect doctest
                4*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 4*5^7 + 4*5^8 + 4*5^9 + 4*5^10 + 4*5^11 + O(5^12)

            TESTS:

            We check the stability of our algorithm::

                sage: for dim in range(3,10):
                ....:     M = matrix(dim, dim, [ R(1) for _ in range(dim^2) ])
                ....:     print M.determinant()
                O(5^20)
                O(5^30)
                O(5^40)
                O(5^50)
                O(5^60)
                O(5^70)
                O(5^80)
            """
            from sage.matrix.matrix_cdv_dense import determinant
            return determinant(M)

    class ElementMethods:
        @abstract_method
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(7); x
                7 + O(7^21)
                sage: x.valuation()
                1
            """

        def denominator(self):
            """
            Return the denominator of this element normalized
            as a power of the uniformizer

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(1/21)
                sage: x.denominator()
                7 + O(7^21)

                sage: x = K(7)
                sage: x.denominator()
                1 + O(7^20)

            Note that the denominator lives in the ring of integers::

                sage: x.denominator().parent()
                7-adic Ring with capped relative precision 20

            An error is raised when the input is indistinguishable from 0::

                sage: x = K(0,5); x
                O(7^5)
                sage: x.denominator()
                Traceback (most recent call last):
                ...
                ValueError: Cannot determine the denominator of an element indistinguishable from 0
            """
            if self == 0:
                raise ValueError("Cannot determine the denominator of an element indistinguishable from 0")
            val = self.valuation()
            R = self.parent().integer_ring()
            if val >= 0:
                return R(1)
            else:
                return R(1) << (-val)
