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

    class ParentMethods:
        def _matrix_smith_form(self, M, transformation):
            """
            Return the Smith normal form of this matrix.

            INPUT:

            - ``transformation`` -- a boolean (default: True)
              Indicates whether the transformation matrices are returned

            NOTE:

            We recall that a Smith decomposition of a matrix `M`
            defined over a complete discrete valuation ring/field
            is a writing of the form `L*M*R = S` where:

            - `L` and `R` are invertible matrices in the ring of
              integers

            - the only non-vanishing entries of `S` are located on
              the diagonal (through `S` might be not a square matrix)

            - if `d_i` denotes the `(i,i)` entry of `D`, then `d_i`
              divides `d_{i+1}` for all `i`.

            The `d_i`'s are uniquely determined provided that they are
            normalized so that they are all either `0` or a power of the 
            distinguished uniformizer of the base ring.
            Normalized this way, the writing `L*M*R = S` is called the
            Smith normal form of `M`.

            EXAMPLES::

                sage: A = Zp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

                sage: S, L, R = M.smith_form()
                sage: S
                [ ...1     0]
                [    0 ...10]
                sage: L
                [...222222223          ...]
                [...444444444         ...2]
                sage: R
                [         ...1 ...2222222214]
                [            0          ...1]

            If not needed, it is possible to do not compute the
            transformations matrices L and R as follows::

                sage: M.smith_form(transformation=False)
                [ ...1     0]
                [    0 ...10]

            This method works for rectangular matrices as well::

                sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
                sage: S, L, R = M.smith_form()
                sage: S
                [ ...1     0]
                [    0 ...10]
                [    0     0]
                sage: L
                [...222222223          ...          ...]
                [...444444444         ...2          ...]
                [...444444443         ...1         ...1]
                sage: R
                [         ...1 ...2222222214]
                [            0          ...1]

            An error is raised if the precision on the entries is
            not enough to determine the Smith normal form::

                sage: M = matrix(A, 2, 2, [1, 1, 1, 1])
                sage: M.smith_form()
                Traceback (most recent call last):
                ...
                PrecisionError: Not enough precision to compute Smith normal form

            TESTS::

                sage: M = random_matrix(A, 4)
                sage: S, L, R = M.smith_form()
                sage: L*M*R == S
                True

            We check that Smith decomposition works over various rings::

                sage: from sage.rings.padics.precision_error import PrecisionError
                sage: ring1 = ZpCA(5,15)
                sage: ring2 = Zq(5^3,names='a')
                sage: ring3 = Zp(5).extension(x^2-5, names='pi')
                sage: ring4 = PowerSeriesRing(GF(5), name='t')
                sage: for A in [ ring1, ring2, ring3, ring4 ]:
                ....:     for _ in range(10):
                ....:         M = random_matrix(A,4)
                ....:         try:
                ....:             S, L, R = M.smith_form()
                ....:         except PrecisionError:
                ....:             continue
                ....:         if L*M*R != S: raise RuntimeError
            """
            from sage.matrix.matrix_cdv_dense import smith_normal_form
            return smith_normal_form(M, transformation)


    class ElementMethods:
        @abstract_method
        def precision_absolute(self):
            """
            Return the absolute precision of this element.

            EXAMPLES::

                sage: R = Zp(7)
                sage: x = R(7); x
                7 + O(7^21)
                sage: x.precision_absolute()
                21
            """

        @abstract_method
        def precision_relative(self):
            """
            Return the relative precision of this element.

            EXAMPLES::

                sage: R = Zp(7)
                sage: x = R(7); x
                7 + O(7^21)
                sage: x.precision_relative()
                20
            """

        @abstract_method
        def lift_to_maximal_precision(self):
            """
            Lift this element to the maximal precision
            allowed by the parent.

            EXAMPLES::

                sage: R = Zp(7,prec=20)
                sage: x = R(7,5); x
                7 + O(7^5)
                sage: x.lift_to_maximal_precision()
                7 + O(7^21)
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
        def _matrix_smith_form(self, M, transformation):
            """
            Return the Smith normal form of this matrix.

            INPUT:

            - ``transformation`` -- a boolean (default: True)
              Indicates whether the transformation matrices are returned

            NOTE:

            We recall that a Smith decomposition of a matrix `M`
            defined over a complete discrete valuation ring/field
            is a writing of the form `L*M*R = S` where:

            - `L` and `R` are invertible matrices in the ring of
              integers

            - the only non-vanishing entries of `S` are located on
              the diagonal (through `S` might be not a square matrix)

            - if `d_i` denotes the `(i,i)` entry of `D`, then `d_i`
              divides `d_{i+1}` for all `i`.

            The `d_i`'s are uniquely determined provided that they are
            normalized so that they are all either `0` or a power of the 
            distinguished uniformizer of the base ring.
            Normalized this way, the writing `L*M*R = S` is called the
            Smith normal form of `M`.

            EXAMPLES::

                sage: A = Qp(5, prec=10, print_mode="digits")
                sage: M = matrix(A, 2, 2, [2, 7, 1, 6])

                sage: S, L, R = M.smith_form()
                sage: S
                [ ...1     0]
                [    0 ...10]
                sage: L
                [...222222223          ...]
                [...444444444         ...2]
                sage: R
                [         ...1 ...2222222214]
                [            0          ...1]

            If not needed, it is possible to do not compute the
            transformations matrices L and R as follows::

                sage: M.smith_form(transformation=False)
                [ ...1     0]
                [    0 ...10]

            This method works for rectangular matrices as well::

                sage: M = matrix(A, 3, 2, [2, 7, 1, 6, 3, 8])
                sage: S, L, R = M.smith_form()
                sage: S
                [ ...1     0]
                [    0 ...10]
                [    0     0]
                sage: L
                [...222222223          ...          ...]
                [...444444444         ...2          ...]
                [...444444443         ...1         ...1]
                sage: R
                [         ...1 ...2222222214]
                [            0          ...1]

            An error is raised if the precision on the entries is
            not enough to determine the Smith normal form::

                sage: M = matrix(A, 2, 2, [1, 1, 1, 1])
                sage: M.smith_form()
                Traceback (most recent call last):
                ...
                PrecisionError: Not enough precision to compute Smith normal form

            TESTS::

                sage: M = random_matrix(A, 4)
                sage: S, L, R = M.smith_form()
                sage: L*M*R == S
                True

            We check that Smith decomposition works over various rings::

                sage: from sage.rings.padics.precision_error import PrecisionError
                sage: ring1 = Qp(7,10)
                sage: ring2 = Qq(7^2,names='a')
                sage: ring3 = Qp(7).extension(x^3-7, names='pi')
                sage: ring4 = LaurentSeriesRing(GF(7), name='t')
                sage: for A in [ ring1, ring2, ring4 ]:  # ring3 causes troubles (see ticket #23464)
                ....:     for _ in range(10):
                ....:         M = random_matrix(A,4)
                ....:         try:
                ....:             S, L, R = M.smith_form()
                ....:         except PrecisionError:
                ....:             continue
                ....:         if L*M*R != S: raise RuntimeError
            """
            from sage.matrix.matrix_cdv_dense import smith_normal_form
            return smith_normal_form(M, transformation)


    class ElementMethods:
        @abstract_method
        def precision_absolute(self):
            """
            Return the absolute precision of this element.

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(7); x
                7 + O(7^21)
                sage: x.precision_absolute()
                21
            """

        @abstract_method
        def precision_relative(self):
            """
            Return the relative precision of this element.

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(7); x
                7 + O(7^21)
                sage: x.precision_relative()
                20
            """


