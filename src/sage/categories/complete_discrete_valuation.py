r"""
This module implements the two following categories :

 -  Complete Discrete Valuation Rings (CDVR)

 -  Complete Discrete Valuation Fields (CDVF)
"""
#**************************************************************************
#  Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#**************************************************************************


from sage.misc.abstract_method import abstract_method

from sage.categories.category_singleton import Category_singleton
from discrete_valuation import DiscreteValuationRings, DiscreteValuationFields
#from sage.misc.cachefunc import cached_method

class CompleteDiscreteValuationRings(Category_singleton):
    """
    The category of complete discrete valuation rings
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationRings().super_categories()
            [Category of discrete valuation rings]
        """
        return [DiscreteValuationRings()]

    class ParentMethods:
        def is_cdvr(self):
            """
            Return True if this ring is complete DVR.

            EXAMPLES::

                sage: Zp(7).is_cdvr()
                True

                sage: K.<u> = QQ[[]]
                sage: K.is_cdvr()
                True
            """
            return True

    class ElementMethods:
        def precision_absolute(self):
            """
            Return the absolute precision of this element.

            EXAMPLES:

                sage: R = Zp(7)
                sage: x = R(7); x
                7 + O(7^21)
                sage: x.precision_absolute()
                21
            """
            pass

        def precision_relative(self):
            """
            Return the relative precision of this element.

            EXAMPLES:

                sage: R = Zp(7)
                sage: x = R(7); x
                7 + O(7^21)
                sage: x.precision_relative()
                20
            """
            pass


class CompleteDiscreteValuationFields(Category_singleton):
    """
    The category of complete discrete valuation fields
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationFields().super_categories()
            [Category of discrete valuation fields]
        """
        return [DiscreteValuationFields()]

    class ParentMethods:
        def is_cdvf(self):
            """
            Return True if this ring is a field equipped with a discrete valuation
            for which it is complete.

            EXAMPLES::

                sage: Qp(7).is_cdvf()
                True
            """
            return True

    class ElementMethods:
        def precision_absolute(self):
            """
            Return the absolute precision of this element.

            EXAMPLES:

                sage: K = Qp(7)
                sage: x = K(7); x
                7 + O(7^21)
                sage: x.precision_absolute()
                21
            """
            pass

        def precision_relative(self):
            """
            Return the relative precision of this element.

            EXAMPLES:

                sage: K = Qp(7)
                sage: x = K(7); x
                7 + O(7^21)
                sage: x.precision_relative()
                20
            """
            pass
