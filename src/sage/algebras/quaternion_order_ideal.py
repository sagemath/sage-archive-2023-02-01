"""
Quaternion ideal

AUTHOR:

- David Kohel (2005-09)
"""

#*****************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module import FreeModule
from sage.algebras.algebra_order_ideal import AlgebraOrderLeftIdeal, AlgebraOrderRightIdeal, AlgebraOrderTwoSidedIdeal
from sage.algebras.quaternion_order import QuaternionOrder_generic
from sage.algebras.quaternion_order_element import QuaternionOrderElement

def QuaternionOrderLeftIdeal(O, gens):
    if not isinstance(O, QuaternionOrder_generic):
        raise TypeError, "Argument O (= %s) must be a quaternion order."%O
    if not isinstance(gens, (list,tuple)):
        raise TypeError, "Argument gens (= %s) must be a list."%gens
    return QuaternionOrderLeftIdeal_generic(O,gens)

def QuaternionOrderRightIdeal(O, gens):
    if not isinstance(O, QuaternionOrder_generic):
        raise TypeError, "Argument O (= %s) must be a quaternion order."%O
    if not isinstance(gens, (list,tuple)):
        raise TypeError, "Argument gens (= %s) must be a list."%gens
    return QuaternionOrderRightIdeal_generic(O, gens)

def QuaternionOrderTwoSidedIdeal(O, gens):
    if not isinstance(O, QuaternionOrder_generic):
        raise TypeError, "Argument O (= %s) must be a quaternion order."%O
    if not isinstance(gens, (list,tuple)):
        raise TypeError, "Argument gens (= %s) must be a list."%gens
    return QuaternionOrderTwoSidedIdeal_generic(O, gens)

class QuaternionOrderLeftIdeal_generic(AlgebraOrderLeftIdeal):
    """
    A left ideal over an order in a quaternion algebra.
    """
    def __init__(self, O, gens):
        """

        """
        AlgebraOrderLeftIdeal.__init__(self, O, gens)

class QuaternionOrderRightIdeal_generic(AlgebraOrderRightIdeal):
    """
    A right ideal over an order in a quaternion algebra.
    """
    def __init__(self, O, gens):
        """

        """
        AlgebraOrderRightIdeal.__init__(self, O, gens)

class QuaternionOrderTwoSidedIdeal_generic(AlgebraOrderTwoSidedIdeal):
    """
    A two-sided ideal over an order in a quaternion algebra.
    """
    def __init__(self, O, gens = []):
        """

        """
        AlgebraOrderTwoSidedIdeal.__init__(self, O, gens)
