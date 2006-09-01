"""
Algebras
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Algebra base classes
from algebra import Algebra, is_Algebra
from commutative_algebra import CommutativeAlgebra, is_CommutativeAlgebra

# Ring element base classes
from algebra_element import AlgebraElement, is_AlgebraElement
from commutative_algebra_element import CommutativeAlgebraElement, is_CommutativeAlgebraElement

from free_algebra import FreeAlgebra
from free_algebra_quotient import FreeAlgebraQuotient
from quaternion_algebra import (QuaternionAlgebra, QuaternionAlgebraWithInnerProduct,
     QuaternionAlgebraWithGramMatrix, QuaternionAlgebraWithDiscriminants,
     hilbert_symbol, fundamental_discriminant)
from quaternion_order import QuaternionOrderWithBasis, QuaternionDefiningOrder
from quaternion_order_ideal import QuaternionOrderLeftIdeal, QuaternionOrderRightIdeal, QuaternionOrderTwoSidedIdeal


def is_R_algebra(Q, R):
    # TODO: do something nontrivial when morphisms are defined.
    return True




