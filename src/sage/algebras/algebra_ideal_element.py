"""
Algebra ideal elements

AUTHOR: David Kohel, 2005-09
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

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.algebras.algebra_element import AlgebraElement

class AlgebraIdealElement(AlgebraElement):
    """
    An element of an ideal in an algebra, stored as an element
    of the algebra in which it embedds.
    """
    def __init__(self, A, x):
        """
        Create the element x of the quaternion order A.
        """
        AlgebraElement.__init__(self, A)
        self.__algebra_element = x

