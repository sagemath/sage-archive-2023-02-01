"""
Algebra order ideal elements
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

from sage.algebras.algebra_order_element import AlgebraOrderElement

class AlgebraOrderIdealElement(AlgebraOrderElement):
    """
    An element of an ideal for an order in an algebra, stored as
    an element of the ambient algebra in which it embeds.
    """
    def __init__(self, O, x, check=True):
        """
        Create the element x of the quaternion order O.
        """
        AlgebraOrderElement.__init__(self, O, x, check=check)
