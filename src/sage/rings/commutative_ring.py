"""
Abstract base class for commutative rings
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

from sage.rings.ring import CommutativeRing

def is_CommutativeRing(R):
    """
    Check to see if ``R`` is a :class:`CommutativeRing`.

    EXAMPLES::

        sage: sage.rings.commutative_ring.is_CommutativeRing(ZZ)
        True
    """
    return isinstance(R, CommutativeRing)
