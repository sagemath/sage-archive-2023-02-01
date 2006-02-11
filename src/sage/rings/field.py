"""
Abstract base class for fields
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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

from sage.ext.ring import Field, is_Field

def is_PrimeField(R):
    from finite_field import is_FiniteField
    from rational_field import is_RationalField
    if is_RationalField(R):
        return True
    if is_FiniteField(R):
        return R.degree() == 1
    return False
