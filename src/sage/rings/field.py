"""
Abstract base class for fields
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

from sage.rings.ring import Field

def is_PrimeField(R):
    """
    Determine if R is a field that is equal to its own prime subfield.

    INPUT:

    - R - a ring or field

    OUTPUT:

    - True - if R is `\QQ` or a finite field `GF(p)` for p prime.
    - False - otherwise

    EXAMPLES::

        sage: sage.rings.field.is_PrimeField(QQ)
        True

    ::

        sage: sage.rings.field.is_PrimeField(GF(7))
        True

    ::

        sage: sage.rings.field.is_PrimeField(GF(7^2,'t'))
        False
    """
    from finite_rings.constructor import is_FiniteField
    from rational_field import is_RationalField

    if is_RationalField(R):
        return True
    if is_FiniteField(R):
        return R.degree() == 1
    return False
