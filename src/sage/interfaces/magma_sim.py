"""
Make MAGMA-style commands available in Sage via aliases.
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

from sage.all import *

Factorization = factor
Factorisation = factor
CharacteristicPolynomial = charpoly

RMatrixSpace = MatrixSpace
def MatrixAlgebra(R, d):
    """
    Return the full matrix algebra of degree d over R.
    """
    return MatrixSpace(R, d)

true  = True
false = False

