"""
Finite Fields.
"""

#*****************************************************************************
#       Copyright (C) 2010 David Roe <roed@math.harvard.edu>
#                          William Stein <wstein@gmail.com>
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

from constructor import FiniteField, is_FiniteField, is_PrimeFiniteField
from conway_polynomials import conway_polynomial, exists_conway_polynomial
GF = FiniteField

from element_base import FinitePolyExtElement as FiniteFieldElement # for backward compatibility; is this needed?
from element_base import is_FiniteFieldElement

