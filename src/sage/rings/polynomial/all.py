"""
Polynomials
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

# Quotient of polynomial ring
from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing, is_PolynomialQuotientRing
from sage.rings.polynomial.polynomial_quotient_ring_element import PolynomialQuotientRingElement

# Univariate Polynomial Rings
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import polygen, polygens, is_PolynomialRing
from sage.rings.polynomial.polynomial_element import is_Polynomial, Polynomial

# Multivariate Polynomial Rings
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.polynomial.term_order import TermOrder
from sage.rings.polynomial.multi_polynomial_element import degree_lowest_rational_function, is_MPolynomial

# Generic convolution
from sage.rings.polynomial.convolution import convolution

# Boolean Polynomial Rings
from sage.rings.polynomial.polynomial_ring_constructor import BooleanPolynomialRing_constructor as BooleanPolynomialRing

from sage.rings.polynomial.multi_polynomial_ideal import is_MPolynomialIdeal

# Laurent Polynomial Rings
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing, is_LaurentPolynomialRing

# Infinite Polynomial Rings
from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing

# Evaluation of cyclotomic polynomials
from sage.rings.polynomial.cyclotomic import cyclotomic_value
