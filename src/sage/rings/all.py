"""
Rings
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

# Ring base classes
from ring import Ring, is_Ring
from commutative_ring import CommutativeRing, is_CommutativeRing
from integral_domain import IntegralDomain, is_IntegralDomain
from dedekind_domain import DedekindDomain, is_DedekindDomain
from principal_ideal_domain import PrincipalIdealDomain, is_PrincipalIdealDomain
from euclidean_domain import EuclideanDomain, is_EuclideanDomain
from field import Field, is_Field, is_PrimeField

from commutative_algebra_element import CommutativeAlgebraElement, is_CommutativeAlgebraElement

# Ring element base classes
from ring_element import RingElement, is_RingElement
from commutative_ring_element import CommutativeRingElement, is_CommutativeRingElement
from integral_domain_element import IntegralDomainElement, is_IntegralDomainElement
from dedekind_domain_element import DedekindDomainElement, is_DedekindDomainElement
from principal_ideal_domain_element import PrincipalIdealDomainElement, is_PrincipalIdealDomainElement
from euclidean_domain_element import EuclideanDomainElement, is_EuclideanDomainElement
from field_element import FieldElement, is_FieldElement


# Ideals
from ideal import Ideal, is_Ideal

# Quotient
from quotient_ring import QuotientRing

# Infinities
from infinity import infinity, Infinity, is_Infinite, InfinityRing, unsigned_infinity, UnsignedInfinityRing
from extended_rational_field import ExtendedRationalField
from extended_integer_ring import ExtendedIntegerRing

# Rational integers.
from integer_ring import IntegerRing, ZZ, crt_basis
from integer import Integer, is_Integer

# Rational numbers
from rational_field import RationalField, QQ, is_RationalField
from rational import Rational
Rationals = RationalField

# Integers modulo n.
from integer_mod_ring import IntegerModRing, Zmod, is_IntegerModRing
from integer_mod import IntegerMod, Mod, mod, is_IntegerMod
Integers = IntegerModRing

# Finite fields
from finite_field import (FiniteField, is_FiniteField, is_PrimeFiniteField,
                          conway_polynomial, exists_conway_polynomial)
GF = FiniteField

from finite_field_element import FiniteFieldElement, is_FiniteFieldElement

# Number field
from number_field.all import *

# p-adic field
from padics.all import *
from padics.padic_printing import _printer_defaults as padic_printing

# Real numbers
from real_mpfr import (RealField_constructor as RealField, is_RealField, is_RealNumber, RR,
                       create_RealNumber as RealNumber)   # this is used by the preparser to wrap real literals -- very important.
Reals = RealField

from real_double import RealDoubleField, RDF, RealDoubleElement, is_RealDoubleElement

from real_lazy import RealLazyField, RLF, ComplexLazyField, CLF

# Quad double
from real_rqdf import RealQuadDoubleField, RQDF, QuadDoubleElement

# Polynomial Rings and Polynomial Quotient Rings
from polynomial.all import *

# Algebraic numbers
from qqbar import (AlgebraicRealField, is_AlgebraicRealField, AA,
                   AlgebraicReal, is_AlgebraicReal,
                   AlgebraicField, is_AlgebraicField, QQbar,
                   AlgebraicNumber, is_AlgebraicNumber,
                   number_field_elements_from_algebraics)

# Intervals
from real_mpfi import (RealIntervalField, is_RealIntervalField,
                       is_RealIntervalFieldElement, RIF,
                       RealInterval)

# Complex numbers
from complex_field import ComplexField, is_ComplexField
from complex_number import (is_ComplexNumber, create_ComplexNumber as ComplexNumber)
Complexes = ComplexField
from complex_interval_field import ComplexIntervalField, is_ComplexIntervalField
from complex_interval import (is_ComplexIntervalFieldElement, create_ComplexIntervalFieldElement as ComplexIntervalFieldElement)

from complex_double import ComplexDoubleField, ComplexDoubleElement, CDF, is_ComplexDoubleElement

# Power series ring in one variable
from power_series_ring import PowerSeriesRing, is_PowerSeriesRing
from power_series_ring_element import PowerSeries, is_PowerSeries

# Laurent series ring in one variable
from laurent_series_ring import LaurentSeriesRing, is_LaurentSeriesRing
from laurent_series_ring_element import LaurentSeries, is_LaurentSeries

# Pseudo-ring of PARI objects.
from pari_ring import PariRing, Pari

# Big-oh notation
from big_oh import O

# Fraction field
from fraction_field import FractionField, is_FractionField
Frac = FractionField
from fraction_field_element import is_FractionFieldElement

# continued fractions
from contfrac import continued_fraction, CFF, ContinuedFractionField

# Arithmetic
from arith import *
from fast_arith import prime_range

from bernoulli_mod_p import bernoulli_mod_p, bernoulli_mod_p_single

from morphism import is_RingHomomorphism

from homset import is_RingHomset

from monomials import monomials

#from fast_polynomial.compiled_polynomial import compiled_polynomial

CC = ComplexField()
CIF = ComplexIntervalField()

# i = I = QuadraticField(-1, 'I').gen()
I = CC.gen()

from residue_field import ResidueField


from misc import composite_field

import tests
