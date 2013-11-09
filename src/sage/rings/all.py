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
from ring import Ring
from commutative_ring import CommutativeRing
from integral_domain import IntegralDomain
from dedekind_domain import DedekindDomain
from principal_ideal_domain import PrincipalIdealDomain
from euclidean_domain import EuclideanDomain
from field import Field

from commutative_algebra_element import CommutativeAlgebraElement

# Ring element base classes
from ring_element import RingElement
from commutative_ring_element import CommutativeRingElement
from integral_domain_element import IntegralDomainElement
from dedekind_domain_element import DedekindDomainElement
from principal_ideal_domain_element import PrincipalIdealDomainElement
from euclidean_domain_element import EuclideanDomainElement
from field_element import FieldElement


# Ideals
from ideal import Ideal

# Quotient
from quotient_ring import QuotientRing

# Infinities
from infinity import infinity, Infinity, InfinityRing, unsigned_infinity, UnsignedInfinityRing

# Rational integers.
from integer_ring import IntegerRing, ZZ, crt_basis
from integer import Integer

# Rational numbers
from rational_field import RationalField, QQ
from rational import Rational
Rationals = RationalField

# Integers modulo n.
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing, Zmod
from sage.rings.finite_rings.integer_mod import IntegerMod, Mod, mod
Integers = IntegerModRing

# Finite fields
from finite_rings.all import *

# Number field
from number_field.all import *

# Function field
from function_field.all import *

# p-adic field
from padics.all import *
from padics.padic_printing import _printer_defaults as padic_printing

# Semirings
from semirings.all import *

# Real numbers
from real_mpfr import (RealField, RR,
                       create_RealNumber as RealNumber)   # this is used by the preparser to wrap real literals -- very important.
Reals = RealField

from real_double import RealDoubleField, RDF, RealDoubleElement

from real_lazy import RealLazyField, RLF, ComplexLazyField, CLF

# Polynomial Rings and Polynomial Quotient Rings
from polynomial.all import *


# Algebraic numbers
from qqbar import (AlgebraicRealField, AA,
                   AlgebraicReal,
                   AlgebraicField, QQbar,
                   AlgebraicNumber,
                   number_field_elements_from_algebraics)

# Intervals
from real_mpfi import (RealIntervalField,
                       RIF,
                       RealInterval)

# Complex numbers
from complex_field import ComplexField
from complex_number import (create_ComplexNumber as ComplexNumber)
Complexes = ComplexField
from complex_interval_field import ComplexIntervalField
from complex_interval import (create_ComplexIntervalFieldElement as ComplexIntervalFieldElement)

from complex_double import ComplexDoubleField, ComplexDoubleElement, CDF

from complex_mpc import MPComplexField

# Power series rings
from power_series_ring import PowerSeriesRing
from power_series_ring_element import PowerSeries

# Laurent series ring in one variable
from laurent_series_ring import LaurentSeriesRing
from laurent_series_ring_element import LaurentSeries

# Pseudo-ring of PARI objects.
from pari_ring import PariRing, Pari

# Big-oh notation
from big_oh import O

# Fraction field
from fraction_field import FractionField
Frac = FractionField

# continued fractions
from contfrac import continued_fraction, CFF, ContinuedFractionField

# Arithmetic
from arith import *
from fast_arith import prime_range

from bernoulli_mod_p import bernoulli_mod_p, bernoulli_mod_p_single

from monomials import monomials

#from fast_polynomial.compiled_polynomial import compiled_polynomial

CC = ComplexField()
CIF = ComplexIntervalField()

# i = I = QuadraticField(-1, 'I').gen()
I = CC.gen()

from residue_field import ResidueField


from misc import composite_field

import tests

# Universal Cyclotomic Field
from sage.rings.universal_cyclotomic_field.all import *

from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.invariant_theory', 'invariant_theory')
