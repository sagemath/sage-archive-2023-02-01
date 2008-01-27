#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from sage.rings.ring import Ring
from sage.rings.integer import Integer

from sage.algebras.algebra import Algebra

import sage.combinat.partition
import sage.combinat.skew_partition
import sage.structure.parent_gens
import sage.libs.symmetrica.all as symmetrica
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra, CombinatorialAlgebraElement
from sage.matrix.constructor import matrix

from sage.rings.integer_ring import IntegerRing
from sage.rings.rational_field import RationalField

from sage.misc.misc import repr_lincomb
from sage.algebras.algebra_element import AlgebraElement

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

import operator
import hall_littlewood

ZZ = IntegerRing()
QQ = RationalField()

translate = {'monomial':'MONOMIAL', 'homogeneous':'HOMSYM', 'power':'POWSYM', 'elementary':'ELMSYM', 'schur':'SCHUR'}

conversion_functions = {}

def init():
    #Set up the conversion functions
    for other_basis in translate:
        for basis in translate:
            try:
                conversion_functions[(other_basis, basis)] = eval('symmetrica.t_' + translate[other_basis] + '_' +  translate[basis])
            except AttributeError:
                pass


init()


