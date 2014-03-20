"""
Lattice constructor
"""

#*****************************************************************************
#
#   Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2014 Martin Albrecht <martinralbecht@googlemailc.om>
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

from integer_lattice import IntegerLattice

def Lattice(basis, lll_reduce=True):
    """Construct a new lattice from ``basis``.

    INPUT:
      - ``basis``
      - a list of vectors or
      - a matrix over the integers
      - an element of an absolute order

    - ``lll_reduce`` -- (default: ``True``) run LLL reduction on the basis
      on construction

    EXAMPLES:

    We construct a lattice from a list of rows::
    
        sage: Lattice([[1,0,3],[0,2,1], [0,2,7]])
        Lattice of degree 3 and rank 3 over Integer Ring
        Basis matrix:
        [-2  0  0]
        [ 0  2  1]
        [ 1 -2  2]

    Sage includes a generator for hard lattices from cryptography::
    
        sage: A = sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True)
        sage: Lattice(A)
        Lattice of degree 10 and rank 10 over Integer Ring
        Basis matrix:
        [ 0  1  2  0  1  2 -1  0 -1 -1]
        [ 0  1  0 -3  0  0  0  0  3 -1]
        [ 1  1 -1  0 -3  0  0  1  2 -2]
        [-1  2 -1 -1  2 -2  1 -1  0 -1]
        [ 1  0 -4  2  0  1 -2 -1  0  0]
        [ 2  3  0  1  1  0 -2  3  0  0]
        [-2 -3 -2  0  0  1 -1  1  3 -2]
        [-3  0 -1  0 -2 -1 -2  1 -1  1]
        [ 1  4 -1  1  2  2  1  0  3  1]
        [-1 -1  0 -3 -1  2  2  3 -1  0]

    You can also construct the lattice directly::
    
        sage: sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True, lattice=True)
        Lattice of degree 10 and rank 10 over Integer Ring
        Basis matrix:
        [ 0  1  2  0  1  2 -1  0 -1 -1]
        [ 0  1  0 -3  0  0  0  0  3 -1]
        [ 1  1 -1  0 -3  0  0  1  2 -2]
        [-1  2 -1 -1  2 -2  1 -1  0 -1]
        [ 1  0 -4  2  0  1 -2 -1  0  0]
        [ 2  3  0  1  1  0 -2  3  0  0]
        [-2 -3 -2  0  0  1 -1  1  3 -2]
        [-3  0 -1  0 -2 -1 -2  1 -1  1]
        [ 1  4 -1  1  2  2  1  0  3  1]
        [-1 -1  0 -3 -1  2  2  3 -1  0]

    We construct an ideal lattice from an element of an absolute order::
    
        sage: K.<a>  = CyclotomicField(17)
        sage: O = K.ring_of_integers()
        sage: f = O.random_element(); f
        -a^15 - a^12 - a^10 - 8*a^9 - a^8 - 4*a^7 + 3*a^6 + a^5 + 2*a^4 + 8*a^3 - a^2 + a + 1

        sage: Lattice(f)
        Lattice of degree 16 and rank 16 over Integer Ring
        Basis matrix:
        [ 1  1 -1  8  2  1  3 -4 -1 -8 -1  0 -1  0  0 -1]
        [-1  0  1  1 -1  8  2  1  3 -4 -1 -8 -1  0 -1  0]
        [ 1  1  0  1  2  2  0  9  3  2  4 -3  0 -7  0  1]
        [ 1  0  1  1  0  1  2  2  0  9  3  2  4 -3  0 -7]
        [ 2 -5 -2 -9 -2 -1 -2 -1 -1 -2 -1  0  0 -2  7  1]
        [ 1  4  0 -5 -3  0  3  5 -2  2  0 -7  4  0 -6 -2]
        [-7  4  0 -6 -2  0  1  4  0 -5 -3  0  3  5 -2  2]
        [-1  0  0 -1  0  1  1 -1  8  2  1  3 -4 -1 -8 -1]
        [-1 -1 -2  4  1  9  1 -1  0  1 -8 -1 -1 -4  3  2]
        [-1 -2  6 -6  8  1 -3  5  3  1  1  0 -2  4  3  2]
        [ 4  8  2  7 -3  2 -1  2  0 -4  0  3  6  3  0 -4]
        [ 0 -1  0  1  1 -1  8  2  1  3 -4 -1 -8 -1  0 -1]
        [ 2 -7 -1  0 -2  5  2  9  2  1  2  1  1  2  1  0]
        [-1  7 -5  9  2 -2  6  4  2  2  1 -1  5  4  3  1]
        [ 3  1 -6  0  3  1 -5  2  2  8 -4  4 -3  2 -6 -7]
        [ 2  6 -4  4  0 -1  7  0 -6  3  9  1 -3 -1  4  3]

    Sage also interfaces with fpLLL's lattice generator::

        sage: from sage.libs.fplll.fplll import gen_simdioph
        sage: gen_simdioph(8, 20, 10)
        [   1024  829556  161099   11567  521155  769480  639201  689979]
        [      0 1048576       0       0       0       0       0       0]
        [      0       0 1048576       0       0       0       0       0]
        [      0       0       0 1048576       0       0       0       0]
        [      0       0       0       0 1048576       0       0       0]
        [      0       0       0       0       0 1048576       0       0]
        [      0       0       0       0       0       0 1048576       0]
        [      0       0       0       0       0       0       0 1048576]

    """
    return IntegerLattice(basis, lll_reduce=lll_reduce)
        