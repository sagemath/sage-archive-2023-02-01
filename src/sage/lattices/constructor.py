"""
Lattice constructor

AUTHORS:

- Martin Albrecht (2014-03): initial version
"""

#*****************************************************************************
#  Copyright (C) 2014 Martin Albrecht <martinralbecht@googlemail.com>
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
from sage.rings.number_field.number_field_element import OrderElement_absolute
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.modules.free_module import FreeModule_ambient_pid

def RealLattice(basis, lll_reduce=True):
    """
    Construct a new lattice from ``basis``.

    INPUT:
      - ``basis``
      - a list of vectors or
      - a matrix over the integers
      - an element of an absolute order

    - ``lll_reduce`` -- (default: ``True``) run LLL reduction on the basis
      on construction

    EXAMPLES:

    We construct a lattice from a list of rows::
    
        sage: RealLattice([[1,0,3],[0,2,1], [0,2,7]])
        Lattice of degree 3 and rank 3 over Integer Ring
        Basis matrix:
        [-2  0  0]
        [ 0  2  1]
        [ 1 -2  2]

    Sage includes a generator for hard lattices from cryptography::
    
        sage: A = sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True)
        sage: RealLattice(A)
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

        sage: RealLattice(f)
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

    We construct `\ZZ^n`::

        sage: RealLattice(ZZ^10)
        Lattice of degree 10 and rank 10 over Integer Ring
        Basis matrix:
        [1 0 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 0 1]

    
    Sage also interfaces with fpLLL's lattice generator::

        sage: from sage.libs.fplll.fplll import gen_simdioph
        sage: RealLattice(gen_simdioph(8, 20, 10), lll_reduce=False)
        Lattice of degree 8 and rank 8 over Integer Ring
        Basis matrix:
        [   1024  829556  161099   11567  521155  769480  639201  689979]
        [      0 1048576       0       0       0       0       0       0]
        [      0       0 1048576       0       0       0       0       0]
        [      0       0       0 1048576       0       0       0       0]
        [      0       0       0       0 1048576       0       0       0]
        [      0       0       0       0       0 1048576       0       0]
        [      0       0       0       0       0       0 1048576       0]
        [      0       0       0       0       0       0       0 1048576]

    For now, however, only integer lattices are implemented::

        sage: RealLattice(random_matrix(QQ, 10, 10))
        Traceback (most recent call last):
        ...
        NotImplementedError: Only integer lattices are currently implemented.

        sage: RealLattice(random_matrix(RR, 10, 10))
        Traceback (most recent call last):
        ...
        NotImplementedError: Only integer lattices are currently implemented.

    """

    if isinstance(basis, OrderElement_absolute):
        basis = basis.matrix()
    elif isinstance(basis, FreeModule_ambient_pid):
        basis = basis.basis_matrix()
        
    try:
        basis = matrix(ZZ, basis)
    except TypeError:
        raise NotImplementedError("Only integer lattices are currently implemented.")
    
    return IntegerLattice(basis, lll_reduce=lll_reduce)
        