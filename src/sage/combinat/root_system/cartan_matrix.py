"""
Cartan matrices

AUTHORS:

- Travis Scrimshaw (2012-04-22): Nicolas M. Thiery moved matrix creation to
    :class:`CartanType` to prepare cartan_matrix() for deprecation.
"""
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
from cartan_type import CartanType

def cartan_matrix(t):
    """
    Returns the Cartan matrix corresponding to type t.

    EXAMPLES::

        sage: cartan_matrix(['A', 4])
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [ 0  0 -1  2]
        sage: cartan_matrix(['B', 6])
        [ 2 -1  0  0  0  0]
        [-1  2 -1  0  0  0]
        [ 0 -1  2 -1  0  0]
        [ 0  0 -1  2 -1  0]
        [ 0  0  0 -1  2 -1]
        [ 0  0  0  0 -2  2]
        sage: cartan_matrix(['C', 4])
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -1  2 -2]
        [ 0  0 -1  2]
        sage: cartan_matrix(['D', 6])
        [ 2 -1  0  0  0  0]
        [-1  2 -1  0  0  0]
        [ 0 -1  2 -1  0  0]
        [ 0  0 -1  2 -1 -1]
        [ 0  0  0 -1  2  0]
        [ 0  0  0 -1  0  2]
        sage: cartan_matrix(['E',6])
        [ 2  0 -1  0  0  0]
        [ 0  2  0 -1  0  0]
        [-1  0  2 -1  0  0]
        [ 0 -1 -1  2 -1  0]
        [ 0  0  0 -1  2 -1]
        [ 0  0  0  0 -1  2]
        sage: cartan_matrix(['E',7])
        [ 2  0 -1  0  0  0  0]
        [ 0  2  0 -1  0  0  0]
        [-1  0  2 -1  0  0  0]
        [ 0 -1 -1  2 -1  0  0]
        [ 0  0  0 -1  2 -1  0]
        [ 0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0 -1  2]
        sage: cartan_matrix(['E', 8])
        [ 2  0 -1  0  0  0  0  0]
        [ 0  2  0 -1  0  0  0  0]
        [-1  0  2 -1  0  0  0  0]
        [ 0 -1 -1  2 -1  0  0  0]
        [ 0  0  0 -1  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0  0 -1  2]
        sage: cartan_matrix(['F', 4])
        [ 2 -1  0  0]
        [-1  2 -1  0]
        [ 0 -2  2 -1]
        [ 0  0 -1  2]

    This is different from MuPAD-Combinat, due to different node
    convention?

    ::

        sage: cartan_matrix(['G', 2])
        [ 2 -3]
        [-1  2]
        sage: cartan_matrix(['A',1,1])
        [ 2 -2]
        [-2  2]
        sage: cartan_matrix(['A', 3, 1])
        [ 2 -1  0 -1]
        [-1  2 -1  0]
        [ 0 -1  2 -1]
        [-1  0 -1  2]
        sage: cartan_matrix(['B', 3, 1])
        [ 2  0 -1  0]
        [ 0  2 -1  0]
        [-1 -1  2 -1]
        [ 0  0 -2  2]
        sage: cartan_matrix(['C', 3, 1])
        [ 2 -1  0  0]
        [-2  2 -1  0]
        [ 0 -1  2 -2]
        [ 0  0 -1  2]
        sage: cartan_matrix(['D', 4, 1])
        [ 2  0 -1  0  0]
        [ 0  2 -1  0  0]
        [-1 -1  2 -1 -1]
        [ 0  0 -1  2  0]
        [ 0  0 -1  0  2]
        sage: cartan_matrix(['E', 6, 1])
        [ 2  0 -1  0  0  0  0]
        [ 0  2  0 -1  0  0  0]
        [-1  0  2  0 -1  0  0]
        [ 0 -1  0  2 -1  0  0]
        [ 0  0 -1 -1  2 -1  0]
        [ 0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0 -1  2]
        sage: cartan_matrix(['E', 7, 1])
        [ 2 -1  0  0  0  0  0  0]
        [-1  2  0 -1  0  0  0  0]
        [ 0  0  2  0 -1  0  0  0]
        [ 0 -1  0  2 -1  0  0  0]
        [ 0  0 -1 -1  2 -1  0  0]
        [ 0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0 -1  2 -1]
        [ 0  0  0  0  0  0 -1  2]
        sage: cartan_matrix(['E', 8, 1])
        [ 2  0  0  0  0  0  0  0 -1]
        [ 0  2  0 -1  0  0  0  0  0]
        [ 0  0  2  0 -1  0  0  0  0]
        [ 0 -1  0  2 -1  0  0  0  0]
        [ 0  0 -1 -1  2 -1  0  0  0]
        [ 0  0  0  0 -1  2 -1  0  0]
        [ 0  0  0  0  0 -1  2 -1  0]
        [ 0  0  0  0  0  0 -1  2 -1]
        [-1  0  0  0  0  0  0 -1  2]
        sage: cartan_matrix(['F', 4, 1])
        [ 2 -1  0  0  0]
        [-1  2 -1  0  0]
        [ 0 -1  2 -1  0]
        [ 0  0 -2  2 -1]
        [ 0  0  0 -1  2]
        sage: cartan_matrix(['G', 2, 1])
        [ 2  0 -1]
        [ 0  2 -3]
        [-1 -1  2]

    .. note::

        This function is likely to be deprecated in favor of
        ``CartanType(...).cartan_matrix()``, to avoid polluting the
        global namespace.
    """
    return CartanType(t).cartan_matrix()
