"""
Miscellaneous module-related functions.

AUTHORS:
    -- William Stein (2007-11-18)
"""

####################################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
####################################################################################

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ

def gram_schmidt(B):
    """
    Return the Gram-Schmidt orthogonalization of the entries in the list
    B of vectors, along with the matrix mu of Gram-Schmidt coefficients.

    Note that the output vectors need not have unit length. We do this
    to avoid having to extract square roots.

    EXAMPLES:
        sage: B = [vector([1,2,1/5]), vector([1,2,3]), vector([-1,0,0])]
        sage: from sage.modules.misc import gram_schmidt
        sage: G, mu = gram_schmidt(B)
        sage: G
        [(1, 2, 1/5), (-1/9, -2/9, 25/9), (-4/5, 2/5, 0)]
        sage: G[0] * G[1]
        0
        sage: G[0] * G[2]
        0
        sage: G[1] * G[2]
        0
        sage: mu
        [      0       0       0]
        [   10/9       0       0]
        [-25/126    1/70       0]
        sage: a = matrix([])
        sage: a.gram_schmidt()
        ([], [])
        sage: a = matrix([[],[],[],[]])
        sage: a.gram_schmidt()
         ([], [])
    """
    if len(B) == 0 or len(B[0]) == 0:
        return B, matrix(ZZ,0,0,[])
    n = len(B)
    Bstar = [B[0]]
    K = B[0].base_ring().fraction_field()
    mu = matrix(K, n, n)
    for i in range(1,n):
        for j in range(i):
            mu[i,j] = B[i].dot_product(Bstar[j]) / (Bstar[j].dot_product(Bstar[j]))
        Bstar.append(B[i] - sum(mu[i,j]*Bstar[j] for j in range(i)))
    return Bstar, mu

