"""
Decoding methods for linear error-correcting codes.
Methods implemented:

 * nearest neighbor
 * syndrome

AUTHOR:
    -- David Joyner (2009-02-01): initial version

TODO:
  Add lots more methods!
"""
#*****************************************************************************
#       Copyright (C) 2009 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.coding.linear_code import hamming_weight, LinearCode

def syndrome(C, v):
    """
    The vector v represents a received word, so should
    be in the same ambient space V as C. Returns the
    elements in V (including v) which belong to the
    syndrome of v (ie, the coset v+C, sorted by weight).

    EXAMPLES:
        sage: C = HammingCode(2,GF(3)); C
        Linear code of length 4, dimension 2 over Finite Field of size 3
        sage: V = VectorSpace(GF(3), 4)
        sage: v = V([0, 2, 0, 1])
        sage: from sage.coding.decoder import syndrome
        sage: syndrome(C, v)
        [(2, 0, 0, 0),
         (0, 2, 0, 1),
         (0, 0, 2, 2),
         (0, 1, 1, 0),
         (1, 2, 2, 0),
         (1, 0, 1, 1),
         (1, 1, 0, 2),
         (2, 2, 1, 2),
         (2, 1, 2, 1)]

    """
    V = C.ambient_space()
    if not(type(v)==list):
        v = v.list()
    v = V(v)
    coset = [[c+v,hamming_weight(c+v)] for c in C]
    coset.sort(lambda x,y: x[1]-y[1])
    return [x[0] for x in coset]

def coset_leader(C, v):
    """
    The vector v represents a received word, so should
    be in the same ambient space V as C. Returns an
    element of the syndrome of v of lowest weight.

    EXAMPLES:
        sage: C = HammingCode(2,GF(3)); C
        Linear code of length 4, dimension 2 over Finite Field of size 3
        sage: V = VectorSpace(GF(3), 4)
        sage: v = V([0, 2, 0, 1])
        sage: from sage.coding.decoder import coset_leader
        sage: coset_leader(C, v)
        ((2, 0, 0, 0), 1)
        sage: coset_leader(C, v)[0]-v in C
        True

    """
    coset = [[c+v, hamming_weight(c+v)] for c in C]
    wts = [x[1] for x in coset]
    min_wt = min(wts)
    s = C[0]                # initializing
    w = hamming_weight(v)   # initializing
    for x in coset:
        if x[1]==min_wt:
            w = x[1]
            s = x[0]
            break
    return s,w

def decode(C, v, method="syndrome"):
    """
    The vector v represents a received word, so should
    be in the same ambient space V as C. Returns an
    element in C which is closest to v in the Hamming
    metric.

    Methods implemented include "nearest neighbor" (essentially
    a brute force search) and "syndrome".

    EXAMPLES:
        sage: C = HammingCode(2,GF(3))
        sage: V = VectorSpace(GF(3), 4)
        sage: v = V([0, 2, 0, 1])
        sage: v in C
        False
        sage: from sage.coding.decoder import decode
        sage: c = decode(C, v);c
        (1, 2, 0, 1)
        sage: c in C
        True
        sage: c = decode(C, v, method="nearest neighbor");c
        (1, 2, 0, 1)
        sage: C = HammingCode(3,GF(3)); C
        Linear code of length 13, dimension 10 over Finite Field of size 3
        sage: V = VectorSpace(GF(3), 13)
        sage: v = V([2]+[0]*12)
        sage: decode(C, v)
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        sage: decode(C, v, method="nearest neighbor")
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    """
    V = C.ambient_space()
    if not(type(v)==list):
        v = v.list()
    v = V(v)
    if method=="nearest neighbor":
        diffs = [[c-v,hamming_weight(c-v)] for c in C]
        diffs.sort(lambda x,y:  x[1]-y[1])
        return diffs[0][0]+v
    if method=="syndrome":
        return -V(syndrome(C, v)[0])+v
