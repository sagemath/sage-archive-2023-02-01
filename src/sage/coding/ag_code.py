r"""nodoctest -- totally broken right now
Algebraic-Geometric Codes

AUTHOR:
    -- David Joyner (2006-01-26): written
    -- William Stein (2006-01-23) -- inclusion in Sage
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner <wdj@usna.edu>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import copy

import linear_code

from sage.matrix.all import MatrixSpace

def ag_code(C, D, E):
    r"""
    INPUT:
        C -- a plane curve over a finite field F of prime order
        D -- a divisor on C
        E -- P1 + ... + Pn, another divisor, a sum of distinct
             points on X whose support is disjoint from that
             of D and with n > deg(D) > 0.

    OUTPUT:
        The linear code defined by C, D, and E as in Stichtenoth.

    Calls Singular's \code{BrillNoether functions}.  Will break on some
    singular curves or if the field size is too large, etc.; when
    this happens a ??? exception is raised.

    EXAMPLES:
        sage: x,y,z = ProjectiveSpace(2, GF(17), names = 'xyz').gens()
        sage: C   = Curve(y^2*z^7 - x^9 - x*z^8)
        sage: pts = C.rational_points(sorted=False)
        sage: D   = C.divisor([(3,pts[0]), (-1,pts[1]), (10,pts[5])]); D
        3*(0 : 0 : 1) - (0 : 1 : 0) + 10*(11 : 0 : 1)
        sage: E   = add([C.divisor(p) for p in [pts[2],pts[3],pts[4]]]); E
        (1 : 11 : 1) + (1 : 6 : 1) + (10 : 0 : 1)
        sage: V = ag_code(C, D,E)
        sage: V
        Linear code of length 2, dimension 2 over Finite Field of size 17
        sage: V.basis()
        [(12, 12), (11, 11), (12, 12), (2, 2), (2, 2), (3, 14), (13, 4), (7, 7), (7, 10)]

        sage: P2 = ProjectiveSpace(2, GF(11), names = ['x','y','z'])
        sage: x, y, z = P2.coordinate_ring().gens()
        sage: f = y^8 + x^8 - z^8
        sage: C = Curve(f)
        sage: pts = C.rational_points()
        sage: D = C.divisor(pts[0])*10 - C.divisor(pts[1]) + C.divisor(pts[5])*10
        sage: E = add([C.divisor(pts[i]) for i in [2,3,4,6,7,8,9,10,11]])
        sage: V = ag_code(C, D,E); V
        Linear code of length 8, dimension 2 over Finite Field of size 11
        sage: V.basis()
        [(6, 6), (11, 11), (13, 13), (13, 13), (15, 15), (8, 9), (3, 14), (7, 7), (9, 8)]

        sage: P2 = ProjectiveSpace(2, GF(11), names = ['x','y','z'])
        sage: x, y, z = P2.coordinate_ring().gens()
        sage: f = x^3*y - y^3*z + z^3*x
        sage: C = Curve(f)
        sage: pts = C.rational_points(sorted=False)
        sage: D = C.divisor(pts[0])*3 - C.divisor(pts[1]) + C.divisor(pts[5])*5
        sage: len(C.riemann_roch_basis(D))
        5
        sage: E = add([C.divisor(pts[i]) for i in [2,3,4,6,7,8,9,10,11]])
        sage: V = ag_code(C, D,E); V
        Linear code of length 8, dimension 5 over Finite Field of size 11
        sage: V.basis()
        [(1, 8, 4, 1, 10, 1, 5, 8), (2, 8, 1, 7, 3, 5, 8, 2), (1, 7, 4, 7, 5, 7, 4, 9), (0, 7, 0, 9, 5, 8, 6, 3), (6, 10, 1, 4, 6, 4, 8, 0)]
        sage: V.minimum_distance()
        3

    AUTHOR: David Joyner (2006-01)
    """
    F = C.base_ring()
    one = F(1)
    B = C.riemann_roch_basis(D)
    k = len(B)
    if k == 0:
        A = MatrixSpace(F, 0, 0)(0)
        return sage.coding.all.LinearCode(A)

    P = C.rational_points(algorithm="bn", sorted=False)  # algorithm="bn" does not work here ....

    N = len(P)
    pts = copy.copy(E.support())
    for p in pts:
        for i in range(k):
            if (B[i].denominator())(F(p[0]), F(p[1]), F(p[2])) == 0:
                pts.remove(p)
                break

    MS = MatrixSpace(F, len(B), len(pts))
    pts.sort()
    G = [[B[i](F(p[0]), F(p[1]), F(p[2])) for p in pts] for i in range(k)]
    G = MS(G)
    return linear_code.LinearCode(G)
