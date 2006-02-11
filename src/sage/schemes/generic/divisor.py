"""
Divisors

AUTHORS:
   -- William Stein
   -- David Kohel
   -- David Joyner

EXAMPLES:
    sage: x,y,z = ProjectiveSpace(2, GF(5), names='xyz').gens()
    sage: C = Curve(y^2*z^7 - x^9 - x*z^8)
    sage: pts = C.rational_points(); pts
    [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1), (3 : 1 : 1), (3 : 4 : 1)]
    sage: D = C.divisor(pts[0])*3 - C.divisor(pts[1]) + C.divisor(pts[5])*10; D
    3*(0 : 0 : 1) - (0 : 1 : 0) + 10*(3 : 4 : 1)
    sage: D[1][0]
    -1
    sage: D[1][1]
    (0 : 1 : 0)
    sage: C.divisor([(3, pts[0]), (-1, pts[1]), (10,pts[5])])
    3*(0 : 0 : 1) - (0 : 1 : 0) + 10*(3 : 4 : 1)
"""

#*******************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#  Copyright (C) 2005 William Stein
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*******************************************************************************


from sage.structure.all import FormalSum

from sage.groups.group import AbelianGroup

from sage.rings.all import Z

from sage.ext.search import search

class Divisor_generic(FormalSum):
    def scheme(self):
        """
        Return the scheme that this divisor is on.

        EXAMPLES:
            sage: x,y = AffineSpace(2, GF(5), names='xyz').gens()
            sage: C = Curve(y^2 - x^9 - x)
            sage: pts = C.rational_points(); pts
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: D = C.divisor(pts[0])*3 - C.divisor(pts[1]); D
            3*(0, 0) - (2, 2)
            sage: D.scheme()
            Closed subscheme of Affine Space of dimension 2 over
            Finite Field of size 5 defined by:
              xyz_1^2 + 4*xyz_0 + 4*xyz_0^9
        """
        return self.parent().scheme()

class Divisor_curve_points(Divisor_generic):
    r"""
    For any curve $C$, use \code{C.divisor(v)} to construct a divisor
    on $C$.  Here $v$ can be either
    \begin{itemize}
       \item a rational point on $C$
       \item a list of rational points
       \item a list of 2-tuples $(c,P)$, where $c$ is
             an integer and $P$ is a rational point.
    \end{itemize}

    TODO: Divisors shouldn't be restricted to rational points.  The
    problem is that the divisor group is the formal sum of the group
    of points on the curve, and there's no implemented notion of point
    on $E/K$ that has coordinates in $L$.   This is what should
    be implemented, by adding an appropriate class to
    \code{schemes/generic/morphism.py}.

    EXAMPLES:
        sage: E = EllipticCurve([0, 0, 1, -1, 0])
        sage: P = E(0,0)
        sage: 10*P
        (161/16 : -2065/64 : 1)
        sage: D = E.divisor(P)
        sage: D
        (0 : 0 : 1)
        sage: 10*D
        10*(0 : 0 : 1)
        sage: E.divisor([P, P])
        2*(0 : 0 : 1)
        sage: E.divisor([(3,P), (-4,5*P)])
        3*(0 : 0 : 1) - 4*(1/4 : -5/8 : 1)
    """
    def __init__(self, v, check=True, reduce=True, parent=None):
        """
        INPUT:
            v -- a list of pairs (c, P), where c is an integer
                 and P is a point on a curve.  The P's must
                 all lie on the same curve.

        To create the 0 divisor use [(0, P)], so as to give
        the curve.
        """
        if not isinstance(v, (list, tuple)):
            v = [(1,v)]

        if len(v) < 1:
            raise ValueError, "v (=%s) must have length at least 1"%v

        if not (isinstance(v[0], tuple) and len(v[0]) == 2):
            C = v[0].scheme()
        else:
            C = v[0][1].scheme()

        if check:
            w = []
            for t in v:
                if isinstance(t, tuple) and len(t) == 2:
                    w.append((Z(t[0]), C(t[1])))
                else:
                    w.append((Z(1), C(t)))
            v = w

        if parent is None:
            parent = DivisorGroup(C)

        Divisor_generic.__init__(self, v, check=False, reduce=True,
                                 parent = parent)


    def support(self):
        """
        Return the support of this divisor, which is the set of points
        that occur in this divisor with nonzero coefficients.

        EXAMPLES:
            sage: x,y = AffineSpace(2, GF(5), names='xyz').gens()
            sage: C = Curve(y^2 - x^9 - x)
            sage: pts = C.rational_points(); pts
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: D = C.divisor(pts[0])*3 - C.divisor(pts[1]); D
            3*(0, 0) - (2, 2)
            sage: D.support()
            [(0, 0), (2, 2)]
        """
        try:
            return self.__support
        except AttributeError:
            self.__support = [self[i][1] for i in range(len(self))]
            return self.__support

    def coeff(self, P):
        """
        Return the coefficient of a given point P in this divisor.

        EXAMPLES:
            sage: x,y = AffineSpace(2, GF(5), names='xyz').gens()
            sage: C = Curve(y^2 - x^9 - x)
            sage: pts = C.rational_points(); pts
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: D = C.divisor(pts[0])*3 - C.divisor(pts[1]); D
            3*(0, 0) - (2, 2)
            sage: D.coeff(pts[0])
            3
            sage: D.coeff(pts[1])
            -1
        """
   	if not(P in self.support()):
       	    return Z(0)
        t, i = search(self.support(), P)
        assert t
        return self[i][0]


class DivisorGroup(AbelianGroup):
    def __init__(self, scheme):
        self.__scheme = scheme

    def _repr_(self):
        return "Group of Divisors on %s"%self.__scheme

    def __cmp__(self, right):
        if not isinstance(right, DivisorGroup):
            return -1
        return cmp(self.__scheme, right.__scheme)

    def scheme(self):
        return self.__scheme

class DivisorGroup_curve_points(DivisorGroup):
    def __call__(self, v):
        return Divisor_curve_points(v)

