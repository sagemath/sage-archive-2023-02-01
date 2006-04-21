"""
Divisors on schemes

AUTHORS:
   -- William Stein
   -- David Kohel
   -- David Joyner

EXAMPLES:
    sage: x,y,z = ProjectiveSpace(2, GF(5), names='x,y,z').gens()
    sage: C = Curve(y^2*z^7 - x^9 - x*z^8)
    sage: pts = C.rational_points(); pts
    [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1), (3 : 1 : 1), (3 : 4 : 1)]
    sage: D = C.divisor(pts[0])*3 - C.divisor(pts[1]) + C.divisor(pts[5])*10; D
    -(z, x) + 3*(y, x) + 10*(z + y, 2*z + x)
    sage: D[1][0]
    3
    sage: D[1][1]
    Ideal (y, x) of Polynomial Ring in x, y, z over Finite Field of size 5
    sage: C.divisor([(3, pts[0]), (-1, pts[1]), (10,pts[5])])
    -(z, x) + 3*(y, x) + 10*(z + y, 2*z + x)
"""
#*******************************************************************************
#  Copyright (C) 2005 David Kohel <kohel@maths.usyd.edu.au>
#  Copyright (C) 2005 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

import sage.misc.misc #as repr_lincomb

from sage.structure.all import FormalSum

from sage.rings.all import Z

from projective_space import is_ProjectiveSpace

from affine_space import is_AffineSpace

from morphism import is_SchemeMorphism

import divisor_group

from sage.ext.search import search


def CurvePointToIdeal(C,P):
    A = C.ambient_space()
    R = A.coordinate_ring()
    n = A.ngens()
    x = A.gens()
    polys = [ ]
    m = n-1
    while m > 0 and P[m] == 0:
        m += -1
    if is_ProjectiveSpace(A):
        a_m = P[m]
        x_m = x[m]
        for i in range(m):
            ai = P[i]
            if ai == 0:
                polys.append(x[i])
            else:
                polys.append(a_m*x[i]-ai*x_m)
    elif is_AffineSpace(A):
        for i in range(m+1):
            ai = P[i]
            if ai == 0:
                polys.append(x[i])
            else:
                polys.append(x[i]-ai)
    for i in range(m+1,n):
        polys.append(x[i])
    return R.ideal(polys)

def is_Divisor(Div):
    return isinstance(Div, Divisor_generic)

def is_DivisorGroup(Div):
    return isinstance(Div, DivisorGroup)

class Divisor_generic(FormalSum):
    def scheme(self):
        """
        Return the scheme that this divisor is on.

        EXAMPLES:
            sage: x, y = AffineSpace(2, GF(5), names='x,y').gens()
            sage: C = Curve(y^2 - x^9 - x)
            sage: pts = C.rational_points(); pts
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: D = C.divisor(pts[0])*3 - C.divisor(pts[1]); D
            3*(y, x) - (3 + y, 3 + x)
            sage: D.scheme()
            Affine Curve over Finite Field of size 5 defined by y^2 + 4*x + 4*x^9
        """
        return self.parent().scheme()

class Divisor_curve(Divisor_generic):
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
        (y, x)
        sage: 10*D
        10*(y, x)
        sage: E.divisor([P, P])
        2*(y, x)
        sage: E.divisor([(3,P), (-4,5*P)])
        3*(y, x) - 4*(-1/4*z + x, 5/8*z + y)
    """
    def __init__(self, v, check=True, reduce=True, parent=None):
        """
        INPUT:
            v -- a list of pairs (c, P), where c is an integer
                 and P is a point on a curve.  The P's must
                 all lie on the same curve.

        To create the 0 divisor use [(0, P)], so as to give
        the curve.

        TODO: Include an extension field in the definition of the
        divisor group.
        """
        if not isinstance(v, (list, tuple)):
            v = [(1,v)]

        if parent is None:
            if len(v) > 0:
                t = v[0]
                if isinstance(t, tuple) and len(t) == 2:
                    try:
                        C = t[1].scheme()
                    except TypeError:
                        raise TypeError, \
                              "Argument v (= %s) must consist of multiplicities and points on a scheme."
                else:
                    try:
                        C = t.scheme()
                    except TypeError:
                        raise TypeError, \
                              "Argument v (= %s) must consist of multiplicities and points on a scheme."
                parent = divisor_group.DivisorGroup_curve(C)
            else:
                raise TypeError, \
                      "Argument v (= %s) must consist of multiplicities and points on a scheme."

        if len(v) < 1:
            raise ValueError, "v (=%s) must have length at least 1"%v
        know_points = False
        if check:
            w = []
            points = []
            know_points = True
            for t in v:
                if isinstance(t, tuple) and len(t) == 2:
                    n = Z(t[0])
                    I = t[1]
                    points.append((n,I))
                else:
                    n = Z(1)
                    I = t
                if is_SchemeMorphism(I):
                    I = CurvePointToIdeal(C,I)
                else:
                    know_points = False
                w.append((n,I))
            v = w
        Divisor_generic.__init__(
            self, v, check=False, reduce=True, parent=parent)

        if know_points:
            self.__points = points

    def _repr_(self):
        ideals = [ z[1] for z in self ]
        coeffs = [ z[0] for z in self ]
        polys = [ tuple(I.gens()) for I in ideals ]
        return sage.misc.misc.repr_lincomb(polys, coeffs)

    def support(self):
        """
        Return the support of this divisor, which is the set of points
        that occur in this divisor with nonzero coefficients.

        EXAMPLES:
            sage: x,y = AffineSpace(2, GF(5), names='xy').gens()
            sage: C = Curve(y^2 - x^9 - x)
            sage: pts = C.rational_points(); pts
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: D = C.divisor([(3,pts[0]), (-1, pts[1])]); D
            3*(y, x) - (3 + y, 3 + x)
            sage: D.support()
            [(0, 0), (2, 2)]
        """
        try:
            return self.__support
        except AttributeError:
            try:
                self.__support = [s[1] for s in self.__points]
                return self.__support
            except AttributeError:
                raise NotImplementedError

    def coeff(self, P):
        """
        Return the coefficient of a given point P in this divisor.

        EXAMPLES:
            sage: x,y = AffineSpace(2, GF(5), names='xy').gens()
            sage: C = Curve(y^2 - x^9 - x)
            sage: pts = C.rational_points(); pts
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: D = C.divisor([(3,pts[0]), (-1,pts[1])]); D
            3*(y, x) - (3 + y, 3 + x)
            sage: D.coeff(pts[0])
            3
            sage: D.coeff(pts[1])
            -1
        """
        P = self.scheme()(P)
   	if not(P in self.support()):
       	    return Z(0)
        t, i = search(self.support(), P)
        assert t
        try:
            return self.__points[i][0]
        except AttributeError:
                raise NotImplementedError

