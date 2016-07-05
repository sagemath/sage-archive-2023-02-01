"""
Projective curves.

EXAMPLES:

We can construct curves in either a projective plane::

    sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
    sage: C = Curve([y*z^2 - x^3], P); C
    Projective Plane Curve over Rational Field defined by -x^3 + y*z^2

or in higher dimensional projective spaces::

    sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
    sage: C = Curve([y*w^3 - x^4, z*w^3 - x^4], P); C
    Projective Curve over Rational Field defined by -x^4 + y*w^3, -x^4 + z*w^3

AUTHORS:

- William Stein (2005-11-13)

- David Joyner (2005-11-13)

- David Kohel (2006-01)

- Moritz Minzlaff (2010-11)
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.fields import Fields
from sage.categories.homset import Hom
from sage.interfaces.all import singular
from sage.matrix.constructor import matrix
from sage.misc.all import add, sage_eval
from sage.rings.all import degree_lowest_rational_function
from sage.rings.qqbar import QQbar
from sage.schemes.affine.affine_space import AffineSpace

from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective
from sage.schemes.projective.projective_space import is_ProjectiveSpace

from curve import Curve_generic

class ProjectiveCurve(Curve_generic, AlgebraicScheme_subscheme_projective):
    def _repr_type(self):
        r"""
        Return a string representation of the type of this curve.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([y^3 - z^3 - w^3, z*x^3 - y^4])
            sage: C._repr_type()
            'Projective'
        """
        return "Projective"

    def __init__(self, A, X):
        r"""
        Initialization function.

        EXMAPLES::

            sage: P.<x,y,z,w,u> = ProjectiveSpace(GF(7), 4)
            sage: C = Curve([y*u^2 - x^3, z*u^2 - x^3, w*u^2 - x^3, y^3 - x^3], P); C
            Projective Curve over Finite Field of size 7 defined by -x^3 + y*u^2,
            -x^3 + z*u^2, -x^3 + w*u^2, -x^3 + y^3

        ::

            sage: K.<u> = CyclotomicField(11)
            sage: P.<x,y,z,w> = ProjectiveSpace(K, 3)
            sage: C = Curve([y*w - u*z^2 - x^2, x*w - 3*u^2*z*w], P); C
            Projective Curve over Cyclotomic Field of order 11 and degree 10 defined
            by -x^2 + (-u)*z^2 + y*w, x*w + (-3*u^2)*z*w
        """
        if not is_ProjectiveSpace(A):
            raise TypeError("A (=%s) must be a projective space"%A)
        Curve_generic.__init__(self, A, X)
        d = self.dimension()
        if d != 1:
            raise ValueError("defining equations (=%s) define a scheme of dimension %s != 1"%(X,d))

    def affine_patch(self, i, AA=None):
        r"""
        Return the i-th affine patch of this projective curve.

        INPUT:

        - ``i`` -- affine coordinate chart of the projective ambient space of this curve to compute affine patch
          with respect to.

        - ``AA`` -- (default: None) ambient affine space, this is constructed if it is not given.

        OUTPUT:

        - a curve in affine space.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(CC, 3)
            sage: C = Curve([y*z - x^2, w^2 - x*y], P)
            sage: C.affine_patch(0)
            Affine Curve over Complex Field with 53 bits of precision defined by
            x0*x1 - 1.00000000000000, x2^2 - x0

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve(x^3 - x^2*y + y^3 - x^2*z, P)
            sage: C.affine_patch(1)
            Affine Plane Curve over Rational Field defined by x0^3 - x0^2*x1 - x0^2 + 1

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: P.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([u^2 - v^2], P)
            sage: C.affine_patch(1, A).ambient_space() == A
            True
        """
        from constructor import Curve
        return Curve(AlgebraicScheme_subscheme_projective.affine_patch(self, i, AA))

    def multiplicity(self, P):
        r"""
        Return the multiplicity of this projective curve at the point ``P``.

        This is computed as the corresponding multiplicity of an affine patch of this curve that
        contains the point. This curve must be defined over a field. An error is returned if ``P``
        not a point on this curve.

        INPUT:

        - ``P`` -- a point in the ambient space of this curve.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^4 - x^3*z - x^2*z^2], P)
            sage: Q = P([0,0,1])
            sage: C.multiplicity(Q)
            2

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(RR, 3)
            sage: C = Curve([y^8 - x^2*z*w^5, w^2 - 2*y^2 - x*z], P)
            sage: Q1 = P([-1,-1,1,1])
            sage: C.multiplicity(Q1)
            1
            sage: Q2 = P([1,0,0,0])
            sage: C.multiplicity(Q2)
            7
            sage: Q3 = P([0,0,1,0])
            sage: C.multiplicity(Q3)
            8

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(GF(29), 3)
            sage: C = Curve([y^17 - x^5*w^4*z^8, x*y - z^2], P)
            sage: Q = P([3,0,0,1])
            sage: C.multiplicity(Q)
            8

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = P.curve([y^2*z^5 - x^7])
            sage: Q = P([-1,-1,1])
            sage: C.multiplicity(Q)
            Traceback (most recent call last):
            ...
            TypeError: (=(-1 : -1 : 1)) is not a point on (=Projective Plane Curve
            over Rational Field defined by -x^7 + y^2*z^5)
        """
        if not self.base_ring() in Fields():
            raise TypeError("curve must be defined over a field")

        # Check whether P is a point on this curve
        try:
            P = self(P)
        except TypeError:
            raise TypeError("(=%s) is not a point on (=%s)"%(P,self))

        # Find an affine chart of the ambient space of self that contains P
        i = 0
        while(P[i] == 0):
            i = i + 1
        C = self.affine_patch(i)
        Q = list(P)
        t = Q.pop(i)
        Q = [1/t*Q[j] for j in range(self.ambient_space().dimension_relative())]
        return C.multiplicity(C.ambient_space()(Q))

    def is_complete_intersection(self):
        r"""
        Return whether this projective curve is or is not a complete intersection.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([x*y - z*w, x^2 - y*w, y^2*w - x*z*w], P)
            sage: C.is_complete_intersection()
            False

        ::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([y*w - x^2, z*w^2 - x^3], P)
            sage: C.is_complete_intersection()
            True

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = Curve([z^2 - y*w, y*z - x*w, y^2 - x*z], P)
            sage: C.is_complete_intersection()
            False
        """
        singular.lib("sing.lib")
        I = singular.simplify(self.defining_ideal(), 10)
        L = singular.is_ci(I).sage()
        return len(self.ambient_space().gens()) - len(I.sage().gens()) == L[-1]

class ProjectivePlaneCurve(ProjectiveCurve):
    def __init__(self, A, f):
        r"""
        Initialization function.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
            sage: C = Curve([y*z - x^2 - QQbar.gen()*z^2], P); C
            Projective Plane Curve over Algebraic Field defined by
            -x^2 + y*z + (-I)*z^2

        ::

            sage: P.<x,y,z> = ProjectiveSpace(GF(5^2, 'v'), 2)
            sage: C = Curve([y^2*z - x*z^2 - z^3], P); C
            Projective Plane Curve over Finite Field in v of size 5^2 defined by y^2*z - x*z^2 - z^3
        """
        if not (is_ProjectiveSpace(A) and A.dimension != 2):
            raise TypeError("Argument A (= %s) must be a projective plane."%A)
        Curve_generic.__init__(self, A, [f])

    def _repr_type(self):
        r"""
        Return a string representation of the type of this curve.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y*z^3 - 5/7*x^4 + 4*x^3*z - 9*z^4], P)
            sage: C._repr_type()
            'Projective Plane'
        """
        return "Projective Plane"

    def arithmetic_genus(self):
        r"""
        Return the arithmetic genus of this curve.

        This is the arithmetic genus `g_a(C)` as defined in
        Hartshorne. If the curve has degree `d` then this is simply
        `(d-1)(d-2)/2`. It need *not* equal the geometric genus
        (the genus of the normalization of the curve).

        EXAMPLE::

            sage: x,y,z = PolynomialRing(GF(5), 3, 'xyz').gens()
            sage: C = Curve(y^2*z^7 - x^9 - x*z^8); C
            Projective Plane Curve over Finite Field of size 5 defined by -x^9 + y^2*z^7 - x*z^8
            sage: C.arithmetic_genus()
            28
            sage: C.genus()
            4
        """
        d = self.defining_polynomial().total_degree()
        return int((d-1)*(d-2)/2)

    def divisor_of_function(self, r):
        """
        Return the divisor of a function on a curve.

        INPUT: r is a rational function on X

        OUTPUT:


        -  ``list`` - The divisor of r represented as a list of
           coefficients and points. (TODO: This will change to a more
           structural output in the future.)


        EXAMPLES::

            sage: FF = FiniteField(5)
            sage: P2 = ProjectiveSpace(2, FF, names = ['x','y','z'])
            sage: R = P2.coordinate_ring()
            sage: x, y, z = R.gens()
            sage: f = y^2*z^7 - x^9 - x*z^8
            sage: C = Curve(f)
            sage: K = FractionField(R)
            sage: r = 1/x
            sage: C.divisor_of_function(r)     # todo: not implemented  !!!!
            [[-1, (0, 0, 1)]]
            sage: r = 1/x^3
            sage: C.divisor_of_function(r)     # todo: not implemented  !!!!
            [[-3, (0, 0, 1)]]
        """
        F = self.base_ring()
        f = self.defining_polynomial()
        x, y, z = f.parent().gens()
        pnts = self.rational_points()
        divf = []
        for P in pnts:
            if P[2] != F(0):
                # What is the '5' in this line and the 'r()' in the next???
                lcs = self.local_coordinates(P,5)
                ldg = degree_lowest_rational_function(r(lcs[0],lcs[1]),z)
                if ldg[0] != 0:
                    divf.append([ldg[0],P])
        return divf


    def local_coordinates(self, pt, n):
        r"""
        Return local coordinates to precision n at the given point.

            Behaviour is flaky - some choices of `n` are worst that
            others.


        INPUT:


        -  ``pt`` - an F-rational point on X which is not a
           point of ramification for the projection (x,y) - x.

        -  ``n`` - the number of terms desired


        OUTPUT: x = x0 + t y = y0 + power series in t

        EXAMPLES::

            sage: FF = FiniteField(5)
            sage: P2 = ProjectiveSpace(2, FF, names = ['x','y','z'])
            sage: x, y, z = P2.coordinate_ring().gens()
            sage: C = Curve(y^2*z^7-x^9-x*z^8)
            sage: pt = C([2,3,1])
            sage: C.local_coordinates(pt,9)     # todo: not implemented  !!!!
                  [2 + t, 3 + 3*t^2 + t^3 + 3*t^4 + 3*t^6 + 3*t^7 + t^8 + 2*t^9 + 3*t^11 + 3*t^12]
        """

        f = self.defining_polynomial()
        R = f.parent()
        F = self.base_ring()
        p = F.characteristic()
        x0 = F(pt[0])
        y0 = F(pt[1])
        astr = ["a"+str(i) for i in range(1,2*n)]
        x,y = R.gens()
        R0 = PolynomialRing(F,2*n+2,names = [str(x),str(y),"t"]+astr)
        vars0 = R0.gens()
        t = vars0[2]
        yt = y0*t**0 + add([vars0[i]*t**(i-2) for i in range(3,2*n+2)])
        xt = x0+t
        ft = f(xt,yt)
        S = singular
        S.eval('ring s = '+str(p)+','+str(R0.gens())+',lp;')
        S.eval('poly f = '+str(ft))
        cmd = 'matrix c = coeffs ('+str(ft)+',t)'
        S.eval(cmd)
        N = int(S.eval('size(c)'))
        b = ["c["+str(i)+",1]," for i in range(2,N/2-4)]
        b = ''.join(b)
        b = b[:len(b)-1] #to cut off the trailing comma
        cmd = 'ideal I = '+b
        S.eval(cmd)
        c = S.eval('slimgb(I)')
        d = c.split("=")
        d = d[1:]
        d[len(d)-1] += "\n"
        e = [x[:x.index("\n")] for x in d]
        vals = []
        for x in e:
            for y in vars0:
                if str(y) in x:
                    if len(x.replace(str(y),"")) != 0:
                        i = x.find("-")
                        if i>0:
                            vals.append([eval(x[1:i]),x[:i],F(eval(x[i+1:]))])
                        i = x.find("+")
                        if i>0:
                            vals.append([eval(x[1:i]),x[:i],-F(eval(x[i+1:]))])
                    else:
                        vals.append([eval(str(y)[1:]),str(y),F(0)])
        vals.sort()
        k = len(vals)
        v = [x0+t,y0+add([vals[i][2]*t**(i+1) for i in range(k)])]
        return v

    def plot(self, *args, **kwds):
        """
        Plot the real points of an affine patch of this projective
        plane curve.


        INPUT:


        -  ``self`` - an affine plane curve

        -  ``patch`` - (optional) the affine patch to be plotted; if not
           specified, the patch corresponding to the last projective
           coordinate being nonzero

        -  ``*args`` - optional tuples (variable, minimum, maximum) for
           plotting dimensions

        -  ``**kwds`` - optional keyword arguments passed on to
           ``implicit_plot``


        EXAMPLES:

        A cuspidal curve::

            sage: R.<x, y, z> = QQ[]
            sage: C = Curve(x^3 - y^2*z)
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive

        The other affine patches of the same curve::

            sage: C.plot(patch=0)
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=1)
            Graphics object consisting of 1 graphics primitive

        An elliptic curve::

            sage: E = EllipticCurve('101a')
            sage: C = Curve(E)
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=0)
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=1)
            Graphics object consisting of 1 graphics primitive

        A hyperelliptic curve::

            sage: P.<x> = QQ[]
            sage: f = 4*x^5 - 30*x^3 + 45*x - 22
            sage: C = HyperellipticCurve(f)
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=0)
            Graphics object consisting of 1 graphics primitive
            sage: C.plot(patch=1)
            Graphics object consisting of 1 graphics primitive
        """
        # if user hasn't specified a favourite affine patch, take the
        # one avoiding "infinity", i.e. the one corresponding to the
        # last projective coordinate being nonzero
        patch = kwds.pop('patch', self.ngens() - 1)
        from constructor import Curve
        C = Curve(self.affine_patch(patch))
        return C.plot(*args, **kwds)

    def is_singular(self, P=None):
        r"""
        Return whether this curve is singular or not, or if a point ``P`` is provided,
        whether ``P`` is a singular point of this curve.

        INPUT:

        - ``P`` -- (default: None) a point on this curve.

        OUTPUT:

        - Boolean. If no point ``P`` is provided, returns True of False depending on whether
          this curve is singular or not. If a point ``P`` is provided, returns True or False
          depending on whether ``P`` is or is not a singular point of this curve.

        EXAMPLES:

        Over `\QQ`::

            sage: F = QQ
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^3-Y^2*Z)
            sage: C.is_singular()
            True

        Over a finite field::

            sage: F = GF(19)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^3+Y^3+Z^3)
            sage: C.is_singular()
            False
            sage: D = Curve(X^4-X*Z^3)
            sage: D.is_singular()
            True
            sage: E = Curve(X^5+19*Y^5+Z^5)
            sage: E.is_singular()
            True
            sage: E = Curve(X^5+9*Y^5+Z^5)
            sage: E.is_singular()
            False

        Over `\CC`::

            sage: F = CC
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X)
            sage: C.is_singular()
            False
            sage: D = Curve(Y^2*Z-X^3)
            sage: D.is_singular()
            True
            sage: E = Curve(Y^2*Z-X^3+Z^3)
            sage: E.is_singular()
            False

        Showing that ticket #12187 is fixed::

            sage: F.<X,Y,Z> = GF(2)[]
            sage: G = Curve(X^2+Y*Z)
            sage: G.is_singular()
            False

        ::

            sage: P.<x,y,z> = ProjectiveSpace(CC, 2)
            sage: C = Curve([y^4 - x^3*z], P)
            sage: Q = P([0,0,1])
            sage: C.is_singular()
            True
        """
        if P is None:
            poly = self.defining_polynomial()
            return poly.parent().ideal(poly.gradient()+[poly]).dimension() > 0
        else:
            return not self.is_smooth(P)

    def tangents(self, P):
            r"""
            Return the tangents of this projective plane curve at the point ``P``.

            These are found by homogenizing the tangents of an affine patch of this curve
            containing ``P``. The point ``P`` must be a point on this curve.

            INPUT:

            - ``P`` -- a point on this curve.

            OUTPUT:

            - a list of polynomials in the coordinate ring of the ambient space of this curve.

            EXAMPLES::

                sage: set_verbose(-1)
                sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
                sage: C = Curve([x^3*y + 2*x^2*y^2 + x*y^3 + x^3*z + 7*x^2*y*z + 14*x*y^2*z + 9*y^3*z], P)
                sage: Q = P([0,0,1])
                sage: C.tangents(Q)
                [x + 4.147899035704788?*y, x + (1.426050482147607? + 0.3689894074818041?*I)*y,
                x + (1.426050482147607? - 0.3689894074818041?*I)*y]

            ::

                sage: P.<x,y,z> = ProjectiveSpace(QQ,2)
                sage: C = P.curve([x^2*y^3*z^4 - y^6*z^3 - 4*x^2*y^4*z^3 - 4*x^4*y^2*z^3 + 3*y^7*z^2 +\
                10*x^2*y^5*z^2 + 9*x^4*y^3*z^2 + 5*x^6*y*z^2 - 3*y^8*z - 9*x^2*y^6*z - 11*x^4*y^4*z -\
                7*x^6*y^2*z - 2*x^8*z + y^9 + 2*x^2*y^7 + 3*x^4*y^5 + 4*x^6*y^3 + 2*x^8*y])
                sage: Q = P([0,1,1])
                sage: C.tangents(Q)
                [-y + z, 3*x^2 - y^2 + 2*y*z - z^2]

            ::

                sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
                sage: C = P.curve([z^3*x + y^4 - x^2*z^2])
                sage: Q = P([1,1,1])
                sage: C.tangents(Q)
                Traceback (most recent call last):
                ...
                TypeError: (=(1 : 1 : 1)) is not a point on (=Projective Plane Curve
                over Rational Field defined by y^4 - x^2*z^2 + x*z^3)
            """
            # Check whether P is a point on this curve
            try:
                P = self(P)
            except TypeError:
                raise TypeError("(=%s) is not a point on (=%s)"%(P,self))

            # Find an affine chart of the ambient space of self that contains P
            i = 0
            while(P[i] == 0):
                i = i + 1
            C = self.affine_patch(i)
            Q = list(P)
            t = Q.pop(i)
            L = C.tangents(C.ambient_space()([1/t*Q[j] for j in range(self.ambient_space().dimension_relative())]))
            R = self.ambient_space().coordinate_ring()
            H = Hom(C.ambient_space().coordinate_ring(), R)
            G = list(R.gens())
            x = G.pop(i)
            phi = H(G)
            return [phi(g).homogenize(x) for g in L]

    def is_ordinary_singularity(self, P):
        r"""
        Return whether the singular point ``P`` of this projective plane curve is an ordinary singularity.

        The point ``P`` is an ordinary singularity of this curve if it is a singular point, and
        if the tangents of this curve at ``P`` are distinct.

        INPUT:

        - ``P`` -- a point on this curve.

        OUTPUT:

        - Boolean. True or False depending on whether ``P`` is or is not an ordinary singularity of this
          curve, respectively. An error is raised if ``P`` is not a singular point of this curve.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^2*z^3 - x^5], P)
            sage: Q = P([0,0,1])
            sage: C.is_ordinary_singularity(Q)
            False

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 - 3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = P.curve([x^2*y^3*z^4 - y^6*z^3 - 4*x^2*y^4*z^3 - 4*x^4*y^2*z^3 + 3*y^7*z^2 + 10*x^2*y^5*z^2\
            + 9*x^4*y^3*z^2 + 5*x^6*y*z^2 - 3*y^8*z - 9*x^2*y^6*z - 11*x^4*y^4*z - 7*x^6*y^2*z - 2*x^8*z + y^9 +\
            2*x^2*y^7 + 3*x^4*y^5 + 4*x^6*y^3 + 2*x^8*y])
            sage: Q = P([0,1,1])
            sage: C.is_ordinary_singularity(Q)
            True

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = P.curve([z^5 - y^5 + x^5 + x*y^2*z^2])
            sage: Q = P([0,1,1])
            sage: C.is_ordinary_singularity(Q)
            Traceback (most recent call last):
            ...
            TypeError: (=(0 : 1 : 1)) is not a singular point of (=Projective Plane
            Curve over Rational Field defined by x^5 - y^5 + x*y^2*z^2 + z^5)
        """
        r = self.multiplicity(P)
        if r < 2:
            raise TypeError("(=%s) is not a singular point of (=%s)"%(P,self))

        T = self.tangents(P)

        # when there is a tangent of higher multiplicity
        if len(T) < r:
            return False

        # otherwise they are distinct
        return True

    def quadratic_transformation(self):
        r"""
        Return the quadratic transformation of this curve created from applying the standard Cremona
        transformation.

        The standard Cremona transformation is the birational automorphism of `\mathbb{P}^{2}` defined
        `(x : y : z)\mapsto (yz : xz : xy)`. The map is not defined at the points `(0 : 0 : 1)`, `(0 : 1 : 0)`,
        and `(1 : 0 : 0)`. The transformed curve is created by applying the map to this curve, and then dividing
        out by any common powers of `x,y,z` in the resulting defining polynomial.

        OUTPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve into its ambient space, and the
          quadratic transform of this curve. A restriction of the map defines a birational map between the curves.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve(x^3*y - z^4 - z^2*x^2, P)
            sage: C.quadratic_transformation()
            (Scheme morphism:
               From: Projective Plane Curve over Rational Field defined by x^3*y - x^2*z^2 - z^4
               To:   Projective Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to
                     (y*z : x*z : x*y),
             Projective Plane Curve over Rational Field defined by -x^3*y - x*y*z^2
            + z^4)
        """
        PP = self.ambient_space()
        R = PP.coordinate_ring()
        L = R.gens()
        coords = [L[1]*L[2], L[0]*L[2], L[0]*L[1]]
        H = Hom(self, PP)
        phi = H(coords)
        G = self.defining_polynomials()[0](coords)
        # divide out by any common powers of the generators of R
        degs = [G.degree()]*len(L)
        for F in G.monomials():
            for i in range(len(L)):
                if F.degree(L[i]) < degs[i]:
                    degs[i] = F.degree(L[i])
        T = []
        for item in G.dict().items():
            tup = tuple([item[0][i] - degs[i] for i in range(len(L))])
            T.append(tuple([tup, item[1]]))
        G = R(dict(T))
        from constructor import Curve
        return tuple([phi, Curve(G, PP)])

    def excellent_position(self, Q):
        r"""
        Return a transformation of this curve into one in excellent position with respect to the point ``Q``.

        Here excellent position is defined as in [Fulton89]_. A curve `C` of degree `d` containing the point
        `(0 : 0 : 1)` with multiplicity `r` is said to be in excellent position if none of the coordinate lines
        are tangent to `C` at any of the fundamental points `(1 : 0 : 0)`, `(0 : 1 : 0)`, and `(0 : 0 : 1)`, and
        if the two coordinate lines containing `(0 : 0 : 1)` intersect `C` transversally in `d - r` distinct
        non-fundamental points, and if the other coordinate line intersects `C` transversally at `d` distinct,
        non-fundamental points.

        INPUT:

        - ``Q`` -- a point on this curve.

        OUPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve into its ambient space, and
          the transformed curve that is the image of that morphism.

        EXAMPLES::

            sage: set_verbose(-1)
            sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
            sage: C = Curve([x*y - z^2], P)
            sage: Q = P([1,1,1])
            sage: C.excellent_position(Q) # long time (4 seconds)
            (Scheme morphism:
               From: Projective Plane Curve over Algebraic Field defined by x*y - z^2
               To:   Projective Space of dimension 2 over Algebraic Field
               Defn: Defined on coordinates by sending (x : y : z) to
                     (2*x - y - z : (-3)*x + y + 2*z : -x + y + z),
             Projective Plane Curve over Algebraic Field defined by (-5)*x^2 +
            (-5)*x*y - y^2 + (-4)*x*z + (-3)*y*z)

        REFERENCES:

        ..  [Fulton89] \W. Fulton. Algebraic curves: an introduction to algebraic geometry. Addison-Wesley,
            Redwood City CA (1989).
        """
        PP = self.ambient_space()
        # currently only implemented for when base ring is QQbar
        if not self.base_ring() == QQbar:
            raise NotImplementedError("base ring must be QQbar")
        # check that Q is on this curve
        try:
            Q = self(Q)
        except TypeError:
            raise TypeError("(=%s) must be a point on this curve"%Q)
        # first check that this curve is irreducible
        poly = self.defining_polynomial()
        coeff = poly.coefficients()
        from sage.rings.qqbar import number_field_elements_from_algebraics
        temp = number_field_elements_from_algebraics(coeff)
        polyK = PP.coordinate_ring().change_ring(base_ring=temp[0])
        poly = polyK(dict([(k,temp[1][coeff.index(v)]) for (k,v) in poly.dict().items()]))
        if not polyK.ideal([poly]).is_prime():
            raise TypeError("this curve must be irreducible")
        r = self.multiplicity(Q)
        # first move Q to (0 : 0 : 1), (1 : 0 : 0), or (0 : 1 : 0)
        i = 0
        while(Q[i] == 0):
            i = i + 1
        coords = [PP.gens()[j] + Q[j]/Q[i]*PP.gens()[i] for j in range(3)]
        coords[i] = PP.gens()[i]
        accoords = [PP.gens()[j] - Q[j]/Q[i]*PP.gens()[i] for j in range(3)] # coords used in map construction
        accoords[i] = PP.gens()[i]
        C = PP.curve(self.defining_polynomial()(coords))
        d = C.defining_polynomial().degree()
        P = [0]*3
        P[i] = 1
        P = C(P)
        l = [0,1,2]
        l.pop(i)
        u = PP.gens()[l[0]]
        v = PP.gens()[l[1]]
        # find the line to map to the line x
        a = 0
        good = False
        while (not good):
            a = a + 1
            L = PP.curve(u + a*v)
            # check that L and the tangents of C at P are distinct
            distinct = True
            for T in C.tangents(P):
                Lpt = PP([L.defining_polynomial().coefficient({PP.gens()[j]:1}) for j in range(3)])
                Tpt = PP([T.coefficient({PP.gens()[j]:1}) for j in range(3)])
                if Lpt == Tpt:
                    distinct = False
                    break
            if not distinct:
                continue
            # check that all intersections of C with L are transverse
            # and that there are the right number of intersections, (d - r) other than P
            pts = C.intersection(L).rational_points()
            X = C.intersection(L)
            pts = X.rational_points()
            pts.remove(X(P))
            if len(pts) != d - r:
                continue
            transverse = True
            for pt in pts:
                if not C.is_transverse(L, pt):
                    transverse = False
                    break
            if not transverse:
                continue
            good = True
        Lx = PP.curve(u + a*v)
        # find the line to map to the line y
        good = False
        b = a
        while (not good):
            b = b + 1
            L = PP.curve(u + b*v)
            # check that L and the tangents of C at P are distinct
            distinct = True
            for T in C.tangents(P):
                Lpt = PP([L.defining_polynomial().coefficient({PP.gens()[j]:1}) for j in range(3)])
                Tpt = PP([T.coefficient({PP.gens()[j]:1}) for j in range(3)])
                if Lpt == Tpt:
                    distinct = False
                    break
            if not distinct:
                continue
            # check that all intersections of C with L are transverse
            # and that there are the right number of intersections (d - r) other than P
            X = C.intersection(L)
            pts = X.rational_points()
            pts.remove(X(P))
            if len(pts) != d - r:
                continue
            transverse = True
            for pt in pts:
                if not C.is_transverse(L, pt):
                    transverse = False
                    break
            if not transverse:
                continue
            good = True
        Ly = PP.curve(u + b*v)
        # find the line to map to the line z
        # pick one point on Lx, and one on Ly, each not on this curve, and such that
        # the line between them does not contain P
        # note that since this curve is irreducible, it does not contain either
        # Lx, Ly, so only finitely many points of those lines will lie on this curve
        good = False
        t = 0
        while(not good):
            # construct points on Lx incrementally
            t = t + 1
            Px = [0]*3
            Px[l[1]] = t
            Px[l[0]] = -a*t
            Px[i] = 0
            Px = PP(Px)
            try:
                C(Px)
                continue
            except TypeError:
                pass
            good = True
        good = False
        t = 0
        while(not good):
            # construct points on Ly incrementally
            t = t + 1
            Py = [0]*3
            Py[l[1]] = t
            Py[l[0]] = -b*t
            Py[i] = 1
            Py = PP(Py)
            try:
                C(Py)
                continue
            except TypeError:
                pass
            Lz = PP.line([Px, Py])
            pts = C.intersection(Lz).rational_points()
            if len(pts) != d:
                continue
            transverse = True
            for pt in pts:
                if not C.is_transverse(Lz, pt):
                    transverse = False
                    break
            if not transverse:
                continue
            good = True
        # now create a change of coordinates sending Lx to the line x, Ly to the line y, and Lz to the line z
        # by construction, P, Px, Py are linearly independent so the following matrix is invertible
        M = matrix([[Py[j], Px[j], P[j]] for j in range(3)])
        # M defines a change of coordinates sending (1:0:0) to Py, (0:1:0) to Px, (0:0:1) to P; the
        # inverse of the transformation we want, used to create defining polynomial
        coords = [sum([M.row(j)[k]*PP.gens()[k] for k in range(3)]) for j in range(3)]
        # coords for map
        M = M.inverse()
        accoords2 = [sum([M.row(j)[k]*PP.gens()[k] for k in range(3)]) for j in range(3)]
        H = Hom(self, PP)
        phi = H([f(accoords) for f in accoords2])
        return tuple([phi, PP.curve(C.defining_polynomial()(coords))])

    def ordinary_model(self):
        r"""
        Return an ordinary plane curve model of this curve.

        OUPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve into its ambient space, and a
          plane curve with only ordinary singularities. A restriction of this map defines a birational map between
          the curves.

        EXAMPLES::

            sage: set_verbose(-1)
            sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
            sage: C = Curve([y^2*z - x^3], P)
            sage: C.ordinary_model() # long time (5 seconds)
            (Scheme morphism:
               From: Projective Plane Curve over Algebraic Field defined by -x^3 + y^2*z
               To:   Projective Space of dimension 2 over Algebraic Field
               Defn: Defined on coordinates by sending (x : y : z) to
                     (x^2 + 3*x*y + 2*y^2 + x*z + 2*y*z : -x^2 + (-2)*x*y - y^2 -
            x*z - y*z : -x^2 + (-3)*x*y + (-2)*y^2),
             Projective Plane Curve over Algebraic Field defined by x^3*y +
            2*x^2*y^2 + x*y^3 + x^3*z + 7*x^2*y*z + 14*x*y^2*z + 9*y^3*z)

        ::

            sage: set_verbose(-1)
            sage: P.<x,y,z> = ProjectiveSpace(QQbar, 2)
            sage: C = Curve([y^2*z^3 - x^5], P)
            sage: C.ordinary_model() # long time (8 seconds)
            (Scheme morphism:
               From: Projective Plane Curve over Algebraic Field defined by -x^5 + y^2*z^3
               To:   Projective Space of dimension 2 over Algebraic Field
               Defn: Defined on coordinates by sending (x : y : z) to
                     (x^2 + 3*x*y + 2*y^2 + x*z + 2*y*z : -x^2 + (-2)*x*y - y^2 -
            x*z - y*z : -x^2 + (-3)*x*y + (-2)*y^2),
             Projective Plane Curve over Algebraic Field defined by x^5*y^3 +
            2*x^4*y^4 + x^3*y^5 + 3*x^4*y^3*z + 6*x^3*y^4*z + 3*x^2*y^5*z +
            3*x^3*y^3*z^2 + 6*x^2*y^4*z^2 + 3*x*y^5*z^2 + x^5*z^3 + 10*x^4*y*z^3 +
            40*x^3*y^2*z^3 + 81*x^2*y^3*z^3 + 82*x*y^4*z^3 + 33*y^5*z^3)
        """
        C = self
        H = Hom(C, C.ambient_space())
        coords = [C.ambient_space().gens()[i] for i in range(3)]
        resolved = False
        while not resolved:
            resolved = True
            for Q in C.singular_points():
                if not C.is_ordinary_singularity(Q):
                    temp_exc = C.excellent_position(Q)
                    temp_qua = temp_exc[1].quadratic_transformation()
                    C = temp_qua[1]
                    coords = [f(temp_exc[0].defining_polynomials()) for f in temp_qua[0].defining_polynomials()]
                    repeat = False
                    break
        phi = H(coords)
        return tuple([phi, C])

    def is_transverse(self, C, P):
        r"""
        Return whether the intersection of this curve with the curve ``C`` at the point ``P`` is transverse.

        The intersection at ``P`` is transverse if ``P`` is a nonsingular point of both curves, and if the
        tangents of the curves at ``P`` are distinct.

        INPUT:

        - ``C`` -- a curve in the ambient space of this curve.

        - ``P`` -- a point in the intersection of both curves.

        OUPUT: Boolean.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([x^2 - y^2], P)
            sage: D = Curve([x - y], P)
            sage: Q = P([1,1,0])
            sage: C.is_transverse(D, Q)
            False

        ::

            sage: K = QuadraticField(-1)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = Curve([y^2*z - K.0*x^3], P)
            sage: D = Curve([z*x + y^2], P)
            sage: Q = P([0,0,1])
            sage: C.is_transverse(D, Q)
            False

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([x^2 - 2*y^2 - 2*z^2], P)
            sage: D = Curve([y - z], P)
            sage: Q = P([2,1,1])
            sage: C.is_transverse(D, Q)
            True
        """
        if not self.intersects_at(C, P):
            raise TypeError("(=%s) must be a point in the intersection of (=%s) and this curve"%(P,C))
        if self.is_singular(P) or C.is_singular(P):
            return False

        # there is only one tangent at a nonsingular point of a plane curve
        return not self.tangents(P)[0] == C.tangents(P)[0]

class ProjectivePlaneCurve_finite_field(ProjectivePlaneCurve):
    def rational_points_iterator(self):
        r"""
        Return a generator object for the rational points on this curve.

        INPUT:

        - ``self`` -- a projective curve

        OUTPUT:

        A generator of all the rational points on the curve defined over its base field.

        EXAMPLE::

            sage: F = GF(37)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^7+Y*X*Z^5*55+Y^7*12)
            sage: len(list(C.rational_points_iterator()))
            37

        ::

            sage: F = GF(2)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X*Y*Z)
            sage: a = C.rational_points_iterator()
            sage: next(a)
            (1 : 0 : 0)
            sage: next(a)
            (0 : 1 : 0)
            sage: next(a)
            (1 : 1 : 0)
            sage: next(a)
            (0 : 0 : 1)
            sage: next(a)
            (1 : 0 : 1)
            sage: next(a)
            (0 : 1 : 1)
            sage: next(a)
            Traceback (most recent call last):
            ...
            StopIteration

        ::

            sage: F = GF(3^2,'a')
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^3+5*Y^2*Z-33*X*Y*X)
            sage: b = C.rational_points_iterator()
            sage: next(b)
            (0 : 1 : 0)
            sage: next(b)
            (0 : 0 : 1)
            sage: next(b)
            (2*a + 2 : a : 1)
            sage: next(b)
            (2 : a + 1 : 1)
            sage: next(b)
            (a + 1 : 2*a + 1 : 1)
            sage: next(b)
            (1 : 2 : 1)
            sage: next(b)
            (2*a + 2 : 2*a : 1)
            sage: next(b)
            (2 : 2*a + 2 : 1)
            sage: next(b)
            (a + 1 : a + 2 : 1)
            sage: next(b)
            (1 : 1 : 1)
            sage: next(b)
            Traceback (most recent call last):
            ...
            StopIteration

        """
        g = self.defining_polynomial()
        K = g.parent().base_ring()
        from sage.rings.polynomial.all import PolynomialRing
        R = PolynomialRing(K,'X')
        X = R.gen()
        one = K.one()
        zero = K.zero()

        # the point with  Z = 0 = Y
        try:
            t = self.point([one,zero,zero])
            yield(t)
        except TypeError:
            pass

        # points with Z = 0, Y = 1
        g10 = R(g(X,one,zero))
        if g10.is_zero():
            for x in K:
                yield(self.point([x,one,zero]))
        else:
            for x in g10.roots(multiplicities=False):
                yield(self.point([x,one,zero]))

        # points with Z = 1
        for y in K:
            gy1 = R(g(X,y,one))
            if gy1.is_zero():
                for x in K:
                    yield(self.point([x,y,one]))
            else:
                for x in gy1.roots(multiplicities=False):
                    yield(self.point([x,y,one]))

    def rational_points(self, algorithm="enum", sort=True):
        r"""
        Return the rational points on this curve computed via enumeration.


        INPUT:

        - ``algorithm`` (string, default: 'enum') -- the algorithm to
          use.  Currently this is ignored.

        - ``sort`` (boolean, default ``True``) -- whether the output
          points should be sorted.  If False, the order of the output
          is non-deterministic.

        OUTPUT:

        A list of all the rational points on the curve defined over
        its base field, possibly sorted.

        .. note::

           This is a slow Python-level implementation.


        EXAMPLES::

            sage: F = GF(7)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^3+Y^3-Z^3)
            sage: C.rational_points()
            [(0 : 1 : 1), (0 : 2 : 1), (0 : 4 : 1), (1 : 0 : 1), (2 : 0 : 1), (3 : 1 : 0), (4 : 0 : 1), (5 : 1 : 0), (6 : 1 : 0)]


        ::

            sage: F = GF(1237)
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^7+7*Y^6*Z+Z^4*X^2*Y*89)
            sage: len(C.rational_points())
            1237

        ::

            sage: F = GF(2^6,'a')
            sage: P2.<X,Y,Z> = ProjectiveSpace(F,2)
            sage: C = Curve(X^5+11*X*Y*Z^3 + X^2*Y^3 - 13*Y^2*Z^3)
            sage: len(C.rational_points())
            104

        ::

            sage: R.<x,y,z> = GF(2)[]
            sage: f = x^3*y + y^3*z + x*z^3
            sage: C = Curve(f); pts = C.rational_points()
            sage: pts
            [(0 : 0 : 1), (0 : 1 : 0), (1 : 0 : 0)]

        """
        points = list(self.rational_points_iterator())
        if sort:
            points.sort()
        return points

class ProjectivePlaneCurve_prime_finite_field(ProjectivePlaneCurve_finite_field):
    def _points_via_singular(self, sort=True):
        r"""
        Return all rational points on this curve, computed using Singular's
        Brill-Noether implementation.

        INPUT:


        -  ``sort`` - bool (default: True), if True return the
           point list sorted. If False, returns the points in the order
           computed by Singular.


        EXAMPLE::

            sage: x, y, z = PolynomialRing(GF(5), 3, 'xyz').gens()
            sage: f = y^2*z^7 - x^9 - x*z^8
            sage: C = Curve(f); C
            Projective Plane Curve over Finite Field of size 5 defined by
            -x^9 + y^2*z^7 - x*z^8
            sage: C._points_via_singular()
            [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1),
             (3 : 1 : 1), (3 : 4 : 1)]
            sage: C._points_via_singular(sort=False)     #random
            [(0 : 1 : 0), (3 : 1 : 1), (3 : 4 : 1), (2 : 2 : 1),
             (0 : 0 : 1), (2 : 3 : 1)]


        .. note::

           The Brill-Noether package does not always work (i.e., the
           'bn' algorithm. When it fails a RuntimeError exception is
           raised.
        """
        f = self.defining_polynomial()._singular_()
        singular = f.parent()
        singular.lib('brnoeth')
        try:
            X1 = f.Adj_div()
        except (TypeError, RuntimeError) as s:
            raise RuntimeError(str(s) + "\n\n ** Unable to use the\
                                          Brill-Noether Singular package to\
                                          compute all points (see above).")

        X2 = singular.NSplaces(1, X1)
        R = X2[5][1][1]
        singular.set_ring(R)

        # We use sage_flattened_str_list since iterating through
        # the entire list through the sage/singular interface directly
        # would involve hundreds of calls to singular, and timing issues with
        # the expect interface could crop up.  Also, this is vastly
        # faster (and more robust).
        v = singular('POINTS').sage_flattened_str_list()
        pnts = [self(int(v[3*i]), int(v[3*i+1]), int(v[3*i+2]))
                for i in range(len(v)//3)]
        # singular always dehomogenizes with respect to the last variable
        # so if this variable divides the curve equation, we need to add
        # points at infinity
        F = self.defining_polynomial()
        z = F.parent().gens()[-1]
        if z.divides(F):
            pnts += [self(1,a,0) for a in self.base_ring()]
            pnts += [self(0,1,0)]
        # remove multiple points
        pnts = list(set(pnts))
        if sort:
            pnts.sort()
        return pnts

    def riemann_roch_basis(self, D):
        r"""
        Return a basis for the Riemann-Roch space corresponding to
        `D`.

        This uses Singular's Brill-Noether implementation.

        INPUT:

        -  ``D`` - a divisor

        OUTPUT:

        A list of function field elements that form a basis of the Riemann-Roch space

        EXAMPLE::

            sage: R.<x,y,z> = GF(2)[]
            sage: f = x^3*y + y^3*z + x*z^3
            sage: C = Curve(f); pts = C.rational_points()
            sage: D = C.divisor([ (4, pts[0]), (4, pts[2]) ])
            sage: C.riemann_roch_basis(D)
            [x/y, 1, z/y, z^2/y^2, z/x, z^2/(x*y)]

        ::

            sage: R.<x,y,z> = GF(5)[]
            sage: f = x^7 + y^7 + z^7
            sage: C = Curve(f); pts = C.rational_points()
            sage: D = C.divisor([ (3, pts[0]), (-1,pts[1]), (10, pts[5]) ])
            sage: C.riemann_roch_basis(D)
            [(-2*x + y)/(x + y), (-x + z)/(x + y)]


        .. NOTE::

            Currently this only works over prime field and divisors
            supported on rational points.
        """
        f = self.defining_polynomial()._singular_()
        singular = f.parent()
        singular.lib('brnoeth')
        try:
            X1 = f.Adj_div()
        except (TypeError, RuntimeError) as s:
            raise RuntimeError(str(s) + "\n\n ** Unable to use the Brill-Noether Singular package to compute all points (see above).")
        X2 = singular.NSplaces(1, X1)
        # retrieve list of all computed closed points (possibly of degree >1)
        v = X2[3].sage_flattened_str_list()    # We use sage_flattened_str_list since iterating through
                                               # the entire list through the sage/singular interface directly
                                               # would involve hundreds of calls to singular, and timing issues with
                                               # the expect interface could crop up.  Also, this is vastly
                                               # faster (and more robust).
        v = [ v[i].partition(',') for i in range(len(v)) ]
        pnts = [ ( int(v[i][0]), int(v[i][2])-1 ) for i in range(len(v))]
        # retrieve coordinates of rational points
        R = X2[5][1][1]
        singular.set_ring(R)
        v = singular('POINTS').sage_flattened_str_list()
        coords = [self(int(v[3*i]), int(v[3*i+1]), int(v[3*i+2])) for i in range(len(v)//3)]
        # build correct representation of D for singular
        Dsupport = D.support()
        Dcoeffs = []
        for x in pnts:
            if x[0] == 1:
                Dcoeffs.append(D.coefficient(coords[x[1]]))
            else:
                Dcoeffs.append(0)
        Dstr = str(tuple(Dcoeffs))
        G = singular(','.join([str(x) for x in Dcoeffs]), type='intvec')
        # call singular's brill noether routine and return
        T = X2[1][2]
        T.set_ring()
        LG = G.BrillNoether(X2)
        LG = [X.split(',\n') for X in LG.sage_structured_str_list()]
        x,y,z = self.ambient_space().coordinate_ring().gens()
        vars = {'x':x, 'y':y, 'z':z}
        V = [(sage_eval(a, vars)/sage_eval(b, vars)) for a, b in LG]
        return V

    def rational_points(self, algorithm="enum", sort=True):
        r"""
        INPUT:


        -  ``algorithm`` - string:

        -  ``'enum'`` - straightforward enumeration

        -  ``'bn'`` - via Singular's brnoeth package.


        EXAMPLE::

            sage: x, y, z = PolynomialRing(GF(5), 3, 'xyz').gens()
            sage: f = y^2*z^7 - x^9 - x*z^8
            sage: C = Curve(f); C
            Projective Plane Curve over Finite Field of size 5 defined by
            -x^9 + y^2*z^7 - x*z^8
            sage: C.rational_points()
            [(0 : 0 : 1), (0 : 1 : 0), (2 : 2 : 1), (2 : 3 : 1),
             (3 : 1 : 1), (3 : 4 : 1)]
            sage: C = Curve(x - y + z)
            sage: C.rational_points()
            [(0 : 1 : 1), (1 : 1 : 0), (1 : 2 : 1), (2 : 3 : 1),
             (3 : 4 : 1), (4 : 0 : 1)]
            sage: C = Curve(x*z+z^2)
            sage: C.rational_points('all')
            [(0 : 1 : 0), (1 : 0 : 0), (1 : 1 : 0), (2 : 1 : 0),
             (3 : 1 : 0), (4 : 0 : 1), (4 : 1 : 0), (4 : 1 : 1),
             (4 : 2 : 1), (4 : 3 : 1), (4 : 4 : 1)]

        .. note::

           The Brill-Noether package does not always work (i.e., the
           'bn' algorithm. When it fails a RuntimeError exception is
           raised.
        """
        if algorithm == "enum":

            return ProjectivePlaneCurve_finite_field.rational_points(self,
                                                                algorithm="enum",
                                                                sort=sort)

        elif algorithm == "bn":

            return self._points_via_singular(sort=sort)

        elif algorithm == "all":

            S_enum = self.rational_points(algorithm = "enum")
            S_bn = self.rational_points(algorithm = "bn")
            if S_enum != S_bn:
                raise RuntimeError("Bug in rational_points -- different\
                                     algorithms give different answers for\
                                     curve %s!"%self)
            return S_enum

        else:

            raise ValueError("No algorithm '%s' known"%algorithm)

def Hasse_bounds(q, genus=1):
    r"""
    Return the Hasse-Weil bounds for the cardinality of a nonsingular
    curve defined over `\GF{q}` of given ``genus``.

    INPUT:

    - ``q`` (int) -- a prime power

    - ``genus`` (int, default 1) -- a non-negative integer,

    OUTPUT:

    (tuple)  The Hasse bounds (lb,ub) for the cardinality of a curve of
    genus ``genus`` defined over `\GF{q}`.

    EXAMPLES::

        sage: Hasse_bounds(2)
        (1, 5)
        sage: Hasse_bounds(next_prime(10^30))
        (999999999999998000000000000058, 1000000000000002000000000000058)
    """
    if genus==1:
        rq = (4*q).isqrt()
    else:
        rq = (4*(genus**2)*q).isqrt()
    return (q+1-rq,q+1+rq)

# Fix pickles from changing class names and plane_curves folder name
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.schemes.plane_curves.projective_curve',
                           'ProjectiveCurve_generic', ProjectivePlaneCurve)
