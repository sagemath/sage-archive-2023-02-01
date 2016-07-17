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
from sage.categories.number_fields import NumberFields
from sage.categories.homset import Hom
from sage.interfaces.all import singular
from sage.matrix.constructor import matrix
from sage.misc.all import add, sage_eval
from sage.rings.all import degree_lowest_rational_function
from sage.rings.qqbar import (number_field_elements_from_algebraics,
                              QQbar)
from sage.rings.rational_field import is_RationalField
from sage.schemes.affine.affine_space import AffineSpace

from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective
from sage.schemes.projective.projective_space import (is_ProjectiveSpace,
                                                      ProjectiveSpace)

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

    def arithmetic_genus(self):
        r"""
        Return the arithmetic genus of this projective curve.

        This is the arithmetic genus `g_a(C)` as defined in [Hartshorne]_. If `P` is the
        Hilbert polynomial of the defining ideal of this curve, then the arithmetic genus
        of this curve is `1 - P(0)`. This curve must be irreducible.

        OUTPUT: Integer.

        EXAMPLES::

            sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
            sage: C = P.curve([w*z - x^2, w^2 + y^2 + z^2])
            sage: C.arithmetic_genus()
            1

        ::

            sage: P.<x,y,z,w,t> = ProjectiveSpace(GF(7), 4)
            sage: C = P.curve([t^3 - x*y*w, x^3 + y^3 + z^3, z - w])
            sage: C.arithmetic_genus()
            10
        """
        if not self.is_irreducible():
            raise TypeError("this curve must be irreducible")
        return 1 - self.defining_ideal().hilbert_polynomial()(0)

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
        Return the arithmetic genus of this projective curve.

        This is the arithmetic genus `g_a(C)` as defined in [Hartshorne]_. For a projective
        plane curve of degree `d`, this is simply `(d-1)(d-2)/2`. It need *not* equal
        the geometric genus (the genus of the normalization of the curve). This curve must be
        irreducible.

        OUTPUT: Integer.

        EXAMPLE::

            sage: x,y,z = PolynomialRing(GF(5), 3, 'xyz').gens()
            sage: C = Curve(y^2*z^7 - x^9 - x*z^8); C
            Projective Plane Curve over Finite Field of size 5 defined by -x^9 + y^2*z^7 - x*z^8
            sage: C.arithmetic_genus()
            28
            sage: C.genus()
            4

        ::

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^3*x - x^2*y*z - 7*z^4])
            sage: C.arithmetic_genus()
            3

        REFERENCES:

        ..  [Hartshorne] \R. Hartshorne. Algebraic Geometry. Springer-Verlag, New York, 1977.
        """
        if not self.is_irreducible():
            raise TypeError("this curve must be irreducible")
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

    def tangents(self, P, factor=True):
            r"""
            Return the tangents of this projective plane curve at the point ``P``.

            These are found by homogenizing the tangents of an affine patch of this curve
            containing ``P``. The point ``P`` must be a point on this curve.

            INPUT:

            - ``P`` -- a point on this curve.

            - ``factor`` -- (default: True) whether to attempt computing the polynomials of the individual tangent
              lines over the base field of this curve, or to just return the polynomial corresponding to the union
              of the tangent lines (which requires fewer computations).

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
                sage: C.tangents(Q, factor=False)
                [6*x^3 + 42*x^2*y + 84*x*y^2 + 54*y^3]

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
            PP = self.ambient_space()
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
            L = C.tangents(C.ambient_space()([1/t*Q[j] for j in range(PP.dimension_relative())]), factor)
            R = PP.coordinate_ring()
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

        # Find an affine chart of the ambient space of self that contains P
        i = 0
        while(P[i] == 0):
            i = i + 1
        C = self.affine_patch(i)
        Q = list(P)
        t = Q.pop(i)
        Q = [1/t*Q[j] for j in range(self.ambient_space().dimension_relative())]
        return C.is_ordinary_singularity(C.ambient_space()(Q))

    def quadratic_transformation(self):
        r"""
        Return the quadratic transformation of this curve created from applying the standard Cremona
        transformation.

        The standard Cremona transformation is the birational automorphism of `\mathbb{P}^{2}` defined
        `(x : y : z)\mapsto (yz : xz : xy)`. The map is not defined at the points `(0 : 0 : 1)`, `(0 : 1 : 0)`,
        and `(1 : 0 : 0)`. The transformed curve is created by applying the map to this curve, and then dividing
        out by any common powers of the variables `x,y,z` in the resulting defining polynomial.

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
        G = self.defining_polynomial()(coords)
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

            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([x*y - z^2], P)
            sage: Q = P([1,1,1])
            sage: C.excellent_position(Q)
            (Scheme morphism:
               From: Projective Plane Curve over Rational Field defined by x*y - z^2
               To:   Projective Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to
                     (-x + 1/2*y + 1/2*z : -1/2*y + 1/2*z : x + 1/2*y - 1/2*z),
             Projective Plane Curve over Rational Field defined by -x^2 - 3*x*y -
            4*y^2 - x*z - 3*y*z)

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 - 3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = P.curve([z^2*y^3*x^4 - y^6*x^3 - 4*z^2*y^4*x^3 - 4*z^4*y^2*x^3 + 3*y^7*x^2 + 10*z^2*y^5*x^2\
            + 9*z^4*y^3*x^2 + 5*z^6*y*x^2 - 3*y^8*x - 9*z^2*y^6*x - 11*z^4*y^4*x - 7*z^6*y^2*x - 2*z^8*x + y^9 +\
            2*z^2*y^7 + 3*z^4*y^5 + 4*z^6*y^3 + 2*z^8*y])
            sage: Q = P([1,0,0])
            sage: C.excellent_position(Q)
            (Scheme morphism:
               From: Projective Plane Curve over Number Field in b with defining
            polynomial a^2 - 3 defined by -x^3*y^6 + 3*x^2*y^7 - 3*x*y^8 + y^9 +
            x^4*y^3*z^2 - 4*x^3*y^4*z^2 + 10*x^2*y^5*z^2 - 9*x*y^6*z^2 + 2*y^7*z^2 -
            4*x^3*y^2*z^4 + 9*x^2*y^3*z^4 - 11*x*y^4*z^4 + 3*y^5*z^4 + 5*x^2*y*z^6 -
            7*x*y^2*z^6 + 4*y^3*z^6 - 2*x*z^8 + 2*y*z^8
               To:   Projective Space of dimension 2 over Number Field in b with
            defining polynomial a^2 - 3
               Defn: Defined on coordinates by sending (x : y : z) to
                     (1/4*y + 1/2*z : -1/4*y + 1/2*z : x + 1/4*y - 1/2*z),
             Projective Plane Curve over Number Field in b with defining polynomial
            a^2 - 3 defined by 900*x^9 - 7410*x^8*y + 29282*x^7*y^2 - 69710*x^6*y^3
            + 110818*x^5*y^4 - 123178*x^4*y^5 + 96550*x^3*y^6 - 52570*x^2*y^7 +
            18194*x*y^8 - 3388*y^9 - 1550*x^8*z + 9892*x^7*y*z - 30756*x^6*y^2*z +
            58692*x^5*y^3*z - 75600*x^4*y^4*z + 67916*x^3*y^5*z - 42364*x^2*y^6*z +
            16844*x*y^7*z - 3586*y^8*z + 786*x^7*z^2 - 3958*x^6*y*z^2 +
            9746*x^5*y^2*z^2 - 14694*x^4*y^3*z^2 + 15174*x^3*y^4*z^2 -
            10802*x^2*y^5*z^2 + 5014*x*y^6*z^2 - 1266*y^7*z^2 - 144*x^6*z^3 +
            512*x^5*y*z^3 - 912*x^4*y^2*z^3 + 1024*x^3*y^3*z^3 - 816*x^2*y^4*z^3 +
            512*x*y^5*z^3 - 176*y^6*z^3 + 8*x^5*z^4 - 8*x^4*y*z^4 - 16*x^3*y^2*z^4 +
            16*x^2*y^3*z^4 + 8*x*y^4*z^4 - 8*y^5*z^4)

        REFERENCES:

        ..  [Fulton89] \W. Fulton. Algebraic curves: an introduction to algebraic geometry. Addison-Wesley,
            Redwood City CA (1989).
        """
        PP = self.ambient_space()
        # check that Q is on this curve
        try:
            Q = self(Q)
        except TypeError:
            raise TypeError("(=%s) must be a point on this curve"%Q)
        r = self.multiplicity(Q)
        d = self.defining_polynomial().degree()
        # first move Q to (0 : 0 : 1), (1 : 0 : 0), or (0 : 1 : 0)
        # this makes it easier to construct the main transformation
        i = 0
        while Q[i] == 0:
            i = i + 1
        coords = [PP.gens()[j] + Q[j]/Q[i]*PP.gens()[i] for j in range(3)]
        coords[i] = PP.gens()[i]
        accoords = [PP.gens()[j] - Q[j]/Q[i]*PP.gens()[i] for j in range(3)] # coords used in map construction
        accoords[i] = PP.gens()[i]
        baseC = PP.curve(self.defining_polynomial()(coords))
        P = [0]*3
        P[i] = 1
        P = PP(P)
        l = [0,1,2]
        l.pop(i)
        # choose points forming a triangle with one vertex at P to map to the coordinate triangle
        good = False
        a = 0
        while not good:
            a = a + 1
            # find points to map to (1 : 0 : 0) and (0 : 1 : 0), not on the curve
            Px = [0]*3
            Px[l[0]] = a
            Px[l[1]] = 1
            Py = [0]*3
            Py[l[0]] = -a
            Py[l[1]] = 1
            Py[i] = 1
            try:
                Px = baseC(Px)
                Py = baseC(Py)
                continue
            except TypeError:
                pass
            # by construction, P, Px, Py are linearly independent so the following matrix is invertible
            M = matrix([[Px[j], Py[j], P[j]] for j in range(3)])
            # M defines a change of coordinates sending (1 : 0 : 0) to Py, (0 : 1 : 0) to Px, (0 : 0 : 1) to P; the
            # inverse of the transformation we want, used to create the new defining polynomial
            coords = [sum([M.row(j)[k]*PP.gens()[k] for k in range(3)]) for j in range(3)]
            C = PP.curve(baseC.defining_polynomial()(coords))
            # check tangents at (0 : 0 : 1)
            T = C.tangents(PP([0,0,1]), factor=False)[0]
            if all([g.degree(PP.gens()[0]) > 0 for g in T.monomials()]):
                continue
            if all([g.degree(PP.gens()[1]) > 0 for g in T.monomials()]):
                continue
            # check that the other intersections of C with the exceptional lines are correct
            need_continue = False
            for j in range(3):
                poly = C.defining_polynomial().subs({PP.gens()[j]: 0})
                # this is a homogeneous polynomial in the other two variables
                # and so should factor completely into homogeneous linear factors
                # each corresponding to an intersection point where the jth coord is 0.
                # check if there are enough roots, up to multiplicity (that is, that PP.gens()[j]
                # doesn't divide the defining polynomial of C)
                if poly.degree() != d:
                    need_continue = True
                    break
                # if j != 2, then there should be d - r multiplicity 1 roots,
                # besides the root corresponding to (0 : 0 : 1)
                # if j == 2, then all roots should have multiplicity 1
                npoly = poly
                if j != 2:
                    # since (0 : 0 : 1) has multiplicity r, divide out by the highest
                    # shared power of the corresponding variable before doing the resultant computations
                    if j == 0:
                        while PP.gens()[1].divides(npoly):
                            npoly = PP.coordinate_ring()(npoly/PP.gens()[1])
                    else:
                        while PP.gens()[0].divides(npoly):
                            npoly = PP.coordinate_ring()(npoly/PP.gens()[0])
                    # check the degree again
                    if npoly.degree() != d - r:
                        need_continue = True
                        break
                    # check that npoly isn't a constant now
                    if npoly.degree() > 0:
                        t = 0
                        while npoly.degree(PP.gens()[t]) == 0:
                            t = t + 1
                        if npoly.resultant(npoly.derivative(PP.gens()[t]), PP.gens()[t]) == 0:
                            need_continue = True
                            break
                else:
                    t = 0
                    while npoly.degree(PP.gens()[t]) == 0:
                        t = t + 1
                    if poly.resultant(poly.derivative(PP.gens()[t]), PP.gens()[t]) == 0:
                        need_continue = True
                        break
                # check that intersections with the line PP.gens()[j] are transverse.
                # at a simple point P of the curve, the tangent at that point is
                # given by F_x(P)*x + F_y(P)*y + F_z(P)*z where F is the defining polynomial
                # of the curve
                tmp_l = [0,1,2]
                tmp_l.pop(j)
                poly1 = npoly.derivative(PP.gens()[tmp_l[0]])
                poly2 = npoly.derivative(PP.gens()[tmp_l[1]])
                if poly1.degree() > 0 or poly2.degree() > 0:
                    t = 0
                    while poly1.degree(PP.gens()[t]) == 0 and poly2.degree(PP.gens()[t]) == 0:
                        t = t + 1
                    # maybe a stricter check than necessary
                    if poly1.resultant(poly2, PP.gens()[t]) == 0:
                        need_continue = True
                        break
            if need_continue:
                continue
            good = True
            # coords for map
            M = M.inverse()
            accoords2 = [sum([M.row(j)[k]*PP.gens()[k] for k in range(3)]) for j in range(3)]
            H = Hom(self, PP)
            phi = H([f(accoords) for f in accoords2])
        return tuple([phi, C])

    def ordinary_model(self, sing=None):
        r"""
        Return an ordinary plane curve model of this curve.

        Currently only implemented over number fields. If not all of the coordinates of the non-ordinary
        singularities of this curve are contained in its base field, then the curve returned will be
        defined over an extension.

        INPUT:

        - ``sing`` -- (default: None) the set of singular points of this curve. If not given, this is constructed.
          For higher degree curves construction can be expensive. This curve must be irreducible.

        OUPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve into its ambient space, and a
          plane curve with only ordinary singularities. A restriction of this map defines a birational map between
          the curves.

        EXAMPLES::

            sage: set_verbose(-1)
            sage: K = QuadraticField(3)
            sage: P.<x,y,z> = ProjectiveSpace(K, 2)
            sage: C = Curve([x^5 - K.0*y*z^4], P)
            sage: C.ordinary_model()
            (Scheme morphism:
               From: Projective Plane Curve over Number Field in a with defining
            polynomial x^2 - 3 defined by x^5 + (-a)*y*z^4
               To:   Projective Space of dimension 2 over Number Field in a with
            defining polynomial x^2 - 3
               Defn: Defined on coordinates by sending (x : y : z) to
                     (-1/4*x^2 - 1/2*x*y + 1/2*x*z + 1/2*y*z - 1/4*z^2 : 1/4*x^2 +
            1/2*x*y + 1/2*y*z - 1/4*z^2 : -1/4*x^2 + 1/4*z^2),
             Projective Plane Curve over Number Field in a with defining polynomial
            x^2 - 3 defined by (-a)*x^5*y + (-4*a)*x^4*y^2 + (-6*a)*x^3*y^3 +
            (-4*a)*x^2*y^4 + (-a)*x*y^5 + (-a - 1)*x^5*z + (-4*a + 5)*x^4*y*z +
            (-6*a - 10)*x^3*y^2*z + (-4*a + 10)*x^2*y^3*z + (-a - 5)*x*y^4*z +
            y^5*z)

        ::

            sage: set_verbose(-1)
            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y^2*z^2 - x^4 - x^3*z], P)
            sage: C.ordinary_model([P([0,0,1]), P([0,1,0])]) # long time (1 second)
            (Scheme morphism:
               From: Projective Plane Curve over Rational Field defined by -x^4 -
            x^3*z + y^2*z^2
               To:   Projective Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x : y : z) to
                     (-3/64*x^4 + 9/64*x^2*y^2 - 3/32*x*y^3 - 1/16*x^3*z -
            1/8*x^2*y*z + 1/4*x*y^2*z - 1/16*y^3*z - 1/8*x*y*z^2 + 1/16*y^2*z^2 :
            -1/64*x^4 + 3/64*x^2*y^2 - 1/32*x*y^3 + 1/16*x*y^2*z - 1/16*y^3*z +
            1/16*y^2*z^2 : 3/64*x^4 - 3/32*x^3*y + 3/64*x^2*y^2 + 1/16*x^3*z -
            3/16*x^2*y*z + 1/8*x*y^2*z - 1/8*x*y*z^2 + 1/16*y^2*z^2),
             Projective Plane Curve over Rational Field defined by 4*x^6*y^3 -
            24*x^5*y^4 + 36*x^4*y^5 + 8*x^6*y^2*z - 40*x^5*y^3*z + 24*x^4*y^4*z +
            72*x^3*y^5*z - 4*x^6*y*z^2 + 8*x^5*y^2*z^2 - 56*x^4*y^3*z^2 +
            104*x^3*y^4*z^2 + 44*x^2*y^5*z^2 + 8*x^6*z^3 - 16*x^5*y*z^3 -
            24*x^4*y^2*z^3 + 40*x^3*y^3*z^3 + 48*x^2*y^4*z^3 + 8*x*y^5*z^3 -
            8*x^5*z^4 + 36*x^4*y*z^4 - 56*x^3*y^2*z^4 + 20*x^2*y^3*z^4 +
            40*x*y^4*z^4 - 16*y^5*z^4)

        ::

            sage: set_verbose(-1)
            sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([(x^2 + y^2 - y*z - 2*z^2)*(y*z - x^2 + 2*z^2) + y^4], P)
            sage: C.ordinary_model() # long time (3 seconds)
            (Scheme morphism:
               From: Projective Plane Curve over Number Field in a with defining
            polynomial y^2 - 2 defined by -x^4 - x^2*y^2 + y^4 + 2*x^2*y*z + y^3*z +
            4*x^2*z^2 + y^2*z^2 - 4*y*z^3 - 4*z^4
               To:   Projective Space of dimension 2 over Number Field in a with
            defining polynomial y^2 - 2
               Defn: Defined on coordinates by sending (x : y : z) to
                     ((-5/128*a - 5/128)*x^4 + (-5/32*a + 5/32)*x^3*y + (-1/16*a +
            3/32)*x^2*y^2 + (1/16*a - 1/16)*x*y^3 + (1/32*a - 1/32)*y^4 - 1/32*x^3*z
            + (3/16*a - 5/8)*x^2*y*z + (1/8*a - 5/16)*x*y^2*z + (1/8*a +
            5/32)*x^2*z^2 + (-3/16*a + 5/16)*x*y*z^2 + (-3/16*a - 1/16)*y^2*z^2 +
            1/16*x*z^3 + (1/4*a + 1/4)*y*z^3 + (-3/32*a - 5/32)*z^4 : (-5/128*a -
            5/128)*x^4 + (5/32*a)*x^3*y + (3/32*a + 3/32)*x^2*y^2 + (-1/16*a)*x*y^3
            + (-1/32*a - 1/32)*y^4 - 1/32*x^3*z + (-11/32*a)*x^2*y*z + (1/8*a +
            5/16)*x*y^2*z + (3/16*a + 1/4)*y^3*z + (1/8*a + 5/32)*x^2*z^2 + (-1/16*a
            - 3/8)*x*y*z^2 + (-3/8*a - 9/16)*y^2*z^2 + 1/16*x*z^3 + (5/16*a +
            1/2)*y*z^3 + (-3/32*a - 5/32)*z^4 : (1/64*a + 3/128)*x^4 + (-1/32*a -
            1/32)*x^3*y + (3/32*a - 9/32)*x^2*y^2 + (1/16*a - 3/16)*x*y^3 - 1/32*y^4
            + (3/32*a + 1/8)*x^2*y*z + (-1/8*a + 1/8)*x*y^2*z + (-1/16*a)*y^3*z +
            (-1/16*a - 3/32)*x^2*z^2 + (1/16*a + 1/16)*x*y*z^2 + (3/16*a +
            3/16)*y^2*z^2 + (-3/16*a - 1/4)*y*z^3 + (1/16*a + 3/32)*z^4),
             Projective Plane Curve over Number Field in a with defining polynomial
            y^2 - 2 defined by (30*a - 59)*x^6*y^4 + (148*a - 138)*x^5*y^5 + (78*a -
            187)*x^4*y^6 + (30*a - 141)*x^6*y^3*z + (432*a - 297)*x^5*y^4*z + (422*a
            - 743)*x^4*y^5*z + (140*a - 355)*x^3*y^6*z + (-25*a - 149)*x^6*y^2*z^2 +
            (413*a - 203)*x^5*y^3*z^2 + (673*a - 817)*x^4*y^4*z^2 + (395*a -
            905)*x^3*y^5*z^2 + (80*a - 230)*x^2*y^6*z^2 + (-30*a - 71)*x^6*y*z^3 +
            (127*a - 87)*x^5*y^2*z^3 + (422*a - 273)*x^4*y^3*z^3 + (350*a -
            665)*x^3*y^4*z^3 + (140*a - 380)*x^2*y^5*z^3 + (15*a - 60)*x*y^6*z^3 +
            (-5*a - 9)*x^6*z^4 + (-2*a - 33)*x^5*y*z^4 + (93*a - 7)*x^4*y^2*z^4 +
            (95*a - 135)*x^3*y^3*z^4 + (60*a - 145)*x^2*y^4*z^4 + (15*a -
            50)*x*y^5*z^4 - 5*y^6*z^4)
        """
        if not self.base_ring() in NumberFields():
            raise NotImplementedError("the base ring of this curve must be a number field")
        if not self.is_irreducible():
            raise TypeError("this curve must be irreducible")
        C_orig = self
        C = self
        PP = C.ambient_space()
        coords = [C.ambient_space().gens()[i] for i in range(3)]
        if sing is None:
            pts = C.singular_subscheme().change_ring(C.base_ring().embeddings(QQbar)[0]).rational_points()
            # if any points have coordinates not in the base field of this curve,
            # extend the base field and repeat
            L = [pt[j] for j in range(3) for pt in pts]
            for el in L:
                temp = number_field_elements_from_algebraics(el)
                if is_RationalField(temp[0]):
                    continue
                try:
                    tempF = C.base_ring().extension(temp[0].defining_polynomial(), C.base_ring().variable_name()\
                                                    + 'z')
                    K = tempF.absolute_field(temp[0].variable_name())
                    PP = PP.change_ring(K)
                    psi = C.base_ring().embeddings(K)[0]
                    C = PP.curve(PP.coordinate_ring()(C.defining_polynomial().map_coefficients(psi)))
                    C_orig = PP.curve(PP.coordinate_ring()(C_orig.defining_polynomial().map_coefficients(psi)))
                except ValueError:
                    pass
            pts = C.singular_points()
        else:
            pts = list(sing)
        while len(pts) > 0:
            for i in range(len(pts) - 1, -1, -1):
                try:
                    if C.is_ordinary_singularity(pts[i]):
                        pts.pop(i)
                except TypeError:
                    pts.pop(i)
            if len(pts) > 0:
                temp_exc = C.excellent_position(pts[0])
                temp_qua = temp_exc[1].quadratic_transformation()
                C = temp_qua[1]
                tmp_coords = [temp_exc[0].defining_polynomials()[i](coords) for i in range(3)]
                coords = [temp_qua[0].defining_polynomials()[i](tmp_coords) for i in range(3)]
                # transform the old points
                for i in range(len(pts) - 1, -1, -1):
                    # find image if it is a point the composition map is defined on
                    try:
                        temp_pt = temp_qua[0](temp_exc[0](pts[i]))
                        pts.pop(i)
                        if not PP(list(temp_pt)) in [PP(list(tpt)) for tpt in pts]:
                            pts.append(temp_pt)
                    except (TypeError, ValueError):
                        pass
                # add points from the intersection of C and the line z
                PPline = ProjectiveSpace(PP.base_ring(), 1)
                # make sure the conversion happens in the right order
                ringH = Hom(PP.coordinate_ring(), PPline.coordinate_ring())
                psi = ringH(list(PPline.gens()) + [0])
                X = PPline.subscheme([psi(f) for f in C.singular_subscheme().defining_polynomials()])
                qqbar_pts = X.change_ring(X.base_ring().embeddings(QQbar)[0]).rational_points()
                L = [pt[j] for j in range(2) for pt in qqbar_pts]
                for el in L:
                    temp = number_field_elements_from_algebraics(el)
                    if is_RationalField(temp[0]):
                        continue
                    try:
                        tempF = C.base_ring().extension(temp[0].defining_polynomial(), C.base_ring().variable_name()\
                                                        + 'z')
                        K = tempF.absolute_field(temp[0].variable_name())
                        PP = PP.change_ring(K)
                        psi = C.base_ring().embeddings(K)[0]
                        C = PP.curve(PP.coordinate_ring()(C.defining_polynomial().map_coefficients(psi)))
                        C_orig = PP.curve(PP.coordinate_ring()(C_orig.defining_polynomial().map_coefficients(psi)))
                    except ValueError:
                        pass
                # convert points
                psi = PPline.base_ring().embeddings(PP.base_ring())[0]
                PPline.change_ring(PP.base_ring())
                X = X.change_ring(psi)
                pts = [PP(pt.change_ring(psi)) for pt in pts]
                coords = [PP.coordinate_ring()(f.map_coefficients(psi)) for f in coords]
                newpts = [PP(list(pt) + [0]) for pt in X.rational_points()]
                for pt in newpts:
                    if not PP(list(pt)) in [PP(list(tpt)) for tpt in pts]:
                        pts.append(pt)
        H = Hom(C_orig, PP)
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
