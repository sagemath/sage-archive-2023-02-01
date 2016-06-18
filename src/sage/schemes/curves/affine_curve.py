"""
Affine curves.

EXAMPLES:

We can construct curves in either an affine plane::

    sage: A.<x,y> = AffineSpace(QQ, 2)
    sage: C = Curve([y - x^2], A); C
    Affine Plane Curve over Rational Field defined by -x^2 + y

or in higher dimensional affine space::

    sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
    sage: C = Curve([y - x^2, z - w^3, w - y^4], A); C
    Affine Curve over Rational Field defined by -x^2 + y, -w^3 + z, -y^4 + w

AUTHORS:

- William Stein (2005-11-13)

- David Joyner (2005-11-13)

- David Kohel (2006-01)
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

from sage.misc.all import add

from sage.rings.all import degree_lowest_rational_function

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.schemes.affine.affine_space import is_AffineSpace

from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_affine

from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.projective.projective_space import ProjectiveSpace

from curve import Curve_generic

class AffineCurve(Curve_generic, AlgebraicScheme_subscheme_affine):
    def _repr_type(self):
        r"""
        Return a string representation of the type of this curve.

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = Curve([x - y, z - w, w - x], A)
            sage: C._repr_type()
            'Affine'
        """
        return "Affine"

    def __init__(self, A, X):
        r"""
        Initialization function.

        EXAMPLES::

            sage: R.<v> = QQ[]
            sage: K.<u> = NumberField(v^2 + 3)
            sage: A.<x,y,z> = AffineSpace(K, 3)
            sage: C = Curve([z - u*x^2, y^2], A); C
            Affine Curve over Number Field in u with defining polynomial v^2 + 3
            defined by (-u)*x^2 + z, y^2

        ::

            sage: A.<x,y,z> = AffineSpace(GF(7), 3)
            sage: C = Curve([x^2 - z, z - 8*x], A); C
            Affine Curve over Finite Field of size 7 defined by x^2 - z, -x + z
        """
        if not is_AffineSpace(A):
            raise TypeError("A (=%s) must be an affine space"%A)
        Curve_generic.__init__(self, A, X)
        d = self.dimension()
        if d != 1:
            raise ValueError("defining equations (=%s) define a scheme of dimension %s != 1"%(X,d))

    def projective_closure(self, i=0, PP=None):
        r"""
        Return the projective closure of this affine curve.

        INPUT:

        - ``i`` -- (default: 0) the index of the affine coordinate chart of the projective space that the affine
          ambient space of this curve embeds into.

        - ``PP`` -- (default: None) ambient projective space to compute the projective closure in. This is
          constructed if it is not given.

        OUTPUT:

        - a curve in projective space.

        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: C = Curve([y-x^2,z-x^3], A)
            sage: C.projective_closure()
            Projective Curve over Rational Field defined by x1^2 - x0*x2,
            x1*x2 - x0*x3, x2^2 - x1*x3

        ::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: C = Curve([y - x^2, z - x^3], A)
            sage: C.projective_closure()
            Projective Curve over Rational Field defined by
            x1^2 - x0*x2, x1*x2 - x0*x3, x2^2 - x1*x3

        ::

            sage: A.<x,y> = AffineSpace(CC, 2)
            sage: C = Curve(y - x^3 + x - 1, A)
            sage: C.projective_closure(1)
            Projective Plane Curve over Complex Field with 53 bits of precision defined by
            x0^3 - x0*x1^2 + x1^3 - x1^2*x2

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: P.<u,v,w> = ProjectiveSpace(QQ, 2)
            sage: C = Curve([y - x^2], A)
            sage: C.projective_closure(1, P).ambient_space() == P
            True
        """
        from constructor import Curve
        return Curve(AlgebraicScheme_subscheme_affine.projective_closure(self, i, PP))

    def projection(self, indices, newvariables=True):
        r"""
        Return the projection of this curve onto the coordinates specified by ``indices``.

        INPUT:

        - ``indices`` -- a list or tuple of distinct integers specifying the indices of the coordinates to use
          in the projection. Can also be a list or tuple consisting of variables of the coordinate ring of the
          ambient space of this curve. If integers are used to specify the coordinates, 0 denotes the first
          coordinate.

        - ``newvariables`` -- Boolean, (default: True). If True, uses new variables when creating the ambient
          space to project into. If False, the variables of the coordinate ring of the ambient space of this curve
          are used.

        OUTPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve to affine space of dimension
          equal to the number of coordinates specified in ``indices``, and the affine subscheme that is the image
          of that morphism. If the image is a curve, the second element of the tuple will be a curve.

        EXAMPLES::

            sage: A.<x,y,z> = AffineSpace(QQ, 3)
            sage: C = Curve([y^7 - x^2 + x^3 - 2*z, z^2 - x^7 - y^2], A)
            sage: C.projection([0,1])
            (Scheme morphism:
               From: Affine Curve over Rational Field defined by y^7 + x^3 - x^2 -
            2*z, -x^7 - y^2 + z^2
               To:   Affine Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x, y, z) to
                     (x, y),
             Affine Plane Curve over Rational Field defined by x1^14 + 2*x0^3*x1^7 -
            2*x0^2*x1^7 - 4*x0^7 + x0^6 - 2*x0^5 + x0^4 - 4*x1^2)

        ::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = Curve([x - 2, y - 3, z - 1], A)
            sage: C.projection([0,1,2])
            (Scheme morphism:
               From: Affine Curve over Rational Field defined by x - 2, y - 3, z - 1
               To:   Affine Space of dimension 3 over Rational Field
               Defn: Defined on coordinates by sending (x, y, z, w) to
                     (x, y, z),
             Closed subscheme of Affine Space of dimension 3 over Rational Field
            defined by:
               x2 - 1,
               x1 - 3,
               x0 - 2)

        ::

            sage: A.<x,y,z,w,u> = AffineSpace(GF(11), 5)
            sage: C = Curve([x^3 - 5*y*z + u^2, x - y^2 + 3*z^2, w^2 + 2*u^3*y, y - u^2 + z*x], A)
            sage: C.projection([1,2,4], False)
            (Scheme morphism:
               From: Affine Curve over Finite Field of size 11 defined by x^3 -
            5*y*z + u^2, -y^2 + 3*z^2 + x, 2*y*u^3 + w^2, x*z - u^2 + y
               To:   Affine Space of dimension 3 over Finite Field of size 11
               Defn: Defined on coordinates by sending (x, y, z, w, u) to
                     (y, z, u),
             Affine Curve over Finite Field of size 11 defined by y^2*z - 3*z^3 -
            u^2 + y, u^6 - 5*y*z^4 + z^3*u^2 - 3*y*u^4 + 3*y^2*u^2 - y^3, y^2*u^4 -
            3*z^2*u^4 - 2*y^3*u^2 - 5*y*z^2*u^2 + y^4 - 5*y*z^3 + 2*z^4 + z^2*u^2 -
            3*z*u^2 + 3*y*z, y^4*u^2 + 2*z^4*u^2 - y^5 - 2*y*z^4 + 5*z*u^4 + y*z*u^2
            - 5*y*z^2 + 4*z^3 + z*u^2 + 5*u^2 - 5*y, y^6 - 5*z^6 - 5*z^3*u^2 +
            5*y*z^3 + 2*u^4 - 4*y*u^2 + 2*y^2 - 5*y*z + u^2)

        ::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = A.curve([x*y - z^3, x*z - w^3, w^2 - x^3])
            sage: C.projection([y,z], False)
            (Scheme morphism:
               From: Affine Curve over Rational Field defined by -z^3 + x*y, -w^3 +
            x*z, -x^3 + w^2
               To:   Affine Space of dimension 2 over Rational Field
               Defn: Defined on coordinates by sending (x, y, z, w) to
                     (y, z),
             Affine Plane Curve over Rational Field defined by z^23 - y^7*z^4)
        """
        AA = self.ambient_space()
        n = AA.dimension_relative()
        if n == 2:
            raise TypeError("this curve is already a plane curve")
        if self.base_ring() not in Fields():
            raise TypeError("this curve must be defined over a field")
        if len(indices) < 2 or len(indices) > AA.dimension_relative() or len(set(indices)) < len(indices):
            raise ValueError("(=%s) must be a list or tuple of distinct indices or variables"%indices)
        indices = list(indices)
        if all([f in AA.gens() for f in indices]):
            indices = [AA.gens().index(f) for f in indices]
            indices.sort()
        else:
            indices = [int(i) for i in indices] # type checking
            indices.sort()
            if indices[0] < 0 or indices[-1] > n - 1:
                raise ValueError("index values must be between 0 and one minus the dimension of the ambient space " \
                                 "of this curve")
        # construct the projection map
        if newvariables:
            AA2 = AffineSpace(self.base_ring(), len(indices))
        else:
            AA2 = AffineSpace(self.base_ring(), len(indices), [AA.gens()[i] for i in indices])
        H = Hom(self, AA2)
        psi = H([AA.gens()[i] for i in indices])
        # compute the image via elimination
        removecoords = list(AA.gens())
        for i in range(len(indices) - 1, -1, -1):
            removecoords.pop(indices[i])
        J = self.defining_ideal().elimination_ideal(removecoords)
        K = Hom(AA.coordinate_ring(), AA2.coordinate_ring())
        l = [0]*(n)
        for i in range(len(indices)):
            l[indices[i]] = AA2.gens()[i]
        phi = K(l)
        G = [phi(f) for f in J.gens()]
        from constructor import Curve
        try:
            C = Curve(G, AA2)
        except (TypeError, ValueError):
            C = AA2.subscheme(G)
        return tuple([psi, C])

    def plane_projection(self, newvariables=True):
        r"""
        Return a projection of this curve into an affine plane so that the image of the projection is
        a plane curve.

        INPUT:

        - ``newvariables`` -- Boolean, (default: True). If True, uses new variables when creating the ambient
          space to project into. If False, the variables of the coordinate ring of the ambient space of this curve
          are used.

        OUTPUT:

        - a tuple consisting of two elements: a scheme morphism from this curve into an affine plane, and the plane
          curve that defines the image of that morphism.

        EXAMPLES::

            sage: A.<x,y,z,w> = AffineSpace(QQ, 4)
            sage: C = Curve([x^2 - y*z*w, z^3 - w, w + x*y - 3*z^3], A)
            sage: C.plane_projection()
            (Scheme morphism:
              From: Affine Curve over Rational Field defined by -y*z*w + x^2, z^3 -
            w, -3*z^3 + x*y + w
              To:   Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y, z, w) to
                    (x, y), Affine Plane Curve over Rational Field defined by
            x0^2*x1^7 - 16*x0^4)

        ::

            sage: R.<a> = QQ[]
            sage: K.<b> = NumberField(a^2 + 2)
            sage: A.<x,y,z> = AffineSpace(K, 3)
            sage: C = A.curve([x - b, y - 2])
            sage: C.plane_projection(False)
            (Scheme morphism:
               From: Affine Curve over Number Field in b with defining polynomial
            a^2 + 2 defined by x + (-b), y - 2
               To:   Affine Space of dimension 2 over Number Field in b with
            defining polynomial a^2 + 2
               Defn: Defined on coordinates by sending (x, y, z) to
                     (x, z),
             Affine Plane Curve over Number Field in b with defining polynomial a^2
            + 2 defined by x + (-b))
        """
        n = self.ambient_space().dimension_relative()
        # finds a projection that will have a plane curve as its image
        # the following iterates over all pairs (i,j) with j > i to test all
        # possible projections
        for i in range(0,n-1):
            for j in range(i + 1, n):
                L = self.projection([i,j], newvariables)
                if isinstance(L[1], Curve_generic):
                    return L

class AffinePlaneCurve(AffineCurve):
    def __init__(self, A, f):
        r"""
        Initialization function.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([x^3 - y^2], A); C
            Affine Plane Curve over Rational Field defined by x^3 - y^2

        ::

            sage: A.<x,y> = AffineSpace(CC, 2)
            sage: C = Curve([y^2 + x^2], A); C
            Affine Plane Curve over Complex Field with 53 bits of precision defined
            by x^2 + y^2
        """
        if not (is_AffineSpace(A) and A.dimension != 2):
            raise TypeError("Argument A (= %s) must be an affine plane."%A)
        Curve_generic.__init__(self, A, [f])

    def _repr_type(self):
        r"""
        Return a string representation of the type of this curve.

        EXAMPLES::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: C = Curve([y - 7/2*x^5 + x - 3], A)
            sage: C._repr_type()
            'Affine Plane'
        """
        return "Affine Plane"

    def divisor_of_function(self, r):
        """
        Return the divisor of a function on a curve.

        INPUT: r is a rational function on X

        OUTPUT:


        -  ``list`` - The divisor of r represented as a list of
           coefficients and points. (TODO: This will change to a more
           structural output in the future.)


        EXAMPLES::

            sage: F = GF(5)
            sage: P2 = AffineSpace(2, F, names = 'xy')
            sage: R = P2.coordinate_ring()
            sage: x, y = R.gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f)
            sage: K = FractionField(R)
            sage: r = 1/x
            sage: C.divisor_of_function(r)     # todo: not implemented (broken)
                  [[-1, (0, 0, 1)]]
            sage: r = 1/x^3
            sage: C.divisor_of_function(r)     # todo: not implemented (broken)
                  [[-3, (0, 0, 1)]]
        """
        F = self.base_ring()
        f = self.defining_polynomial()
        pts = self.places_on_curve()
        numpts = len(pts)
        R = f.parent()
        x,y = R.gens()
        R0 = PolynomialRing(F,3,names = [str(x),str(y),"t"])
        vars0 = R0.gens()
        t = vars0[2]
        divf = []
        for pt0 in pts:
            if pt0[2] != F(0):
                lcs = self.local_coordinates(pt0,5)
                yt = lcs[1]
                xt = lcs[0]
                ldg = degree_lowest_rational_function(r(xt,yt),t)
                if ldg[0] != 0:
                    divf.append([ldg[0],pt0])
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

            sage: F = GF(5)
            sage: pt = (2,3)
            sage: R = PolynomialRing(F,2, names = ['x','y'])
            sage: x,y = R.gens()
            sage: f = y^2-x^9-x
            sage: C = Curve(f)
            sage: C.local_coordinates(pt, 9)
            [t + 2, -2*t^12 - 2*t^11 + 2*t^9 + t^8 - 2*t^7 - 2*t^6 - 2*t^4 + t^3 - 2*t^2 - 2]
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
        yt = y0*t**0+add([vars0[i]*t**(i-2) for i in range(3,2*n+2)])
        xt = x0+t
        ft = f(xt,yt)
        S = singular
        S.eval('ring s = '+str(p)+','+str(R0.gens())+',lp;')
        S.eval('poly f = '+str(ft) + ';')
        c = S('coeffs(%s, t)'%ft)
        N = int(c.size())
        b = ["%s[%s,1],"%(c.name(), i) for i in range(2,N//2-4)]
        b = ''.join(b)
        b = b[:len(b)-1] # to cut off the trailing comma
        cmd = 'ideal I = '+b
        S.eval(cmd)
        S.eval('short=0')    # print using *'s and ^'s.
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
        Plot the real points on this affine plane curve.

        INPUT:


        -  ``self`` - an affine plane curve

        -  ``*args`` - optional tuples (variable, minimum, maximum) for
           plotting dimensions

        -  ``**kwds`` - optional keyword arguments passed on to
           ``implicit_plot``


        EXAMPLES:

        A cuspidal curve::

            sage: R.<x, y> = QQ[]
            sage: C = Curve(x^3 - y^2)
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive

        A 5-nodal curve of degree 11.  This example also illustrates
        some of the optional arguments::

            sage: R.<x, y> = ZZ[]
            sage: C = Curve(32*x^2 - 2097152*y^11 + 1441792*y^9 - 360448*y^7 + 39424*y^5 - 1760*y^3 + 22*y - 1)
            sage: C.plot((x, -1, 1), (y, -1, 1), plot_points=400)
            Graphics object consisting of 1 graphics primitive

        A line over `\mathbf{RR}`::

            sage: R.<x, y> = RR[]
            sage: C = Curve(R(y - sqrt(2)*x))
            sage: C.plot()
            Graphics object consisting of 1 graphics primitive
        """
        I = self.defining_ideal()
        return I.plot(*args, **kwds)

class AffinePlaneCurve_finite_field(AffinePlaneCurve):
    def rational_points(self, algorithm="enum"):
        r"""
        Return sorted list of all rational points on this curve.

        Use *very* naive point enumeration to find all rational points on
        this curve over a finite field.

        EXAMPLE::

            sage: A.<x,y> = AffineSpace(2,GF(9,'a'))
            sage: C = Curve(x^2 + y^2 - 1)
            sage: C
            Affine Plane Curve over Finite Field in a of size 3^2 defined by x^2 + y^2 - 1
            sage: C.rational_points()
            [(0, 1), (0, 2), (1, 0), (2, 0), (a + 1, a + 1), (a + 1, 2*a + 2), (2*a + 2, a + 1), (2*a + 2, 2*a + 2)]
        """
        f = self.defining_polynomial()
        R = f.parent()
        K = R.base_ring()
        points = []
        for x in K:
            for y in K:
                if f(x,y) == 0:
                    points.append(self((x,y)))
        points.sort()
        return points


class AffinePlaneCurve_prime_finite_field(AffinePlaneCurve_finite_field):
    # CHECK WHAT ASSUMPTIONS ARE MADE REGARDING AFFINE VS. PROJECTIVE MODELS!!!
    # THIS IS VERY DIRTY STILL -- NO DATASTRUCTURES FOR DIVISORS.

    def riemann_roch_basis(self,D):
        """
        Interfaces with Singular's BrillNoether command.

        INPUT:


        -  ``self`` - a plane curve defined by a polynomial eqn f(x,y)
           = 0 over a prime finite field F = GF(p) in 2 variables x,y
           representing a curve X: f(x,y) = 0 having n F-rational
           points (see the Sage function places_on_curve)

        -  ``D`` - an n-tuple of integers
           `(d1, ..., dn)` representing the divisor
           `Div = d1*P1+...+dn*Pn`, where
           `X(F) = \{P1,...,Pn\}`.
           *The ordering is that dictated by places_on_curve.*


        OUTPUT: basis of L(Div)

        EXAMPLE::

            sage: R = PolynomialRing(GF(5),2,names = ["x","y"])
            sage: x, y = R.gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f)
            sage: D = [6,0,0,0,0,0]
            sage: C.riemann_roch_basis(D)
            [1, (y^2*z^4 - x*z^5)/x^6, (y^2*z^5 - x*z^6)/x^7, (y^2*z^6 - x*z^7)/x^8]
        """
        f = self.defining_polynomial()
        R = f.parent()
        F = self.base_ring()
        p = F.characteristic()
        Dstr = str(tuple(D))
        G = singular(','.join([str(x) for x in D]), type='intvec')

        singular.LIB('brnoeth.lib')

        S = singular.ring(p, R.gens(), 'lp')
        fsing = singular(str(f))

        X = fsing.Adj_div()
        P = singular.NSplaces(1, X)
        T = P[1][2]
        T.set_ring()

        LG = G.BrillNoether(P)

        dim = len(LG)
        basis = [(LG[i][1], LG[i][2]) for i in range(1,dim+1)]
        x, y, z = PolynomialRing(F, 3, names = ["x","y","z"]).gens()
        V = []
        for g in basis:
            T.set_ring()  # necessary...
            V.append(eval(g[0].sage_polystring())/eval(g[1].sage_polystring()))
        return V


    def rational_points(self, algorithm="enum"):
        r"""
        Return sorted list of all rational points on this curve.

        INPUT:


        -  ``algorithm`` - string:

           +  ``'enum'`` - straightforward enumeration

           +  ``'bn'`` - via Singular's Brill-Noether package.

           +  ``'all'`` - use all implemented algorithms and
              verify that they give the same answer, then return it


        .. note::

           The Brill-Noether package does not always work. When it
           fails a RuntimeError exception is raised.

        EXAMPLE::

            sage: x, y = (GF(5)['x,y']).gens()
            sage: f = y^2 - x^9 - x
            sage: C = Curve(f); C
            Affine Plane Curve over Finite Field of size 5 defined by -x^9 + y^2 - x
            sage: C.rational_points(algorithm='bn')
            [(0, 0), (2, 2), (2, 3), (3, 1), (3, 4)]
            sage: C = Curve(x - y + 1)
            sage: C.rational_points()
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)]

        We compare Brill-Noether and enumeration::

            sage: x, y = (GF(17)['x,y']).gens()
            sage: C = Curve(x^2 + y^5 + x*y - 19)
            sage: v = C.rational_points(algorithm='bn')
            sage: w = C.rational_points(algorithm='enum')
            sage: len(v)
            20
            sage: v == w
            True
        """
        if algorithm == "enum":

            return AffinePlaneCurve_finite_field.rational_points(self, algorithm="enum")

        elif algorithm == "bn":
            f = self.defining_polynomial()._singular_()
            singular = f.parent()
            singular.lib('brnoeth')
            try:
                X1 = f.Adj_div()
            except (TypeError, RuntimeError) as s:
                raise RuntimeError(str(s) + "\n\n ** Unable to use the Brill-Noether Singular package to compute all points (see above).")

            X2 = singular.NSplaces(1, X1)
            R = X2[5][1][1]
            singular.set_ring(R)

            # We use sage_flattened_str_list since iterating through
            # the entire list through the sage/singular interface directly
            # would involve hundreds of calls to singular, and timing issues
            # with the expect interface could crop up.  Also, this is vastly
            # faster (and more robust).
            v = singular('POINTS').sage_flattened_str_list()
            pnts = [self(int(v[3*i]), int(v[3*i+1])) for i in range(len(v)//3) if int(v[3*i+2])!=0]
            # remove multiple points
            pnts = sorted(set(pnts))
            return pnts

        elif algorithm == "all":

            S_enum = self.rational_points(algorithm = "enum")
            S_bn = self.rational_points(algorithm = "bn")
            if S_enum != S_bn:
                raise RuntimeError("Bug in rational_points -- different algorithms give different answers for curve %s!"%self)
            return S_enum

        else:
            raise ValueError("No algorithm '%s' known"%algorithm)
