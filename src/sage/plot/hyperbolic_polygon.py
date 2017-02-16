"""
Polygons and triangles in hyperbolic geometry

AUTHORS:

- Hartmut Monien (2011-08)
- Vincent Delecroix (2014-11)
"""
#*****************************************************************************
#       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>,
#                     2014 Vincent Delecroix <20100.delecroix@gmail.com>,
#                     2015 Stefan Kraemer <skraemer@th.physik.uni-bonn.de>
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
from __future__ import print_function

from sage.plot.bezier_path import BezierPath
from sage.plot.graphics import Graphics
from sage.plot.circle import circle
from sage.plot.arc import arc
from sage.plot.misc import options, rename_keyword
from sage.rings.all import CC
from sage.symbolic.constants import pi
from sage.functions.trig import cos, sin, tan

class HyperbolicPolygon(BezierPath):
    """
    Primitive class for hyberbolic polygon type.

    See ``hyperbolic_polygon?`` for information about plotting a hyperbolic
    polygon in the complex plane.

    INPUT:

    - ``pts`` -- coordinates of the polygon (as complex numbers)

    - ``options`` -- dict of valid plot options to pass to constructor

    EXAMPLES:

    Note that constructions should use :func:`hyperbolic_polygon` or
    :func:`hyperbolic_triangle`::

         sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
         sage: print(HyperbolicPolygon([0, 1/2, I], "UHP", {}))
         Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
    """
    def __init__(self, pts, model, options):
        """
        Initialize HyperbolicPolygon.

        EXAMPLES::

            sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
            sage: print(HyperbolicPolygon([0, 1/2, I], "UHP", {}))
            Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
        """
        pts = [CC(_) for _ in pts]
        self.path = []
        if model == "UHP":
            self._UHP_hyperbolic_arc(pts[0], pts[1], True)
            for i in range(1, len(pts) - 1):
                self._UHP_hyperbolic_arc(pts[i], pts[i + 1])
            self._UHP_hyperbolic_arc(pts[-1], pts[0])
        else:
            if model == "PD":
                self._PD_hyperbolic_arc(pts[0], pts[1],True)
                for i in range(1, len(pts) - 1):
                    self._PD_hyperbolic_arc(pts[i], pts[i + 1])
                self._PD_hyperbolic_arc(pts[-1], pts[0])
        self._pts = pts
        BezierPath.__init__(self, self.path, options)
    def _repr_(self):
        """
        String representation of HyperbolicPolygon.

        TESTS::

            sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
            sage: HyperbolicPolygon([0, 1/2, I], "UHP", {})._repr_()
            'Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)'
        """
        return "Hyperbolic polygon ({})".format(", ".join(map(str, self._pts)))
    def _UHP_hyperbolic_arc(self, z0, z3, first=False):
        """
        Function to construct Bezier path as an approximation to
        the hyperbolic arc between the complex numbers z0 and z3 in the
        hyperbolic plane.
        """
        z0, z3 = (CC(z0), CC(z3))
        p = (abs(z0)*abs(z0)-abs(z3)*abs(z3))/(z0-z3).real()/2
        r = abs(z0 - p)

        if abs(z3 - z0)/r < 0.1:
            self.path.append([(z0.real(),z0.imag()), (z3.real(),z3.imag())])
            return

        if z0.imag() == 0 and z3.imag() == 0:
            p = (z0.real()+z3.real()) / 2
            zm = CC(p, r)
            self._UHP_hyperbolic_arc(z0, zm, first)
            self._UHP_hyperbolic_arc(zm, z3)
            return
        else:
            zm = ((z0+z3)/2-p)/abs((z0+z3)/2-p)*r + p
            t = (8*zm-4*(z0+z3)).imag()/3/(z3-z0).real()
            z1 = z0 + t*CC(z0.imag(), (p-z0.real()))
            z2 = z3 - t*CC(z3.imag(), (p-z3.real()))
        if first:
            self.path.append([(z0.real(), z0.imag()),
                              (z1.real(), z1.imag()),
                              (z2.real(), z2.imag()),
                              (z3.real(), z3.imag())]);
            first = False
        else:
            self.path.append([(z1.real(), z1.imag()),
                              (z2.real(), z2.imag()),
                              (z3.real(), z3.imag())]);
    def _bezier_path(self, arc0, z0, first):
        """
        Returns the corresponding bezier path
        """
        ma = arc0._matplotlib_arc()
        transform = ma.get_transform().get_matrix()
        cA, cC, cE = transform[0]
        cB, cD, cF = transform[1]
        points = []
        for u in ma._path.vertices:
            x, y = list(u)
            points += [(cA * x + cC * y + cE, cB * x + cD * y + cF)]
        if abs(CC(points[0]) - z0)>0.00001:
            points.reverse() #order is important
        if first:
            self.path.append(points[0: 4])
            N = 4
        else:
            N = 0
            #the first point is equal to the last of the previous arc
            del points[0]
        #Complete the last tuple of control points
        tail = self.path[-1]
        ltail = len(tail)
        while ltail < 3:
            self.path[-1].append(points[N])
            ltail +=1
            N +=1
        #Add new triplets
        while N < len(points):
            self.path.append(points[N: N + 3])
            N += 3
        return
    def _CenterAndRadiusOrthogonalcirclegiven2points(self,z1,z2):
        """
        Calculate center and radius of an orthogonal circle to the
        unit disc through z1, z2 if they are ideal points or
        through z1, z2 and their inverses.
        """
        z1,z2 = (CC(z1),CC(z2))
        if abs(z1.abs()-1)<0.000001 and abs(z2.abs()-1)<0.000001:
            #both z1,z2 are ideal points
            phi1 = z1.arg()
            phi2 = z2.arg()
            deltaphi = 0.5 * abs(phi1-phi2)
            phi = 0.5 * (phi1+phi2)
            extR = sin(deltaphi) * tan(deltaphi)
            R = cos(deltaphi) + extR
            #center is at R*exp(i*(min(phi1,phi2)+deltaphi)
            argcenter = min(phi1,phi2)+deltaphi
            c = CC(R*cos(argcenter),R*sin(argcenter))
            r = abs(tan(deltaphi))
        else: 
            if z1.abs()<1:
                z3=(1/z1).conjugate()
            else:
                if z2.abs()<1:
                    z3=(1/z2).conjugate()
            k1=z1-z3
            k2=z3-z2
            s=(k1*k2.conjugate()).real()/(k1*k2.conjugate()).imag()
            c=(z1+z2)/2+CC(0,1)*s*(z1-z2)/2
            r=(z1-c).abs()
        return c,r
    def _PD_hyperbolic_arc(self, z0, z3, first=False):
        """
        Function to construct an hyperbolic arc between the complez numbers z0
        and z3 in the Poincare Disc model
        """
        z0, z3 = (CC(z0), CC(z3))
        phi0 = z0.arg()
        phi3 = z3.arg()
        if abs(phi0 - phi3) == 0 or abs(phi0 - phi3) == pi:
            #The points lie in a geodesic of the first kind
            self.path.append([(z0.real(),z0.imag()), (z3.real(),z3.imag())])
            return
        else:
            #The points lie in a geodesic of the second kind
            T=self._CenterAndRadiusOrthogonalcirclegiven2points(z0,z3)
            a1=(z0-T[0]).arg()
            a2=(z3-T[0]).arg()
            if (T[0]).real()>0:
                if a1<0:
                    a1=a1+2*pi
                if a2<0:
                    a2=a2+2*pi
            pic = arc((T[0].real(),T[0].imag()),T[1],sector=(a1,a2))[0]
            #Transform arc into a bezier path
            self._bezier_path(pic,z0,first)
            return


@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=False, thickness=1, rgbcolor="blue", zorder=2,
         linestyle='solid')
def hyperbolic_polygon(pts, model="UHP", **options):
    r"""
    Return a hyperbolic polygon in the hyperbolic plane with vertices ``pts``.

    Type ``?hyperbolic_polygon`` to see all options.

    INPUT:

    - ``pts`` -- a list or tuple of complex numbers

    OPTIONS:

    
    - ``model`` -- default: ``UHP`` Model used for hyperbolic plane
    
    - ``alpha`` -- default: 1

    - ``fill`` -- default: ``False``

    - ``thickness`` -- default: 1

    - ``rgbcolor`` -- default: ``'blue'``

    - ``linestyle`` -- (default: ``'solid'``) The style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or ``'--'``,
      ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Show a hyperbolic polygon with coordinates `-1`, `3i`, `2+2i`, `1+i`::

        sage: hyperbolic_polygon([-1,3*I,2+2*I,1+I])
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_polygon([-1,3*I,2+2*I,1+I])
        sphinx_plot(P)

    With more options::

        sage: hyperbolic_polygon([-1,3*I,2+2*I,1+I], fill=True, color='red')
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        sphinx_plot(hyperbolic_polygon([-1,3*I,2+2*I,1+I], fill=True, color='red'))

    Poincare disc model is supported via the parameter ``model``.
    Show a hyperbolic polygon in the Poincare disc model with coordinates
    `1`, `i`, `-1`, `-i`::

        sage: hyperbolic_polygon([1,I,-1,-I], model="PD", color='green')
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        sphinx_plot(hyperbolic_polygon([1,I,-1,-I], model="PD", color='green'))

    With more options::

        sage: hyperbolic_polygon([1,I,-1,-I], model="PD", color='green', fill=True, linestyle="-")
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        P = hyperbolic_polygon([1,I,-1,-I], model="PD", color='green', fill=True, linestyle="-")
        sphinx_plot(P)
    """
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    g.add_primitive(HyperbolicPolygon(pts, model, options))
    if model=="PD":
        g  = g + circle((0,0),1, rgbcolor='black')
    g.set_aspect_ratio(1)
    return g


def hyperbolic_triangle(a, b, c, model="UHP", **options):
    """
    Return a hyperbolic triangle in the hyperbolic plane with vertices ``(a,b,c)``.

    Type ``?hyperbolic_polygon`` to see all options.

    INPUT:

    - ``a, b, c`` -- complex numbers in the upper half complex plane

    OPTIONS:

    - ``alpha`` -- default: 1

    - ``fill`` -- default: ``False``

    - ``thickness`` -- default: 1

    - ``rgbcolor`` -- default: ``'blue'``

    - ``linestyle`` - (default: ``'solid'``) The style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or ``'--'``,
      ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Show a hyperbolic triangle with coordinates `0, 1/2+i\sqrt{3}/2` and
    `-1/2+i\sqrt{3}/2`::

         sage: hyperbolic_triangle(0, -1/2+I*sqrt(3)/2, 1/2+I*sqrt(3)/2)
         Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_triangle(0, 0.5*(-1+I*sqrt(3)), 0.5*(1+I*sqrt(3)))
        sphinx_plot(P)

    A hyperbolic triangle with coordinates `0, 1` and `2+i` and a dashed line::

         sage: hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red', linestyle='--')
         Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red', linestyle='--')
        sphinx_plot(P)

    It can also plot a hyperbolic triangle in the Poincare Disc model::

        sage: z1 = CC((cos(pi/3),sin(pi/3)))
        sage: z2 = CC((0.6*cos(3*pi/4),0.6*sin(3*pi/4)))
        sage: z3 = 1
        sage: hyperbolic_triangle(z1, z2, z3, model="PD", color="red")
        Graphics object consisting of 2 graphics primitives

    .. PLOT::
 
        z1 = CC((cos(pi/3),sin(pi/3)))
        z2 = CC((0.6*cos(3*pi/4),0.6*sin(3*pi/4)))
        z3 = 1
        P = hyperbolic_triangle(z1, z2, z3, model="PD", color="red")
        sphinx_plot(P)

    ::

        sage: hyperbolic_triangle(0.3+0.3*I, 0.8*I, -0.5-0.5*I, model="PD", color='magenta')
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        P = hyperbolic_triangle(0.3+0.3*I, 0.8*I, -0.5-0.5*I, model="PD", color='magenta')
        sphinx_plot(P)

    """
    return hyperbolic_polygon((a, b, c), model, **options)
