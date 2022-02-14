"""
Polygons and triangles in hyperbolic geometry

AUTHORS:

- Hartmut Monien (2011-08)
- Vincent Delecroix (2014-11)
"""
# *****************************************************************************
#       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>,
#                     2014 Vincent Delecroix <20100.delecroix@gmail.com>,
#                     2015 Stefan Kraemer <skraemer@th.physik.uni-bonn.de>
#                     2021 Javier Honrubia <jhonrubia6@alumno.uned.es
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
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.plot.bezier_path import BezierPath
from sage.plot.circle import circle
from sage.misc.decorators import options, rename_keyword
from sage.rings.cc import CC
from sage.plot.hyperbolic_arc import HyperbolicArcCore


class HyperbolicPolygon(HyperbolicArcCore):
    """
    Primitive class for hyperbolic polygon type.

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
            sage: HP = HyperbolicPolygon([0, 1/2, I], "UHP", {})
            sage: TestSuite(HP).run(skip ="_test_pickling")
        """
        pts = [CC(_) for _ in pts]
        self.path = []
        if model == "UHP":
            if(pts[0].is_infinity()):
                # Check for more than one Infinite vertex
                for i in range(1, len(pts)):
                    if(pts[i].is_infinity()):
                        raise ValueError("No more than one infinite vertex allowed")
            else:
                # If any Infinity vertex exist it must be the first
                for i in range(1, len(pts)):
                    if(pts[i].is_infinity()):
                        pts = pts[i:] + pts[0:i]
            self._UHP_hyperbolic_arc(pts[0], pts[1], True)
            for i in range(1, len(pts) - 1):
                self._UHP_hyperbolic_arc(pts[i], pts[i + 1])
            self._UHP_hyperbolic_arc(pts[-1], pts[0])
        elif model == "PD":
            self._PD_hyperbolic_arc(pts[0], pts[1], True)
            for i in range(1, len(pts) - 1):
                self._PD_hyperbolic_arc(pts[i], pts[i + 1])
            self._PD_hyperbolic_arc(pts[-1], pts[0])
        elif model == "KM":
            self._KM_hyperbolic_arc(pts[0], pts[1], True)
            for i in range(1, len(pts) - 1):
                self._KM_hyperbolic_arc(pts[i], pts[i + 1])
            self._KM_hyperbolic_arc(pts[-1], pts[0])
        else:
            raise ValueError("%s is not a valid model for Hyperbolic plane" % model)
        self._pts = pts
        BezierPath.__init__(self, self.path, options)

    def _repr_(self):
        """
        String representation of HyperbolicPolygon.

        TESTS::

            sage: from sage.plot.hyperbolic_polygon import HyperbolicPolygon
            sage: HyperbolicPolygon([0, 1/2, I], "UHP", {})
            Hyperbolic polygon (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
        """
        return "Hyperbolic polygon ({})".format(", ".join(map(str, self._pts)))


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

    - ``linestyle`` -- (default: ``'solid'``) the style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or
      ``'--'``, ``':'``, ``'-'``, ``'-.'``, respectively

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

    With vertex at Infinity::

        sage: hyperbolic_polygon([-1,0,1,Infinity], color='green'))
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        from sage.rings.infinity import infinity
        sphinx_plot(hyperbolic_polygon([-1,0,1,infinity], color='green'))

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

    Klein model is also supported via the paraeter ``model``.
    Show a hyperbolic polygon in the Klein model with coordinates
    `1`, `e^{i\pi/3}` , `e^{i2\pi/3}` , `-1` , `e^{i4\pi/3}` , `e^{i5\pi/3}`::

        sage: p1=1
        sage: p2=(cos(pi/3),sin(pi/3))
        sage: p3=(cos(2*pi/3),sin(2*pi/3))
        sage: p4=-1
        sage: p5=(cos(4*pi/3),sin(4*pi/3))
        sage: p6=(cos(5*pi/3),sin(5*pi/3))
        hyperbolic_polygon([p1,p2,p3,p4,p5,p6],model="KM", fill=True, color='purple')
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        p1=1
        p2=(cos(pi/3),sin(pi/3))
        p3=(cos(2*pi/3),sin(2*pi/3))
        p4=-1
        p5=(cos(4*pi/3),sin(4*pi/3))
        p6=(cos(5*pi/3),sin(5*pi/3))
        P = hyperbolic_polygon([p1,p2,p3,p4,p5,p6],model="KM", fill=True, color='purple')
        sphinx_plot(P)

    Hiperboloid model is supported partially, (neither fill nor color option yet) via the paraeter ``model``.
    Show a hyperbolic polygon in the Klein model with coordinates
    `(3,3,sqrt(19)),(3,-3,sqrt(19)),(-3,-3,sqrt(19)),(-3,3,sqrt(19))`::

        sage: hyperbolic_polygon([(3,3,sqrt(19)),(3,-3,sqrt(19)),(-3,-3,sqrt(19)),(-3,3,sqrt(19))],model = "HM")
        Launched html viewer for Graphics3d Object

    .. PLOT::

        P = hyperbolic_polygon([(3,3,sqrt(19)),(3,-3,sqrt(19)),(-3,-3,sqrt(19)),(-3,3,sqrt(19))],model = "HM")
        sphinx_plot(P)

    """
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    if model != "HM":
        g.add_primitive(HyperbolicPolygon(pts, model, options))
    else:
        from sage.geometry.hyperbolic_space.hyperbolic_interface import HyperbolicPlane
        from sage.plot.plot3d.shapes2 import line3d

        # First side
        HM = HyperbolicPlane().HM()
        geodesic = HM.get_geodesic(pts[0], pts[1])
        points = geodesic._plot_vertices()

        # n - 2 sides
        for j in range(1, len(pts) - 1):
            geodesic = HM.get_geodesic(pts[j], pts[j+1])
            points = points + geodesic._plot_vertices()

        # last side
        geodesic = HM.get_geodesic(pts[-1], pts[0])
        points = points + geodesic._plot_vertices()

        g = line3d(points)  # No options passed in this version
      
    if model == "PD" or model == "KM":
        g = g + circle((0, 0), 1, rgbcolor='black')
        g.set_aspect_ratio(1)
    return g


def hyperbolic_triangle(a, b, c, model="UHP", **options):
    r"""
    Return a hyperbolic triangle in the hyperbolic plane with
    vertices ``(a,b,c)``.

    Type ``?hyperbolic_polygon`` to see all options.

    INPUT:

    - ``a, b, c`` -- complex numbers in the upper half complex plane

    OPTIONS:

    - ``alpha`` -- default: 1

    - ``fill`` -- default: ``False``

    - ``thickness`` -- default: 1

    - ``rgbcolor`` -- default: ``'blue'``

    - ``linestyle`` -- (default: ``'solid'``) the style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or
      ``'--'``, ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Show a hyperbolic triangle with coordinates `0`, `1/2+i\sqrt{3}/2` and
    `-1/2+i\sqrt{3}/2`::

         sage: hyperbolic_triangle(0, -1/2+I*sqrt(3)/2, 1/2+I*sqrt(3)/2)
         Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_triangle(0, 0.5*(-1+I*sqrt(3)), 0.5*(1+I*sqrt(3)))
        sphinx_plot(P)

    A hyperbolic triangle with coordinates `0`, `1` and `2+i` and a dashed line::

         sage: hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red', linestyle='--')
         Graphics object consisting of 1 graphics primitive

    .. PLOT::

        P = hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red', linestyle='--')
        sphinx_plot(P)

    A hyperbolic triangle with a vertex at Infinity::

        sage: hyperbolic_triangle(-5,Infinity,5)
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        from sage.rings.infinity import infinity
        sphinx_plot(hyperbolic_triangle(-5,infinity,5))


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
