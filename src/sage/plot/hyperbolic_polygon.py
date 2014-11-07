"""
Triangles in hyperbolic geometry

AUTHORS:

- Hartmut Monien (2011 - 08)
"""
#*****************************************************************************
#       Copyright (C) 2011 Hartmut Monien <monien@th.physik.uni-bonn.de>,
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
from sage.plot.bezier_path import BezierPath
from sage.plot.colors import to_mpl_color
from sage.plot.misc import options, rename_keyword
from sage.rings.all import CC

class HyperbolicTriangle(BezierPath):
    """
    Primitive class for hyberbolic triangle type. See ``hyperbolic_triangle?``
    for information about plotting a hyperbolic triangle in the complex plane.

    INPUT:

    - ``a,b,c`` - coordinates of the hyperbolic triangle in the upper
      complex plane

    - ``options`` - dict of valid plot options to pass to constructor

    EXAMPLES:

    Note that constructions should use ``hyperbolic_triangle``::

         sage: from sage.plot.hyperbolic_triangle import HyperbolicTriangle
         sage: print HyperbolicTriangle(0, 1/2, I, {})
         Hyperbolic triangle (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
    """
    def __init__(self, A, B, C, options):
        """
        Initialize HyperbolicTriangle:

        Examples::

            sage: from sage.plot.hyperbolic_triangle import HyperbolicTriangle
            sage: print HyperbolicTriangle(0, 1/2, I, {})
            Hyperbolic triangle (0.000000000000000, 0.500000000000000, 1.00000000000000*I)
        """
        A, B, C = (CC(A), CC(B), CC(C))
        self.path = []
        self._hyperbolic_arc(A, B, True);
        self._hyperbolic_arc(B, C);
        self._hyperbolic_arc(C, A);
        BezierPath.__init__(self, self.path, options)
        self.A, self.B, self.C = (A, B, C)

    def _repr_(self):
        """
        String representation of HyperbolicArc.

        TESTS::

            sage: from sage.plot.hyperbolic_triangle import HyperbolicTriangle
            sage: HyperbolicTriangle(0, 1/2, I,{})._repr_()
            'Hyperbolic triangle (0.000000000000000, 0.500000000000000, 1.00000000000000*I)'
        """
        return "Hyperbolic triangle (%s, %s, %s)" % (self.A, self.B, self.C)

    def _hyperbolic_arc(self, z0, z3, first=False):
        """
        Function to construct Bezier path as an approximation to
        the hyperbolic arc between the complex numbers z0 and z3 in the
        hyperbolic plane.
        """
        if (z0-z3).real() == 0:
            self.path.append([(z0.real(),z0.imag()), (z3.real(),z3.imag())])
            return
        z0, z3 = (CC(z0), CC(z3))
        if z0.imag() == 0 and z3.imag() == 0:
            p = (z0.real()+z3.real())/2
            r = abs(z0-p)
            zm = CC(p, r)
            self._hyperbolic_arc(z0, zm, first)
            self._hyperbolic_arc(zm, z3)
            return
        else:
            p = (abs(z0)*abs(z0)-abs(z3)*abs(z3))/(z0-z3).real()/2
            r = abs(z0-p)
            zm = ((z0+z3)/2-p)/abs((z0+z3)/2-p)*r+p
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

@rename_keyword(color='rgbcolor')
@options(alpha=1, fill=False, thickness=1, rgbcolor="blue", zorder=2, linestyle='solid')

def hyperbolic_triangle(a, b, c, **options):
    """
    Return a hyperbolic triangle in the complex hyperbolic plane with points
    (a, b, c). Type ``?hyperbolic_triangle`` to see all options.

    INPUT:

    - ``a, b, c`` - complex numbers in the upper half complex plane

    OPTIONS:

    - ``alpha`` - default: 1

    - ``fill`` - default: False

    - ``thickness`` - default: 1

    - ``rgbcolor`` - default: 'blue'

    - ``linestyle`` - (default: ``'solid'``) The style of the line, which is
      one of ``'dashed'``, ``'dotted'``, ``'solid'``, ``'dashdot'``, or ``'--'``,
      ``':'``, ``'-'``, ``'-.'``, respectively.

    EXAMPLES:

    Show a hyperbolic triangle with coordinates 0, `1/2+i\sqrt{3}/2` and
    `-1/2+i\sqrt{3}/2`::

         sage: hyperbolic_triangle(0, -1/2+I*sqrt(3)/2, 1/2+I*sqrt(3)/2)
         Graphics object consisting of 1 graphics primitive

    A hyperbolic triangle with coordinates 0, 1 and 2+i and a dashed line::

         sage: hyperbolic_triangle(0, 1, 2+i, fill=true, rgbcolor='red', linestyle='--')
         Graphics object consisting of 1 graphics primitive
    """
    from sage.plot.all import Graphics
    g = Graphics()
    g._set_extra_kwds(g._extract_kwds_for_show(options))
    g.add_primitive(HyperbolicTriangle(a, b, c, options))
    g.set_aspect_ratio(1)
    return g
