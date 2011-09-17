r"""
Triangulations of a point configuration

A point configuration is a finite set of points in Euclidean space or,
more generally, in projective space. A triangulation is a simplicial
decomposition of the convex hull of a given point configuration such
that all vertices of the simplices end up lying on points of the
configuration. That is, there are no new vertices apart from the
initial points.

Note that points that are not vertices of the convex hull need not be
used in the triangulation. A triangulation that does make use of all
points of the configuration is called fine, and you can restrict
yourself to such triangulations if you want. See
:class:`PointConfiguration` and
:meth:`~PointConfiguration.restrict_to_fine_triangulations` for
more details.

Finding a single triangulation and listing all connected
triangulations is implemented natively in this package. However, for
more advanced options [TOPCOM]_ needs to be installed. You can find an
experimental spkg at http://trac.sagemath.org/sage_trac/ticket/8169

NOTE:

TOPCOM and the internal algorithms tend to enumerate triangulations in
a different order. This is why we always explicitly specify the engine
as ``engine='TOPCOM'`` or ``engine='internal'`` in the doctests. In
your own applications, you do not need to specify the engine. By
default, TOPCOM is used if it is available and the internal algorithms
are used otherwise.

EXAMPLES:

First, we select the internal implementation for enumerating
triangulations::

   sage: PointConfiguration.set_engine('internal')   # to make doctests independent of TOPCOM

A 2-dimensional point configuration::

    sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
    sage: p
    A point configuration in QQ^2 consisting of 5 points. The
    triangulations of this point configuration are assumed to
    be connected, not necessarily fine, not necessarily regular.
    sage: t = p.triangulate()  # a single triangulation
    sage: t
    (<1,3,4>, <2,3,4>)
    sage: len(t)
    2
    sage: t[0]
    (1, 3, 4)
    sage: t[1]
    (2, 3, 4)
    sage: list(t)
    [(1, 3, 4), (2, 3, 4)]
    sage: t.plot(axes=False)
    sage: list( p.triangulations() )
    [(<1,3,4>, <2,3,4>),
     (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
     (<1,2,3>, <1,2,4>),
     (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)]
    sage: p_fine = p.restrict_to_fine_triangulations()
    sage: p_fine
    A point configuration in QQ^2 consisting of 5 points. The
    triangulations of this point configuration are assumed to
    be connected, fine, not necessarily regular.
    sage: list( p_fine.triangulations() )
    [(<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
     (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)]

A 3-dimensional point configuration::

    sage: p = [[0,-1,-1],[0,0,1],[0,1,0], [1,-1,-1],[1,0,1],[1,1,0]]
    sage: points = PointConfiguration(p)
    sage: triang = points.triangulate()
    sage: triang.plot(axes=False)

The standard example of a non-regular triangulation::

    sage: p = PointConfiguration([[-1,-5/9],[0,10/9],[1,-5/9],[-2,-10/9],[0,20/9],[2,-10/9]])
    sage: regular = p.restrict_to_regular_triangulations(True).triangulations_list()      # optional - TOPCOM
    sage: nonregular = p.restrict_to_regular_triangulations(False).triangulations_list()  # optional - TOPCOM
    sage: len(regular)     # optional - TOPCOM
    16
    sage: len(nonregular)  # optional - TOPCOM
    2
    sage: nonregular[0].plot(aspect_ratio=1, axes=False)   # optional - TOPCOM

Note that the points need not be in general position. That is, the
points may lie in a hyperplane and the linear dependencies will be
removed before passing the data to TOPCOM which cannot handle it::

    sage: points = [[0,0,0,1],[0,3,0,1],[3,0,0,1],[0,0,1,1],[0,3,1,1],[3,0,1,1],[1,1,2,1]]
    sage: points = [ p+[1,2,3] for p in points ]
    sage: pc = PointConfiguration(points)
    sage: pc.ambient_dim()
    7
    sage: pc.dim()
    3
    sage: pc.triangulate()
    (<0,1,2,3>, <1,2,3,4>, <2,3,4,5>, <3,4,5,6>)
    sage: len( pc.triangulations_list() )
    26

REFERENCES:

    .. [TOPCOM]
       J. Rambau,
       TOPCOM <http://www.rambau.wm.uni-bayreuth.de/TOPCOM/>.

    .. [GKZ]
       Gel'fand, I. M.; Kapranov, M. M.; and Zelevinsky, A. V.
       "Discriminants, Resultants and Multidimensional Determinants" Birkhauser 1994.

    .. [PUNTOS]
       Jesus A. De Loera
       http://www.math.ucdavis.edu/~deloera/RECENT_WORK/puntos2000

AUTHORS:

    - Volker Braun: initial version, 2010

    - Josh Whitney: added functionality for computing
      volumes and secondary polytopes of PointConfigurations

    - Marshall Hampton: improved documentation and doctest coverage

    - Volker Braun: rewrite using Parent/Element and catgories. Added
      a Point class. More doctests. Less zombies.

    - Volker Braun: Cythonized parts of it, added a C++ implementation
      of the bistellar flip algorithm to enumerate all connected
      triangulations.
"""

########################################################################
# Note: The doctests that make use of TOPCOM are
#       marked # optional - TOPCOM
#       If you have it installed, run doctests as
#
#   sage -tp 4 -long -optional TOPCOM sage/geometry/triangulation/
########################################################################


########################################################################
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2010 Josh Whitney <josh.r.whitney@gmail.com>
#       Copyright (C) 2010 Marshall Hampton <hamptonio@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element

from sage.combinat.combination import Combinations
from sage.rings.all import QQ, ZZ
from sage.matrix.constructor import matrix
from sage.modules.all import vector
from sage.groups.perm_gps.permgroup import PermutationGroup

from copy import copy
import sys
import pexpect


from sage.geometry.triangulation.base import \
    PointConfiguration_base, Point, ConnectedTriangulationsIterator


########################################################################
def triangulation_render_2d(triangulation, **kwds):
    r"""
    Return a graphical representation of a 2-d triangulation.

    INPUT:

    - ``triangulation`` -- a :class:`Triangulation`.

    - ``**kwds`` -- keywords that are passed on to the graphics primitives.

    OUTPUT:

    A 2-d graphics object.

    EXAMPLES::

        sage: points = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
        sage: triang = points.triangulate()
        sage: triang.plot(axes=False)   # indirect doctest
    """
    from sage.plot.all import point2d, line2d, arrow, polygon2d
    points = [ point.reduced_affine() for point in triangulation.point_configuration() ]
    coord = [ [p[0], p[1]] for p in points ]
    plot_points = sum([ point2d(p,
                                zorder=2, pointsize=10, **kwds)
                        for p in coord ])

    tmp_lines = []
    for t in triangulation:
        if len(t)>=2:
            tmp_lines.append([t[0], t[1]])
        if len(t)>=3:
            tmp_lines.append([t[0], t[2]])
            tmp_lines.append([t[1], t[2]])
    all_lines = []
    interior_lines = []
    for l in tmp_lines:
        if l not in all_lines:
            all_lines.append(l)
        else:
            interior_lines.append(l)
    exterior_lines = [ l for l in all_lines if not l in interior_lines ]

    plot_interior_lines = sum([ line2d([ coord[l[0]], coord[l[1]] ],
                                       zorder=1, rgbcolor=(0,1,0), **kwds)
                                for l in interior_lines ])
    plot_exterior_lines = sum([ line2d([ coord[l[0]], coord[l[1]] ],
                                       zorder=1, rgbcolor=(0,0,1), **kwds)
                                for l in exterior_lines ])

    plot_triangs = sum([ polygon2d([coord[t[0]], coord[t[1]], coord[t[2]]],
                                   zorder=0, rgbcolor=(0.8, 1, 0.8), **kwds)
                         for t in triangulation if len(t)>=3 ])

    return \
        plot_points + \
        plot_interior_lines + plot_exterior_lines + \
        plot_triangs




def triangulation_render_3d(triangulation, **kwds):
    r"""
    Return a graphical representation of a 3-d triangulation.

    INPUT:

    - ``triangulation`` -- a :class:`Triangulation`.

    - ``**kwds`` -- keywords that are  passed on to the graphics primitives.

    OUTPUT:

    A 3-d graphics object.

    EXAMPLES::

        sage: p = [[0,-1,-1],[0,0,1],[0,1,0], [1,-1,-1],[1,0,1],[1,1,0]]
        sage: points = PointConfiguration(p)
        sage: triang = points.triangulate()
        sage: triang.plot(axes=False)     # indirect doctest
    """
    from sage.plot.plot3d.all import point3d, line3d, arrow3d, polygon3d
    points = [ point.reduced_affine() for point in triangulation.point_configuration() ]
    coord = [ [p[0], p[1], p[2] ] for p in points ]
    plot_points = sum([ point3d(p, size=15,
                                **kwds)
                        for p in coord ])

    tmp_lines = []
    for t in triangulation:
        if len(t)>=2:
            tmp_lines.append([t[0], t[1]])
        if len(t)>=3:
            tmp_lines.append([t[0], t[2]])
            tmp_lines.append([t[1], t[2]])
        if len(t)>=4:
            tmp_lines.append([t[0], t[3]])
            tmp_lines.append([t[1], t[3]])
            tmp_lines.append([t[2], t[3]])
    all_lines = []
    interior_lines = []
    for l in tmp_lines:
        if l not in all_lines:
            all_lines.append(l)
        else:
            interior_lines.append(l)
    exterior_lines = [ l for l in all_lines if not l in interior_lines ]

    from sage.plot.plot3d.texture import Texture
    line_int = Texture(color='darkblue', ambient=1, diffuse=0)
    line_ext = Texture(color='green', ambient=1, diffuse=0)
    triang_int = Texture(opacity=0.3, specular=0, shininess=0, diffuse=0, ambient=1, color='yellow')
    triang_ext = Texture(opacity=0.6, specular=0, shininess=0, diffuse=0, ambient=1, color='green')

    plot_interior_lines = sum([ line3d([ coord[l[0]], coord[l[1]] ],
                                       thickness=2, texture=line_int, **kwds)
                                for l in interior_lines ])
    plot_exterior_lines = sum([ line3d([ coord[l[0]], coord[l[1]] ],
                                       thickness=3, texture=line_ext, **kwds)
                                for l in exterior_lines ])

    tmp_triangs = []
    for t in triangulation:
        if len(t)>=3:
            tmp_triangs.append([t[0], t[1], t[2]])
        if len(t)>=4:
            tmp_triangs.append([t[0], t[1], t[3]])
            tmp_triangs.append([t[0], t[2], t[3]])
            tmp_triangs.append([t[1], t[2], t[3]])
    all_triangs = []
    interior_triangs = []
    for l in tmp_triangs:
        if l not in all_triangs:
            all_triangs.append(l)
        else:
            interior_triangs.append(l)
    exterior_triangs = [ l for l in all_triangs if not l in interior_triangs ]

    plot_interior_triangs = \
        sum([ polygon3d([coord[t[0]], coord[t[1]], coord[t[2]]],
                        texture = triang_int, **kwds)
              for t in interior_triangs ])
    plot_exterior_triangs = \
        sum([ polygon3d([coord[t[0]], coord[t[1]], coord[t[2]]],
                        texture = triang_ext, **kwds)
              for t in exterior_triangs ])

    return \
        plot_points + \
        plot_interior_lines + plot_exterior_lines + \
        plot_interior_triangs + plot_exterior_triangs





########################################################################
class Triangulation(Element):
    """
    A triangulation of a :class:`PointConfiguration`.

    .. WARNING::

        You should never create :class:`Triangulation` objects
        manually. See
        :meth:`~sage.geometry.triangulate.PointConfiguration.triangulate`
        and
        :meth:`~sage.geometry.triangulate.PointConfiguration.triangulations`
        to triangulate point configurations.
    """

    def __init__(self, triangulation, parent, check=True):
        """
        The constructor of a ``Triangulation`` object. Note that an
        internal reference to the underlying ``PointConfiguration`` is
        kept.

        INPUT:

        - ``parent`` -- a :class:`PointConfiguration`

        - ``triangulation`` -- an iterable of integers or iterable of
          iterables (e.g. a list of lists). In the first case, the
          integers specify simplices via
          :meth:`PointConfiguration.simplex_to_int`. In the second
          case, the point indices of the maximal simplices of the
          triangulation.

        - ``check`` -- boolean. Whether to perform checks that the
          triangulation is, indeed, a triangulation of the point
          configuration.

        NOTE:

        Passing ``check=False`` allows you to create triangulations of
        subsets of the points of the configuration, see
        :meth:`PointConfiguration.bistellar_flips`.

        EXAMPLES::

            sage: p = [[0,1],[0,0],[1,0]]
            sage: points = PointConfiguration(p)
            sage: from sage.geometry.triangulation.point_configuration import Triangulation
            sage: Triangulation([(0,1,2)], points)
            (<0,1,2>)
            sage: Triangulation([1], points)
            (<0,1,2>)
        """
        Element.__init__(self, parent=parent)
        self._point_configuration = parent

        try:
            triangulation = tuple(sorted( tuple(sorted(t)) for t in triangulation))
        except TypeError:
            triangulation = tuple( self.point_configuration().int_to_simplex(i)
                                   for i in triangulation )
        assert not check or all( len(t)==self.point_configuration().dim()+1
                                 for t in triangulation)
        self._triangulation = triangulation


    def point_configuration(self):
        """
        Returns the point configuration underlying the triangulation.

        EXAMPLES::

            sage: pconfig = PointConfiguration([[0,0],[0,1],[1,0]])
            sage: pconfig
            A point configuration in QQ^2 consisting of 3 points. The
            triangulations of this point configuration are assumed to
            be connected, not necessarily fine, not necessarily regular.
            sage: triangulation = pconfig.triangulate()
            sage: triangulation
            (<0,1,2>)
            sage: triangulation.point_configuration()
            A point configuration in QQ^2 consisting of 3 points. The
            triangulations of this point configuration are assumed to
            be connected, not necessarily fine, not necessarily regular.
            sage: pconfig == triangulation.point_configuration()
            True
        """
        return self._point_configuration


    def __cmp__(self, right):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        - 0 if ``right`` is the same triangulation as ``self``, 1 or
          -1 otherwise.

        TESTS::

            sage: pc = PointConfiguration([[0,0],[0,1],[1,0]])
            sage: t1 = pc.triangulate()
            sage: from sage.geometry.triangulation.point_configuration import Triangulation
            sage: t2 = Triangulation([[2,1,0]], pc)
            sage: t1 is t2
            False
            sage: cmp(t1, t2)
            0
            sage: t1 == t2    # indirect doctest
            True
            sage: abs( cmp(t1, Triangulation(((0,1),(1,2)), pc, check=False) ))
            1
            sage: abs( cmp(t2, "not a triangulation") )
            1
        """
        left = self
        c = cmp(isinstance(left,Triangulation), isinstance(right,Triangulation))
        if c: return c
        c = cmp(left.point_configuration(), right.point_configuration())
        if c: return c
        return cmp(left._triangulation, right._triangulation)


    def __iter__(self):
        """
        Iterate through the simplices of the triangulation.

        EXAMPLES::

            sage: pc = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: triangulation = pc.triangulate()
            sage: iter = triangulation.__iter__()
            sage: iter.next()
            (1, 3, 4)
            sage: iter.next()
            (2, 3, 4)
            sage: iter.next()
            Traceback (most recent call last):
            ...
            StopIteration
        """
        for p in self._triangulation:
            yield p


    def __getitem__(self, i):
        """
        Access the point indices of the i-th simplex of the triangulation.

        INPUT:

        - ``i`` -- integer. The index of a simplex.

        OUTPUT:

        A tuple of integers. The vertex indices of the i-th simplex.

        EXAMPLES::

            sage: pc = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: triangulation = pc.triangulate()
            sage: triangulation[1]
            (2, 3, 4)
        """
        return self._triangulation[i]


    def __len__(self):
        """
        Returns the length of the triangulation.

        TESTS::

            sage: pc = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: triangulation = pc.triangulations().next()
            sage: triangulation.__len__()
            2
            sage: len(triangulation)    # equivalent
            2
        """
        return len(self._triangulation)


    def _repr_(self):
        r"""
        Return a string representation.

        TESTS::

            sage: pc = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1],[2,2]])
            sage: t = pc.triangulations()
            sage: t.next()._repr_()
            '(<1,4,5>, <2,4,5>)'
        """
        #s = 'A triangulation'
        #s += ' in QQ^'+str(self.point_configuration().ambient_dim())
        #s += ' consisting of '+str(len(self))+' simplices.'
        s = '('
        s += ', '.join([ '<'+','.join(map(str,t))+'>' for t in self._triangulation])
        s += ')'
        return s


    def plot(self, **kwds):
        r"""
        Produce a graphical representation of the triangulation.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: triangulation = p.triangulate()
            sage: triangulation
            (<1,3,4>, <2,3,4>)
            sage: triangulation.plot(axes=False)
        """
        dim = self.point_configuration().dim()

        if dim == 2:
            return triangulation_render_2d(self, **kwds)

        if dim == 3:
            return triangulation_render_3d(self, **kwds)

        raise NotImplementedError, \
            'Plotting '+str(dim)+'-dimensional triangulations not implemented!'


    def gkz_phi(self):
        r"""
        Calculate the GKZ phi vector of the triangulation.

        OUTPUT:

        Vector -- the phi vector of self.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[1,0],[2,1],[1,2],[0,1]])
            sage: p.triangulate().gkz_phi()
            (1, 3, 4, 2, 5)

        NOTE:

        For a definition of the phi vector, see [GKZ]_ page 220 equation 1.4.
        """
        vec = vector(ZZ, self.point_configuration().n_points())
        for simplex in self:
            vol = self.point_configuration().volume(simplex)
            for i in simplex:
                vec[i] = vec[i] + vol
        return vec


    def enumerate_simplices(self):
        r"""
        Return the enumerated simplices.

        OUTPUT:

        A tuple of integers that uniquely specifies the triangulation.

        EXAMPLES::

            sage: pc = PointConfiguration(matrix([
            ...      [ 0, 0, 0, 0, 0, 2, 4,-1, 1, 1, 0, 0, 1, 0],
            ...      [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0],
            ...      [ 0, 2, 0, 0, 0, 0,-1, 0, 1, 0, 1, 0, 0, 1],
            ...      [ 0, 1, 1, 0, 0, 1, 0,-2, 1, 0, 0,-1, 1, 1],
            ...      [ 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0]
            ...   ]).columns())
            sage: triangulation = pc.lexicographic_triangulation()
            sage: triangulation.enumerate_simplices()
            (1678, 1688, 1769, 1779, 1895, 1905, 2112, 2143, 2234, 2360, 2555, 2580,
             2610, 2626, 2650, 2652, 2654, 2661, 2663, 2667, 2685, 2755, 2757, 2759,
             2766, 2768, 2772, 2811, 2881, 2883, 2885, 2892, 2894, 2898)

        You can recreate the triangulation from this list by passing
        it to the constructor::

            sage: from sage.geometry.triangulation.point_configuration import Triangulation
            sage: Triangulation([1678, 1688, 1769, 1779, 1895, 1905, 2112, 2143,
            ...    2234, 2360, 2555, 2580, 2610, 2626, 2650, 2652, 2654, 2661, 2663,
            ...    2667, 2685, 2755, 2757, 2759, 2766, 2768, 2772, 2811, 2881, 2883,
            ...    2885, 2892, 2894, 2898], pc)
            (<1,3,4,7,10,13>, <1,3,4,8,10,13>, <1,3,6,7,10,13>, <1,3,6,8,10,13>,
             <1,4,6,7,10,13>, <1,4,6,8,10,13>, <2,3,4,6,7,12>, <2,3,4,7,12,13>,
             <2,3,6,7,12,13>, <2,4,6,7,12,13>, <3,4,5,6,9,12>, <3,4,5,8,9,12>,
             <3,4,6,7,11,12>, <3,4,6,9,11,12>, <3,4,7,10,11,13>, <3,4,7,11,12,13>,
             <3,4,8,9,10,12>, <3,4,8,10,12,13>, <3,4,9,10,11,12>, <3,4,10,11,12,13>,
             <3,5,6,8,9,12>, <3,6,7,10,11,13>, <3,6,7,11,12,13>, <3,6,8,9,10,12>,
             <3,6,8,10,12,13>, <3,6,9,10,11,12>, <3,6,10,11,12,13>, <4,5,6,8,9,12>,
             <4,6,7,10,11,13>, <4,6,7,11,12,13>, <4,6,8,9,10,12>, <4,6,8,10,12,13>,
             <4,6,9,10,11,12>, <4,6,10,11,12,13>)
        """
        pc = self._point_configuration
        return tuple( pc.simplex_to_int(t) for t in self )


    def fan(self, origin=None):
        r"""
        Construct the fan of cones over the simplices of the triangulation.

        INPUT:

        - ``origin`` -- ``None`` (default) or coordinates of a
          point. The common apex of all cones of the fan. If ``None``,
          the triangulation must be a star triangulation and the
          distinguished central point is used as the origin.

        OUTPUT:

        A :class:`sage.geometry.fan.RationalPolyhedralFan`. The
        coordinates of the points are shifted so that the apex of
        the fan is the origin of the coordinate system.

        .. note:: If the set of cones over the simplices is not a fan, a
            suitable exception is raised.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0), (1,0), (0,1), (-1,-1)], star=0, fine=True)
            sage: pc.set_engine('internal')   # to make doctests independent of TOPCOM
            sage: triangulation = pc.triangulate()
            sage: fan = triangulation.fan(); fan
            Rational polyhedral fan in 2-d lattice N
            sage: fan.is_equivalent( toric_varieties.P2().fan() )
            True

        Toric diagrams (the `\ZZ_5` hyperconifold)::

            sage: vertices=[(0, 1, 0), (0, 3, 1), (0, 2, 3), (0, 0, 2)]
            sage: interior=[(0, 1, 1), (0, 1, 2), (0, 2, 1), (0, 2, 2)]
            sage: points = vertices+interior
            sage: pc = PointConfiguration(points, fine=True)
            sage: pc.set_engine('internal')   # to make doctests independent of TOPCOM
            sage: triangulation = pc.triangulate()
            sage: fan = triangulation.fan( (-1,0,0) )
            sage: fan
            Rational polyhedral fan in 3-d lattice N
            sage: fan.rays()
            (N(1, 1, 0), N(1, 3, 1), N(1, 2, 3), N(1, 0, 2),
             N(1, 1, 1), N(1, 1, 2), N(1, 2, 1), N(1, 2, 2))
        """
        from sage.geometry.fan import Fan
        if origin is None:
            origin = self.point_configuration().star_center()
        R = self.base_ring()
        origin = vector(R, origin)
        points = self.point_configuration().points()
        return Fan(self, (vector(R, p) - origin for p in points))



########################################################################
class PointConfiguration(UniqueRepresentation, PointConfiguration_base):
    """
    A collection of points in Euclidean (or projective) space.

    This is the parent class for the triangulations of the point
    configuration. There are a few options to specifically select what
    kind of triangulations are admissible.

    INPUT:

    The constructor accepts the following arguments:

    - ``points`` -- the points. Technically, any iterable of iterables
      will do. In particular, a :class:`PointConfiguration` can be passed.

    - ``projective`` -- boolean (default: ``False``). Whether the
      point coordinates should be interpreted as projective (``True``)
      or affine (``False``) coordinates. If necessary, points are
      projectivized by setting the last homogeneous coordinate to one
      and/or affine patches are chosen internally.

    - ``connected`` -- boolean (default: ``True``). Whether the
      triangulations should be connected to the regular triangulations
      via bistellar flips. These are much easier to compute than all
      triangulations.

    - ``fine`` -- boolean (default: ``False``). Whether the
      triangulations must be fine, that is, make use of all points of
      the configuration.

    - ``regular`` -- boolean or ``None`` (default: ``None``). Whether
      the triangulations must be regular. A regular triangulation is
      one that is induced by a piecewise-linear convex support
      function. In other words, the shadows of the faces of a
      polyhedron in one higher dimension.

        * ``True``: Only regular triangulations.

        * ``False``: Only non-regular triangulations.

        * ``None`` (default): Both kinds of triangulation.

    - ``star`` -- either ``None`` or a point. Whether the
      triangulations must be star. A triangulation is star if all
      maximal simplices contain a common point. The central point can
      be specified by its index (an integer) in the given points or by
      its coordinates (anything iterable.)

    EXAMPLES::

        sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
        sage: p
        A point configuration in QQ^2 consisting of 5 points. The
        triangulations of this point configuration are assumed to
        be connected, not necessarily fine, not necessarily regular.
        sage: p.triangulate()  # a single triangulation
        (<1,3,4>, <2,3,4>)
    """


    # we cache the output of _have_TOPCOM() in this class variable
    _have_TOPCOM_cached = None

    # whether to use TOPCOM. Will be set to True or False during
    # initialization. All implementations should check this boolean
    # variable to decide whether to call TOPCOM or not
    _use_TOPCOM = None


    @classmethod
    def _have_TOPCOM(cls):
        r"""
        Return whether TOPCOM is installed.

        EXAMPLES::

            sage: PointConfiguration._have_TOPCOM()    # optional - TOPCOM
            True
        """
        if PointConfiguration._have_TOPCOM_cached != None:
            return PointConfiguration._have_TOPCOM_cached

        try:
            out = PointConfiguration._TOPCOM_exec('points2placingtriang',
                                                  '[[0,1],[1,1]]', verbose=False).next()
            PointConfiguration._have_TOPCOM_cached = True
            assert out=='{{0,1}}',\
                'TOPCOM ran but did not produce the correct output!'
        except pexpect.ExceptionPexpect:
            PointConfiguration._have_TOPCOM_cached = False

        PointConfiguration.set_engine('auto')
        return PointConfiguration._have_TOPCOM_cached


    @staticmethod
    def __classcall__(cls, points, projective=False, connected=True, fine=False, regular=None, star=None):
        r"""
        Normalize the constructor arguments to be unique keys.

        EXAMPLES::

            sage: pc1 = PointConfiguration([[1,2],[2,3],[3,4]], connected=True)
            sage: pc2 = PointConfiguration(((1,2),(2,3),(3,4)), regular=None)
            sage: pc1 is pc2   # indirect doctest
            True
        """
        if isinstance(points, PointConfiguration):
            points = tuple( p.projective() for p in points )
            projective = True
        elif projective:
            points = tuple( tuple(p) for p in points )
        else:
            points = tuple( tuple(p)+(1,) for p in points )
        if star!=None and star not in ZZ:
            star_point = tuple(star)
            if len(star_point)<len(points[0]):
                star_point = tuple(star)+(1,)
            star = points.index(star_point)
        return super(PointConfiguration, cls)\
            .__classcall__(cls, points, connected, fine, regular, star)


    def __init__(self, points, connected, fine, regular, star):
        """
        Initialize a :class:`PointConfiguration` object.

        EXAMPLES::

            sage: p = PointConfiguration([[0,4],[2,3],[3,2],[4,0],[3,-2],[2,-3],[0,-4],[-2,-3],[-3,-2],[-4,0],[-3,2],[-2,3]])
            sage: len(p.triangulations_list())    # long time
            16796

        TESTS::

            sage: TestSuite(p).run()
        """
        # first, test if we have TOPCOM and set up class variables accordingly
        PointConfiguration._have_TOPCOM()

        assert connected in [True, False], 'Unknown value: connected='+str(connected)
        self._connected = connected
        if connected!=True and not PointConfiguration._have_TOPCOM():
            raise ValueError, 'You must install TOPCOM to find non-connected triangulations.'

        assert fine in [True, False], 'Unknown value: fine='+str(fine)
        self._fine = fine

        assert regular in [True, False, None], 'Unknown value: regular='+str(regular)
        self._regular = regular
        if regular!=None and not PointConfiguration._have_TOPCOM():
           raise ValueError, 'You must install TOPCOM to test for regularity.'

        assert star==None or star in ZZ, 'Unknown value: fine='+str(star)
        self._star = star

        PointConfiguration_base.__init__(self, points)


    @classmethod
    def set_engine(cls, engine='auto'):
        r"""
        Set the engine used to compute triangulations.

        INPUT:

        - ``engine`` -- either 'auto' (default), 'internal', or
          'TOPCOM'. The latter two instruct this package to always use
          its own triangulation algorithms or TOPCOM's algorithms,
          respectively. By default ('auto'), TOPCOM is used if it is
          available and internal routines otherwise.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: p.set_engine('internal')   # to make doctests independent of TOPCOM
            sage: p.triangulate()
            (<1,3,4>, <2,3,4>)
            sage: p.set_engine('TOPCOM')   # optional - TOPCOM
            sage: p.triangulate()          # optional - TOPCOM
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: p.set_engine('internal') # optional - TOPCOM
        """
        if engine not in ['auto', 'TOPCOM', 'internal']:
            raise ValueError, 'Unknown value for "engine": '+str(engine)

        have_TOPCOM = PointConfiguration._have_TOPCOM()
        PointConfiguration._use_TOPCOM = \
            (engine=='TOPCOM') or (engine=='auto' and have_TOPCOM)


    def star_center(self):
        r"""
        Return the center used for star triangulations.

        .. seealso:: :meth:`restrict_to_star_triangulations`.

        OUTPUT:

        A :class:`sage.geometry.triangulation.base.Point` if a
        distinguished star central point has been fixed.
        ``ValueError`` exception is raised otherwise.

        EXAMPLES::

            sage: pc = PointConfiguration([(1,0),(-1,0),(0,1),(0,2)], star=(0,1)); pc
            A point configuration in QQ^2 consisting of 4 points. The
            triangulations of this point configuration are assumed to be
            connected, not necessarily fine, not necessarily regular, and
            star with center P(0, 1).
            sage: pc.star_center()
            P(0, 1)

            sage: pc_nostar = pc.restrict_to_star_triangulations(None)
            sage: pc_nostar
            A point configuration in QQ^2 consisting of 4 points. The
            triangulations of this point configuration are assumed to be
            connected, not necessarily fine, not necessarily regular.
            sage: pc_nostar.star_center()
            Traceback (most recent call last):
            ...
            ValueError: The point configuration has no star center defined.
        """
        if self._star is None:
            raise ValueError, 'The point configuration has no star center defined.'
        else:
            return self[self._star]


    def __reduce__(self):
        r"""
        Override __reduce__ to correctly pickle/unpickle.

        TESTS::

            sage: p = PointConfiguration([[0, 1], [0, 0], [1, 0], [1,1]])
            sage: loads(p.dumps()) is p
            True
        """
        points = tuple( p.projective() for p in self )
        return (PointConfiguration, (points, True, self._connected, self._fine, self._regular, self._star))


    def an_element(self):
        """
        Synonymous for :meth:`triangulate`.

        TESTS::

            sage: p = PointConfiguration([[0, 1], [0, 0], [1, 0], [1,1]])
            sage: p.an_element()
            (<0,1,3>, <1,2,3>)
        """
        return self.triangulate()


    def _element_constructor_(self, e):
        """
        Construct a triangulation.

        TESTS::

            sage: p = PointConfiguration([[0, 1], [0, 0], [1, 0], [1,1]])
            sage: p._element_constructor_([ (0,1,2), (2,3,0) ])
            (<0,1,2>, <0,2,3>)
        """
        return self.element_class(e, parent=self)


    Element = Triangulation


    def __iter__(self):
        """
        Iterate through the points of the point configuration.

        OUTPUT:

        Returns projective coordinates of the points. See also the
        ``PointConfiguration.points()`` method, which returns affine
        coordinates.

        EXAMPLES::

            sage: p = PointConfiguration([[1,1], [2,2], [3,3]]);
            sage: list(p)     # indirect doctest
            [P(1, 1), P(2, 2), P(3, 3)]
            sage: [ p[i] for i in range(0,p.n_points()) ]
            [P(1, 1), P(2, 2), P(3, 3)]
            sage: list(p.points())
            [P(1, 1), P(2, 2), P(3, 3)]
            sage: [ p.point(i) for i in range(0,p.n_points()) ]
            [P(1, 1), P(2, 2), P(3, 3)]
        """
        for p in self.points():
            yield p


    def _repr_(self):
        r"""
        Return a string representation.

        TESTS::

            sage: p = PointConfiguration([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [-1, -1, 1], [1, 1, -1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1], [0, 0, 0]])
            sage: p._repr_()
            'A point configuration in QQ^3 consisting of 9 points. The triangulations of this point configuration are assumed to be connected, not necessarily fine, not necessarily regular.'
        """
        s = 'A point configuration'
        s += ' in QQ^'+str(self.ambient_dim())
        if len(self)==1:
            s += ' consisting of '+str(len(self))+' point. '
        else:
            s += ' consisting of '+str(len(self))+' points. '

        s += 'The triangulations of this point configuration are assumed to be'

        if self._connected:
            s += ' connected,'
        else:
            s += ' not necessarily connected,'

        if self._fine:
            s += ' fine,'
        else:
            s += ' not necessarily fine,'

        if self._regular==True:
            s += ' regular'
        elif self._regular==False:
            s += ' irregular'
        else:
            s += ' not necessarily regular'

        if self._star==None:
            s += '.'
        else:
            s += ', and star with center '+str(self.star_center())+'.'
        return s


    def _TOPCOM_points(self):
        r"""
        Convert the list of input points to a string that can be fed
        to TOPCOM.

        TESTS::

            sage: p = PointConfiguration([[1,1,1], [-1,1,1], [1,-1,1], [-1,-1,1], [1,1,-1]])
            sage: p._TOPCOM_points()
            '[[0,0,0,1],[-2,0,0,1],[0,-2,0,1],[-2,-2,0,1],[0,0,-2,1]]'
        """
        s = '['
        s += ','.join([
                '[' + ','.join(map(str,p.reduced_projective())) + ']'
                for p in self ])
        s += ']'
        return s


    @classmethod
    def _TOPCOM_exec(cls, executable, input_string, verbose=True):
        r"""
        Run TOPCOM.

        INPUT:

        - ``executable`` -- string. The name of the executable.

        - ``input_string`` -- string. Will be piped into the running
          executable's stdin.

        - ``verbose`` -- boolean. Whether to print out the TOPCOM
          interaction.

        TESTS::

            sage: p = PointConfiguration([[1,1,1], [-1,1,1], [1,-1,1], [-1,-1,1], [1,1,-1]])
            sage: out = p._TOPCOM_exec('points2placingtriang', '[[0,0,0,1],[-2,0,0,1],[0,-2,0,1],[-2,-2,0,1],[0,0,-2,1]]', verbose=True)
            sage: list(out)       # optional - TOPCOM
            #### TOPCOM input ####
            # points2placingtriang
            # [[0,0,0,1],[-2,0,0,1],[0,-2,0,1],[-2,-2,0,1],[0,0,-2,1]]
            #### TOPCOM output ####
            # {{0,1,2,4},{1,2,3,4}}
            #######################
            ['{{0,1,2,4},{1,2,3,4}}']
        """
        timeout = 600
        proc = pexpect.spawn(executable, timeout=timeout)
        proc.expect('Evaluating Commandline Options \.\.\.')
        proc.expect('\.\.\. done\.')
        proc.setecho(0)
        assert proc.readline().strip() == ''

        if verbose:
            print "#### TOPCOM input ####"
            print "# " + executable
            print "# " + input_string
            sys.stdout.flush()

        proc.send(input_string)
        proc.send('X\nX\n')

        if verbose:
            print "#### TOPCOM output ####"
            sys.stdout.flush()

        while True:
            try:
                line = proc.readline().strip()
            except pexpect.TIMEOUT:
                if verbose:
                    print '# Still runnnig '+str(executable)
                continue
            if len(line)==0: # EOF
                break;
            if verbose:
                print "# " + line
                sys.stdout.flush()

            try:
                yield line.strip()
            except GeneratorExit:
                proc.close(force=True)
                raise StopIteration

        if verbose:
            print "#######################"
            sys.stdout.flush()


    def _TOPCOM_communicate(self, executable, verbose=True):
        r"""
        Execute TOPCOM and parse the output into a :class:`Triangulation`.

        TESTS::

            sage: p = PointConfiguration([[1,1,1], [-1,1,1], [1,-1,1], [-1,-1,1], [1,1,-1]])
            sage: out = p._TOPCOM_communicate('points2placingtriang', verbose=True)
            sage: list(out)       # optional - TOPCOM
            #### TOPCOM input ####
            # points2placingtriang
            # [[0,0,0,1],[-2,0,0,1],[0,-2,0,1],[-2,-2,0,1],[0,0,-2,1]]
            #### TOPCOM output ####
            # {{0,1,2,4},{1,2,3,4}}
            #######################
            [(<0,1,2,4>, <1,2,3,4>)]
        """
        for line in self._TOPCOM_exec(executable,
                                      self._TOPCOM_points(), verbose):
            triangulation = line[ line.find('{{')+2 : line.rfind('}}') ]
            triangulation = triangulation.split('},{')
            triangulation = [ [ QQ(t) for t in triangle.split(',') ]
                              for triangle in triangulation ]

            if self._star!=None:
                o = self._star
                if not all( t.count(o)>0 for t in triangulation):
                    continue

            yield self(triangulation)


    def _TOPCOM_triangulations(self, verbose=True):
        r"""
        Returns all triangulations satisfying the restrictions imposed.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: iter = p._TOPCOM_triangulations(verbose=True)
            sage: iter.next()     # optional - TOPCOM
            #### TOPCOM input ####
            # points2triangs
            # [[0,0,1],[0,1,1],[1,0,1],[1,1,1],[-1,-1,1]]
            #### TOPCOM output ####
            # T[1]:=[5,3:{{0,1,2},{1,2,3},{0,2,4},{0,1,4}}];
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
        """
        command = 'points2'

        if not self._connected:
            command += 'all'

        if self._fine:
            command += 'fine'

        command += 'triangs'

        if self._regular==True:
            command += ' --regular'
        if self._regular==False:
            command += ' --nonregular'

        for t in self._TOPCOM_communicate(command, verbose):
            yield t


    def _TOPCOM_triangulate(self, verbose=True):
        r"""
        Return one (in no particular order) triangulation subject
        to all restrictions imposed previously.

        INPUT:

        - ``verbose`` -- boolean. Whether to print out the TOPCOM
          interaction.

        OUTPUT:

        A :class:`sage.geometry.triangulation.Triangulation`
        satisfying all restrictions imposed. Raises a ``ValueError``
        if no such triangulation exists.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: p.set_engine('TOPCOM')                 # optional - TOPCOM
            sage: p._TOPCOM_triangulate(verbose=False)   # optional - TOPCOM
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: list( p.triangulate() )                # optional - TOPCOM
            [(0, 1, 2), (0, 1, 4), (0, 2, 4), (1, 2, 3)]
            sage: p.set_engine('internal')               # optional - TOPCOM
        """
        assert self._regular!=False, \
            'When asked for a single triangulation TOPCOM ' + \
            'always returns a regular triangulation.'

        command = "points2"
        if self._fine:
            command += "finetriang"
        else:
            command += "placingtriang"

        return self._TOPCOM_communicate(command, verbose).next()


    def restrict_to_regular_triangulations(self, regular=True):
        """
        Restrict to regular triangulations.

        NOTE:

        Regularity testing requires the optional TOPCOM package.

        INPUT:

        - ``regular`` -- ``True``, ``False``, or ``None``. Whether to
          restrict to regular triangulations, irregular
          triangulations, or lift any restrictions on regularity.

        OUTPUT:

        A new :class:`PointConfiguration` with the same points, but
        whose triangulations will all be regular as specified. See
        :class:`PointConfiguration` for details.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: p
            A point configuration in QQ^2 consisting of 5 points. The
            triangulations of this point configuration are assumed to
            be connected, not necessarily fine, not necessarily regular.
            sage: len(p.triangulations_list())
            4
            sage: p_regular = p.restrict_to_regular_triangulations() # optional - TOPCOM
            sage: len(p_regular.triangulations_list())               # optional - TOPCOM
            4
            sage: p == p_regular.restrict_to_regular_triangulations(regular=None) # optional - TOPCOM
            True
        """
        return PointConfiguration(self,
                                  connected=self._connected,
                                  fine=self._fine,
                                  regular=regular,
                                  star=self._star)


    def restrict_to_connected_triangulations(self, connected=True):
        """
        Restrict to connected triangulations.

        NOTE:

        Finding non-connected triangulations requires the optional
        TOPCOM package.

        INPUT:

        - ``connected`` -- boolean. Whether to restrict to
          triangulations that are connected by bistellar flips to the
          regular triangulations.

        OUTPUT:

        A new :class:`PointConfiguration` with the same points, but
        whose triangulations will all be in the connected
        component. See :class:`PointConfiguration` for details.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: p
            A point configuration in QQ^2 consisting of 5 points. The
            triangulations of this point configuration are assumed to
            be connected, not necessarily fine, not necessarily regular.
            sage: len(p.triangulations_list())
            4
            sage: p_all = p.restrict_to_connected_triangulations(connected=False)  # optional - TOPCOM
            sage: len(p_all.triangulations_list())                                 # optional - TOPCOM
            4
            sage: p == p_all.restrict_to_connected_triangulations(connected=True)  # optional - TOPCOM
            True
        """
        return PointConfiguration(self,
                                  connected=connected,
                                  fine=self._fine,
                                  regular=self._regular,
                                  star=self._star)


    def restrict_to_fine_triangulations(self, fine=True):
        """
        Restrict to fine triangulations.

        INPUT:

        - ``fine`` -- boolean. Whether to restrict to fine triangulations.

        OUTPUT:

        A new :class:`PointConfiguration` with the same points, but
        whose triangulations will all be fine. See
        :class:`PointConfiguration` for details.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: p
            A point configuration in QQ^2 consisting of 5 points. The
            triangulations of this point configuration are assumed to
            be connected, not necessarily fine, not necessarily regular.
            sage: len(p.triangulations_list())
            4
            sage: p_fine = p.restrict_to_fine_triangulations()
            sage: len(p.triangulations_list())
            4
            sage: p == p_fine.restrict_to_fine_triangulations(fine=False)
            True
        """
        return PointConfiguration(self,
                                  connected=self._connected,
                                  fine=fine,
                                  regular=self._regular,
                                  star=self._star)


    def restrict_to_star_triangulations(self, star):
        """
        Restrict to star triangulations with the given point as the
        center.

        INPUT:

        - ``origin`` -- ``None`` or an integer or the coordinates of a
          point. An integer denotes the index of the central point. If
          ``None`` is passed, any restriction on the starshape will be
          removed.

        OUTPUT:

        A new :class:`PointConfiguration` with the same points, but
        whose triangulations will all be star. See
        :class:`PointConfiguration` for details.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: len(list( p.triangulations() ))
            4
            sage: p_star =  p.restrict_to_star_triangulations(0)
            sage: p_star is p.restrict_to_star_triangulations((0,0))
            True
            sage: p_star.triangulations_list()
            [(<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)]
            sage: p_newstar = p_star.restrict_to_star_triangulations(1)  # pick different origin
            sage: p_newstar.triangulations_list()
            [(<1,2,3>, <1,2,4>)]
            sage: p == p_star.restrict_to_star_triangulations(star=None)
            True
        """
        return PointConfiguration(self,
                                  connected=self._connected,
                                  fine=self._fine,
                                  regular=self._regular,
                                  star=star)


    def triangulations(self, verbose=False):
        r"""
        Returns all triangulations.

        - ``verbose`` -- boolean (default: ``False``). Whether to
          print out the TOPCOM interaction, if any.

        OUTPUT:

        A generator for the triangulations satisfying all the
        restrictions imposed. Each triangulation is returned as a
        :class:`Triangulation` object.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: iter = p.triangulations()
            sage: iter.next()
            (<1,3,4>, <2,3,4>)
            sage: iter.next()
            (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)
            sage: iter.next()
            (<1,2,3>, <1,2,4>)
            sage: iter.next()
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: p.triangulations_list()
            [(<1,3,4>, <2,3,4>),
             (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
             (<1,2,3>, <1,2,4>),
             (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)]
            sage: p_fine = p.restrict_to_fine_triangulations()
            sage: p_fine.triangulations_list()
            [(<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
             (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)]

         Note that we explicitly asked the internal algorithm to
         compute the triangulations. Using TOPCOM, we obtain the same
         triangulations but in a different order::

            sage: p.set_engine('TOPCOM')                       # optional - TOPCOM
            sage: iter = p.triangulations()                    # optional - TOPCOM
            sage: iter.next()                                  # optional - TOPCOM
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: iter.next()                                  # optional - TOPCOM
            (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)
            sage: iter.next()                                  # optional - TOPCOM
            (<1,2,3>, <1,2,4>)
            sage: iter.next()                                  # optional - TOPCOM
            (<1,3,4>, <2,3,4>)
            sage: p.triangulations_list()                      # optional - TOPCOM
            [(<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>),
             (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>),
             (<1,2,3>, <1,2,4>),
             (<1,3,4>, <2,3,4>)]
            sage: p_fine = p.restrict_to_fine_triangulations() # optional - TOPCOM
            sage: p_fine.set_engine('TOPCOM')                  # optional - TOPCOM
            sage: p_fine.triangulations_list()                 # optional - TOPCOM
            [(<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>),
             (<0,1,3>, <0,1,4>, <0,2,3>, <0,2,4>)]
            sage: p.set_engine('internal')                     # optional - TOPCOM
        """
        if self._use_TOPCOM:
            for triangulation in self._TOPCOM_triangulations(verbose):
                yield triangulation
        else:
            if (self._connected!=True):
                raise ValueError, 'Need TOPCOM to find disconnected triangulations.'
            if (self._regular!=None):
                raise ValueError, 'Need TOPCOM to test for regularity.'
            ci = ConnectedTriangulationsIterator(self, star=self._star, fine=self._fine)
            for encoded_triangulation in ci:
                yield self(encoded_triangulation)


    def triangulations_list(self, verbose=False):
        r"""
        Return all triangulations.

        INPUT:

        - ``verbose`` -- boolean. Whether to print out the TOPCOM
          interaction, if any.

        OUTPUT:

        A list of triangulations (see :class:`Triangulation`)
        satisfying all restrictions imposed previously.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1]])
            sage: p.triangulations_list()
            [(<0,1,2>, <1,2,3>), (<0,1,3>, <0,2,3>)]
            sage: map(list, p.triangulations_list() )
            [[(0, 1, 2), (1, 2, 3)], [(0, 1, 3), (0, 2, 3)]]
            sage: p.set_engine('TOPCOM')       # optional - TOPCOM
            sage: p.triangulations_list()      # optional - TOPCOM
            [(<0,1,2>, <1,2,3>), (<0,1,3>, <0,2,3>)]
            sage: p.set_engine('internal')     # optional - TOPCOM
        """
        return list(self.triangulations(verbose))


    def triangulate(self, verbose=False):
        r"""
        Return one (in no particular order) triangulation.

        INPUT:

        - ``verbose`` -- boolean. Whether to print out the TOPCOM
          interaction, if any.

        OUTPUT:

        A :class:`Triangulation` satisfying all restrictions
        imposed. Raises a ``ValueError`` if no such triangulation
        exists.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: p.triangulate()
            (<1,3,4>, <2,3,4>)
            sage: list( p.triangulate() )
            [(1, 3, 4), (2, 3, 4)]

        Using TOPCOM yields a different, but equally good, triangulation::

            sage: p.set_engine('TOPCOM')           # optional - TOPCOM
            sage: p.triangulate()                  # optional - TOPCOM
            (<0,1,2>, <0,1,4>, <0,2,4>, <1,2,3>)
            sage: list( p.triangulate() )          # optional - TOPCOM
            [(0, 1, 2), (0, 1, 4), (0, 2, 4), (1, 2, 3)]
            sage: p.set_engine('internal')         # optional - TOPCOM
        """
        if self._use_TOPCOM and self._regular!=False:
            try:
                return self._TOPCOM_triangulate(verbose)
            except StopIteration:
                # either topcom did not return a triangulation or we filtered it out
                pass

        if self._connected and not self._fine and self._regular!=False and self._star==None:
            return self.lexicographic_triangulation()

        try:
            return self.triangulations(verbose).next()
        except StopIteration:
            # there is no triangulation
            pass
        raise ValueError, 'No triangulation with the required properties.'


    def convex_hull(self):
        """
        Return the convex hull of the point configuration.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: p.convex_hull()
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices.
        """
        try:
            return self._polyhedron
        except AttributeError:
            pass

        from sage.geometry.polyhedra import Polyhedron
        pts = [ p.reduced_affine() for p in self.points() ];
        self._polyhedron = Polyhedron(vertices=pts);
        return self._polyhedron


    def restricted_automorphism_group(self):
        r"""
        Return the restricted automorphism group.

        First, let the linear automorphism group be the subgroup of
        the Euclidean group `E(d) = GL(d,\RR) \ltimes \RR^d`
        preserving the `d`-dimensional point configuration. The
        Euclidean group acts in the usual way `\vec{x}\mapsto
        A\vec{x}+b` on the ambient space.

        The restricted automorphism group is the subgroup of the
        linear automorphism group generated by permutations of
        points. See [BSS]_ for more details and a description of the
        algorithm.

        OUTPUT:

        A
        :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        that is isomorphic to the restricted automorphism group is
        returned.

        Note that in Sage, permutation groups always act on positive
        integers while lists etc. are indexed by nonnegative
        integers. The indexing of the permutation group is chosen to
        be shifted by ``+1``. That is, the transposition ``(i,j)`` in
        the permutation group corresponds to exchange of `self[i-1]``
        and `self[j-1]`.

        EXAMPLES::

            sage: pyramid = PointConfiguration([[1,0,0],[0,1,1],[0,1,-1],[0,-1,-1],[0,-1,1]])
            sage: pyramid.restricted_automorphism_group()
            Permutation Group with generators [(3,5), (2,3)(4,5), (2,4)]
            sage: DihedralGroup(4).is_isomorphic(_)
            True

        The square with an off-center point in the middle. Note thath
        the middle point breaks the restricted automorphism group
        `D_4` of the convex hull::

            sage: square = PointConfiguration([(3/4,3/4),(1,1),(1,-1),(-1,-1),(-1,1)])
            sage: square.restricted_automorphism_group()
            Permutation Group with generators [(3,5)]
            sage: DihedralGroup(1).is_isomorphic(_)
            True
        """
        if '_restricted_automorphism_group' in self.__dict__:
            return self._restricted_automorphism_group

        v_list = [ vector(p.projective()) for p in self ]
        Qinv = sum( v.column() * v.row() for v in v_list ).inverse()

        # construct the graph
        from sage.graphs.graph import Graph
        G = Graph(dense=True)
        for i in range(0,len(v_list)):
            for j in range(i+1,len(v_list)):
                v_i = v_list[i]
                v_j = v_list[j]
                G.add_edge(i,j, v_i * Qinv * v_j)

        group, node_dict = G.automorphism_group(edge_labels=True, translation=True)

        # Relabel the permutation group
        perm_to_vertex = dict( (i,v+1) for v,i in node_dict.items() )
        group = PermutationGroup([ [ tuple([ perm_to_vertex[i] for i in cycle ])
                                     for cycle in generator.cycle_tuples() ]
                                   for generator in group.gens() ])

        self._restricted_automorphism_group = group
        return group


    def face_codimension(self, point):
        r"""
        Return the smallest `d\in\mathbb{Z}` such that ``point`` is
        contained in the interior of a codimension-`d` face.

        EXAMPLES::

            sage: triangle = PointConfiguration([[0,0], [1,-1], [1,0], [1,1]]);
            sage: triangle.point(2)
            P(1, 0)
            sage: triangle.face_codimension(2)
            1
            sage: triangle.face_codimension( [1,0] )
            1

        This also works for degenerate cases like the tip of the
        pyramid over a square (which saturates four inequalities)::

            sage: pyramid = PointConfiguration([[1,0,0],[0,1,1],[0,1,-1],[0,-1,-1],[0,-1,1]])
            sage: pyramid.face_codimension(0)
            3
        """
        try:
            p = vector(self.point(point).reduced_affine())
        except TypeError:
            p = vector(point);

        inequalities = []
        for ieq in self.convex_hull().inequality_generator():
            if (ieq.A()*p + ieq.b() == 0):
                inequalities += [ ieq.vector() ];
        return matrix(inequalities).rank();


    def face_interior(self, dim=None, codim=None):
        """
        Return points by the codimension of the containing face in the convex hull.

        EXAMPLES::

            sage: triangle = PointConfiguration([[-1,0], [0,0], [1,-1], [1,0], [1,1]]);
            sage: triangle.face_interior()
            ((1,), (3,), (0, 2, 4))
            sage: triangle.face_interior(dim=0)    # the vertices of the convex hull
            (0, 2, 4)
            sage: triangle.face_interior(codim=1)  # interior of facets
            (3,)
        """
        assert not (dim!=None and codim!=None), "You cannot specify both dim and codim."

        if (dim!=None):
            return self.face_interior()[self.convex_hull().dim()-dim]
        if (codim!=None):
            return self.face_interior()[codim]

        try:
            return self._face_interior
        except AttributeError:
            pass

        d = [ self.face_codimension(i) for i in range(0,self.n_points()) ]

        return tuple( tuple(filter( lambda i: d[i]==codim, range(0,self.n_points())) )
                      for codim in range(0,self.dim()+1) )


    def exclude_points(self, point_idx_list):
        """
        Return a new point configuration with the given points
        removed.

        INPUT:

        - ``point_idx_list`` -- a list of integers. The indices of
          points to exclude.

        OUTPUT:

        A new :class:`PointConfiguration` with the given points
        removed.

        EXAMPLES::

            sage: p = PointConfiguration([[-1,0], [0,0], [1,-1], [1,0], [1,1]]);
            sage: list(p)
            [P(-1, 0), P(0, 0), P(1, -1), P(1, 0), P(1, 1)]
            sage: q = p.exclude_points([3])
            sage: list(q)
            [P(-1, 0), P(0, 0), P(1, -1), P(1, 1)]
            sage: p.exclude_points( p.face_interior(codim=1) ).points()
            (P(-1, 0), P(0, 0), P(1, -1), P(1, 1))
        """
        points = [ self.point(i) for i in range(0,self.n_points())
                   if not i in point_idx_list ]
        return PointConfiguration(points,
                                  projective=False,
                                  connected=self._connected,
                                  fine=self._fine,
                                  regular=self._regular,
                                  star=self._star)


    def volume(self, simplex=None):
        """
        Find n! times the n-volume of a simplex of dimension n.

        INPUT:

        - ``simplex`` (optional argument) -- a simplex from a
          triangulation T specified as a list of point indices.

        OUTPUT:

        * If a simplex was passed as an argument: n!*(volume of the simplex simp).

        * Without argument: n!*(the total volume of the convex hull).

        EXAMPLES:

        The volume of the standard simplex should always be 1::

            sage: p = PointConfiguration([[0,0],[1,0],[0,1],[1,1]])
            sage: p.volume( [0,1,2] )
            1
            sage: simplex = p.triangulate()[0]  # first simplex of triangulation
            sage: p.volume(simplex)
            1

        The square can be triangulated into two minimal simplices, so
        in the "integral" normalization ts volume equals two::

            sage: p.volume()
            2

        NOTES:

        We need to have n!*(volume of the simplex) to ensure that
        the volume is an integer.  Essentially, this normalizes things so that
        the volume of the standard n-simplex is 1.  See [GKZ]_ page 182.
        """
        if (simplex==None):
            return sum([ self.volume(s) for s in self.triangulate() ])

        #Form a matrix whose columns are the points of simplex
        #with the first point of simplex shifted to the origin.
        v = [ vector(self.point(i).reduced_affine()) for i in simplex ]
        m = matrix([ v_i - v[0] for v_i in v[1:] ]).transpose()
        return abs(m.det())


    def secondary_polytope(self):
        r"""
        Calculate the secondary polytope of the point configuration.

        For a definition of the secondary polytope, see [GKZ]_ page 220
        Definition 1.6.

        Note that if you restricted the admissible triangulations of
        the point configuration then the output will be the
        corresponding face of the whole secondary polytope.

        OUTPUT:

        The secondary polytope of the point configuration as an
        instance of
        :class:`sage.geometry.polyhedra.Polyhedron`.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[1,0],[2,1],[1,2],[0,1]])
            sage: poly = p.secondary_polytope()
            sage: matrix(poly.vertices()).transpose()
            [1 3 1 5 3]
            [3 1 5 1 4]
            [4 5 2 4 2]
            [2 2 4 4 5]
            [5 4 3 1 1]
            sage: poly.Vrepresentation()
            [A vertex at (1, 3, 4, 2, 5),
             A vertex at (3, 1, 5, 2, 4),
             A vertex at (1, 5, 2, 4, 3),
             A vertex at (5, 1, 4, 4, 1),
             A vertex at (3, 4, 2, 5, 1)]
            sage: poly.Hrepresentation()
            [An equation (1/2, 1, 1, 0, 0) x - 15/2 == 0,
             An equation (-1, -1, 0, 1, 0) x + 2 == 0,
             An equation (3/2, 1, 0, 0, 1) x - 19/2 == 0,
             An inequality (-1, -2, 0, 0, 0) x + 11 >= 0,
             An inequality (1, 0, 0, 0, 0) x - 1 >= 0,
             An inequality (1, 1, 0, 0, 0) x - 4 >= 0,
             An inequality (0, 1, 0, 0, 0) x - 1 >= 0,
             An inequality (-3/2, -1, 0, 0, 0) x + 17/2 >= 0]
        """
        from sage.geometry.polyhedra import Polyhedron
        #TODO: once restriction to regular triangulations is fixed,
        #change the next line to only take the regular triangulations,
        #since they are the vertices of the secondary polytope anyway.
        l = self.triangulations_list()
        return Polyhedron(vertices = [x.gkz_phi() for x in l])


    def circuits_support(self):
        r"""
        A generator for the supports of the circuits of the point configuration.

        See :meth:`circuits` for details.

        OUTPUT:

        A generator for the supports `C_-\cup C_+` (returned as a
        Python tuple) for all circuits of the point configuration.

        EXAMPLES::

            sage: p = PointConfiguration([(0,0),(+1,0),(-1,0),(0,+1),(0,-1)])
            sage: list( p.circuits_support() )
            [(0, 3, 4), (0, 1, 2), (1, 2, 3, 4)]
        """
        n = len(self)
        U = [ self[i].reduced_projective() for i in range(0,n) ]

        # the index set of U
        I = set(range(0,n))
        # The (indices of) known independent elements of U
        independent_k = [ (i,) for i in range(0,n) ]
        supports_k = []

        supports = ()   # supports of circuits
        for k in range(2, self.dim()+3):

            # possibly linear dependent subsets
            supports_knext = set()
            possible_dependency = set()
            for indep in independent_k:
                indep_plus_one = [ tuple(sorted(indep+(i,))) for i in (I-set(indep)) ]
                possible_dependency.update(indep_plus_one)
            for supp in supports_k:
                supp_plus_one = [ tuple(sorted(supp+(i,))) for i in (I-set(supp)) ]
                possible_dependency.difference_update(supp_plus_one)
                supports_knext.update(supp_plus_one)

            # remember supports and independents for the next k-iteration
            supports_k = list(supports_knext)
            independent_k = []
            for idx in possible_dependency:
                rk = matrix([ U[i] for i in idx ]).rank()
                if rk==k:
                    independent_k.append(idx)
                else:
                    supports_k.append(idx)
                    yield idx
        assert independent_k==[]  # there are no independent (self.dim()+3)-tuples


    def circuits(self):
        r"""
        Return the circuits of the point configuration.

        Roughly, a circuit is a minimal linearly dependent subset of
        the points. That is, a circuit is a partition

        .. MATH::

            \{ 0, 1, \dots, n-1 \} = C_+ \cup C_0 \cup C_-

        such that there is an (unique up to an overall normalization) affine
        relation

        .. MATH::

            \sum_{i\in C_+}  \alpha_i \vec{p}_i =
            \sum_{j\in C_-}  \alpha_j \vec{p}_j

        with all positive (or all negative) coefficients, where
        `\vec{p}_i=(p_1,\dots,p_k,1)` are the projective coordinates
        of the `i`-th point.

        OUTPUT:

        The list of (unsigned) circuits as triples `(C_+, C_0,
        C_-)`. The swapped circuit `(C_-, C_0, C_+)` is not returned
        separately.

        EXAMPLES::

            sage: p = PointConfiguration([(0,0),(+1,0),(-1,0),(0,+1),(0,-1)])
            sage: p.circuits()
            (((0,), (1, 2), (3, 4)), ((0,), (3, 4), (1, 2)), ((1, 2), (0,), (3, 4)))


        TESTS::

            sage: U=matrix([
            ...      [ 0, 0, 0, 0, 0, 2, 4,-1, 1, 1, 0, 0, 1, 0],
            ...      [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0],
            ...      [ 0, 2, 0, 0, 0, 0,-1, 0, 1, 0, 1, 0, 0, 1],
            ...      [ 0, 1, 1, 0, 0, 1, 0,-2, 1, 0, 0,-1, 1, 1],
            ...      [ 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0]
            ...   ])
            sage: p = PointConfiguration(U.columns())
            sage: len( p.circuits() )    # long time
            218
        """
        try:
            return self._circuits
        except AttributeError:
            pass

        n = len(self)
        U = [ self[i].reduced_projective() for i in range(0,n) ]

        Circuits = ()
        for support in self.circuits_support():
            m = matrix([ U[i] for i in support ]).transpose()
            ker = m.right_kernel().basis()[0]
            assert len(ker)==len(support)
            Cplus  = [ support[i] for i in range(0,len(support)) if ker[i]>0 ]
            Cminus = [ support[i] for i in range(0,len(support)) if ker[i]<0 ]
            Czero  = set( range(0,n) ).difference(support)
            Circuits += ( (tuple(Cplus), tuple(Czero), tuple(Cminus)), )
        self._circuits = Circuits
        return Circuits


    def positive_circuits(self, *negative):
        r"""
        Returns the positive part of circuits with fixed negative part.

        A circuit is a pair `(C_+, C_-)`, each consisting of a subset
        (actually, an ordered tuple) of point indices.

        INPUT:

        - ``*negative`` -- integer. The indices of points.

        OUTPUT:

        A tuple of all circuits with `C_- = ` ``negative``.

        EXAMPLE::

            sage: p = PointConfiguration([(1,0,0),(0,1,0),(0,0,1),(-2,0,-1),(-2,-1,0),(-3,-1,-1),(1,1,1),(-1,0,0),(0,0,0)])
            sage: p.positive_circuits(8)
            ((0, 7), (0, 1, 4), (0, 2, 3), (0, 5, 6), (0, 1, 2, 5), (0, 3, 4, 6))
            sage: p.positive_circuits(0,5,6)
            ((8,),)
        """
        pos = ()
        negative = tuple(sorted(negative))
        for circuit in self.circuits():
            Cpos = circuit[0]
            Cneg = circuit[2]
            if Cpos == negative:
                pos += ( Cneg, )
            elif Cneg == negative:
                pos += ( Cpos, )
        return pos


    def bistellar_flips(self):
        r"""
        Return the bistellar flips.

        OUTPUT:

        The bistellar flips as a tuple. Each flip is a pair
        `(T_+,T_-)` where `T_+` and `T_-` are partial triangulations
        of the point configuration.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0),(1,0),(0,1),(1,1)])
            sage: pc.bistellar_flips()
            (((<0,1,3>, <0,2,3>), (<0,1,2>, <1,2,3>)),)
            sage: Tpos, Tneg = pc.bistellar_flips()[0]
            sage: Tpos.plot(axes=False)
            sage: Tneg.plot(axes=False)

        The 3d analog::

            sage: pc = PointConfiguration([(0,0,0),(0,2,0),(0,0,2),(-1,0,0),(1,1,1)])
            sage: pc.bistellar_flips()
            (((<0,1,2,3>, <0,1,2,4>), (<0,1,3,4>, <0,2,3,4>, <1,2,3,4>)),)

        A 2d flip on the base of the pyramid over a square::

            sage: pc = PointConfiguration([(0,0,0),(0,2,0),(0,0,2),(0,2,2),(1,1,1)])
            sage: pc.bistellar_flips()
            (((<0,1,3>, <0,2,3>), (<0,1,2>, <1,2,3>)),)
            sage: Tpos, Tneg = pc.bistellar_flips()[0]
            sage: Tpos.plot(axes=False)
        """
        flips = []
        for C in self.circuits():
            Cpos = list(C[0])
            Cneg = list(C[2])
            support = sorted(Cpos+Cneg)
            Tpos = [ Cpos+Cneg[0:i]+Cneg[i+1:len(Cneg)] for i in range(0,len(Cneg)) ]
            Tneg = [ Cneg+Cpos[0:i]+Cpos[i+1:len(Cpos)] for i in range(0,len(Cpos)) ]
            flips.append( (self.element_class(Tpos, parent=self, check=False),
                           self.element_class(Tneg, parent=self, check=False)) )
        return tuple(flips)


    def lexicographic_triangulation(self):
        r"""
        Return the lexicographic triangulation.

        The algorithm was taken from [PUNTOS]_.

        EXAMPLES::

            sage: p = PointConfiguration([(0,0),(+1,0),(-1,0),(0,+1),(0,-1)])
            sage: p.lexicographic_triangulation()
            (<1,3,4>, <2,3,4>)

        TESTS::

            sage: U=matrix([
            ...      [ 0, 0, 0, 0, 0, 2, 4,-1, 1, 1, 0, 0, 1, 0],
            ...      [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0],
            ...      [ 0, 2, 0, 0, 0, 0,-1, 0, 1, 0, 1, 0, 0, 1],
            ...      [ 0, 1, 1, 0, 0, 1, 0,-2, 1, 0, 0,-1, 1, 1],
            ...      [ 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0]
            ...   ])
            sage: pc = PointConfiguration(U.columns())
            sage: pc.lexicographic_triangulation()
            (<1,3,4,7,10,13>, <1,3,4,8,10,13>, <1,3,6,7,10,13>, <1,3,6,8,10,13>,
             <1,4,6,7,10,13>, <1,4,6,8,10,13>, <2,3,4,6,7,12>, <2,3,4,7,12,13>,
             <2,3,6,7,12,13>, <2,4,6,7,12,13>, <3,4,5,6,9,12>, <3,4,5,8,9,12>,
             <3,4,6,7,11,12>, <3,4,6,9,11,12>, <3,4,7,10,11,13>, <3,4,7,11,12,13>,
             <3,4,8,9,10,12>, <3,4,8,10,12,13>, <3,4,9,10,11,12>, <3,4,10,11,12,13>,
             <3,5,6,8,9,12>, <3,6,7,10,11,13>, <3,6,7,11,12,13>, <3,6,8,9,10,12>,
             <3,6,8,10,12,13>, <3,6,9,10,11,12>, <3,6,10,11,12,13>, <4,5,6,8,9,12>,
             <4,6,7,10,11,13>, <4,6,7,11,12,13>, <4,6,8,9,10,12>, <4,6,8,10,12,13>,
             <4,6,9,10,11,12>, <4,6,10,11,12,13>)
            sage: len(_)
            34
        """
        lex_supp = set()
        for circuit in self.circuits():
            Cplus = circuit[0]
            Cminus = circuit[2]
            s0 = min(Cplus + Cminus)
            if s0 in Cplus:
                lex_supp.add(Cplus)
            else:
                lex_supp.add(Cminus)

        lex_supp = sorted(lex_supp, key=lambda x:-len(x))
        basepts = copy(lex_supp)
        for i in range(0,len(lex_supp)-1):
            for j in range(i+1,len(lex_supp)):
                if set(lex_supp[j]).issubset(set(lex_supp[i])):
                    try:
                        basepts.remove(lex_supp[i])
                    except ValueError:
                        pass

        basepts = [ (len(b),)+b for b in basepts ]     # decorate
        basepts = sorted(basepts)                      # sort
        basepts = [ b[1:] for b in basepts ]           # undecorate

        def make_cotriang(basepts):
            if len(basepts)==0:
                return [frozenset()]
            triangulation = set()
            for tail in make_cotriang(basepts[1:]):
                for head in basepts[0]:
                    triangulation.update([ frozenset([head]).union(tail) ])

            nonminimal = set()
            for rel in Combinations(triangulation,2):
                if rel[0].issubset(rel[1]): nonminimal.update([rel[1]])
                if rel[1].issubset(rel[0]): nonminimal.update([rel[0]])
            triangulation.difference_update(nonminimal)

            triangulation = [ [len(t)]+sorted(t) for t in triangulation ] # decorate
            triangulation = sorted(triangulation)                         # sort
            triangulation = [ frozenset(t[1:]) for t in triangulation ]   # undecorate

            return triangulation

        triangulation = make_cotriang(basepts)
        I = frozenset(range(0,self.n_points()))
        triangulation = [ tuple(I.difference(t)) for t in triangulation ]

        return self(triangulation)


