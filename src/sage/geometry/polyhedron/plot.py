"""
Functions for plotting polyhedra
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.rings.all import QQ, ZZ, RDF
from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix, identity_matrix
from sage.misc.functional import norm
from sage.structure.sequence import Sequence

from sage.plot.all import point2d, line2d, arrow, polygon2d
from sage.plot.plot3d.all import point3d, line3d, arrow3d, polygon3d
from sage.graphs.graph import Graph

from base import is_Polyhedron
from constructor import Polyhedron



#############################################################
def render_2d(projection, **kwds):
    """
    Return 2d rendering of the projection of a polyhedron into
    2-dimensional ambient space.

    EXAMPLES::

        sage: p1 = Polyhedron(vertices=[[1,1]], rays=[[1,1]])
        sage: q1 = p1.projection()
        sage: p2 = Polyhedron(vertices=[[1,0], [0,1], [0,0]])
        sage: q2 = p2.projection()
        sage: p3 = Polyhedron(vertices=[[1,2]])
        sage: q3 = p3.projection()
        sage: p4 = Polyhedron(vertices=[[2,0]], rays=[[1,-1]], lines=[[1,1]])
        sage: q4 = p4.projection()
        sage: q1.show() + q2.show() + q3.show() + q4.show()
        sage: from sage.geometry.polyhedron.plot import render_2d
        sage: q = render_2d(p1.projection())
        sage: q._Graphics__objects
        [Point set defined by 1 point(s),
         Arrow from (1.0,1.0) to (2.0,2.0),
         Polygon defined by 3 points]
    """
    if is_Polyhedron(projection):
        projection = Projection(projection)
    return \
        projection.render_points_2d(zorder=2, pointsize=10, **kwds) + \
        projection.render_outline_2d(zorder=1, **kwds) + \
        projection.render_fill_2d(zorder=0, rgbcolor=(0,1,0), **kwds)


def render_3d(projection, **kwds):
    """
    Return 3d rendering of a polyhedron projected into
    3-dimensional ambient space.

    NOTE:

    This method, ``render_3d``, is used in the ``show()``
    method of a polyhedron if it is in 3 dimensions.

    EXAMPLES::

        sage: p1 = Polyhedron(vertices=[[1,1,1]], rays=[[1,1,1]])
        sage: p2 = Polyhedron(vertices=[[2,0,0], [0,2,0], [0,0,2]])
        sage: p3 = Polyhedron(vertices=[[1,0,0], [0,1,0], [0,0,1]], rays=[[-1,-1,-1]])
        sage: p1.projection().show() + p2.projection().show() + p3.projection().show()

    It correctly handles various degenerate cases::

        sage: Polyhedron(lines=[[1,0,0],[0,1,0],[0,0,1]]).show()                              # whole space
        sage: Polyhedron(vertices=[[1,1,1]], rays=[[1,0,0]], lines=[[0,1,0],[0,0,1]]).show()  # half space
        sage: Polyhedron(vertices=[[1,1,1]], lines=[[0,1,0],[0,0,1]]).show()                  # R^2 in R^3
        sage: Polyhedron(rays=[[0,1,0],[0,0,1]], lines=[[1,0,0]]).show()                      # quadrant wedge in R^2
        sage: Polyhedron(rays=[[0,1,0]], lines=[[1,0,0]]).show()                              # upper half plane in R^3
        sage: Polyhedron(lines=[[1,0,0]]).show()                                              # R^1 in R^2
        sage: Polyhedron(rays=[[0,1,0]]).show()                                               # Half-line in R^3
        sage: Polyhedron(vertices=[[1,1,1]]).show()                                           # point in R^3
    """
    if is_Polyhedron(projection):
        projection = Projection(projection)
    return \
        projection.render_vertices_3d(width=3, color='green', **kwds) +\
        projection.render_wireframe_3d(width=3, color='green', **kwds) + \
        projection.render_solid_3d(**kwds)


def render_4d(polyhedron, **kwds):
    """
    Return a 3d rendering of the Schlegel projection of a 4d
    polyhedron projected into 3-dimensional space.

    NOTES:

    The ``show()`` method of ``Polyhedron()`` uses this to draw itself
    if the ambient dimension is 4.

    INPUT:

    - ``polyhedron`` -- A
      :mod:`~sage.geometry.polyhedron.constructor.Polyhedron` object.

    - ``kwds`` -- plot keywords. Passing
      ``projection_direction=<list>`` sets the projetion direction of
      the Schlegel projection. If it is not given, the center of a
      facet is used.

    EXAMPLES::

        sage: poly = polytopes.twenty_four_cell()
        sage: poly
        A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 24 vertices
        sage: poly.show()
        sage: poly.show(projection_direction=[2,5,11,17])
        sage: type( poly.show() )
        <class 'sage.plot.plot3d.base.Graphics3dGroup'>

    TESTS::

        sage: from sage.geometry.polyhedron.plot import render_4d
        sage: p = polytopes.n_cube(4)
        sage: q = render_4d(p)
        sage: tach_str = q.tachyon()
        sage: tach_str.count('FCylinder')
        32
    """
    projection_direction = None
    try:
        projection_direction = kwds.pop('projection_direction')
    except KeyError:
        for ineq in polyhedron.inequality_generator():
            center = [v() for v in ineq.incident() if v.is_vertex()]
            center = sum(center) / len(center)
            if not center.is_zero():
                projection_direction = center
                break
    projection_3d = Projection(polyhedron).schlegel(projection_direction)
    return render_3d(projection_3d, **kwds)



#############################################################
def cyclic_sort_vertices_2d(Vlist):
    """
    Return the vertices/rays in cyclic order if possible.

    NOTES:

    This works if and only if each vertex/ray is adjacent to exactly
    two others. For example, any 2-dimensional polyhedron satisfies
    this.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
        sage: square = Polyhedron([[1,0],[-1,0],[0,1],[0,-1]])
        sage: vertices = [v for v in square.vertex_generator()]
        sage: vertices
        [A vertex at (-1, 0),
         A vertex at (0, -1),
         A vertex at (0, 1),
         A vertex at (1, 0)]
        sage: cyclic_sort_vertices_2d(vertices)
        [A vertex at (1, 0),
         A vertex at (0, -1),
         A vertex at (-1, 0),
         A vertex at (0, 1)]
    """
    if len(Vlist)==0: return Vlist

    adjacency_matrix = Vlist[0].polyhedron().vertex_adjacency_matrix()
    result = [ Vlist.pop() ]
    while len(Vlist)>0:
        for i in range(len(Vlist)):
            if adjacency_matrix[Vlist[i].index(), result[-1].index()] == 1:
                result.append( Vlist.pop(i) )
                break;
        else:
            raise ValueError
    return result




#########################################################################
def projection_func_identity(x):
    """
    The identity projection.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import projection_func_identity
        sage: projection_func_identity((1,2,3))
        [1, 2, 3]
    """
    return list(x)



class ProjectionFuncStereographic():
    """
    The stereographic (or perspective) projection.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
        sage: cube = polytopes.n_cube(3).vertices()
        sage: proj = ProjectionFuncStereographic([1.2, 3.4, 5.6])
        sage: ppoints = [proj(vector(x)) for x in cube]
        sage: ppoints[0]
        (-0.0486511..., 0.0859565...)
    """
    def __init__(self, projection_point):
        """
        Create a stereographic projection function.

        INPUT:

        - ``projection_point`` -- a list of coordinates in the
          appropriate dimension, which is the point projected from.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
            sage: proj = ProjectionFuncStereographic([1.0,1.0])
            sage: proj.__init__([1.0,1.0])
            sage: proj.house
            [-0.7071067811...  0.7071067811...]
            [ 0.7071067811...  0.7071067811...]
            sage: TestSuite(proj).run(skip='_test_pickling')
        """
        self.projection_point = vector(projection_point)
        self.dim = self.projection_point.degree()

        pproj = vector(RDF,self.projection_point)
        self.psize = norm(pproj)
        if (self.psize).is_zero():
            raise ValueError, "projection direction must be a non-zero vector."
        v = vector(RDF, [0.0]*(self.dim-1) + [self.psize]) - pproj
        polediff = matrix(RDF,v).transpose()
        denom = RDF((polediff.transpose()*polediff)[0][0])
        if denom.is_zero():
            self.house = identity_matrix(RDF,self.dim)
        else:
            self.house = identity_matrix(RDF,self.dim) \
            - 2*polediff*polediff.transpose()/denom   # Householder reflector

    def __call__(self, x):
        """
        Action of the stereographic projection.

        INPUT:

        - ``x`` -- a vector or anything convertible to a vector.

        OUTPUT:

        First reflects ``x`` with a Householder reflection which takes
        the projection point to ``(0,...,0,self.psize)`` where
        ``psize`` is the length of the projection point, and then
        dilates by ``1/(zdiff)`` where ``zdiff`` is the difference
        between the last coordinate of ``x`` and ``psize``.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
            sage: proj = ProjectionFuncStereographic([1.0,1.0])
            sage: proj.__call__(vector([1,2]))
            (-1.0)
            sage: proj = ProjectionFuncStereographic([2.0,1.0])
            sage: proj.__call__(vector([1,2]))
            (3.0)
            sage: proj = ProjectionFuncStereographic([0,0,2])
            sage: proj.__call__(vector([0,0,1]))
            (0.0, 0.0)
            sage: proj.__call__(vector([1,0,0]))
            (0.5, 0.0)
        """
        img = self.house * x
        denom = self.psize-img[self.dim-1]
        if denom.is_zero():
            raise ValueError, 'Point cannot coincide with ' \
                'coordinate singularity at ' + repr(x)
        return vector(RDF, [img[i]/denom for i in range(self.dim-1)])


class ProjectionFuncSchlegel():
    """
    The Schlegel projection from the given input point.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import ProjectionFuncSchlegel
        sage: proj = ProjectionFuncSchlegel([2,2,2])
        sage: proj(vector([1.1,1.1,1.11]))[0]
        0.0302...
    """
    def __init__(self, projection_direction, height = 1.1):
        """
        Initializes the projection.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncSchlegel
            sage: proj = ProjectionFuncSchlegel([2,2,2])
            sage: proj.__init__([2,2,2])
            sage: proj(vector([1.1,1.1,1.11]))[0]
            0.0302...
            sage: TestSuite(proj).run(skip='_test_pickling')
        """
        self.projection_dir = vector(RDF, projection_direction)
        if norm(self.projection_dir).is_zero():
            raise ValueError, "projection direction must be a non-zero vector."
        self.dim = self.projection_dir.degree()
        spcenter = height * self.projection_dir/norm(self.projection_dir)
        self.height = height
        v = vector(RDF, [0.0]*(self.dim-1) + [self.height]) - spcenter
        polediff = matrix(RDF,v).transpose()
        denom = (polediff.transpose()*polediff)[0][0]
        if denom.is_zero():
            self.house = identity_matrix(RDF,self.dim)
        else:
            self.house = identity_matrix(RDF,self.dim) \
            - 2*polediff*polediff.transpose()/denom #Householder reflector

    def __call__(self, x):
        """
        Apply the projection to a vector.

        - ``x`` -- a vector or anything convertible to a vector.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncSchlegel
            sage: proj = ProjectionFuncSchlegel([2,2,2])
            sage: proj.__call__([1,2,3])
            (0.56162854..., 2.09602626...)
        """
        v = vector(RDF,x)
        if v.is_zero():
            raise ValueError, "The origin must not be a vertex."
        v = v/norm(v)         # normalize vertices to unit sphere
        v = self.house*v      # reflect so self.projection_dir is at "north pole"
        denom = self.height-v[self.dim-1]
        if denom.is_zero():
            raise ValueError, 'Point cannot coincide with ' \
                'coordinate singularity at ' + repr(x)
        return vector(RDF, [ v[i]/denom for i in range(self.dim-1) ])



#########################################################################
class Projection(SageObject):
    """
    The projection of a :class:`Polyhedron`.

    This class keeps track of the necessary data to plot the input
    polyhedron.
    """

    def __init__(self, polyhedron, proj=projection_func_identity):
        """
        Initialize the projection of a Polyhedron() object.

        INPUT:

          - ``polyhedron`` -- a ``Polyhedron()`` object

          - ``proj`` -- a projection function for the points

        NOTES:

        Once initialized, the polyhedral data is fixed. However, the
        projection can be changed later on.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: Projection(p)
            The projection of a polyhedron into 3 dimensions
            sage: def pr_12(x): return [x[1],x[2]]
            sage: Projection(p, pr_12)
            The projection of a polyhedron into 2 dimensions
            sage: Projection(p,  lambda x: [x[1],x[2]] )   # another way of doing the same projection
            The projection of a polyhedron into 2 dimensions
            sage: _.show()   # plot of the projected icosahedron in 2d
            sage: proj = Projection(p)
            sage: proj.stereographic([1,2,3])
            The projection of a polyhedron into 2 dimensions
            sage: proj.show()
            sage: TestSuite(proj).run(skip='_test_pickling')
        """
        self.coords = Sequence([])
        self.points = Sequence([])
        self.lines  = Sequence([])
        self.arrows = Sequence([])
        self.polygons = Sequence([])
        self.polyhedron_ambient_dim = polyhedron.ambient_dim()

        if polyhedron.ambient_dim() == 2:
            self._init_from_2d(polyhedron)
        elif polyhedron.ambient_dim() == 3:
            self._init_from_3d(polyhedron)
        else:
            self._init_points(polyhedron)
            self._init_lines_arrows(polyhedron)

        self.coords.set_immutable()
        self.points.set_immutable()
        self.lines.set_immutable()
        self.arrows.set_immutable()
        self.polygons.set_immutable()

        self.__call__(proj)


    def _repr_(self):
        """
        Return a string describing the projection.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: proj = Projection(p)
            sage: print proj._repr_()
            The projection of a polyhedron into 3 dimensions
        """
        s = 'The projection of a polyhedron into ' + \
            repr(self.dimension) + ' dimensions'
        return s


    def __call__(self, proj=projection_func_identity):
        """
        Apply a projection.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: pproj = Projection(p)
            sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
            sage: pproj_stereo = pproj.__call__(proj = ProjectionFuncStereographic([1,2,3]))
            sage: pproj_stereo.polygons[0]
            [10, 4, 6]
        """
        self.transformed_coords = \
            Sequence([proj(p) for p in self.coords])
        self._init_dimension()
        return self


    def identity(self):
        """
        Return the identity projection of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: pproj = Projection(p)
            sage: ppid = pproj.identity()
            sage: ppid.dimension
            3
        """
        return self.__call__(projection_func_identity)


    def stereographic(self, projection_point=None):
        r"""
        Rteurn the stereographic projection.

        INPUT:

        - ``projection_point`` - The projection point. This must be
          distinct from the polyhedron's vertices. Default is `(1,0,\dots,0)`

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import Projection
            sage: proj = Projection(polytopes.buckyball())  #long time
            sage: proj                                      #long time
            The projection of a polyhedron into 3 dimensions
            sage: proj.stereographic([5,2,3]).show()        #long time
            sage: Projection( polytopes.twenty_four_cell() ).stereographic([2,0,0,0])
            The projection of a polyhedron into 3 dimensions
        """
        if projection_point == None:
            projection_point = [1] + [0]*(self.polyhedron_ambient_dim-1)
        return self.__call__(ProjectionFuncStereographic(projection_point))


    def schlegel(self, projection_direction=None, height = 1.1):
        """
        Return the Schlegel projection.

        The vertices are normalized to the unit sphere, and
        stereographically projected from a point slightly outside of
        the sphere.

        INPUT:

          - ``projection_direction`` - The direction of the Schlegel
            projection. By default, the vector consisting of the first
            n primes is chosen.

        EXAMPLES::

            sage: cube4 = polytopes.n_cube(4)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: Projection(cube4).schlegel([1,0,0,0])
            The projection of a polyhedron into 3 dimensions
            sage: _.show()
        """
        if projection_direction == None:
            for poly in self.polygons:
                center = sum([self.coords[i] for i in poly]) / len(poly)
                print center, "\n"
                if not center.is_zero():
                    projection_direction = center
                    break
        if projection_direction == None:
            projection_direction = primes_first_n(self.polyhedron_ambient_dim)
        return self.__call__(ProjectionFuncSchlegel(projection_direction, height = height))


    def coord_index_of(self, v):
        """
        Convert a coordinate vector to its internal index.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: proj = p.projection()
            sage: proj.coord_index_of(vector((1,1,1)))
            7
        """
        try:
            return self.coords.index(v)
        except ValueError:
            self.coords.append(v)
            return len(self.coords)-1


    def coord_indices_of(self, v_list):
        """
        Convert list of coordinate vectors to the corresponding list
        of internal indices.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: proj = p.projection()
            sage: proj.coord_indices_of([vector((1,1,1)),vector((1,-1,1))])
            [7, 5]
        """
        return [self.coord_index_of(v) for v in v_list]


    def coordinates_of(self, coord_index_list):
        """
        Given a list of indices, return the projected coordinates.

        EXAMPLES::

            sage: p = polytopes.n_simplex(4).projection()
            sage: p.coordinates_of([1])
            [[0, -81649/100000, 7217/25000, 22361/100000]]
        """
        return [self.transformed_coords[i] for i in coord_index_list]


    def _init_dimension(self):
        """
        Internal function: Initialize from 2d polyhedron. Must always
        be called after a coordinate projection.

        TESTS::

            sage: from sage.geometry.polyhedron.plot import Projection, render_2d
            sage: p = polytopes.n_simplex(2).projection()
            sage: test = p._init_dimension()
            sage: p.show.__doc__ == render_2d.__doc__
            True
        """
        self.dimension = len(self.transformed_coords[0])

        if self.dimension == 2:
            self.show = lambda **kwds: render_2d(self,**kwds)
            self.show.__doc__ = render_2d.__doc__
        elif self.dimension == 3:
            self.show = lambda **kwds: render_3d(self,**kwds)
            self.show.__doc__ = render_3d.__doc__
        else:
            try:
                del self.show
            except AttributeError:
                pass


    def _init_from_2d(self, polyhedron):
        """
        Internal function: Initialize from 2d polyhedron.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0],[0,1],[1,0],[1,1]])
            sage: proj = p.projection()
            sage: [proj.coordinates_of([i]) for i in proj.points]
            [[[0, 0]], [[0, 1]], [[1, 0]], [[1, 1]]]
            sage: proj._init_from_2d
            <bound method Projection._init_from_2d of The projection
            of a polyhedron into 2 dimensions>
        """
        assert polyhedron.ambient_dim() == 2, "Requires polyhedron in 2d"
        self.dimension = 2
        self._init_points(polyhedron)
        self._init_lines_arrows(polyhedron)
        self._init_area_2d(polyhedron)


    def _init_from_3d(self, polyhedron):
        """
        Internal function: Initialize from 3d polyhedron.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,1],[0,1,2],[1,0,3],[1,1,5]])
            sage: proj = p.projection()
            sage: [proj.coordinates_of([i]) for i in proj.points]
            [[[0, 0, 1]], [[0, 1, 2]], [[1, 0, 3]], [[1, 1, 5]]]
            sage: proj._init_from_3d
            <bound method Projection._init_from_3d of The projection
            of a polyhedron into 3 dimensions>
        """
        assert polyhedron.ambient_dim() == 3, "Requires polyhedron in 3d"
        self.dimension = 3
        self._init_points(polyhedron)
        self._init_lines_arrows(polyhedron)
        self._init_solid_3d(polyhedron)


    def _init_points(self, polyhedron):
        """
        Internal function: Initialize points (works in arbitrary
        dimensions).

        TESTS::

            sage: p = polytopes.n_cube(2)
            sage: pp = p.projection()
            sage: del pp.points
            sage: pp.points = Sequence([])
            sage: pp._init_points(p)
            sage: pp.points
            [0, 1, 2, 3]
        """
        for v in polyhedron.vertex_generator():
            self.points.append( self.coord_index_of(v.vector()) )


    def _init_lines_arrows(self, polyhedron):
        """
        Internal function: Initialize compact and non-compact edges
        (works in arbitrary dimensions).

        TESTS::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: pp = p.projection()
            sage: pp.arrows
            [[0, 1], [0, 2], [0, 3], [0, 4]]
            sage: del pp.arrows
            sage: pp.arrows = Sequence([])
            sage: pp._init_lines_arrows(p)
            sage: pp.arrows
            [[0, 1], [0, 2], [0, 3], [0, 4]]
        """
        obj = polyhedron.Vrepresentation()
        for i in range(len(obj)):
            if not obj[i].is_vertex(): continue
            for j in range(len(obj)):
                if polyhedron.vertex_adjacency_matrix()[i,j] == 0: continue
                if i<j and obj[j].is_vertex():
                    l = [obj[i].vector(), obj[j].vector()]
                    self.lines.append( [ self.coord_index_of(l[0]),
                                         self.coord_index_of(l[1]) ] )
                if obj[j].is_ray():
                    l = [obj[i].vector(), obj[i].vector() + obj[j].vector()]
                    self.arrows.append( [ self.coord_index_of(l[0]),
                                          self.coord_index_of(l[1]) ] )
                if obj[j].is_line():
                    l1 = [obj[i].vector(), obj[i].vector() + obj[j].vector()]
                    l2 = [obj[i].vector(), obj[i].vector() - obj[j].vector()]
                    self.arrows.append( [ self.coord_index_of(l1[0]),
                                          self.coord_index_of(l1[1]) ] )
                    self.arrows.append( [ self.coord_index_of(l2[0]),
                                          self.coord_index_of(l2[1]) ] )


    def _init_area_2d(self, polyhedron):
        """
        Internal function: Initialize polygon area for 2d polyhedron.

        TESTS::

            sage: p = polytopes.cyclic_polytope(2,4)
            sage: proj = p.projection()
            sage: proj.polygons = Sequence([])
            sage: proj._init_area_2d(p)
            sage: proj.polygons
            [[3, 0, 1, 2]]
        """
        assert polyhedron.ambient_dim() == 2, "Requires polyhedron in 2d"
        vertices = [v for v in polyhedron.Vrep_generator()]
        vertices = cyclic_sort_vertices_2d(vertices)
        coords = []

        def adjacent_vertices(i):
            n = len(vertices)
            if vertices[(i-1) % n].is_vertex(): yield vertices[(i-1) % n]
            if vertices[(i+1) % n].is_vertex(): yield vertices[(i+1) % n]

        for i in range(len(vertices)):
            v = vertices[i]
            if v.is_vertex():
                coords.append(v())
            if v.is_ray():
                for a in adjacent_vertices(i):
                    coords.append(a() + v())

        if polyhedron.n_lines() == 0:
            self.polygons.append( self.coord_indices_of(coords) )
            return

        polygons = []

        if polyhedron.n_lines() == 1:
            aline = polyhedron.line_generator().next()
            for shift in [aline(), -aline()]:
                for i in range(len(coords)):
                    polygons.append( [ coords[i-1],coords[i],
                                       coords[i]+shift, coords[i-1]+shift ] )

        if polyhedron.n_lines() == 2:
            [line1, line2] = [l for l in polyhedron.lines()]
            assert len(coords)==1, "Can have only a single vertex!"
            v = coords[0]
            l1 = line1()
            l2 = line2()
            polygons = [ [v-l1-l2, v+l1-l2, v+l1+l2, v-l1+l2] ]

        polygons = [ self.coord_indices_of(p) for p in polygons ]
        self.polygons.extend(polygons)



    def _init_solid_3d(self, polyhedron):
        """
        Internal function: Initialize facet polygons for 3d polyhedron.

        TESTS::

            sage: p = polytopes.cyclic_polytope(3,4)
            sage: proj = p.projection()
            sage: proj.polygons = Sequence([])
            sage: proj._init_solid_3d(p)
            sage: proj.polygons
            [[2, 0, 1], [3, 0, 1], [3, 0, 2], [3, 1, 2]]
        """
        assert polyhedron.ambient_dim() == 3, "Requires polyhedron in 3d"

        if polyhedron.dim() <= 1: # empty or 0d or 1d polyhedron => no polygon
            return None

        def defining_equation():  # corresponding to a polygon
            if polyhedron.dim() < 3:
                yield polyhedron.equation_generator().next()
            else:
                for ineq in polyhedron.inequality_generator():
                    yield ineq

        faces = []
        for facet_equation in defining_equation():
            vertices = [v for v in facet_equation.incident()]
            vertices = cyclic_sort_vertices_2d(vertices)
            coords = []

            def adjacent_vertices(i):
                n = len(vertices)
                if vertices[(i-1) % n].is_vertex(): yield vertices[(i-1) % n]
                if vertices[(i+1) % n].is_vertex(): yield vertices[(i+1) % n]

            for i in range(len(vertices)):
                v = vertices[i]
                if v.is_vertex():
                    coords.append(v())
                if v.is_ray():
                    for a in adjacent_vertices(i):
                        coords.append(a() + v())

            faces.append(coords)

        if polyhedron.n_lines() == 0:
            assert len(faces)>0, "no vertices?"
            self.polygons.extend( [self.coord_indices_of(f) for f in faces] )
            return

        # now some special cases if there are lines (dim < ambient_dim)
        polygons = []

        if polyhedron.n_lines()==1:
            assert len(faces)>0, "no vertices?"
            aline = polyhedron.line_generator().next()
            for shift in [aline(), -aline()]:
                for coords in faces:
                    assert len(coords)==2, "There must be two points."
                    polygons.append( [ coords[0],coords[1],
                                       coords[1]+shift, coords[0]+shift ] )

        if polyhedron.n_lines()==2:
            [line1, line2] = [l for l in polyhedron.line_generator()]
            l1 = line1()
            l2 = line2()
            for v in polyhedron.vertex_generator():
                polygons.append( [v()-l1-l2, v()+l1-l2, v()+l1+l2, v()-l1+l2] )

        self.polygons.extend( [self.coord_indices_of(p) for p in polygons] )


    def render_points_2d(self, **kwds):
        """
        Return the points of a polyhedron in 2d.

        EXAMPLES::

            sage: hex = polytopes.regular_polygon(6)
            sage: proj = hex.projection()
            sage: hex_points = proj.render_points_2d()
            sage: hex_points._Graphics__objects
            [Point set defined by 6 point(s)]
        """
        return point2d(self.coordinates_of(self.points), **kwds)


    def render_outline_2d(self, **kwds):
        """
        Return the outline (edges) of a polyhedron in 2d.

        EXAMPLES::

            sage: penta = polytopes.regular_polygon(5)
            sage: outline = penta.projection().render_outline_2d()
            sage: outline._Graphics__objects[0]
            Line defined by 2 points
        """
        wireframe = [];
        for l in self.lines:
            l_coords = self.coordinates_of(l)
            wireframe.append( line2d(l_coords, **kwds) )
        for a in self.arrows:
            a_coords = self.coordinates_of(a)
            wireframe.append( arrow(a_coords[0], a_coords[1], **kwds) )
        return sum(wireframe)


    def render_fill_2d(self, **kwds):
        """
        Return the filled interior (a polygon) of a polyhedron in 2d.

        EXAMPLES::

            sage: cps = [i^3 for i in srange(-2,2,1/5)]
            sage: p = Polyhedron(vertices = [[(t^2-1)/(t^2+1),2*t/(t^2+1)] for t in cps])
            sage: proj = p.projection()
            sage: filled_poly = proj.render_fill_2d()
            sage: filled_poly.axes_width()
            0.8
        """
        poly = [polygon2d(self.coordinates_of(p), **kwds)
                 for p in self.polygons]
        return sum(poly)


    def render_vertices_3d(self, **kwds):
        """
        Return the 3d rendering of the vertices.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(3)
            sage: proj = p.projection()
            sage: verts = proj.render_vertices_3d()
            sage: verts.bounding_box()
            ((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
        """
        return point3d(self.coordinates_of(self.points), **kwds)


    def render_wireframe_3d(self, **kwds):
        r"""
        Return the 3d wireframe rendering.

        EXAMPLES::

            sage: cube = polytopes.n_cube(3)
            sage: cube_proj = cube.projection()
            sage: wire = cube_proj.render_wireframe_3d()
            sage: print wire.tachyon().split('\n')[77]  # for testing
            FCylinder base -1.0 1.0 -1.0 apex -1.0 -1.0 -1.0 rad 0.005 texture...
        """
        wireframe = [];
        for l in self.lines:
            l_coords = self.coordinates_of(l)
            wireframe.append( line3d(l_coords, **kwds))
        for a in self.arrows:
            a_coords = self.coordinates_of(a)
            wireframe.append(arrow3d(a_coords[0], a_coords[1], **kwds))
        return sum(wireframe)


    def render_solid_3d(self, **kwds):
        """
        Return solid 3d rendering of a 3d polytope.

        EXAMPLES::

            sage: p = polytopes.n_cube(3).projection()
            sage: p_solid = p.render_solid_3d(opacity = .7)
            sage: type(p_solid)
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
        """
        return sum([ polygon3d(self.coordinates_of(f), **kwds)
                     for f in self.polygons ])
