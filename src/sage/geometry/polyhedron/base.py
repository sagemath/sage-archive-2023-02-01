r"""
Base class for polyhedra
"""


########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.structure.sage_object import SageObject

from sage.misc.all import union, cached_method, prod
from sage.misc.package import is_package_installed

from sage.rings.all import Integer, QQ, ZZ, primes_first_n
from sage.rings.rational import Rational
from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix, identity_matrix
from sage.functions.other import sqrt, floor, ceil

from sage.plot.all import point2d, line2d, arrow, polygon2d
from sage.plot.plot3d.all import point3d, line3d, arrow3d, polygon3d
from sage.graphs.graph import Graph

from sage.combinat.cartesian_product import CartesianProduct
from sage.groups.perm_gps.permgroup_named import AlternatingGroup

from constructor import Polyhedron
from representation import (
    PolyhedronRepresentation,
    Hrepresentation,
    Inequality, Equation,
    Vrepresentation,
    Vertex, Ray, Line )


#########################################################################
# Notes if you want to implement your own backend:
#
#  * derive from Polyhedron_base
#
#  * you must implement _init_from_Vrepresentation and
#    _init_from_Vrepresentationa
#
#  * You might want to override _init_empty_polyhedron,
#    _init_facet_adjacency_matrix, _init_vertex_adjacency_matrix, and
#    _make_polyhedron_face.
#
#  * You can of course also override any other method for which you
#    have a faster implementation.


#########################################################################
def is_Polyhedron(X):
    """
    Test whether ``X`` is a Polyhedron.

    INPUT:

    - ``X`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: p = polytopes.n_cube(2)
        sage: from sage.geometry.polyhedron.base import is_Polyhedron
        sage: is_Polyhedron(p)
        True
        sage: is_Polyhedron(123456)
        False
    """
    return isinstance(X, Polyhedron_base)


#########################################################################
class Polyhedron_base(SageObject):
    """
    Base class for Polyhedron objects

    INPUT:

    - ``ambient_dim`` -- integer. The dimension of the ambient space.

    - ``Vrep`` -- a list `[vertices, rays, lines]`` or ``None``. The
      V-representation of the polyhedron. If ``None``, the polyhedron
      is determined by the H-representation.

    - ``Hrep`` -- a list `[ieqs, eqns]`` or ``None``. The
      H-representation of the polyhedron. If ``None``, the polyhedron
      is determined by the V-representation.

    Only one of ``Vrep`` or ``Hrep`` can be different from ``None``.

    TESTS::

        sage: p = Polyhedron()
        sage: TestSuite(p).run()
    """

    def __init__(self, ambient_dim, Vrep, Hrep, **kwds):
        """
        Initializes the polyhedron.

        See :class:`Polyhedron_base` for a description of the input
        data.

        TESTS::

            sage: p = Polyhedron()    # indirect doctests
        """
        self._ambient_dim = ambient_dim
        if Vrep is not None:
            vertices, rays, lines = Vrep
            if len(vertices)==0:
                vertices = [[0] * ambient_dim]
            self._init_from_Vrepresentation(ambient_dim, vertices, rays, lines, **kwds)
        elif Hrep is not None:
            ieqs, eqns = Hrep
            self._init_from_Hrepresentation(ambient_dim, ieqs, eqns, **kwds)
        else:
            self._init_empty_polyhedron(ambient_dim)


    def _init_from_Vrepresentation(self, ambient_dim, vertices, rays, lines, **kwds):
        """
        Construct polyhedron from V-representation data.

        INPUT:

        - ``ambient_dim`` -- integer. The dimension of the ambient space.

        - ``vertices`` -- list of point. Each point can be specified
           as any iterable container of
           :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``rays`` -- list of rays. Each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``lines`` -- list of lines. Each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        EXAMPLES::

            sage: p = Polyhedron()
            sage: from sage.geometry.polyhedron.base import Polyhedron_base
            sage: Polyhedron_base._init_from_Vrepresentation(p, 2, [], [], [])
            Traceback (most recent call last):
            ...
            NotImplementedError: A derived class must implement this method.
        """
        raise NotImplementedError('A derived class must implement this method.')


    def _init_from_Hrepresentation(self, ambient_dim, ieqs, eqns, **kwds):
        """
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ambient_dim`` -- integer. The dimension of the ambient space.

        - ``ieqs`` -- list of inequalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``eqns`` -- list of equalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        EXAMPLES::

            sage: p = Polyhedron()
            sage: from sage.geometry.polyhedron.base import Polyhedron_base
            sage: Polyhedron_base._init_from_Hrepresentation(p, 2, [], [])
            Traceback (most recent call last):
            ...
            NotImplementedError: A derived class must implement this method.
        """
        raise NotImplementedError('A derived class must implement this method.')


    def _init_empty_polyhedron(self, ambient_dim):
        """
        Initializes an empty polyhedron.

        INPUT:

        - ``ambient_dim`` -- integer. The dimension of the ambient space.

        TESTS::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in QQ^0
            sage: empty.Vrepresentation()
            ()
            sage: empty.Hrepresentation()
            (An equation -1 == 0,)
            sage: Polyhedron(vertices = [])
            The empty polyhedron in QQ^0
            sage: Polyhedron()._init_empty_polyhedron(0)
        """
        self._Vrepresentation = []
        self._Hrepresentation = []
        Equation(self, [-1] + [0]*ambient_dim);
        self._Vrepresentation = tuple(self._Vrepresentation)
        self._Hrepresentation = tuple(self._Hrepresentation)

        self._V_adjacency_matrix = matrix(ZZ, 0, 0, 0)
        self._V_adjacency_matrix.set_immutable()

        self._H_adjacency_matrix = matrix(ZZ, 1, 1, 0)
        self._H_adjacency_matrix.set_immutable()


    def _init_facet_adjacency_matrix(self):
        """
        Compute the facet adjacency matrix in case it has not been
        computed during initialization.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)])
            sage: '_H_adjacency_matrix' in p.__dict__
            False
            sage: p._init_facet_adjacency_matrix()
            sage: p._H_adjacency_matrix
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        # TODO: This implementation computes the whole face lattice,
        # which is much more information than necessary.
        M = matrix(ZZ, self.n_Hrepresentation(), self.n_Hrepresentation(), 0)
        def set_adjacent(h1,h2):
            if h1 is h2:
                return
            i = h1.index()
            j = h2.index()
            M[i,j]=1
            M[j,i]=1

        face_lattice = self.face_lattice()
        for face in face_lattice:
            Hrep = face.element.ambient_Hrepresentation()
            if len(Hrep) == 2:
                set_adjacent(Hrep[0], Hrep[1])

        self._H_adjacency_matrix = M


    def _init_vertex_adjacency_matrix(self):
        """
        Compute the vertex adjacency matrix in case it has not been
        computed during initialization.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)])
            sage: '_V_adjacency_matrix' in p.__dict__
            False
            sage: p._init_vertex_adjacency_matrix()
            sage: p._V_adjacency_matrix
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        # TODO: This implementation computes the whole face lattice,
        # which is much more information than necessary.
        M = matrix(ZZ, self.n_Vrepresentation(), self.n_Vrepresentation(), 0)
        def set_adjacent(v1,v2):
            if v1 is v2:
                return
            i = v1.index()
            j = v2.index()
            M[i,j]=1
            M[j,i]=1

        face_lattice = self.face_lattice()
        for face in face_lattice:
            Vrep = face.element.ambient_Vrepresentation()
            if len(Vrep) == 2:
                set_adjacent(Vrep[0], Vrep[1])

        for l in self.line_generator():
            for vrep in self.Vrep_generator():
                set_adjacent(l, vrep)
        for r in self.ray_generator():
            for vrep in self.Vrep_generator():
                set_adjacent(r, vrep)

        self._V_adjacency_matrix = M


    def __lt__(self, other):
        """
        Test whether ``self`` is a strict sub-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P < Q   # indirect doctest
            False
            sage: P < P   # indirect doctest
            False
            sage: Q < P   # indirect doctest
            True
        """
        return self._is_subpolyhedron(other) and not other._is_subpolyhedron(self)


    def __le__(self, other):
        """
        Test whether ``self`` is a (not necessarily strict)
        sub-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P <= Q   # indirect doctest
            False
            sage: P <= P   # indirect doctest
            True
            sage: Q <= P   # indirect doctest
            True
        """
        return self._is_subpolyhedron(other)


    def __eq__(self, other):
        """
        Test whether ``self`` is a strict sub-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P == Q   # indirect doctest
            False
            sage: P == P   # indirect doctest
            True
            sage: Q == P   # indirect doctest
            False
        """
        return self._is_subpolyhedron(other) and other._is_subpolyhedron(self)


    def __ne__(self, other):
        """
        Test whether ``self`` is not equal to ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P != Q   # indirect doctest
            True
            sage: P != P   # indirect doctest
            False
            sage: Q != P   # indirect doctest
            True
        """
        return not self.__eq__(other)


    def __gt__(self, other):
        """
        Test whether ``self`` is a strict super-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P > Q   # indirect doctest
            True
            sage: P > P   # indirect doctest
            False
            sage: Q > P   # indirect doctest
            False
        """
        return other._is_subpolyhedron(self) and not self._is_subpolyhedron(other)


    def __ge__(self, other):
        """
        Test whether ``self`` is a (not necessarily strict)
        super-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P >= Q   # indirect doctest
            True
            sage: P >= P   # indirect doctest
            True
            sage: Q >= P   # indirect doctest
            False
        """
        return other._is_subpolyhedron(self)


    def _is_subpolyhedron(self, other):
        """
        Test whether ``self`` is a (not necessarily strict)
        sub-polyhdedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P._is_subpolyhedron(Q)
            False
            sage: Q._is_subpolyhedron(P)
            True
        """
        if not is_Polyhedron(other):
            raise ValueError('Can only compare Polyhedron objects.')
        return all( other_H.contains(self_V)
                    for other_H, self_V in
                    CartesianProduct(other.Hrep_generator(), self.Vrep_generator()) )


    def plot(self, **kwds):
        """
        Return a graphical representation.

        INPUT:

        - ``**kwds`` -- optional keyword parameters.

        See :func:`render_2d`, :func:`render_3d`, :func:`render_4d`
        for a description of available options for different ambient
        space dimensions.

        OUTPUT:

        A graphics object.

        TESTS::

            sage: polytopes.n_cube(2).plot()
        """
        from plot import render_2d, render_3d, render_4d
        render_method = [ None, None, render_2d, render_3d, render_4d ]
        if self.ambient_dim() < len(render_method):
            render = render_method[self.ambient_dim()]
            if render != None:
                return render(self,**kwds)
        raise NotImplementedError('Plotting of '+str(self.ambient_dim())+
                                  '-dimensional polyhedra not implemented')


    show = plot


    def _repr_(self):
        """
        Return a description of the polyhedron.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1]])
            sage: poly_test._repr_()
            'A 2-dimensional polyhedron in QQ^4 defined as the convex hull of 3 vertices'
            sage: grammar_test = Polyhedron(vertices = [[1,1,1,1,1,1]])
            sage: grammar_test._repr_()
            'A 0-dimensional polyhedron in QQ^6 defined as the convex hull of 1 vertex'
        """
        desc = ''
        if self.n_vertices()==0:
            desc += 'The empty polyhedron'
        else:
            desc += 'A ' + repr(self.dim()) + '-dimensional polyhedron'
        desc += ' in '
        if self.field()==QQ: desc += 'QQ'
        else:                desc += 'RDF'
        desc += '^' + repr(self.ambient_dim())

        if self.n_vertices()>0:
            desc += ' defined as the convex hull of '
            desc += repr(self.n_vertices())
            if self.n_vertices()==1: desc += ' vertex'
            else:                    desc += ' vertices'

            if self.n_rays()>0:
                if self.n_lines()>0: desc += ", "
                else:                desc += " and "
                desc += repr(self.n_rays())
                if self.n_rays()==1: desc += ' ray'
                else:                desc += ' rays'

            if self.n_lines()>0:
                if self.n_rays()>0: desc += ", "
                else:               desc += " and "
                desc += repr(self.n_lines())
                if self.n_lines()==1: desc +=' line'
                else:                 desc +=' lines'

        return desc


    def cdd_Hrepresentation(self):
        """
        Write the inequalities/equations data of the polyhedron in
        cdd's H-representation format.

        OUTPUT:

        A string. If you save the output to filename.ine then you can
        run the stand-alone cdd via ``cddr+ filename.ine``

        EXAMPLES::

            sage: p = polytopes.n_cube(2)
            sage: print p.cdd_Hrepresentation()
            H-representation
            begin
             4 3 rational
             1 1 0
             1 0 1
             1 -1 0
             1 0 -1
            end
        """
        from cdd_file_format import cdd_Hrepresentation
        try:
            cdd_type = self._cdd_type
        except AttributeError:
            ring_to_cdd = { QQ:'rational', RDF:'real' }
            cdd_type = ring_to_cdd[self.base_ring()]
        return cdd_Hrepresentation(cdd_type,
                                   list(self.inequality_generator()),
                                   list(self.equation_generator()) )


    def cdd_Vrepresentation(self):
        """
        Write the vertices/rays/lines data of the polyhedron in cdd's
        V-representation format.

        OUTPUT:

        A string. If you save the output to filename.ext then you can
        run the stand-alone cdd via ``cddr+ filename.ext``

        EXAMPLES::

            sage: q = Polyhedron(vertices = [[1,1],[0,0],[1,0],[0,1]])
            sage: print q.cdd_Vrepresentation()
            V-representation
            begin
             4 3 rational
             1 0 0
             1 0 1
             1 1 0
             1 1 1
            end
        """
        from cdd_file_format import cdd_Vrepresentation
        try:
            cdd_type = self._cdd_type
        except AttributeError:
            ring_to_cdd = { QQ:'rational', RDF:'real' }
            cdd_type = ring_to_cdd[self.base_ring()]
        return cdd_Vrepresentation(cdd_type,
                                   list(self.vertex_generator()),
                                   list(self.ray_generator()),
                                   list(self.line_generator()) )


    def n_equations(self):
        """
        Return the number of equations. The representation will
        always be minimal, so the number of equations is the
        codimension of the polyhedron in the ambient space.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0,0],[0,1,0],[0,0,1]])
            sage: p.n_equations()
            1
        """
        try:
            return self._n_equations
        except AttributeError:
            self._n_equations = len(self.equations())
            return self._n_equations


    def n_inequalities(self):
        """
        Return the number of inequalities. The representation will
        always be minimal, so the number of inequalities is the
        number of facets of the polyhedron in the ambient space.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0,0],[0,1,0],[0,0,1]])
            sage: p.n_inequalities()
            3
        """
        try:
            return self._n_inequalities
        except AttributeError:
            self._n_inequalities = 0
            for i in self.inequalities(): self._n_inequalities += 1
            return self._n_inequalities


    def n_facets(self):
        """
        Return the number of facets in the polyhedron.  This is the
        same as the n_inequalities function.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)])
            sage: p.n_facets()
            8
        """
        return self.n_inequalities()


    def n_vertices(self):
        """
        Return the number of vertices. The representation will
        always be minimal.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0],[0,1],[1,1]], rays=[[1,1]])
            sage: p.n_vertices()
            2
        """
        try:
            return self._n_vertices
        except AttributeError:
            self._n_vertices = 0
            for v in self.vertex_generator(): self._n_vertices += 1
            return self._n_vertices


    def n_rays(self):
        """
        Return the number of rays. The representation will
        always be minimal.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0],[0,1]], rays=[[1,1]])
            sage: p.n_rays()
            1
        """
        try:
            return self._n_rays
        except AttributeError:
            self._n_rays = 0
            for r in self.rays(): self._n_rays += 1
            return self._n_rays


    def n_lines(self):
        """
        Return the number of lines. The representation will
        always be minimal.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0]], rays=[[0,1],[0,-1]])
            sage: p.n_lines()
            1
        """
        try:
            return self._n_lines
        except AttributeError:
            self._n_lines = len(self.lines())
            return self._n_lines


    def Hrepresentation(self, index=None):
        """
        Return the objects of the H-representaton. Each entry is
        either an inequality or a equation.

        INPUT:

        - ``index`` -- either an integer or ``None``.

        OUTPUT:

        The optional argument is an index in
        `0...n_Hrepresentations()`. If present, the H-representation
        object at the given index will be returned. Without an
        argument, returns the list of all H-representation objects.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p.Hrepresentation(0)
            An inequality (0, 0, -1) x + 1 >= 0
            sage: p.Hrepresentation(0) == p.Hrepresentation() [0]
            True
        """
        if index==None:
            return self._Hrepresentation
        else:
            return self._Hrepresentation[index]


    def Hrep_generator(self):
        """
        Return an iterator over the objects of the H-representation
        (inequalities or equations).

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p.Hrep_generator().next()
            An inequality (0, 0, -1) x + 1 >= 0
        """
        for H in self.Hrepresentation():
            yield H


    def n_Hrepresentation(self):
        """
        Return the number of objects that make up the
        H-representation of the polyhedron.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: p.n_Hrepresentation()
            16
            sage: p.n_Hrepresentation() == p.n_inequalities() + p.n_equations()
            True
        """
        return len(self.Hrepresentation())


    def Vrepresentation(self, index=None):
        """
        Return the objects of the V-representation. Each entry is
        either a vertex, a ray, or a line.

        INPUT:

        - ``index`` -- either an integer or ``None``.

        OUTPUT:

        The optional argument is an index in
        `0...n_Vrepresentation()`. If present, the V-representation
        object at the given index will be returned. Without an
        argument, returns the list of all V-representation objects.

        EXAMPLES::

            sage: p = polytopes.n_simplex(4)
            sage: p.Vrepresentation(0)
            A vertex at (-7071/10000, 1633/4000, 7217/25000, 22361/100000)
            sage: p.Vrepresentation(0) == p.Vrepresentation() [0]
            True
        """
        if index==None:
            return self._Vrepresentation
        else:
            return self._Vrepresentation[index]


    def n_Vrepresentation(self):
        """
        Return the number of objects that make up the
        V-representation of the polyhedron.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = polytopes.n_simplex(4)
            sage: p.n_Vrepresentation()
            5
            sage: p.n_Vrepresentation() == p.n_vertices() + p.n_rays() + p.n_lines()
            True
        """
        return len(self.Vrepresentation())


    def Vrep_generator(self):
        """
        Returns an iterator over the objects of the V-representation
        (vertices, rays, and lines).

        EXAMPLES::

            sage: p = polytopes.cyclic_polytope(3,4)
            sage: vg = p.Vrep_generator()
            sage: vg.next()
            A vertex at (0, 0, 0)
            sage: vg.next()
            A vertex at (1, 1, 1)
        """
        for V in self.Vrepresentation():
            yield V


    def facial_adjacencies(self):
        """
        Return the list of face indices (i.e. indices of
        H-representation objects) and the indices of faces adjacent to
        them.

        .. NOTE::

            Instead of working with face indices, it is recommended
            that you use the H-representation objects directly (see
            example).

        EXAMPLES::

            sage: p = polytopes.permutahedron(4)
            sage: p.facial_adjacencies()[0:3]
            [[0, [1, 2, 5, 10, 12, 13]], [1, [0, 2, 5, 7, 9, 11]], [2, [0, 1, 10, 11]]]
            sage: f0 = p.Hrepresentation(0)
            sage: f0.index() == 0
            True
            sage: f0_adjacencies = [f0.index(), [n.index() for n in f0.neighbors()]]
            sage: p.facial_adjacencies()[0] == f0_adjacencies
            True
        """
        try:
            return self._facial_adjacencies
        except AttributeError:
            self._facial_adjacencies = \
                [ [ h.index(),
                    [n.index() for n in h.neighbors()]
                  ] for h in self.Hrepresentation() ]
            return self._facial_adjacencies


    def facial_incidences(self):
        """
        Return the face-vertex incidences in the form `[f_i, [v_{i_0}, v_{i_1},\dots ,v_{i_2}]]`.

        .. NOTE::

            Instead of working with face/vertex indices, it is
            recommended that you use the
            H-representation/V-representation objects directly (see
            examples). Or use :meth:`incidence_matrix`.

        OUTPUT:

        The face indices are the indices of the H-representation
        objects, and the vertex indices are the indices of the
        V-representation objects.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[5,0,0],[0,5,0],[5,5,0],[0,0,0],[2,2,5]])
            sage: p.facial_incidences()
            [[0, [0, 1, 3, 4]],
             [1, [0, 1, 2]],
             [2, [0, 2, 3]],
             [3, [2, 3, 4]],
             [4, [1, 2, 4]]]

            sage: f0 = p.Hrepresentation(0)
            sage: f0.index() == 0
            True
            sage: f0_incidences = [f0.index(), [v.index() for v in f0.incident()]]
            sage: p.facial_incidences()[0] == f0_incidences
            True

            sage: p.incidence_matrix().column(0)
            (1, 1, 0, 1, 1)
            sage: p.incidence_matrix().column(1)
            (1, 1, 1, 0, 0)
            sage: p.incidence_matrix().column(2)
            (1, 0, 1, 1, 0)
            sage: p.incidence_matrix().column(3)
            (0, 0, 1, 1, 1)
            sage: p.incidence_matrix().column(4)
            (0, 1, 1, 0, 1)
        """
        try:
            return self._facial_incidences
        except AttributeError:
            self._facial_incidences = \
                [ [ h.index(),
                    [v.index() for v in h.incident()]
                  ] for h in self.Hrepresentation() ]
            return self._facial_incidences


    def vertex_adjacencies(self):
        """
        Return a list of vertex indices and their adjacent vertices.

        .. NOTE::

            Instead of working with vertex indices, you can use the
            V-representation objects directly (see examples).

        OUTPUT:

        The vertex indices are the indices of the V-representation
        objects.

        EXAMPLES::

            sage: permuta3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: permuta3.vertex_adjacencies()[0:3]
            [[0, [1, 2, 6]], [1, [0, 3, 7]], [2, [0, 4, 8]]]
            sage: v0 = permuta3.Vrepresentation(0)
            sage: v0.index() == 0
            True
            sage: list( v0.neighbors() )
            [A vertex at (1, 2, 4, 3), A vertex at (1, 3, 2, 4), A vertex at (2, 1, 3, 4)]
            sage: v0_adjacencies = [v0.index(), [v.index() for v in v0.neighbors()]]
            sage: permuta3.vertex_adjacencies()[0] == v0_adjacencies
            True
        """
        try:
            return self._vertex_adjacencies
        except AttributeError:
            self._vertex_adjacencies = \
                [ [ v.index(),
                    [n.index() for n in v.neighbors()]
                  ] for v in self.Vrepresentation() ]
            return self._vertex_adjacencies


    def vertex_incidences(self):
        """
        Return the vertex-face incidences in the form `[v_i, [f_{i_0}, f_{i_1},\dots ,f_{i_2}]]`.

        .. NOTE::

            Instead of working with face/vertex indices, you can use
            the H-representation/V-representation objects directly
            (see examples).

        EXAMPLES::

            sage: p = polytopes.n_simplex(3)
            sage: p.vertex_incidences()
            [[0, [0, 1, 2]], [1, [0, 1, 3]], [2, [0, 2, 3]], [3, [1, 2, 3]]]
            sage: v0 = p.Vrepresentation(0)
            sage: v0.index() == 0
            True
            sage: p.vertex_incidences()[0] == [ v0.index(), [h.index() for h in v0.incident()] ]
            True
        """
        try:
            return self._vertex_incidences
        except AttributeError:
            self._vertex_incidences = \
                [ [ v.index(),
                    [h.index() for h in v.incident()]
                  ] for v in self.Vrepresentation() ]
            return self._vertex_incidences


    def inequality_generator(self):
        """
        Return  a generator for the defining inequalities of the
        polyhedron.

        OUTPUT:

        A generator of the inequality Hrepresentation objects.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: for v in triangle.inequality_generator(): print(v)
            An inequality (1, 1) x - 1 >= 0
            An inequality (0, -1) x + 1 >= 0
            An inequality (-1, 0) x + 1 >= 0
            sage: [ v for v in triangle.inequality_generator() ]
            [An inequality (1, 1) x - 1 >= 0,
             An inequality (0, -1) x + 1 >= 0,
             An inequality (-1, 0) x + 1 >= 0]
            sage: [ [v.A(), v.b()] for v in triangle.inequality_generator() ]
            [[(1, 1), -1], [(0, -1), 1], [(-1, 0), 1]]
        """
        for H in self.Hrepresentation():
            if H.is_inequality():
                yield H


    def inequalities(self):
        """
        Return a list of inequalities as coefficient lists.

        .. NOTE::

            It is recommended to use :meth:`inequality_generator`
            instead to iterate over the list of :class:`Inequality`
            objects.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[1,0,0],[2,2,2]])
            sage: p.inequalities()[0:3]
            [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
            sage: p3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: ieqs = p3.inequalities()
            sage: ieqs[0]
            [-6, 0, 1, 1, 1]
            sage: ieqs[-1]
            [-3, 0, 1, 0, 1]
            sage: ieqs == [list(x) for x in p3.inequality_generator()]
            True
        """
        try:
            return self._inequalities
        except AttributeError:
            self._ieqs = [list(x) for x in self.inequality_generator()]
            return self._ieqs


    def ieqs(self):
        """
        Deprecated. Alias for inequalities()

        EXAMPLES::

            sage: p3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: p3.ieqs() == p3.inequalities()
            True
        """
        return self.inequalities()


    def equation_generator(self):
        """
        Return a generator for the linear equations satisfied by the
        polyhedron.

        EXAMPLES::

            sage: p = polytopes.regular_polygon(8,base_ring=RDF)
            sage: p3 = Polyhedron(vertices = [x+[0] for x in p.vertices()], base_ring=RDF)
            sage: p3.equation_generator().next()
            An equation (0.0, 0.0, 1.0) x + 0.0 == 0
        """
        for H in self.Hrepresentation():
            if H.is_equation():
                yield H


    def equations(self):
        """
        Return the linear constraints of the polyhedron. As with
        inequalities, each constraint is given as [b -a1 -a2 ... an]
        where for variables x1, x2,..., xn, the polyhedron satisfies
        the equation b = a1*x1 + a2*x2 + ... + an*xn.

        .. NOTE::

            It is recommended to use :meth:`equation_generator()` instead
            to iterate over the list of :class:`Equation` objects.

        EXAMPLES::

            sage: test_p = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1],[3,4,1,2]])
            sage: test_p.equations()
            [[-10, 1, 1, 1, 1]]
        """
        try:
            return self._equations
        except AttributeError:
            self._equations = [list(eq) for eq in self.equation_generator()]
            return self._equations


    def linearities(self):
        """
        Deprecated.  Use equations() instead.
        Returns the linear constraints of the polyhedron. As with
        inequalities, each constraint is given as [b -a1 -a2 ... an]
        where for variables x1, x2,..., xn, the polyhedron satisfies
        the equation b = a1*x1 + a2*x2 + ... + an*xn.

        EXAMPLES::

            sage: test_p = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1],[3,4,1,2]])
            sage: test_p.linearities()
            [[-10, 1, 1, 1, 1]]
            sage: test_p.linearities() == test_p.equations()
            True
        """
        return self.equations()


    def vertices(self):
        """
        Return a list of vertices of the polyhedron.

        .. NOTE::

            It is recommended to use :meth:`vertex_generator` instead to
            iterate over the list of :class:`Vertex` objects.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: triangle.vertices()
            [[0, 1], [1, 0], [1, 1]]
            sage: a_simplex = Polyhedron(ieqs = [
            ...            [0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]
            ...        ], eqns = [[1,-1,-1,-1,-1]])
            sage: a_simplex.vertices()
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
            sage: a_simplex.vertices() == [list(v) for v in a_simplex.vertex_generator()]
            True
        """
        try:
            return self._vertices
        except AttributeError:
            self._vertices = [list(x) for x in self.vertex_generator()]
            return self._vertices


    def vertex_generator(self):
        """
        Return a generator for the vertices of the polyhedron.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: for v in triangle.vertex_generator(): print(v)
            A vertex at (0, 1)
            A vertex at (1, 0)
            A vertex at (1, 1)
            sage: v_gen = triangle.vertex_generator()
            sage: v_gen.next()   # the first vertex
            A vertex at (0, 1)
            sage: v_gen.next()   # the second vertex
            A vertex at (1, 0)
            sage: v_gen.next()   # the third vertex
            A vertex at (1, 1)
            sage: try: v_gen.next()   # there are only three vertices
            ... except StopIteration: print "STOP"
            STOP
            sage: type(v_gen)
            <type 'generator'>
            sage: [ v for v in triangle.vertex_generator() ]
            [A vertex at (0, 1), A vertex at (1, 0), A vertex at (1, 1)]
        """
        for V in self.Vrepresentation():
            if V.is_vertex():
                yield V


    def ray_generator(self):
        """
        Return a generator for the rays of the polyhedron.

        EXAMPLES::

            sage: pi = Polyhedron(ieqs = [[1,1,0],[1,0,1]])
            sage: pir = pi.ray_generator()
            sage: [x.vector() for x in pir]
            [(1, 0), (0, 1)]
        """
        for V in self.Vrepresentation():
            if V.is_ray():
                yield V


    def rays(self):
        """
        Return a list of rays as coefficient lists.

        .. NOTE::

            It is recommended to use :meth:`ray_generator` instead to
            iterate over the list of :class:`Ray` objects.

        OUTPUT:

        A list of rays as lists of coordinates.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0],[1,1,0,0]])
            sage: p.rays()
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            sage: p.rays() == [list(r) for r in p.ray_generator()]
            True
        """
        try:
            return self._rays
        except AttributeError:
            self._rays = [list(x) for x in self.ray_generator()]
            return self._rays


    def line_generator(self):
        """
        Return a generator for the lines of the polyhedron.

        EXAMPLES::

            sage: pr = Polyhedron(rays = [[1,0],[-1,0],[0,1]], vertices = [[-1,-1]])
            sage: pr.line_generator().next().vector()
            (1, 0)
        """
        for V in self.Vrepresentation():
            if V.is_line():
                yield V


    def lines(self):
        """
        Return a list of lines of the polyhedron.  The line data is given
        as a list of coordinates rather than as a Hrepresentation object.

        .. NOTE::

            It is recommended to use :meth:`line_generator` instead to
            iterate over the list of :class:`Line` objects.

        EXAMPLES::

            sage: p = Polyhedron(rays = [[1,0],[-1,0],[0,1],[1,1]], vertices = [[-2,-2],[2,3]])
            sage: p.lines()
            [[1, 0]]
            sage: p.lines() == [list(x) for x in p.line_generator()]
            True
        """
        try:
            return self._lines
        except AttributeError:
            self._lines = [list(x) for x in self.line_generator()]
            return self._lines


    def bounded_edges(self):
        """
        Return the bounded edges (excluding rays and lines).

        OUTPUT:

        A generator for pairs of vertices, one pair per edge.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1,0],[0,1]], rays=[[1,0],[0,1]])
            sage: [ e for e in p.bounded_edges() ]
            [(A vertex at (0, 1), A vertex at (1, 0))]
            sage: for e in p.bounded_edges(): print e
            (A vertex at (0, 1), A vertex at (1, 0))
        """
        obj = self.Vrepresentation()
        edges = []
        for i in range(len(obj)):
            if not obj[i].is_vertex(): continue
            for j in range(i+1,len(obj)):
                if not obj[j].is_vertex(): continue
                if self.vertex_adjacency_matrix()[i,j] == 0: continue
                yield (obj[i], obj[j])


    @cached_method
    def ambient_space(self):
        r"""
        Return the ambient vector space.

        OUTPUT:

        A free module over the base ring of dimension :meth:`ambient_dim`.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.ambient_space()
            Vector space of dimension 4 over Rational Field
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self.base_ring(), self.ambient_dim())

    Vrepresentation_space = ambient_space

    @cached_method
    def Hrepresentation_space(self):
        r"""
        Return the linear space containing the H-representation vectors.

        OUTPUT:

        A free module over the base ring of dimension :meth:`ambient_dim` + 1.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.ambient_space()
            Vector space of dimension 4 over Rational Field
        """
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self.base_ring(), self.ambient_dim()+1)


    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.ambient_dim()
            4
        """
        return self._ambient_dim


    def dim(self):
        """
        Return the dimension of the polyhedron.

        EXAMPLES::

            sage: simplex = Polyhedron(vertices = [[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]])
            sage: simplex.dim()
            3
            sage: simplex.ambient_dim()
            4
       """
        return self.ambient_dim() - self.n_equations()


    def adjacency_matrix(self):
        """
        This is an alias for :meth:`vertex_adjacency_matrix`

        EXAMPLES::

            sage: polytopes.n_cube(3).adjacency_matrix()
            [0 1 1 0 1 0 0 0]
            [1 0 0 1 0 1 0 0]
            [1 0 0 1 0 0 1 0]
            [0 1 1 0 0 0 0 1]
            [1 0 0 0 0 1 1 0]
            [0 1 0 0 1 0 0 1]
            [0 0 1 0 1 0 0 1]
            [0 0 0 1 0 1 1 0]
        """
        return self.vertex_adjacency_matrix()


    def vertex_adjacency_matrix(self):
        """
        Return the binary matrix of vertex adjacencies.

        EXAMPLES::

            sage: polytopes.n_simplex(4).vertex_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]
        """
        if '_V_adjacency_matrix' not in self.__dict__:
            self._init_vertex_adjacency_matrix()
        return self._V_adjacency_matrix;


    def facet_adjacency_matrix(self):
        """
        Return the adjacency matrix for the facets and hyperplanes.

        EXAMPLES::

            sage: polytopes.n_simplex(4).facet_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]
        """
        if '_H_adjacency_matrix' not in self.__dict__:
            self._init_facet_adjacency_matrix()
        return self._H_adjacency_matrix;


    def incidence_matrix(self):
        """
        Return the incidence matrix.

        .. NOTE::

            The columns correspond to inequalities/equations in the
            order :meth:`Hrepresentation`, the rows correspond to
            vertices/rays/lines in the order
            :meth:`Vrepresentation`

        EXAMPLES::

            sage: p = polytopes.cuboctahedron()
            sage: p.incidence_matrix()
            [0 0 1 1 0 1 0 0 0 0 1 0 0 0]
            [0 0 0 1 0 0 1 0 1 0 1 0 0 0]
            [0 0 1 1 1 0 0 1 0 0 0 0 0 0]
            [1 0 0 1 1 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 1 0 0 1 1 1 0 0 0]
            [0 0 1 0 0 1 0 1 0 0 0 1 0 0]
            [1 0 0 0 0 0 1 0 1 0 0 0 1 0]
            [1 0 0 0 1 0 0 1 0 0 0 0 0 1]
            [0 1 0 0 0 1 0 0 0 1 0 1 0 0]
            [0 1 0 0 0 0 0 0 1 1 0 0 1 0]
            [0 1 0 0 0 0 0 1 0 0 0 1 0 1]
            [1 1 0 0 0 0 0 0 0 0 0 0 1 1]
            sage: v = p.Vrepresentation(0)
            sage: v
            A vertex at (-1/2, -1/2, 0)
            sage: h = p.Hrepresentation(2)
            sage: h
            An inequality (1, 1, -1) x + 1 >= 0
            sage: h.eval(v)        # evaluation (1, 1, -1) * (-1/2, -1/2, 0) + 1
            0
            sage: h*v              # same as h.eval(v)
            0
            sage: p.incidence_matrix() [0,2]   # this entry is (v,h)
            1
            sage: h.contains(v)
            True
            sage: p.incidence_matrix() [2,0]   # note: not symmetric
            0
        """
        try:
            return self._incidence_matrix
        except AttributeError:
            self._incidence_matrix = matrix(ZZ, len(self.Vrepresentation()),
                                                len(self.Hrepresentation()), 0)
            for V in self.Vrep_generator():
                for H in self.Hrep_generator():
                    if self._is_zero(H*V):
                        self._incidence_matrix[V.index(),H.index()] = 1

            return self._incidence_matrix


    def base_ring(self):
        """
        Return the base ring.

        OUTPUT:

        Either ``QQ`` (exact arithmetic using gmp, default) or ``RDF``
        (double precision floating-point arithmetic)

        EXAMPLES::

            sage: triangle = Polyhedron(vertices = [[1,0],[0,1],[1,1]])
            sage: triangle.base_ring() == QQ
            True
        """
        return self._base_ring

    field = base_ring


    def coerce_field(self, other):
        """
        Return the common field for both ``self`` and ``other``.

        INPUT:

        - ``other`` -- must be either:

            * another ``Polyhedron`` object

            * `\QQ` or `RDF`

            * a constant that can be coerced to `\QQ` or `RDF`

        OUTPUT:

        Either `\QQ` or `RDF`. Raises ``TypeError`` if ``other`` is not a
        suitable input.

        .. NOTE::

            "Real" numbers in sage are not necessarily elements of
            `RDF`. For example, the literal `1.0` is not.

        EXAMPLES::

            sage: triangle_QQ  = Polyhedron(vertices = [[1,0],[0,1],[1,1]], base_ring=QQ)
            sage: triangle_RDF = Polyhedron(vertices = [[1,0],[0,1],[1,1]], base_ring=RDF)
            sage: triangle_QQ.coerce_field(QQ)
            Rational Field
            sage: triangle_QQ.coerce_field(triangle_RDF)
            Real Double Field
            sage: triangle_RDF.coerce_field(triangle_QQ)
            Real Double Field
            sage: triangle_QQ.coerce_field(RDF)
            Real Double Field
            sage: triangle_QQ.coerce_field(ZZ)
            Rational Field
            sage: triangle_QQ.coerce_field(1/2)
            Rational Field
            sage: triangle_QQ.coerce_field(0.5)
            Real Double Field
        """
        try:
            # other is a Polyhedron object?
            other_field = other.field()
        except AttributeError:

            try:
                # other is a constant?
                other_parent = other.parent()
            except AttributeError:
                other_parent = other

            # other is a field?
            if QQ.coerce_map_from(other_parent) != None:
                other_field = QQ
            elif RDF.coerce_map_from(other_parent) != None:
                other_field = RDF
            else:
                raise TypeError("cannot determine field from %s!" % other)

        assert other_field==QQ or other_field==RDF

        if self.field()==RDF or other_field==RDF:
            return RDF
        else:
            return QQ


    @cached_method
    def center(self):
        """
        Return the average of the vertices.

        OUTPUT:

        The center of the polyhedron. All rays and lines are
        ignored. Raises a ``ZeroDivisionError`` for the empty
        polytope.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p = p + vector([1,0,0])
            sage: p.center()
            (1, 0, 0)
        """
        vertex_sum = vector(ZZ, [0]*self.ambient_dim())
        for v in self.vertex_generator():
            vertex_sum += v.vector()
        return vertex_sum / self.n_vertices()


    @cached_method
    def radius_square(self):
        """
        Return the square of the maximal distance from the
        :meth:`center` to a vertex. All rays and lines are ignored.

        OUTPUT:

        The square of the radius, which is in :meth:`field`.

        EXAMPLES::

            sage: p = polytopes.permutahedron(4, project = False)
            sage: p.radius_square()
            5
        """
        vertices = [ v.vector() - self.center() for v in self.vertex_generator() ]
        return max( v.dot_product(v) for v in vertices )


    def radius(self):
        """
        Return the maximal distance from the center to a vertex. All
        rays and lines are ignored.

        OUTPUT:

        The radius for a rational polyhedron is, in general, not
        rational.  use :meth:`radius_square` if you need a rational
        distance measure.

        EXAMPLES::

            sage: p = polytopes.n_cube(4)
            sage: p.radius()
            2
        """
        return sqrt(self.radius_square())


    def is_compact(self):
        """
        Test for boundedness of the polytope.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: p.is_compact()
            True
            sage: p = Polyhedron(ieqs = [[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,-1,0,0]])
            sage: p.is_compact()
            False
        """
        return self.n_rays()==0 and self.n_lines()==0


    def is_simple(self):
        """
        Test for simplicity of a polytope.

        EXAMPLES::

            sage: p = Polyhedron([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
            sage: p.is_simple()
            True
            sage: p = Polyhedron([[0,0,0],[4,4,0],[4,0,0],[0,4,0],[2,2,2]])
            sage: p.is_simple()
            False

        REFERENCES:

            http://en.wikipedia.org/wiki/Simple_polytope
        """
        if not self.is_compact(): return False

        for v in self.vertex_generator():
            adj = [a for a in v.neighbors()]
            if len(adj) != self.dim():
                return False

        return True

    def gale_transform(self):
        """
        Return the Gale transform of a polytope as described in the
        reference below.

        OUTPUT:

        A list of vectors, the Gale transform.  The dimension is the
        dimension of the affine dependencies of the vertices of the
        polytope.

        EXAMPLES:

        This is from the reference, for a triangular prism::

            sage: p = Polyhedron(vertices = [[0,0],[0,1],[1,0]])
            sage: p2 = p.prism()
            sage: p2.gale_transform()
            [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)]

        REFERENCES:

            Lectures in Geometric Combinatorics, R.R.Thomas, 2006, AMS Press.
        """
        if not self.is_compact(): raise ValueError('Not a polytope.')

        A = matrix(self.n_vertices(),
                   [ [1]+list(x) for x in self.vertex_generator()])
        A = A.transpose()
        A_ker = A.right_kernel()
        return A_ker.basis_matrix().transpose().rows()

    def triangulate(self, engine='auto', connected=True, fine=False, regular=None, star=None):
        r"""
        Returns a triangulation of the polytope.

        INPUT:

        - ``engine`` -- either 'auto' (default), 'internal', or
          'TOPCOM'.  The latter two instruct this package to always
          use its own triangulation algorithms or TOPCOM's algorithms,
          respectively. By default ('auto'), TOPCOM is used if it is
          available and internal routines otherwise.

        The remaining keyword parameters are passed through to the
        :class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`
        constructor:

        - ``connected`` -- boolean (default: ``True``). Whether the
          triangulations should be connected to the regular
          triangulations via bistellar flips. These are much easier to
          compute than all triangulations.

        - ``fine`` -- boolean (default: ``False``). Whether the
          triangulations must be fine, that is, make use of all points
          of the configuration.

        - ``regular`` -- boolean or ``None`` (default:
          ``None``). Whether the triangulations must be regular. A
          regular triangulation is one that is induced by a
          piecewise-linear convex support function. In other words,
          the shadows of the faces of a polyhedron in one higher
          dimension.

          * ``True``: Only regular triangulations.

          * ``False``: Only non-regular triangulations.

          * ``None`` (default): Both kinds of triangulation.

        - ``star`` -- either ``None`` (default) or a point. Whether
          the triangulations must be star. A triangulation is star if
          all maximal simplices contain a common point. The central
          point can be specified by its index (an integer) in the
          given points or by its coordinates (anything iterable.)

        OUTPUT:

        A triangulation of the convex hull of the vertices as a
        :class:`~sage.geometry.triangulation.point_configuration.Triangulation`. The
        indices in the triangulation correspond to the
        :meth:`Vrepresentation` objects.

        EXAMPLES::

            sage: cube = polytopes.n_cube(3)
            sage: triangulation = cube.triangulate(
            ...      engine='internal') # to make doctest independent of TOPCOM
            sage: triangulation
            (<0,1,2,7>, <0,1,4,7>, <0,2,4,7>, <1,2,3,7>, <1,4,5,7>, <2,4,6,7>)
            sage: simplex_indices = triangulation[0]; simplex_indices
            (0, 1, 2, 7)
            sage: simplex_vertices = [ cube.Vrepresentation(i) for i in simplex_indices ]
            sage: simplex_vertices
            [A vertex at (-1, -1, -1), A vertex at (-1, -1, 1),
             A vertex at (-1, 1, -1), A vertex at (1, 1, 1)]
            sage: Polyhedron(simplex_vertices)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
        """
        if not self.is_compact():
            raise NotImplementedError('I can only triangulate compact polytopes.')
        from sage.geometry.triangulation.point_configuration import PointConfiguration
        pc = PointConfiguration((v.vector() for v in self.vertex_generator()),
                                connected=connected, fine=fine, regular=regular, star=star)
        pc.set_engine(engine)
        return pc.triangulate()

    def triangulated_facial_incidences(self):
        """
        Return a list of the form [face_index, [v_i_0,
        v_i_1,...,v_i_{n-1}]] where the face_index refers to the
        original defining inequality.  For a given face, the
        collection of triangles formed by each list of v_i should
        triangulate that face.

        In dimensions greater than 3, this is computed by randomly
        lifting each face up a dimension; this does not always work!
        This should eventually be fixed by using lrs or another
        program that computes triangulations.

        EXAMPLES:

        If the figure is already composed of triangles, then all is well::

            sage: Polyhedron(vertices = [[5,0,0],[0,5,0],[5,5,0],[2,2,5]]
            ...             ).triangulated_facial_incidences()
            doctest:...: DeprecationWarning: This method is
            deprecated. Use triangulate() instead.
            See http://trac.sagemath.org/11634 for details.
            [[0, [0, 1, 2]], [1, [0, 1, 3]], [2, [0, 2, 3]], [3, [1, 2, 3]]]

        Otherwise some faces get split up to triangles::

            sage: Polyhedron(vertices = [[2,0,0],[4,1,0],[0,5,0],[5,5,0],
            ...       [1,1,0],[0,0,1]]).triangulated_facial_incidences()
            doctest:...: DeprecationWarning: This method is
            deprecated. Use triangulate() instead.
            See http://trac.sagemath.org/11634 for details.
            [[0, [1, 2, 5]], [0, [2, 5, 3]], [0, [5, 3, 4]], [1, [0, 1, 2]],
             [2, [0, 2, 3]], [3, [0, 3, 4]], [4, [0, 4, 5]], [5, [0, 1, 5]]]
        """
        from sage.misc.superseded import deprecation
        deprecation(11634, 'This method is deprecated. Use triangulate() instead.')
        try:
            return self._triangulated_facial_incidences
        except AttributeError:
            t_fac_incs = []
            for a_face in self.facial_incidences():
                vert_number = len(a_face[1])
                if vert_number == self.dim():
                    t_fac_incs.append(a_face)
                elif self.dim() >= 4:
                    lifted_verts = []
                    for vert_index in a_face[1]:
                        lifted_verts.append(self.vertices()[vert_index] +
                                            [randint(-vert_index,5000+vert_index + vert_number**2)])
                    temp_poly = Polyhedron(vertices = lifted_verts)
                    for t_face in temp_poly.facial_incidences():
                        if len(t_face[1]) != self.dim():
                            print 'Failed for face: ' + str(a_face)
                            print 'Attempted simplicial face: ' + str(t_face)
                            print 'Attempted lifted vertices: ' + str(lifted_verts)
                            raise RuntimeError, "triangulation failed"
                        normal_fdir = temp_poly.ieqs()[t_face[0]][-1]
                        if normal_fdir >= 0:
                            t_fac_verts = [temp_poly.vertices()[i] for i in t_face[1]]
                            proj_verts = [q[0:self.dim()] for q in t_fac_verts]
                            t_fac_incs.append([a_face[0],
                                               [self.vertices().index(q) for q in proj_verts]])
                else:
                    vs = a_face[1][:]
                    adj = dict([a[0], filter(lambda p: p in a_face[1], a[1])]
                               for a in filter(lambda va: va[0] in a_face[1],
                                               self.vertex_adjacencies()))
                    t = vs[0]
                    vs.remove(t)
                    ts = adj[t]
                    for v in ts:
                        vs.remove(v)
                    t_fac_incs.append([a_face[0], [t] + ts])
                    while vs:
                        t = ts[0]
                        ts = ts[1:]
                        for v in adj[t]:
                            if v in vs:
                                vs.remove(v)
                                ts.append(v)
                                t_fac_incs.append([a_face[0], [t] + ts])
                                break
        self._triangulated_facial_incidences = t_fac_incs
        return t_fac_incs


    def simplicial_complex(self):
        """
        Return a simplicial complex from a triangulation of the polytope.

        Warning: This first triangulates the polytope using
        ``triangulated_facial_incidences``, and this function may fail
        in dimensions greater than 3, although it usually doesn't.

        OUTPUT:

        A simplicial complex.

        EXAMPLES::

            sage: p = polytopes.cuboctahedron()
            sage: sc = p.simplicial_complex()
            doctest:...: DeprecationWarning:
            This method is deprecated. Use triangulate().simplicial_complex() instead.
            See http://trac.sagemath.org/11634 for details.
            doctest:...: DeprecationWarning:
            This method is deprecated. Use triangulate() instead.
            See http://trac.sagemath.org/11634 for details.
            sage: sc
            Simplicial complex with 12 vertices and 20 facets
        """
        from sage.misc.superseded import deprecation
        deprecation(11634, 'This method is deprecated. Use triangulate().simplicial_complex() instead.')
        from sage.homology.simplicial_complex import SimplicialComplex
        return SimplicialComplex([x[1] for x in self.triangulated_facial_incidences()])

    def __add__(self, other):
        """
        The Minkowski sum of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        EXAMPLES::

            sage: four_cube = polytopes.n_cube(4)
            sage: four_simplex = Polyhedron(vertices = [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: unholy_union = four_cube + four_simplex
            sage: unholy_union.dim()
            4
            sage: poly_spam = Polyhedron([[3,4,5,2],[1,0,0,1],[0,0,0,0],[0,4,3,2],[-3,-3,-3,-3]])
            sage: poly_eggs = Polyhedron([[5,4,5,4],[-4,5,-4,5],[4,-5,4,-5],[0,0,0,0]])
            sage: poly_spam_and_eggs = poly_spam + poly_spam + poly_eggs
            sage: poly_spam_and_eggs.n_vertices()
            12
        """
        if is_Polyhedron(other):
            new_vertices = []
            for v1 in self.vertex_generator():
                for v2 in other.vertex_generator():
                    new_vertices.append(list(v1() + v2()))
            new_rays = self.rays() + other.rays()
            new_lines = self.lines() + other.lines()
            other_field = other.field()

        else:  # assume other is a vector and try to add vertices
            displacement = vector(other)
            new_vertices = [list(x() + displacement) for x in self.vertex_generator()]
            new_rays = self.rays()
            new_lines = self.lines()
            other_field = displacement.base_ring()

        return Polyhedron(vertices=new_vertices,
                          rays=new_rays, lines=new_lines,
                          base_ring=self.coerce_field(other_field))


    def __mul__(self, other):
        """
        Multiplication by ``other``.

        INPUT:

        - ``other`` -- A scalar, not necessarily in :meth:`field`, or
          a :class:`Polyhedron`.

        OUTPUT:

        Multiplication by another polyhedron returns the product
        polytope. Multiplication by a scalar returns the polytope
        dilated by that scalar, possibly coerced to the bigger field.

        EXAMPLES::

             sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,6)])
             sage: p.vertex_generator().next()
             A vertex at (2, 4, 8)
             sage: p2 = p*2
             sage: p2.vertex_generator().next()
             A vertex at (4, 8, 16)
        """
        if is_Polyhedron(other):
            new_vertices = [ list(x)+list(y)
                             for x in self.vertex_generator() for y in other.vertex_generator()]
            new_rays = []
            new_rays.extend( [ list(r)+[0]*other.ambient_dim()
                               for r in self.ray_generator() ] )
            new_rays.extend( [ [0]*self.ambient_dim()+list(r)
                               for r in other.ray_generator() ] )
            new_lines = []
            new_lines.extend( [ list(l)+[0]*other.ambient_dim()
                                for l in self.line_generator() ] )
            new_lines.extend( [ [0]*self.ambient_dim()+list(l)
                               for l in other.line_generator() ] )
        else:
            new_vertices = [ list(other*v()) for v in self.vertex_generator()]
            new_rays =  self.rays()
            new_lines = self.lines()

        return Polyhedron(vertices=new_vertices,
                          rays=new_rays, lines=new_lines,
                          base_ring=self.coerce_field(other))


    def __rmul__(self,other):
        """
        Right multiplication.

        See :meth:`__mul__` for details.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,4)])
            sage: p2 = 3*p + p
            sage: p2.vertex_generator().next()
            A vertex at (8, 16, 32)
        """
        return self.__mul__(other)


    def union(self, other):
        """
        Deprecated.  Use ``self.convex_hull(other)`` instead.

        EXAMPLES::

            sage: Polyhedron(vertices=[[0]]).union( Polyhedron(vertices=[[1]]) )
            doctest:...: DeprecationWarning:
            The function union is replaced by convex_hull.
            See http://trac.sagemath.org/11634 for details.
            A 1-dimensional polyhedron in QQ^1 defined as the convex hull of 2 vertices
        """
        from sage.misc.superseded import deprecation
        deprecation(11634, 'The function union is replaced by convex_hull.')
        return self.convex_hull(other)


    def convex_hull(self, other):
        """
        Return the convex hull of the set-theoretic union of the two
        polyhedra.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        The convex hull.

        EXAMPLES::

            sage: a_simplex = polytopes.n_simplex(3)
            sage: verts = a_simplex.vertices()
            sage: verts = [[x[0]*3/5+x[1]*4/5, -x[0]*4/5+x[1]*3/5, x[2]] for x in verts]
            sage: another_simplex = Polyhedron(vertices = verts)
            sage: simplex_union = a_simplex.convex_hull(another_simplex)
            sage: simplex_union.n_vertices()
            7
        """
        hull_vertices = self.vertices() + other.vertices()
        hull_rays = self.rays() + other.rays()
        hull_lines = self.lines() + other.lines()
        hull_field = self.coerce_field(other)
        return Polyhedron(vertices=hull_vertices,
                          rays=hull_rays, lines=hull_lines,
                          base_ring=hull_field)


    def intersection(self, other):
        """
        Return the intersection of one polyhedron with another.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        The intersection.

        EXAMPLES::

            sage: cube = polytopes.n_cube(3)
            sage: oct = polytopes.cross_polytope(3)
            sage: cube_oct = cube.intersection(oct*2)
            sage: len(list( cube_oct.vertex_generator() ))
            12
            sage: cube_oct
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 12 vertices
        """
        new_ieqs = []
        new_ieqs.extend(self.inequalities())
        new_ieqs.extend(other.inequalities())

        new_eqns = []
        new_eqns.extend(self.equations())
        new_eqns.extend(other.equations())

        return Polyhedron(ieqs = new_ieqs, eqns = new_eqns,
                          base_ring=self.coerce_field(other))


    def edge_truncation(self, cut_frac = Integer(1)/3):
        r"""
        Return a new polyhedron formed from two points on each edge
        between two vertices.

        INPUT:

        - ``cut_frac`` -- integer. how deeply to cut into the edge.
            Default is `\frac{1}{3}`.

        OUTPUT:

        A Polyhedron object, truncated as described above.

        EXAMPLES::

            sage: cube = polytopes.n_cube(3)
            sage: trunc_cube = cube.edge_truncation()
            sage: trunc_cube.n_vertices()
            24
            sage: trunc_cube.n_inequalities()
            14
        """
        new_vertices = []
        for e in self.bounded_edges():
            new_vertices.append((1-cut_frac)*e[0]() + cut_frac *e[1]())
            new_vertices.append(cut_frac *e[0]() + (1-cut_frac)*e[1]())

        new_vertices = [list(v) for v in new_vertices]
        new_rays =  self.rays()
        new_lines = self.lines()

        return Polyhedron(vertices=new_vertices, rays=new_rays,
                          lines=new_lines,
                          base_ring=self.coerce_field(cut_frac))


    def _make_polyhedron_face(self, Vindices, Hindices):
        """
        Construct a face of the polyhedron.

        INPUT:

        - ``Vindices`` -- a tuple of integers. The indices of the
          V-represenation objects that span the face.

        - ``Hindices`` -- a tuple of integers. The indices of the
          H-representation objects that hold as equalities on the
          face.

        OUTPUT:

        A new :class:`PolyhedronFace_base` instance. It is not checked
        whether the input data actually defines a face.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: square._make_polyhedron_face((0,2), (1,))
            <0,2>
        """
        return PolyhedronFace_base(self, Vindices, Hindices)


    def face_lattice(self):
        """
        Return the face-lattice poset.

        OUTPUT:

        A :class:`~sage.combinat.posets.posets.FinitePoset`. Elements
        are given as
        :class:`~sage.geometry.polyhedron.PolyhedronFace_base`.

        In the case of a full-dimensional polytope, the faces are
        pairs (vertices, inequalities) of the spanning vertices and
        corresponding saturated inequalities. In general, a face is
        defined by a pair (V-rep. objects, H-rep. objects). The
        V-representation objects span the face, and the corresponding
        H-representation objects are those inequalities and equations
        that are saturated on the face.

        The bottom-most element of the face lattice is the "empty
        face". It contains no V-representation object. All
        H-representation objects are incident.

        The top-most element is the "full face". It is spanned by all
        V-representation objects. The incident H-representation
        objects are all equations and no inequalities.

        In the case of a full-dimensional polytope, the "empty face"
        and the "full face" are the empty set (no vertices, all
        inequalities) and the full polytope (all vertices, no
        inequalities), respectively.

        ALGORITHM:

        For a full-dimensional polytope, the basic algorithm is
        described in
        :func:`~sage.geometry.hasse_diagram.Hasse_diagram_from_incidences`.
        There are three generalizations of [KP2002]_ necessary to deal
        with more general polytopes, corresponding to the extra
        H/V-representation objects:

        * Lines are removed before calling
          :func:`Hasse_diagram_from_incidences`, and then added back
          to each face V-representation except for the "empty face".

        * Equations are removed before calling
          :func:`Hasse_diagram_from_incidences`, and then added back
          to each face H-representation.

        * Rays: Consider the half line as an example. The
          V-representation objects are a point and a ray, which we can
          think of as a point at infinity. However, the point at
          infinity has no inequality associated to it, so there is
          only one H-representation object alltogether. The face
          lattice does not contain the "face at infinity". This means
          that in :func:`Hasse_diagram_from_incidences`, one needs to
          drop faces with V-representations that have no matching
          H-representation. In addition, one needs to ensure that
          every non-empty face contains at least one vertex.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: square.face_lattice()
            Finite poset containing 10 elements
            sage: list(_)
            [<>, <0>, <1>, <2>, <3>, <0,1>, <0,2>, <2,3>, <1,3>, <0,1,2,3>]
            sage: poset_element = _[6]
            sage: a_face = poset_element.element
            sage: a_face
            <0,2>
            sage: a_face.dim()
            1
            sage: set(a_face.ambient_Vrepresentation()) == \
            ...   set([square.Vrepresentation(0), square.Vrepresentation(2)])
            True
            sage: a_face.ambient_Vrepresentation()
            (A vertex at (-1, -1), A vertex at (1, -1))
            sage: a_face.ambient_Hrepresentation()
            (An inequality (0, 1) x + 1 >= 0,)

        A more complicated example::

            sage: c5_10 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5] for i in range(1,11)])
            sage: c5_10_fl = c5_10.face_lattice()
            sage: [len(x) for x in c5_10_fl.level_sets()]
            [1, 10, 45, 100, 105, 42, 1]

        Note that if the polyhedron contains lines then there is a
        dimension gap between the empty face and the first non-empty
        face in the face lattice::

            sage: line = Polyhedron(vertices=[(0,)], lines=[(1,)])
            sage: [ fl.element.dim() for fl in line.face_lattice() ]
            [-1, 1]

        TESTS::

            sage: c5_20 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5]
            ...       for i in range(1,21)])
            sage: c5_20_fl = c5_20.face_lattice() # long time
            sage: [len(x) for x in c5_20_fl.level_sets()] # long time
            [1, 20, 190, 580, 680, 272, 1]
            sage: polytopes.n_cube(2).face_lattice().plot()
            sage: level_sets = polytopes.cross_polytope(2).face_lattice().level_sets()
            sage: print level_sets[0], level_sets[-1]
            [<>] [<0,1,2,3>]

        Various degenerate polyhedra::

            sage: Polyhedron(vertices=[[0,0,0],[1,0,0],[0,1,0]]).face_lattice().level_sets()
            [[<>], [<0>, <1>, <2>], [<0,1>, <0,2>, <1,2>], [<0,1,2>]]
            sage: Polyhedron(vertices=[(1,0,0),(0,1,0)], rays=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<1>, <2>], [<0,1>, <0,2>, <1,2>], [<0,1,2>]]
            sage: Polyhedron(rays=[(1,0,0),(0,1,0)], vertices=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<0>], [<0,1>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(rays=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<0>], [<0,1>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(vertices=[(1,),(0,)]).face_lattice().level_sets()
            [[<>], [<0>, <1>], [<0,1>]]
            sage: Polyhedron(vertices=[(1,0,0),(0,1,0)], lines=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<0,1>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(lines=[(1,0,0)], vertices=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<0,1>]]
            sage: Polyhedron(lines=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<0,1,2>]]
            sage: Polyhedron(lines=[(1,0)], rays=[(0,1)], vertices=[(0,0)])\
            ...       .face_lattice().level_sets()
            [[<>], [<0,1>], [<0,1,2>]]
            sage: Polyhedron(vertices=[(0,)], lines=[(1,)]).face_lattice().level_sets()
            [[<>], [<0,1>]]
            sage: Polyhedron(lines=[(1,0)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<0,1>]]

        REFERENCES:

        ..  [KP2002]

            Volker Kaibel and Marc E. Pfetsch, "Computing the Face
            Lattice of a Polytope from its Vertex-Facet Incidences",
            Computational Geometry: Theory and Applications, Volume
            23, Issue 3 (November 2002), 281-290.  Available at
            http://portal.acm.org/citation.cfm?id=763203 and free of
            charge at http://arxiv.org/abs/math/0106043
        """
        try:
            return self._face_lattice
        except AttributeError:
            pass

        coatom_to_Hindex = [ h.index() for h in self.inequality_generator() ]
        Hindex_to_coatom = [None] * self.n_Hrepresentation()
        for i in range(0,len(coatom_to_Hindex)):
            Hindex_to_coatom[ coatom_to_Hindex[i] ] = i

        atom_to_Vindex = [ v.index() for v in self.Vrep_generator() if not v.is_line() ]
        Vindex_to_atom = [None] * self.n_Vrepresentation()
        for i in range(0,len(atom_to_Vindex)):
                        Vindex_to_atom[ atom_to_Vindex[i] ] = i

        atoms_incidences   = [ tuple([ Hindex_to_coatom[h.index()]
                                       for h in v.incident() if h.is_inequality() ])
                               for v in self.Vrepresentation() if not v.is_line() ]

        coatoms_incidences = [ tuple([ Vindex_to_atom[v.index()]
                                       for v in h.incident() if not v.is_line() ])
                               for h in self.Hrepresentation() if h.is_inequality() ]

        atoms_vertices = [ Vindex_to_atom[v.index()] for v in self.vertex_generator() ]
        equations = [ e.index() for e in self.equation_generator() ]
        lines     = [ l.index() for l in self.line_generator() ]

        def face_constructor(atoms,coatoms):
            if len(atoms)==0:
                Vindices = ()
            else:
                Vindices = tuple(sorted([   atom_to_Vindex[i] for i in   atoms ]+lines))
            Hindices = tuple(sorted([ coatom_to_Hindex[i] for i in coatoms ]+equations))
            return self._make_polyhedron_face(Vindices, Hindices)

        from sage.geometry.hasse_diagram import Hasse_diagram_from_incidences
        self._face_lattice = Hasse_diagram_from_incidences\
            (atoms_incidences, coatoms_incidences,
             face_constructor=face_constructor, required_atoms=atoms_vertices)
        return self._face_lattice


    def f_vector(self):
        r"""
        Return the f-vector.

        OUTPUT:

        Returns a vector whose ``i``-th entry is the number of
        ``i``-dimensional faces of the polytope.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1, 2, 3], [1, 3, 2],
            ...       [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1], [0, 0, 0]])
            sage: p.f_vector()
            (1, 7, 12, 7, 1)
        """
        try:
            return self._f_vector
        except AttributeError:
            self._f_vector = vector(ZZ,[len(x) for x in self.face_lattice().level_sets()])
            return self._f_vector


    def vertex_graph(self):
        """
        Return a graph in which the vertices correspond to vertices
        of the polyhedron, and edges to edges.

        EXAMPLES::

            sage: g3 = polytopes.n_cube(3).vertex_graph()
            sage: len(g3.automorphism_group())
            48
            sage: s4 = polytopes.n_simplex(4).vertex_graph()
            sage: s4.is_eulerian()
            True
        """
        try:
            return self._graph
        except AttributeError:
            self._graph = Graph(self.vertex_adjacency_matrix(), loops=True)
            return self._graph


    graph = vertex_graph


    def polar(self):
        """
        Return the polar (dual) polytope.  The original vertices are
        translated so that their barycenter is at the origin, and then
        the vertices are used as the coefficients in the polar inequalities.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,1],[0,1,0],[1,0,0],[0,0,0],[1,1,1]])
            sage: p
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 5 vertices
            sage: p.polar()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 6 vertices
        """
        assert self.is_compact(), "Not a polytope."

        verts = [list(v() - self.center()) for v in self.vertex_generator()]
        return Polyhedron(ieqs=[[1] + list(v) for v in verts],
                          base_ring=self.base_ring())


    def pyramid(self):
        """
        Returns a polyhedron that is a pyramid over the original.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: egyptian_pyramid = square.pyramid()
            sage: egyptian_pyramid.n_vertices()
            5
            sage: for v in egyptian_pyramid.vertex_generator(): print v
            A vertex at (0, -1, -1)
            A vertex at (0, -1, 1)
            A vertex at (0, 1, -1)
            A vertex at (0, 1, 1)
            A vertex at (1, 0, 0)
        """
        new_verts = \
            [[0] + list(x) for x in self.Vrep_generator()] + \
            [[1] + list(self.center())]

        return Polyhedron(vertices = new_verts, base_ring=self.field())


    def bipyramid(self):
        """
        Return a polyhedron that is a bipyramid over the original.

        EXAMPLES::

            sage: octahedron = polytopes.cross_polytope(3)
            sage: cross_poly_4d = octahedron.bipyramid()
            sage: cross_poly_4d.n_vertices()
            8
            sage: q = [list(v) for v in cross_poly_4d.vertex_generator()]
            sage: q
            [[-1, 0, 0, 0],
             [0, -1, 0, 0],
             [0, 0, -1, 0],
             [0, 0, 0, -1],
             [0, 0, 0, 1],
             [0, 0, 1, 0],
             [0, 1, 0, 0],
             [1, 0, 0, 0]]

        Now check that bipyramids of cross-polytopes are cross-polytopes::

            sage: q2 = [list(v) for v in polytopes.cross_polytope(4).vertex_generator()]
            sage: [v in q2 for v in q]
            [True, True, True, True, True, True, True, True]
        """
        new_verts = \
            [[ 0] + list(x) for x in self.vertex_generator()] + \
            [[ 1] + list(self.center())] + \
            [[-1] + list(self.center())]
        new_rays = [[0] + r for r in self.rays()]
        new_lines = [[0] + list(l) for l in self.lines()]
        return Polyhedron(vertices=new_verts,
                          rays=new_rays, lines=new_lines, base_ring=self.field())


    def prism(self):
        """
        Return a prism of the original polyhedron.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: cube = square.prism()
            sage: cube
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices
            sage: hypercube = cube.prism()
            sage: hypercube.n_vertices()
            16
        """
        new_verts = []
        new_verts.extend( [ [0] + v for v in self.vertices()] )
        new_verts.extend( [ [1] + v for v in self.vertices()] )
        new_rays =        [ [0] + r for r in self.rays()]
        new_lines =       [ [0] + list(l) for l in self.lines()]
        return Polyhedron(vertices=new_verts,
                          rays=new_rays, lines=new_lines, base_ring=self.field())


    def projection(self):
        """
        Return a projection object.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: proj = p.projection()
            sage: proj
            The projection of a polyhedron into 3 dimensions
        """
        from plot import Projection
        self.projection = Projection(self)
        return self.projection


    def render_solid(self, **kwds):
        """
        Return a solid rendering of a 2- or 3-d polytope.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p_solid = p.render_solid(opacity = .7)
            sage: type(p_solid)
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
        """
        proj = self.projection()
        if self.ambient_dim()==3:
            return proj.render_solid_3d(**kwds)
        if self.ambient_dim()==2:
            return proj.render_fill_2d(**kwds)
        raise ValueError, "render_solid is only defined for 2 and 3 dimensional polyhedra."


    def render_wireframe(self, **kwds):
        """
        For polytopes in 2 or 3 dimensions, return the edges
        as a list of lines.

        EXAMPLES::

            sage: p = Polyhedron([[1,2,],[1,1],[0,0]])
            sage: p_wireframe = p.render_wireframe()
            sage: p_wireframe._objects
            [Line defined by 2 points, Line defined by 2 points, Line defined by 2 points]
        """
        proj = self.projection()
        if self.ambient_dim()==3:
            return proj.render_wireframe_3d(**kwds)
        if self.ambient_dim()==2:
            return proj.render_outline_2d(**kwds)
        raise ValueError, "render_wireframe is only defined for 2 and 3 dimensional polyhedra."


    def schlegel_projection(self, projection_dir = None, height = 1.1):
        """
        Returns a projection object whose transformed coordinates are
        a Schlegel projection of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: sch_proj = p.schlegel_projection()
            sage: schlegel_edge_indices = sch_proj.lines
            sage: schlegel_edges = [sch_proj.coordinates_of(x) for x in schlegel_edge_indices]
            sage: len([x for x in schlegel_edges if x[0][0] > 0])
            4
        """
        proj = self.projection()
        if projection_dir == None:
            v = self.vertices()
            f0 = (self.facial_incidences()[0])[1]
            projection_dir = [sum([v[f0[i]][j]/len(f0) for i in range(len(f0))])
                              for j in range(len(v[0]))]
        return proj.schlegel(projection_direction = projection_dir, height = height)


    def lrs_volume(self, verbose = False):
        """
        Computes the volume of a polytope.

        OUTPUT:

        The volume, cast to RDF (although lrs seems to output a
        rational value this must be an approximation in some cases).

        EXAMPLES::

            sage: polytopes.n_cube(3).lrs_volume() #optional, needs lrs package installed
            8.0
            sage: (polytopes.n_cube(3)*2).lrs_volume() #optional, needs lrs package installed
            64.0
            sage: polytopes.twenty_four_cell().lrs_volume() #optional, needs lrs package installed
            2.0

        REFERENCES:

             David Avis's lrs program.
        """
        if is_package_installed('lrs') != True:
            print 'You must install the optional lrs package ' \
                  'for this function to work'
            raise NotImplementedError

        in_str = self.cdd_Vrepresentation()
        in_str += 'volume'
        in_filename = tmp_filename()
        in_file = file(in_filename,'w')
        in_file.write(in_str)
        in_file.close()
        if verbose: print in_str

        lrs_procs = Popen(['lrs',in_filename],
                          stdin = PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = lrs_procs.communicate()
        if verbose:
            print ans
        # FIXME: check err

        for a_line in ans.splitlines():
            if 'Volume=' in a_line:
                volume = a_line.split('Volume=')[1]
                volume = RDF(QQ(volume))
                return volume

        raise ValueError, "lrs did not return a volume"


    def contains(self, point):
        """
        Test whether the polyhedron contains the given ``point``.

        See also :meth:`interior_contains` and
        :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point (an iterable).

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1,1],[1,-1],[0,0]])
            sage: P.contains( [1,0] )
            True
            sage: P.contains( P.center() )  # true for any convex set
            True

        The point need not have coordinates in the same field as the
        polyhedron::

            sage: ray = Polyhedron(vertices=[(0,0)], rays=[(1,0)], base_ring=QQ)
            sage: ray.contains([sqrt(2)/3,0])        # irrational coordinates are ok
            True
            sage: a = var('a')
            sage: ray.contains([a,0])                # a might be negative!
            False
            sage: assume(a>0)
            sage: ray.contains([a,0])
            True
            sage: ray.contains(['hello', 'kitty'])   # no common ring for coordinates
            False

        The empty polyhedron needs extra care, see trac #10238::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in QQ^0
            sage: empty.contains([])
            False
            sage: empty.contains([0])               # not a point in QQ^0
            False
            sage: full = Polyhedron(vertices=[()]); full
            A 0-dimensional polyhedron in QQ^0 defined as the convex hull of 1 vertex
            sage: full.contains([])
            True
            sage: full.contains([0])
            False
        """
        try:
            p = vector(point)
        except TypeError: # point not iterable or no common ring for elements
            if len(point)>0:
                return False
            else:
                p = vector(self.field(), [])

        if len(p)!=self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.contains(p):
                return False
        return True


    def interior_contains(self, point):
        """
        Test whether the interior of the polyhedron contains the
        given ``point``.

        See also :meth:`contains` and
        :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0,0],[1,1],[1,-1]])
            sage: P.contains( [1,0] )
            True
            sage: P.interior_contains( [1,0] )
            False

        If the polyhedron is of strictly smaller dimension than the
        ambient space, its interior is empty::

            sage: P = Polyhedron(vertices=[[0,1],[0,-1]])
            sage: P.contains( [0,0] )
            True
            sage: P.interior_contains( [0,0] )
            False

        The empty polyhedron needs extra care, see trac #10238::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in QQ^0
            sage: empty.interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError: # point not iterable or no common ring for elements
            if len(point)>0:
                return False
            else:
                p = vector(self.field(), [])

        if len(p)!=self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.interior_contains(p):
                return False
        return True


    def relative_interior_contains(self, point):
        """
        Test whether the relative interior of the polyhedron
        contains the given ``point``.

        See also :meth:`contains` and :meth:`interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (-1,0)])
            sage: P.contains( (0,0) )
            True
            sage: P.interior_contains( (0,0) )
            False
            sage: P.relative_interior_contains( (0,0) )
            True
            sage: P.relative_interior_contains( (1,0) )
            False

        The empty polyhedron needs extra care, see trac #10238::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in QQ^0
            sage: empty.relative_interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError: # point not iterable or no common ring for elements
            if len(point)>0:
                return False
            else:
                p = vector(self.field(), [])

        if len(p)!=self.ambient_dim():
            return False

        for eq in self.equation_generator():
            if not eq.contains(p):
                return False

        for ine in self.inequality_generator():
            if not ine.interior_contains(p):
                return False

        return True

    def is_simplex(self):
        r"""
        Return whether the polyhedron is a simplex.

        EXAMPLES::

            sage: Polyhedron([(0,0,0), (1,0,0), (0,1,0)]).is_simplex()
            True
            sage: polytopes.n_simplex(3).is_simplex()
            True
            sage: polytopes.n_cube(3).is_simplex()
            False
        """
        return self.is_compact() and (self.dim()+1==self.n_vertices())

    def is_lattice_polytope(self):
        r"""
        Return whether the polyhedron is a lattice polytope.

        OUTPUT:

        ``True`` if the polyhedron is compact and has only integral
        vertices, ``False`` otherwise.

        EXAMPLES::

            sage: polytopes.cross_polytope(3).is_lattice_polytope()
            True
            sage: polytopes.regular_polygon(5).is_lattice_polytope()
            False
        """
        try:
            return self._is_lattice_polytope
        except AttributeError:
            pass
        self._is_lattice_polytope = self.is_compact() and \
            all(v.is_integral() for v in self.vertex_generator())
        return self._is_lattice_polytope

    def lattice_polytope(self, envelope=False):
        r"""
        Return an encompassing lattice polytope.

        INPUT:

        - ``envelope`` -- boolean (default: ``False``). If the
          polyhedron has non-integral vertices, this option decides
          whether to return a strictly larger lattice polytope or
          raise a ``ValueError``. This option has no effect if the
          polyhedron has already integral vertices.

        OUTPUT:

        A :class:`LatticePolytope
        <sage.geometry.lattice_polytope.LatticePolytopeClass>`. If the
        polyhedron is compact and has integral vertices, the lattice
        polytope equals the polyhedron. If the polyhedron is compact
        but has at least one non-integral vertex, a strictly larger
        lattice polytope is returned.

        If the polyhedron is not compact, a ``NotImplementedError`` is
        raised.

        If the polyhedron is not integral and ``envelope=False``, a
        ``ValueError`` is raised.

        ALGORITHM:

        For each non-integral vertex, a bounding box of integral
        points is added and the convex hull of these integral points
        is returned.

        EXAMPLES:

        First, a polyhedron with integral vertices::

            sage: P = Polyhedron( vertices = [(1, 0), (0, 1), (-1, 0), (0, -1)])
            sage: lp = P.lattice_polytope(); lp
            A lattice polytope: 2-dimensional, 4 vertices.
            sage: lp.vertices()
            [-1  0  0  1]
            [ 0 -1  1  0]

        Here is a polyhedron with non-integral vertices::

            sage: P = Polyhedron( vertices = [(1/2, 1/2), (0, 1), (-1, 0), (0, -1)])
            sage: lp = P.lattice_polytope()
            Traceback (most recent call last):
            ...
            ValueError: Some vertices are not integral. You probably want
            to add the argument "envelope=True" to compute an enveloping
            lattice polytope.
            sage: lp = P.lattice_polytope(True); lp
            A lattice polytope: 2-dimensional, 5 vertices.
            sage: lp.vertices()
            [-1  0  0  1  1]
            [ 0 -1  1  0  1]
        """
        if not self.is_compact():
            raise NotImplementedError, 'Only compact lattice polytopes are allowed.'

        def nonintegral_error():
            raise ValueError, 'Some vertices are not integral. '+\
                'You probably want to add the argument '+\
                '"envelope=True" to compute an enveloping lattice polytope.'

        # try to make use of cached values, if possible
        if envelope:
            try:
                return self._lattice_polytope
            except AttributeError:
                pass
        else:
            try:
                assert self._is_lattice_polytope
                return self._lattice_polytope
            except AttributeError:
                pass
            except AssertionError:
                nonintegral_error()

        # find the integral vertices
        try:
            vertices = matrix(ZZ, self.vertices()).transpose()
            self._is_lattice_polytope = True
        except TypeError:
            self._is_lattice_polytope = False
            if envelope==False: nonintegral_error()
            vertices = []
            for v in self.vertex_generator():
                vbox = [ set([floor(x),ceil(x)]) for x in v ]
                vertices.extend( CartesianProduct(*vbox) )
            vertices = matrix(ZZ, vertices).transpose()

        # construct the (enveloping) lattice polytope
        from sage.geometry.lattice_polytope import LatticePolytope
        self._lattice_polytope = LatticePolytope(vertices)
        return self._lattice_polytope

    def _integral_points_PALP(self):
        r"""
        Return the integral points in the polyhedron using PALP.

        This method is for testing purposes and will eventually be removed.

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)])._integral_points_PALP()
            [(-1, -1), (0, 1), (1, 0), (1, 1), (0, 0)]
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)]).lattice_polytope(True).points()
            [ 0 -1 -1  0  1  1  0]
            [-1  0 -1  1  0  1  0]
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)])._integral_points_PALP()
            [(0, 1), (1, 0), (1, 1), (0, 0)]
        """
        if not self.is_compact():
            raise ValueError, 'Can only enumerate points in a compact polyhedron.'
        lp = self.lattice_polytope(True)
        # remove cached values to get accurate timings
        try:
            del lp._points
            del lp._npoints
        except AttributeError:
            pass
        if self.is_lattice_polytope():
            return lp.points().columns()
        points = filter(lambda p: self.contains(p),
                        lp.points().columns())
        return points

    @cached_method
    def bounding_box(self, integral=False):
        r"""
        Return the coordinates of a rectangular box containing the non-empty polytope.

        INPUT:

        - ``integral`` -- Boolean (default: ``False``). Whether to
          only allow integral coordinates in the bounding box.

        OUTPUT:

        A pair of tuples ``(box_min, box_max)`` where ``box_min`` are
        the coordinates of a point bounding the coordinates of the
        polytope from below and ``box_max`` bounds the coordinates
        from above.

        EXAMPLES::

            sage: Polyhedron([ (1/3,2/3), (2/3, 1/3) ]).bounding_box()
            ((1/3, 1/3), (2/3, 2/3))
            sage: Polyhedron([ (1/3,2/3), (2/3, 1/3) ]).bounding_box(integral=True)
            ((0, 0), (1, 1))
            sage: polytopes.buckyball().bounding_box()
            ((-1059/1309, -1059/1309, -1059/1309), (1059/1309, 1059/1309, 1059/1309))
        """
        box_min = []
        box_max = []
        if self.n_vertices==0:
            raise ValueError('Empty polytope is not allowed')
        for i in range(0,self.ambient_dim()):
            coords = [ v[i] for v in self.Vrep_generator() ]
            max_coord = max(coords)
            min_coord = min(coords)
            if integral:
                box_max.append(ceil(max_coord))
                box_min.append(floor(min_coord))
            else:
                box_max.append(max_coord)
                box_min.append(min_coord)
        return (tuple(box_min), tuple(box_max))

    def integral_points(self, threshold=100000):
        r"""
        Return the integral points in the polyhedron.

        Uses either the naive algorithm (iterate over a rectangular
        bounding box) or triangulation + Smith form.

        INPUT:

        - ``threshold`` -- integer (default: 100000). Use the naive
          algorith as long as the bounding box is smaller than this.

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)]).integral_points()
            ((-1, -1), (0, 0), (0, 1), (1, 0), (1, 1))

            sage: simplex = Polyhedron([(1,2,3), (2,3,7), (-2,-3,-11)])
            sage: simplex.integral_points()
            ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

        The polyhedron need not be full-dimensional::

            sage: simplex = Polyhedron([(1,2,3,5), (2,3,7,5), (-2,-3,-11,5)])
            sage: simplex.integral_points()
            ((-2, -3, -11, 5), (0, 0, -2, 5), (1, 2, 3, 5), (2, 3, 7, 5))

            sage: point = Polyhedron([(2,3,7)])
            sage: point.integral_points()
            ((2, 3, 7),)

            sage: empty = Polyhedron()
            sage: empty.integral_points()
            ()

        Here is a simplex where the naive algorithm of running over
        all points in a rectangular bounding box no longer works fast
        enough::

            sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
            sage: simplex = Polyhedron(v); simplex
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 5 vertices
            sage: len(simplex.integral_points())
            49

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ...        (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = Polyhedron(v)
            sage: pts1 = P.integral_points()                     # Sage's own code
            sage: all(P.contains(p) for p in pts1)
            True
            sage: pts2 = LatticePolytope(v).points().columns()   # PALP
            sage: for p in pts1: p.set_immutable()
            sage: for p in pts2: p.set_immutable()
            sage: set(pts1) == set(pts2)
            True

            sage: timeit('Polyhedron(v).integral_points()')   # random output
            sage: timeit('LatticePolytope(v).points()')       # random output
        """
        if not self.is_compact():
            raise ValueError('Can only enumerate points in a compact polyhedron.')
        if self.n_vertices() == 0:
            return tuple()

        # for small bounding boxes, it is faster to naively iterate over the points of the box
        box_min, box_max = self.bounding_box(integral=True)
        box_points = prod(max_coord-min_coord+1 for min_coord, max_coord in zip(box_min, box_max))
        if  not self.is_lattice_polytope() or \
                (self.is_simplex() and box_points<1000) or \
                box_points<threshold:
            from sage.geometry.integral_points import rectangular_box_points
            return rectangular_box_points(box_min, box_max, self)

        # for more complicate polytopes, triangulate & use smith normal form
        from sage.geometry.integral_points import simplex_points
        if self.is_simplex():
            return simplex_points(self.Vrepresentation())
        triangulation = self.triangulate()
        points = set()
        for simplex in triangulation:
            triang_vertices = [ self.Vrepresentation(i) for i in simplex ]
            new_points = simplex_points(triang_vertices)
            for p in new_points:
                p.set_immutable()
            points.update(new_points)
        # assert all(self.contains(p) for p in points)   # slow
        return tuple(points)

    def combinatorial_automorphism_group(self):
        """
        Computes the combinatorial automorphism group of the vertex
        graph of the polyhedron.

        OUTPUT:

        A
        :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        that is isomorphic to the combinatorial automorphism group is
        returned.

        Note that in Sage, permutation groups always act on positive
        integers while ``self.Vrepresentation()`` is indexed by
        nonnegative integers. The indexing of the permutation group is
        chosen to be shifted by ``+1``. That is, ``i`` in the
        permutation group corresponds to the V-representation object
        ``self.Vrepresentation(i-1)``.

        EXAMPLES::

            sage: quadrangle = Polyhedron(vertices=[(0,0),(1,0),(0,1),(2,3)])
            sage: quadrangle.combinatorial_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)(3,4)]
            sage: quadrangle.restricted_automorphism_group()
            Permutation Group with generators [()]

        Permutations can only exchange vertices with vertices, rays
        with rays, and lines with lines::

            sage: P = Polyhedron(vertices=[(1,0,0), (1,1,0)], rays=[(1,0,0)], lines=[(0,0,1)])
            sage: P.combinatorial_automorphism_group()
            Permutation Group with generators [(3,4)]
        """
        if '_combinatorial_automorphism_group' in self.__dict__:
            return self._combinatorial_automorphism_group

        from sage.groups.perm_gps.permgroup import PermutationGroup

        G = Graph()
        for edge in self.vertex_graph().edges():
            i = edge[0]
            j = edge[1]
            G.add_edge(i, j, (self.Vrepresentation(i).type(), self.Vrepresentation(j).type()) )

        group, node_dict = G.automorphism_group(edge_labels=True, translation=True)

        # Relabel the permutation group
        perm_to_vertex = dict( (i,v+1) for v,i in node_dict.items() )
        group = PermutationGroup([ [ tuple([ perm_to_vertex[i] for i in cycle ])
                                     for cycle in generator.cycle_tuples() ]
                                   for generator in group.gens() ])

        self._combinatorial_automorphism_group = group
        return group


    def _affine_coordinates(self, Vrep_object):
        r"""
        Return affine coordinates for a V-representation object.

        INPUT:

        - ``v`` -- a V-representation object or any iterable
          containing ``self.ambient_dim()`` coordinates. The
          coordinates must specify a point in the affine plane
          containing the polyhedron, or the output will be invalid. No
          checks on the input are performed.

        OUTPUT:

        A ``self.dim()``-dimensional coordinate vector. It contains
        the coordinates of ``v`` in an arbitrary but fixed basis for
        the affine span of the polyhedron.

        EXAMPLES::

            sage: P = Polyhedron(rays=[(1,0,0),(0,1,0)])
            sage: P._affine_coordinates( (-1,-2,-3) )
            (-1, -2)
            sage: [ P._affine_coordinates(v) for v in P.Vrep_generator() ]
            [(0, 0), (0, 1), (1, 0)]
        """
        if '_affine_coordinates_pivots' not in self.__dict__:
            v_list = [ vector(v) for v in self.Vrepresentation() ]
            if len(v_list)>0:
                origin = v_list[0]
                v_list = [ v - origin for v in v_list ]
            coordinates = matrix(v_list)
            self._affine_coordinates_pivots = coordinates.pivots()

        v = list(Vrep_object)
        if len(v) != self.ambient_dim():
            raise ValueError('Incorrect dimension: '+str(v))

        return vector(self.field(), [ v[i] for i in self._affine_coordinates_pivots ])


    def restricted_automorphism_group(self):
        r"""
        Return the restricted automorphism group.

        First, let the linear automorphism group be the subgroup of
        the Euclidean group `E(d) = GL(d,\RR) \ltimes \RR^d`
        preserving the `d`-dimensional polyhedron. The Euclidean group
        acts in the usual way `\vec{x}\mapsto A\vec{x}+b` on the
        ambient space.

        The restricted automorphism group is the subgroup of the linear
        automorphism group generated by permutations of the generators
        of the same type. That is, vertices can only be permuted with
        vertices, ray generators with ray generators, and line
        generators with line generators.

        For example, take the first quadrant

        .. MATH::

            Q = \Big\{ (x,y) \Big| x\geq 0,\; y\geq0 \Big\}
            \subset \QQ^2

        Then the linear automorphism group is

        .. MATH::

            \mathrm{Aut}(Q) =
            \left\{
            \begin{pmatrix}
            a & 0 \\ 0 & b
            \end{pmatrix}
            ,~
            \begin{pmatrix}
            0 & c \\ d & 0
            \end{pmatrix}
            :~
            a, b, c, d \in \QQ_{>0}
            \right\}
            \subset
            GL(2,\QQ)
            \subset
            E(d)

        Note that there are no translations that map the quadrant `Q`
        to itself, so the linear automorphism group is contained in
        the subgroup of rotations of the whole Euclidean group. The
        restricted automorphism group is

        .. MATH::

            \mathrm{Aut}(Q) =
            \left\{
            \begin{pmatrix}
            1 & 0 \\ 0 & 1
            \end{pmatrix}
            ,~
            \begin{pmatrix}
            0 & 1 \\ 1 & 0
            \end{pmatrix}
            \right\}
            \simeq \ZZ_2

        OUTPUT:

        A :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        that is isomorphic to the restricted automorphism group is
        returned.

        Note that in Sage, permutation groups always act on positive
        integers while ``self.Vrepresentation()`` is indexed by
        nonnegative integers. The indexing of the permutation group is
        chosen to be shifted by ``+1``. That is, ``i`` in the
        permutation group corresponds to the V-representation object
        ``self.Vrepresentation(i-1)``.

        REFERENCES:

        ..  [BSS]
            David Bremner, Mathieu Dutour Sikiric, Achill Schuermann:
            Polyhedral representation conversion up to symmetries.
            http://arxiv.org/abs/math/0702239

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3)
            sage: AutP = P.restricted_automorphism_group();  AutP
            Permutation Group with generators [(3,4), (2,3)(4,5), (2,5), (1,2)(5,6), (1,6)]
            sage: P24 = polytopes.twenty_four_cell()
            sage: AutP24 = P24.restricted_automorphism_group()
            sage: PermutationGroup([
            ...     '(3,6)(4,7)(10,11)(14,15)(18,21)(19,22)',
            ...     '(2,3)(7,8)(11,12)(13,14)(17,18)(22,23)',
            ...     '(2,5)(3,10)(6,11)(8,17)(9,13)(12,16)(14,19)(15,22)(20,23)',
            ...     '(2,10)(3,5)(6,12)(7,18)(9,14)(11,16)(13,19)(15,23)(20,22)',
            ...     '(2,11)(3,12)(4,21)(5,6)(9,15)(10,16)(13,22)(14,23)(19,20)',
            ...     '(1,2)(3,4)(6,7)(8,9)(12,13)(16,17)(18,19)(21,22)(23,24)',
            ...     '(1,24)(2,13)(3,14)(5,9)(6,15)(10,19)(11,22)(12,23)(16,20)'
            ...   ]) == AutP24
            True

        Here is the quadrant example mentioned in the beginning::

            sage: P = Polyhedron(rays=[(1,0),(0,1)])
            sage: P.Vrepresentation()
            (A vertex at (0, 0), A ray in the direction (0, 1), A ray in the direction (1, 0))
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(2,3)]

        Also, the polyhedron need not be full-dimensional::

            sage: P = Polyhedron(vertices=[(1,2,3,4,5),(7,8,9,10,11)])
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(1,2)]

        Translations do not change the restricted automorphism
        group. For example, any non-degenerate triangle has the
        dihedral group with 6 elements, `D_6`, as its automorphism
        group::

            sage: initial_points = [vector([1,0]), vector([0,1]), vector([-2,-1])]
            sage: points = initial_points
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - initial_points[0] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - initial_points[1] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - 2*initial_points[1] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]

        Floating-point computations are supported with a simple fuzzy
        zero implementation::

            sage: P = Polyhedron(vertices=[(1.0/3.0,0,0),(0,1.0/3.0,0),(0,0,1.0/3.0)], base_ring=RDF)
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]

        TESTS::

            sage: p = Polyhedron(vertices=[(1,0), (1,1)], rays=[(1,0)])
            sage: p.restricted_automorphism_group()
            Permutation Group with generators [(2,3)]
        """
        if '_restricted_automorphism_group' in self.__dict__:
            return self._restricted_automorphism_group

        from sage.groups.perm_gps.permgroup import PermutationGroup

        if self.field() is QQ:
            def rational_approximation(c):
                return c

        else:  # self.field() is RDF
            c_list = []
            def rational_approximation(c):
                # Implementation detail: Return unique integer if two
                # c-values are the same up to machine precision. But
                # you can think of it as a uniquely-chosen rational
                # approximation.
                for i,x in enumerate(c_list):
                    if self._is_zero(x-c):
                        return i
                c_list.append(c)
                return len(c_list)-1

        # The algorithm identifies the restricted automorphism group
        # with the automorphism group of a edge-colored graph. The
        # nodes of the graph are the V-representation objects. If all
        # V-representation objects are vertices, the edges are
        # labelled by numbers (to be computed below). Roughly
        # speaking, the edge label is the inner product of the
        # coordinate vectors with some orthogonalization thrown in
        # [BSS].
        def edge_label_compact(i,j,c_ij):
            return c_ij

        # In the non-compact case we also label the edges by the type
        # of the V-representation object. This ensures that vertices,
        # rays, and lines are only permuted amongst themselves.
        def edge_label_noncompact(i,j,c_ij):
            return (self.Vrepresentation(i).type(), c_ij, self.Vrepresentation(j).type())

        if self.is_compact():
            edge_label = edge_label_compact
        else:
            edge_label = edge_label_noncompact

        # good coordinates for the V-representation objects
        v_list = []
        for v in self.Vrepresentation():
            v_coords = list(self._affine_coordinates(v))
            if v.is_vertex():
                v_coords = [1]+v_coords
            else:
                v_coords = [0]+v_coords
            v_list.append(vector(v_coords))

        # Finally, construct the graph
        Qinv = sum( v.column() * v.row() for v in v_list ).inverse()
        G = Graph()
        for i in range(0,len(v_list)):
            for j in range(i+1,len(v_list)):
                v_i = v_list[i]
                v_j = v_list[j]
                c_ij = rational_approximation( v_i * Qinv * v_j )
                G.add_edge(i,j, edge_label(i,j,c_ij))

        group, node_dict = G.automorphism_group(edge_labels=True, translation=True)

        # Relabel the permutation group
        perm_to_vertex = dict( (i,v+1) for v,i in node_dict.items() )
        group = PermutationGroup([ [ tuple([ perm_to_vertex[i] for i in cycle ])
                                     for cycle in generator.cycle_tuples() ]
                                   for generator in group.gens() ])

        self._restricted_automorphism_group = group
        return group







#########################################################################
class PolyhedronFace_base(SageObject):
    r"""
    A face of a polyhedron.

    This class is for use in
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_lattice`.

    INPUT:

    No checking is performed whether the H/V-representation indices
    actually determine a face of the polyhedron. You should not
    manually create :class:`PolyhedronFace_base` objects unless you know
    what you are doing.

    OUTPUT:

    A :class:`PolyhedronFace_base`.

    EXAMPLES::

        sage: octahedron = polytopes.cross_polytope(3)
        sage: inequality = octahedron.Hrepresentation(2)
        sage: face_h = tuple([ inequality ])
        sage: face_v = tuple( inequality.incident() )
        sage: face_h_indices = [ h.index() for h in face_h ]
        sage: face_v_indices = [ v.index() for v in face_v ]
        sage: from sage.geometry.polyhedron.base import PolyhedronFace_base
        sage: face = PolyhedronFace_base(octahedron, face_v_indices, face_h_indices)
        sage: face
        <0,1,2>
        sage: face.dim()
        2
        sage: face.ambient_Hrepresentation()
        (An inequality (1, 1, 1) x + 1 >= 0,)
        sage: face.ambient_Vrepresentation()
        (A vertex at (-1, 0, 0), A vertex at (0, -1, 0), A vertex at (0, 0, -1))
    """

    def __init__(self, polyhedron, V_indices, H_indices):
        r"""
        The constructor.

        See :class:`PolyhedronFace_base` for more information.

        INPUT:

        - ``polyhedron`` -- a :class:`Polyhedron`. The ambient
          polyhedron.

        - ``V_indices`` -- list of integers. The indices of the
          face-spanning V-representation objects in the ambient
          polyhedron.

        - ``H_indices`` -- list of integers. The indices of the
          H-representation objects of the ambient polyhedron that are
          saturated on the face.

        TESTS::

            sage: from sage.geometry.polyhedron.base import PolyhedronFace_base
            sage: PolyhedronFace_base(Polyhedron(), [], [])   # indirect doctest
            <>
        """
        self._polyhedron = polyhedron
        self._ambient_Vrepresentation_indices = tuple(V_indices)
        self._ambient_Hrepresentation_indices = tuple(H_indices)
        self._ambient_Vrepresentation = tuple( polyhedron.Vrepresentation(i) for i in V_indices )
        self._ambient_Hrepresentation = tuple( polyhedron.Hrepresentation(i) for i in H_indices )
        # self._Vrepresentation =
        # self._Hrepresentation =


    def ambient_Hrepresentation(self, index=None):
        r"""
        Return the H-representation objects of the ambient polytope
        defining the face.

        INPUT:

        - ``index`` -- optional. Either an integer or ``None``
          (default).

        OUTPUT:

        If the optional argument is not present, a tuple of
        H-representation objects. Each entry is either an inequality
        or an equation.

        If the optional integer ``index`` is specified, the
        ``index``-th element of the tuple is returned.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: for fl in square.face_lattice():
            ...       print fl.element.ambient_Hrepresentation()
            ...
            (An inequality (1, 0) x + 1 >= 0, An inequality (0, 1) x + 1 >= 0,
             An inequality (-1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0)
            (An inequality (1, 0) x + 1 >= 0, An inequality (0, 1) x + 1 >= 0)
            (An inequality (1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0)
            (An inequality (0, 1) x + 1 >= 0, An inequality (-1, 0) x + 1 >= 0)
            (An inequality (-1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0)
            (An inequality (1, 0) x + 1 >= 0,)
            (An inequality (0, 1) x + 1 >= 0,)
            (An inequality (-1, 0) x + 1 >= 0,)
            (An inequality (0, -1) x + 1 >= 0,)
            ()
        """
        if index==None:
            return self._ambient_Hrepresentation
        else:
            return self._ambient_Hrepresentation[index]


    def ambient_Vrepresentation(self, index=None):
        r"""
        Return the V-representation objects of the ambient polytope
        defining the face.

        INPUT:

        - ``index`` -- optional. Either an integer or ``None``
          (default).

        OUTPUT:

        If the optional argument is not present, a tuple of
        V-representation objects. Each entry is either a vertex, a
        ray, or a line.

        If the optional integer ``index`` is specified, the
        ``index``-th element of the tuple is returned.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: for fl in square.face_lattice():
            ...       print fl.element.ambient_Vrepresentation()
            ...
            ()
            (A vertex at (-1, -1),)
            (A vertex at (-1, 1),)
            (A vertex at (1, -1),)
            (A vertex at (1, 1),)
            (A vertex at (-1, -1), A vertex at (-1, 1))
            (A vertex at (-1, -1), A vertex at (1, -1))
            (A vertex at (1, -1), A vertex at (1, 1))
            (A vertex at (-1, 1), A vertex at (1, 1))
            (A vertex at (-1, -1), A vertex at (-1, 1),
             A vertex at (1, -1), A vertex at (1, 1))
        """
        if index==None:
            return self._ambient_Vrepresentation
        else:
            return self._ambient_Vrepresentation[index]


    def n_ambient_Hrepresentation(self):
        """
        Return the number of objects that make up the ambient
        H-representation of the polyhedron.

        See also :meth:`ambient_Hrepresentation`.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: face = p.face_lattice()[10].element
            sage: face
            <0,2>
            sage: face.ambient_Hrepresentation()
            (An inequality (1, -1, 1, -1) x + 1 >= 0,
             An inequality (1, 1, 1, 1) x + 1 >= 0,
             An inequality (1, 1, 1, -1) x + 1 >= 0,
             An inequality (1, -1, 1, 1) x + 1 >= 0)
            sage: face.n_ambient_Hrepresentation()
            4
        """
        return len(self._ambient_Hrepresentation)


    def n_ambient_Vrepresentation(self):
        """
        Return the number of objects that make up the ambient
        V-representation of the polyhedron.

        See also :meth:`ambient_Vrepresentation`.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: face = p.face_lattice()[10].element
            sage: face
            <0,2>
            sage: face.ambient_Vrepresentation()
            (A vertex at (-1, 0, 0, 0), A vertex at (0, 0, -1, 0))
            sage: face.n_ambient_Vrepresentation()
            2
        """
        return len(self._ambient_Vrepresentation)


    def dim(self):
        """
        Return the dimension of the face.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: fl = polytopes.dodecahedron().face_lattice()
            sage: [ x.element.dim() for x in fl ]
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3]
        """
        if '_dim' in self.__dict__:
            return self._dim

        if self.n_ambient_Vrepresentation()==0:
            self._dim = -1
        else:
            origin = vector(self.ambient_Vrepresentation(0))
            v_list = [ vector(v)-origin for v in self.ambient_Vrepresentation() ]
            self._dim = matrix(v_list).rank()
        return self._dim


    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:

        A string listing the V-representation indices of the face.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: a_face = list( square.face_lattice() )[8].element
            sage: a_face.__repr__()
            '<1,3>'
        """
        s = '<'
        s += ','.join([ str(v.index()) for v in self.ambient_Vrepresentation() ])
        s += '>'
        return s

