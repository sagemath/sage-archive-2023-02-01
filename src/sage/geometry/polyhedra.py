"""
Polyhedra

This module contains functions for computing with exact (rational)
or floating-point convex polyhedra in arbitrary dimensions. The bulk
of this functionality is provided through the cddlib library of Komei
Fukuda.

There seems to be some inconsistency in the use of the word polyhedra.
In this module, a polyhedron is a convex (possibly unbounded) set
defined either as
- the intersection of a finite set of half-planes and hyperplanes (H-form),
- the minkowsky sum of the convex hull of a finite set of vertices and
  the positive (aka conic) combination of a finite set of rays (V-form)

The half-planes are also referred to as inequalities, and abbreviated
ieq, and the hyperplanes are also referred to as linearities.

A polytope is defined as a bounded polyhedron.

More information about polyhedra, polytope and triangulation, consistent
with the definitions chosen here can be found in Komei Fukuda's `FAQ in
Polyhedral Computation`_.

Although the ability to specify polyhedra by inequalities and
linearities is supported, much of the code was designed for the case
of polytopes specified by vertices (with rational coordinates).

.. _`FAQ in Polyhedral Computation`: http://www.ifor.math.ethz.ch/~fukuda/polyfaq/polyfaq.html

AUTHOR:
    -- Marshall Hampton: first version and bugfixes, 2008
    -- Arnaud Bergeron: improvements to triangulation and rendering, 2008
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.structure.sage_object import SageObject
from sage.misc.all import SAGE_TMP, tmp_filename
from sage.misc.functional import norm
from random import randint
from sage.rings.all import Integer, QQ, ZZ
from sage.rings.rational import Rational
from sage.rings.real_mpfr import RR
from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.plot.plot3d.shapes2 import point3d
from sage.plot.all import line
from sage.combinat.combinat import permutations
from sage.groups.perm_gps.permgroup_named import AlternatingGroup
from sage.misc.package import is_package_installed
from sage.graphs.graph import Graph

import os
from subprocess import Popen, PIPE

def mink_sum(polyhedra_list, verbose = False):
    """
    Returns the Minkowski sum of a list of polyhedra.
    TODO: currently only works for polytopes, not polyhedra with rays

    EXAMPLES:
        sage: poly_spam = Polyhedron([[3,4,5,2],[1,0,0,1],[0,0,0,0],[0,4,3,2],[-3,-3,-3,-3]])
        sage: poly_eggs = Polyhedron([[5,4,5,4],[-4,5,-4,5],[4,-5,4,-5],[0,0,0,0]])
        sage: poly_spam_and_eggs = mink_sum([poly_spam,poly_spam,poly_eggs])
        sage: poly_spam_and_eggs.n_vertices()
        12
    """
    answer = Polyhedron()
    temp_vertex_list = []
    for vertex1 in polyhedra_list[0].vertices():
        for vertex2 in polyhedra_list[1].vertices():
            temp_vertex_list.append([vertex1[i]+vertex2[i] for i in range(len(vertex1))])
    answer._vertices = temp_vertex_list
    answer._remove_redundant_vertices()
    if verbose:
	print 'Another pair of polytopes computed, sum has ' + str(len(answer._vertices)) + ' vertices'
    if len(polyhedra_list) == 2:
        return answer
    else:
        return mink_sum([answer] + polyhedra_list[2:], verbose = verbose)

class Polyhedron(SageObject):

    def __init__(self, vertices = [], rays = [], ieqs = [], linearities = [], cdd_type = 'rational', verbose = False):
        """
        A class for polyhedral objects.  Architecture is still in flux.
        Do not instantiate with both vertex/ray and ieq/linearity data:
        no consistency check is made.

        EXAMPLES:
            sage: square_from_vertices = Polyhedron(vertices = [[1, 1], [1, -1], [-1, 1], [-1, -1]])
            sage: square_from_ieqs = Polyhedron(ieqs = [[1, 0, 1], [1, 1, 0], [1, 0, -1], [1, -1, 0]])
            sage: square_from_ieqs.vertices()
            [[1, -1], [1, 1], [-1, 1], [-1, -1]]
            sage: square_from_vertices.ieqs()
            [[1, 0, 1], [1, 1, 0], [1, 0, -1], [1, -1, 0]]
	    sage: p = Polyhedron(vertices = [[1.1, 2.2], [3.3, 4.4]], cdd_type = 'real')
	    sage: len(p.ieqs())
	    2

        NOTE:
             Although the option cdd_type = 'real' allows numerical data to be used, its use is
        discouraged - the results can depend upon the tolerance setting of cddlib, whose
        default is rather large since cddlib was designed with rational and integer data in mind.
        """
	self._cdd_type = cdd_type
        self._vertices = vertices
        self._ieqs = ieqs
        self._linearities = linearities
        self._rays = rays
        self._verbose = verbose
	if self._vertices != []:
	    self._remove_redundant_vertices(cdd_type = self._cdd_type)
        else:
            if self._ieqs != []:
                vertices = [x for x in self.vertices()]
                temp = Polyhedron(vertices = vertices)
                self._ieqs = temp.ieqs()


    def _repr_(self):
        """
        Returns a description of the object.

        EXAMPLES:
            sage: poly_test = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1]])
            sage: poly_test._repr_()
            'A Polyhedron with 3 vertices.'
            sage: grammar_test = Polyhedron(vertices = [[1,1,1,1,1,1]])
            sage: grammar_test._repr_()
            'A Polyhedron with 1 vertex.'
        """
        desc = "A Polyhedron"
        # ADD .n_vertices(), etc!!
        if self._vertices != []:
            if len(self._vertices) > 1:
                desc += ' with ' + str(len(self._vertices)) + ' vertices'
            else:
                desc += ' with ' + str(len(self._vertices)) + ' vertex'
        if self._rays != []:
            desc += ' with ' + str(len(self._rays)) + ' rays'
        if self._ieqs != []:
            if len(self._ieqs) > 1:
                desc += ' with ' + str(len(self._ieqs)) + ' supporting hyperplanes'
            else:
                desc += ' with ' + str(len(self._ieqs)) + ' supporting hyperplane'
        if self._linearities != []:
            if len(self._linearities) > 1:
                desc += ' with ' + str(len(self._linearities)) + ' linear constraints'
            else:
                desc += ' with ' + str(len(self._linearities)) + ' linear constraint'
        desc += '.'
        return desc

    def show(self, fill = False, **kwds):
        """
        Returns a wireframe or solid rendering if possible.
        Only works currently in dimensions 2,3, and 4.

        EXAMPLES:
            sage: p = polytopes.n_cube(4)
            sage: p_show = p.show()
            sage: p_show.bounding_box()
            ((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
        """
        if self.ambient_dim() == 2 or self.ambient_dim() == 4:
            return self.render_wireframe()
        if self.ambient_dim() == 3:
            if fill == False:
                return self.render_wireframe()
            else:
                return self.render_solid()
        else:
            raise TypeError,"only polytopes in dimensions 2 to 4 are shown"

    def dim(self):
        """
        Returns the dimension of the polyhedron.

        EXAMPLES:
            sage: standard_simplex = Polyhedron(vertices = [[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]])
            sage: standard_simplex.dim()
            3
        """
        try:
            return self._dim
        except AttributeError:
            if self.vertices() != []:
                self._dim = len(self.vertices()[0]) - len(self.linearities())
            elif self.ieqs() != []:
                self._dim = len(self.ieqs()[0]) - 1 - len(self.linearities())
            else:
                self._dim = 0
            return self._dim

    def ambient_dim(self):
        """
        Returns the ambient dimension of the polyhedron.

        EXAMPLES:
            sage: standard_simplex = Polyhedron(vertices = [[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]])
            sage: standard_simplex.ambient_dim()
            4
        """
        try:
            return self._ambient_dim
        except AttributeError:
            if self.vertices() != []:
                self._ambient_dim = len(self.vertices()[0])
            elif self.ieqs() != []:
                self._ambient_dim = len(self.ieqs()[0]) - 1
            else:
                self._ambient_dim = 0
            return self._ambient_dim

    def is_simple(self):
        """
        Tests for simplicity of a polytope.

        EXAMPLES:
            sage: p = Polyhedron([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
            sage: p._remove_redundant_vertices()
            sage: p.is_simple()
            True
            sage: p = Polyhedron([[0,0,0],[4,4,0],[4,0,0],[0,4,0],[2,2,2]])
            sage: p.is_simple()
            False
        """
        self._remove_redundant_vertices()
        d = self.dim()
        for edge_data in self.vertex_adjacencies():
            if len(edge_data[1]) != d:
                return False
        return True

    def __add__(self, other):
        """
        Addition of two polyhedra is defined as their Minkowski sum.

        EXAMPLES:
            sage: four_cube = polytopes.n_cube(4)
            sage: four_simplex = Polyhedron(vertices = [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: unholy_union = four_cube + four_simplex
            sage: unholy_union.dim()
            4
        """
        if isinstance(other,Polyhedron):
            return mink_sum([self,other])
        else:
            # assume it is a vector and try to add vertices
            try:
                if len(other) != self.dim():
                    raise TypeError, "other (=%s) must be the same dimension"%other
                new_vertices = [list(vector(x) + vector(other)) for x in self.vertices()]
                return Polyhedron(vertices = new_vertices)

            except:
                raise TypeError, "other (=%s) must be a Polyhedron or a vector of dimension "%other + str(self.dim())

    def __mul__(self, other):
        """
        Multiplication by another polyhedron returns the product polytope.
        Multiplication by a scalar returns the polytope dilated by that scalar.

        EXAMPLES:
            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(4)])
            sage: p2 = p*2
            sage: p2.vertices()[1]
            [2, 2, 2]
        """
        if isinstance(other,Polyhedron):
            new_vertices = [x+y for x in self.vertices() for y in other.vertices()]
            return Polyhedron(vertices = new_vertices)
        try:
            new_vertices = [[other*vi for vi in vertex] for vertex in self.vertices()]
            return Polyhedron(vertices = new_vertices)
        except:
            raise TypeError, "other (=%s) must be a Polyhedron or a scalar"%other

    def __rmul__(self,other):
        """
        Multiplication by another polyhedron returns the product polytope.
        Multiplication by a scalar returns the polytope dilated by that scalar.

        EXAMPLES:
            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(4)])
            sage: p2 = 3*p + p
            sage: p2.vertices()[1]
            [4, 4, 4]
        """
        return self.__mul__(other)


    def union(self, other):
        """
        Returns a polyhedron whose vertices are the union of the vertices
        of the two polyhedra.

        EXAMPLES:
            sage: a_simplex = polytopes.n_simplex(3)
            sage: verts = [x for x in a_simplex.vertices()]
            sage: verts = [[x[0]*3/5+x[1]*4/5, -x[0]*4/5+x[1]*3/5, x[2]] for x in verts]
            sage: another_simplex = Polyhedron(vertices = verts)
            sage: simplex_union = a_simplex.union(another_simplex)
            sage: simplex_union.n_vertices()
            7
        """
        return Polyhedron(vertices = self.vertices() + other.vertices())


    def intersection(self, other):
        """
        Returns the intersection of one polyhedron with another.

        EXAMPLES:
            sage: cube = polytopes.n_cube(3)
            sage: oct = polytopes.cross_polytope(3)
            sage: cube_oct = cube.intersection(oct*2)
            sage: cube_oct
	    A Polyhedron with 12 vertices with 14 supporting hyperplanes.
            sage: len(cube_oct.vertices())
            12
            sage: cube_oct
            A Polyhedron with 12 vertices with 14 supporting hyperplanes.
        """
        new_ieqs = self.ieqs() + other.ieqs()
        return Polyhedron(ieqs = new_ieqs)


    def vertices(self, force_from_ieqs = False):
        """
        Returns the vertices of the polyhedron.

        EXAMPLES:
            sage: a_simplex = Polyhedron(ieqs = [[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]], linearities = [[1,-1,-1,-1,-1]])
            sage: a_simplex.vertices()
            [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]
        """
        if (self._vertices == [] and self._ieqs != []) or force_from_ieqs:
            temp_poly_object = ieq_to_vert(self._ieqs, linearities = self._linearities, verbose = self._verbose)
            self._vertices = temp_poly_object._vertices
            self._vertex_incidences = temp_poly_object._vertex_incidences
            self._vertex_adjacencies = temp_poly_object._vertex_adjacencies
            self._rays = temp_poly_object._rays
            self._linearities = temp_poly_object._linearities
            self._linearities_checked = True
            self._checked_rays = True
            self._n_vertices = len(self._vertices)
            return self._vertices
        return self._vertices

    def n_vertices(self):
        """
        Returns the number of vertices in the polyhedron.

        EXAMPLES:
            sage: p = Polyhedron(vertices = [[t,t^2,t^3,randint(0,10)] for t in range(6)])
            sage: p.n_vertices()
            6
        """
        try:
            return self._n_vertices
        except AttributeError:
            self._n_vertices = len(self.vertices())
            return self._n_vertices

    def n_facets(self):
        """
        Returns the number of facets in the polyhedron.

        EXAMPLES:
            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)])
            sage: p.n_facets()
            8
        """
        try:
            return self._n_facets
        except AttributeError:
            self._n_facets = len(self.ieqs())
            return self._n_facets

    def ieqs(self, force_from_vertices = False):
        """
        Return the inequalities (half-planes) defining the polyhedron.

        EXAMPLES:
            sage: permuta4 = Polyhedron(vertices = permutations([1,2,3,4,5]))
            sage: ieqs = permuta4.ieqs()
            sage: ieqs[0]
            [-1, 0, 0, 0, 1, 0]
            sage: ieqs[-1]
            [5, -1, 0, 0, 0, 0]
        """
        if self._vertices != []:
            temp_poly_object = vert_to_ieq(self._vertices, rays = self._rays, cdd_type = self._cdd_type)
            self._linearities = temp_poly_object._linearity
            self._ieqs = temp_poly_object._inequalities
            self._facial_incidences = temp_poly_object._facial_incidences
            self._facial_adjacencies = temp_poly_object._facial_adjacencies
            return self._ieqs
        elif force_from_vertices:
            if self._vertices == [] and self.ieqs() != []:
                temp_poly_object = ieq_to_vert(self._ieqs, linearities = self._linearities, cdd_type = self._cdd_type)
                self._vertices = temp_poly_object.vertices()
                self._linearities = temp_poly_object._linearities
                self._facial_incidences = temp_poly_object._facial_incidences
                self._facial_adjacencies = temp_poly_object._facial_adjacencies
                return self._linearities
        else:
            return self._ieqs

    def vertex_adjacencies(self):
        """
        Returns a list of vertex indices and their adjacent vertices.

        EXAMPLES:
            sage: permuta3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: permuta3.vertex_adjacencies()[0:3]
            [[3, [5, 1, 9]], [5, [3, 4, 11]], [16, [22, 10, 17]]]
        """
        try:
            return self._vertex_adjacencies
        except AttributeError:
            if self._ieqs != []:
                temp_poly_object = ieq_to_vert(self._ieqs, linearities = self._linearities)
                if self._vertices == []:
                    self._vertices = temp_poly_object._vertices
                    self._vertex_incidences = temp_poly_object._vertex_incidences
                    self._vertex_adjacencies = temp_poly_object._vertex_adjacencies
                    return self._vertex_adjacencies
                # avoid changing vertex indexing
                vdict = dict([[temp_poly_object._vertices.index(self._vertices[ind]), ind] for ind in range(len(self._vertices))])
                self._vertex_incidences = [[vdict[vi[0]],vi[1]] for vi in temp_poly_object._vertex_incidences]
                self._vertex_adjacencies = [[vdict[vi[0]],[vdict[q] for q in vi[1]]] for vi in temp_poly_object._vertex_adjacencies]
                return self._vertex_adjacencies
            else:
                if self._vertices == []:
                    return []
                else:
                    dummy = self.ieqs()
                    return self.vertex_adjacencies()

    def facial_adjacencies(self):
        """
        Returns a list of face indices and the indices of faces adjacent to them.

        EXMAPLES:
            sage: permuta3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: permuta3.facial_adjacencies()[0:3]
            [[0, [1, 2, 3, 8, 12, 13, 14]],
            [1, [0, 2, 9, 13, 14]],
            [2, [0, 1, 3, 4, 5, 9, 14]]]
        """
        try:
            return self._facial_adjacencies
        except AttributeError:
            if self._vertices != []:
                temp_poly_object = vert_to_ieq(self._vertices, rays = self._rays)
                self._ieqs = temp_poly_object._ieqs
                self._facial_incidences = temp_poly_object._facial_incidences
                self._facial_adjacencies = temp_poly_object._facial_adjacencies
                return self._facial_adjacencies
            else:
                if self._ieqs == []:
                    return []
                else:
                    dummy = self.vertices()
                    return self.facial_adjacencies()

    def linearities(self):
        """
        Returns the linear constraints of the polyhedron. As with
        inequalities, each constraint is given as [b -a1 -a2 ... an]
        where for variables x1, x2,..., xn, the polyhedron satisfies
        the equation b = a1*x1 + a2*x2 + ... + an*xn.

        EXAMPLES:
            sage: test_p = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1],[3,4,1,2]])
            sage: test_p.linearities()
            [[-10, 1, 1, 1, 1]]
        """
        try:
            if self._linearities_checked == True:
                return self._linearities
        except AttributeError:
            # need to force a check
            dummy_ieqs = self.ieqs(force_from_vertices = True)
            self._linearities_checked = True
            return self._linearities

    def _remove_redundant_vertices(self, cdd_type = 'rational'):
        """
        Removes vertices from the description of the polyhedron which
        do not affect the convex hull.  Now done automatically when a
        polyhedron is defined by vertices.

        EXAMPLES:
            sage: a_triangle = Polyhedron([[0,0,0],[4,0,0],[0,4,0],[1,1,0]])
            sage: a_triangle.vertices()
            [[0, 0, 0], [4, 0, 0], [0, 4, 0]]
            sage: a_triangle._remove_redundant_vertices()
            sage: a_triangle.vertices()
            [[0, 0, 0], [4, 0, 0], [0, 4, 0]]
        """
        try:
            if self._checked_vertices == True:
                pass
        except AttributeError:
            if self.vertices() != []:
                ex_verts = extreme_verts(self.vertices(), cdd_type = cdd_type)
            else:
                return []
            self._vertices = ex_verts
            self._checked_vertices = True
            pass

    def rays(self):
        """
        Returns the rays.  This function does not ensure that redundant
        rays have been removed.

        EXAMPLES:
            sage: positive_orthant = Polyhedron(ieqs = [[0,1,0,0],[0,0,1,0],[0,0,0,1]])
            sage: positive_orthant.rays()
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        """
        if self._rays != []:
            return self._rays
        try:
            if self._checked_rays == True:
                return self._rays
        except AttributeError:
            dummy_verts = self.vertices(force_from_ieqs = True)
            self._checked_rays = True
            return self._rays

    def facial_incidences(self):
        """
        Returns the face-vertex incidences in the form [face_index, [v_i_0, v_i_1,,v_i_2]]

        EXAMPLES:
            sage: p = Polyhedron(vertices = [[5,0,0],[0,5,0],[5,5,0],[0,0,0],[2,2,5]])
            sage: p.facial_incidences()
            [[0, [1, 3, 4]],
            [1, [0, 3, 4]],
            [2, [0, 2, 4]],
            [3, [1, 2, 4]],
            [4, [0, 1, 2, 3]]]
        """
        try:
            return self._facial_incidences
        except:
            dummy_ieqs = self.ieqs(force_from_vertices = True) # force computation of facial incidences
            return self._facial_incidences

    def vertex_adjacency_matrix(self):
        """
        Returns the binary matrix of vertex adjacencies.

        EXAMPLES:
            sage: polytopes.n_simplex(4).vertex_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]
        """
        m = []
        v_adjs = [x for x in self.vertex_adjacencies()]
        v_adjs.sort()
        for va in v_adjs:
            temp = []
            for j in range(self.n_vertices()):
                if j in va[1]:
                    temp.append(1)
                else:
                    temp.append(0)
            m.append(temp)
        return matrix(ZZ,m)

    def graph(self):
        """
        Returns a graph in which the vertices correspond to vertices of the polyhedron,
        and edges to edges.

        EXAMPLES:
            sage: g3 = polytopes.n_cube(3).graph()
            sage: g3.automorphism_group()
            Permutation Group with generators [(2,4)(3,5), (1,2)(5,6),
            (1,8)(2,3)(4,5)(6,7)]
            sage: s4 = polytopes.n_simplex(4).graph()
            sage: s4.is_eulerian()
            True
        """
        return Graph(self.vertex_adjacency_matrix(), loops=True)

    def polar(self):
        """
        Returns the polar (dual) polytope.  The original vertices are
        translated so that their barycenter is at the origin, and then
        the vertices are used as the coefficients in the polar inequalities.

        EXAMPLES:
            sage: p24 = polytopes.twenty_four_cell()
            sage: dual_24 = p24.polar()
            sage: dual_dual_24 = dual_24.polar() #dual of the dual should bring us back to the original in this case
            sage: p24.vertices()[0] in dual_dual_24.vertices()
            True
        """
        old_verts = [x for x in self.vertices()]
        old_verts_center = sum([vector(x) for x in old_verts])/self.n_vertices()
        old_verts = [vector(x) - old_verts_center for x in old_verts]
        old_verts = [list(x) for x in old_verts]
        return Polyhedron(ieqs = [[1] + x for x in old_verts])

    def pyramid(self):
        """
        Returns a polyhedron that is a pyramid over the original.

        EXAMPLES:
            sage: square = polytopes.n_cube(2)
            sage: egyptian_pyramid = square.pyramid()
            sage: egyptian_pyramid.n_vertices()
            5
            sage: egyptian_pyramid.vertices()
            [[1, 1, 0], [1, -1, 0], [-1, 1, 0], [-1, -1, 0], [0, 0, 1]]
        """
        old_verts = [x for x in self.vertices()]
        old_verts_center = sum([vector(x) for x in old_verts])/self.n_vertices()
        new_verts = [x+[0] for x in old_verts] + [list(old_verts_center) + [1]]
        return Polyhedron(vertices = new_verts)

    def prism(self):
        """
        Returns a prism of the original polyhedron.

        EXAMPLES:
            sage: square = polytopes.n_cube(2)
            sage: cube = square.prism()
            sage: cube.ieqs()
            [[0, 0, 0, 1], [1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, -1], [1, 0, -1, 0], [1, -1, 0, 0]]
            sage: hypercube = cube.prism()
            sage: hypercube.n_vertices()
            16
        """
        new_verts = [x + [0] for x in self.vertices()]
        new_verts = new_verts + [x + [1] for x in self.vertices()]
        return Polyhedron(vertices = new_verts)


    def triangulated_facial_incidences(self):
        """
        Returns a list of the form [face_index, [v_i_0, v_i_1,...,v_i_{n-1}]]
        where the face_index refers to the original defining inequality.
        For a given face, the collection of triangles formed by each
        list of v_i should triangulate that face.

        In dimensions greater than 3, this is computed by randomly lifting each
        face up a dimension; this does not always work!  Work is being done on this.

        EXAMPLES:
            If the figure is already composed of triangles, then all is well
            sage: Polyhedron(vertices = [[5,0,0],[0,5,0],[5,5,0],[2,2,5]]).triangulated_facial_incidences()
            [[0, [0, 2, 3]], [1, [0, 1, 2]], [2, [1, 2, 3]], [3, [0, 1, 3]]]

 	    Otherwise some faces get split up to triangles
 	    sage: Polyhedron(vertices = [[2,0,0],[4,1,0],[0,5,0],[5,5,0],[1,1,0],[0,0,1]]).triangulated_facial_incidences()
            [[0, [0, 1, 5]], [1, [0, 4, 5]], [2, [2, 4, 5]], [3, [2, 3, 5]], [4, [1, 3, 5]], [5, [0, 4, 1]], [5, [4, 1, 2]], [5, [1, 2, 3]]]

        """
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
                        lifted_verts.append(self.vertices()[vert_index] + [randint(-vert_index,500+vert_index + vert_number**2)])
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
                            t_fac_incs.append([a_face[0], [self.vertices().index(q) for q in proj_verts]])
                else:
                    vs = a_face[1]
 		    adj = dict([a[0], filter(lambda p: p in a_face[1], a[1])] for a in filter(lambda va: va[0] in a_face[1], self.vertex_adjacencies()))
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

    def render_wireframe(self, rgbcolor = (0,0,1)):
        """
        For polytopes in 2 or 3 dimensions, returns the edges
        as a list of lines.

        EXAMPLES:
            sage: p = Polyhedron([[1,2,],[1,1],[0,0]])
            sage: p_wireframe = p.render_wireframe()
            sage: p_wireframe._Graphics__objects
            [Line defined by 2 points, Line defined by 2 points, Line defined by 2 points]
        """
        if self.dim() > 4:
            print "Dimension is too large for wireframe"
            return NotImplementedError
        if self.dim() < 4:
            if self.ieqs() == self.vertices() == []:
                return Graphics()
            edges = []
            verts = self.vertices()
            for adj in self.vertex_adjacencies():
                for vert in adj[1]:
                    if vert > adj[0]:
                        edges.append([verts[adj[0]],verts[vert]])
            return sum([line(an_edge, rgbcolor = rgbcolor) for an_edge in edges])
        # Now we must be in 4 dimensions, so return Schlegel diagram from the first face.
        v = self.vertices()
        f0 = (self.facial_incidences()[0])[1]
        vcenter = [sum([v[f0[i]][j]/len(f0) for i in range(len(f0))]) for j in range(len(v[0]))] # compute center of face
        spcenter = [vi/norm(vector(vcenter)) for vi in vcenter] # normalize to unit sphere
        spverts = [[vi/norm(vector(vp)) for vi in vp] for vp in v] # normalize vertices to unit sphere
        polediff = matrix(RDF,vector([0.0,0.0,0.0,1.0])-vector(spcenter)).transpose()
        denom = RDF((polediff.transpose()*polediff)[0][0])
        house = matrix(RDF,[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]) - 2*polediff*polediff.transpose()/denom #Householder reflector
        spverts = [house*vector(RDF,x) for x in spverts] # reflect so face center is at "north pole"
        spverts = [list(x) for x in spverts]
        proj_verts = [[x[0]/(1-x[3]),x[1]/(1-x[3]),x[2]/(1-x[3])] for x in spverts] # stereographically project
        hlines = point3d([0,0,0], pointsize = 0.01, rgbcolor = rgbcolor) # workaround due to bug in line not giving a bounding box
        for an_edge in self.vertex_adjacencies():
            for j in an_edge[1]:
                hlines += line([proj_verts[an_edge[0]],proj_verts[j]], rgbcolor = rgbcolor)
        return hlines

    def render_solid(self, rgbcolor = (1,0,0), **kwds):
        """
        Returns solid 3d rendering of a 3d polytope.

        EXAMPLES:
            sage: p = polytopes.n_cube(3)
            sage: p_solid = p.render_solid(opacity = .7)
            sage: type(p_solid)
            <type 'sage.plot.plot3d.index_face_set.IndexFaceSet'>
        """
        if self.ambient_dim() != 3:
            raise TypeError, "render_solid currently only works on polytopes in 3D"
        tri_faces = self.triangulated_facial_incidences()
        from sage.plot.plot3d.index_face_set import IndexFaceSet
        return IndexFaceSet([q[1] for q in tri_faces], self.vertices(), enclosed = True, rgbcolor = rgbcolor, **kwds)

    def lrs_volume(self, verbose = False):
        """
        Computes the volume of a polytope using David Avis's lrs program.

        OUTPUT:
            the volume, cast to RDF (although lrs seems to output a rational
        value this must be an approximation in some cases).

        EXAMPLES:
            sage: polytopes.n_cube(3).lrs_volume() #optional, needs lrs package installed
            8.0
            sage: (polytopes.n_cube(3)*2).lrs_volume() #optional, needs lrs package installed
            64.0
            sage: polytopes.twenty_four_cell().lrs_volume() #optional, needs lrs package installed
            2.0
            """
        if is_package_installed('lrs') != True:
            print 'You must install the optional lrs package for this function to work'
            return False
        vertices = self.vertices()
        cdd_type = self._cdd_type
        v_dim = len(vertices[0])
        v_sage = type(vertices[0][0])
        in_str = 'V-representation\nbegin\n'
        in_str += str(len(vertices)) + ' ' + str(v_dim + 1) + ' ' + cdd_type + '\n'
        for vert in vertices:
            in_str += '1 '
            for numb in vert:
                in_str += str(numb) + ' '
            in_str += '\n'
        in_str += 'end\n'
        in_str += 'volume'
        in_filename = tmp_filename()
        in_file = file(in_filename,'w')
        in_file.write(in_str)
        in_file.close()
        if verbose: print in_str
        lrs_procs = Popen(['lrs',in_filename],stdin = PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = lrs_procs.communicate()
        ans_lines = ans.splitlines()
        if verbose:
            print ans
        for a_line in ans_lines:
            if 'Volume=' in a_line:
                volume = a_line.split('Volume=')[1]
                volume = RDF(QQ(volume))
                break
        return volume



def cdd_convert(a_str,a_type = Rational):
    """
    A little utility function for changing a space-delimited list into
    a python list.

    INPUT:
        a string of space-delimited numbers without newline characters

    EXAMPLES:
        sage: from sage.geometry.polyhedra import cdd_convert
        sage: cdd_convert(' 1 1 0 0')
        [1, 1, 0, 0]
    """
    temp = a_str.split(' ')
    temp = [x for x in temp if x != '']
    return [a_type(x) for x in temp]

def extreme_verts(vertices, cdd_type = 'rational', verbose = False):
    """
    Returns a list of the vertices which define the convex hull of the
    input vertices; that is, this eliminates the redundant vertices.

    INPUT:
        vertices -- a list of vertices, integer or rational or real
        cdd_type -- optional.  Defaults to rational, and the exact
            arithmetic of testcdd1_gmp.  Using real will invoke the
            floating-point program testcdd1 instead.  This is much
            faster but possibly will not give the correct results for
            highly degenerate polytopes in high dimensions.
        verbose -- optional.  Only interesting for debugging purposes,
            if True this may print out extra information.

    EXAMPLES:
        sage: extreme_verts([[1,3],[2,2],[3,1],[0,0]])
        [[1, 3], [3, 1], [0, 0]]
    """
    # first we create the input for cddlib:
    v_dim = len(vertices[0])
    v_sage = type(vertices[0][0])
    in_str = 'V-representation\nbegin\n'
    in_str += str(len(vertices)) + ' ' + str(v_dim + 1) + ' ' + cdd_type + '\n'
    for vert in vertices:
        in_str += '1 '
        for numb in vert:
            in_str += str(numb) + ' '
        in_str += '\n'
    in_str += 'end\n'
    in_filename = tmp_filename()
    in_file = file(in_filename,'w')
    in_file.write(in_str)
    in_file.close()
    if verbose: print in_str
    # now send the input to the appropriate cddlib function
    if cdd_type == 'real':
        cdd_procs = Popen('redcheck',stdin = PIPE, stdout=PIPE, stderr=PIPE)
        cast = RR
    else:
        cdd_procs = Popen('redcheck_gmp',stdin = PIPE, stdout=PIPE, stderr=PIPE)
        cast = Rational
    ans, err = cdd_procs.communicate(input = in_filename)
    # and now parse the output into Sage objects
    ans_lines = ans.splitlines()
    ex_verts = []
    if verbose: print ans
    for index in range(len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('V-representation') != -1:
            v_index = index + 3
    for index in range(v_index,len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('end') != -1: break
        ex_verts.append(cdd_convert(a_line, a_type = cast)[1:])
    return ex_verts

def vert_to_ieq(vertices, rays = [], cdd_type = 'rational', verbose = False):
    """
    Converts a list of vertices to an inequality representation, using cddlib.

    INPUT:
        vertices -- a list of vertices
        rays -- optional, a list of rays
        cdd_type -- optional.  Defaults to rational, and the exact
            arithmetic of testcdd1_gmp.  Using real will invoke the
            floating-point program testcdd1 instead.  This is much
            faster but possibly will not give the correct results for
            highly degenerate polytopes in high dimensions.
        verbose -- optional.  Only interesting for debugging purposes,
            if True this may print out extra information.

    OUTPUT:
        A Polyhedron containing information about the vertices,
        faces, incidences, and adjancencies of the polyhedron.

    EXAMPLES:
        sage: a = vert_to_ieq([[1,2,3],[1,2,2],[1,2,4]])
        sage: a._facial_adjacencies
        [[0, [1, 2, 3]], [1, [0, 2, 3]]]
        sage: a.vertices()
        [[1, 2, 2], [1, 2, 4]]
    """
    # first we create the input for cddlib:
    v_dim = len(vertices[0])
    v_sage = type(vertices[0][0])
    in_str = 'V-representation\nbegin\n'
    in_str += str(len(vertices)+len(rays)) + ' ' + str(v_dim + 1) + ' ' + cdd_type + '\n'
    for ray in rays:
        in_str += '0 '
        for numb in ray:
            in_str += str(numb) + ' '
        in_str += '\n'
    for vert in vertices:
        in_str += '1 '
        for numb in vert:
            in_str += str(numb) + ' '
        in_str += '\n'
    in_str += 'end\n'
    in_filename = tmp_filename()
    in_file = file(in_filename,'w')
    in_file.write(in_str)
    in_file.close()
    if verbose: print in_str
    # now send the input to the appropriate cddlib function
    if cdd_type == 'real':
        cdd_procs = Popen('testcdd1',stdin = PIPE, stdout=PIPE, stderr=PIPE)
        cast = RR
    else:
        cdd_procs = Popen('testcdd1_gmp',stdin = PIPE, stdout=PIPE, stderr=PIPE)
        cast = Rational
    ans, err = cdd_procs.communicate(input = in_filename)
    # and now parse the output into Sage objects
    ans_lines = ans.splitlines()
    linearity_indices = []
    linearities = []
    hieqs = []
    incidences = []
    adjacencies = []
    if verbose: print ans
    for index in range(len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('linearity') != -1:
            linearity_indices = a_line.split('  ')[1]
            linearity_indices = [int(x) for x in linearity_indices.split(' ')]
        if a_line.find('H-representation') != -1:
            hrep_index = index
        if a_line.find('incidence') != -1:
            inc_index = index
        if a_line.find('adjacency') != -1:
            adj_index = index
    try:
        test = hrep_index
    except:
        print ans
    # add the inequalities, skipping linearities
    ieq_index = 0
    for index in range(hrep_index+2,len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('end') != -1: break
        if a_line.find('begin') == -1 and a_line.find('r') == -1 and a_line.find('i') == -1:
            ieq_index += 1
            if linearity_indices.count(ieq_index) == 0:
                hieqs.append(cdd_convert(a_line, a_type = cast))
            else:
                linearities.append(cdd_convert(a_line, a_type = cast))
    # add incidences, skipping linearities
    ieq_index = 0
    for index in range(inc_index+2,len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('end') != -1: break
        if a_line.find(':') != -1:
            ieq_index += 1
            if linearity_indices.count(ieq_index) == 0:
                pre = a_line.split(':')[0]
                pre = cdd_convert(pre)[0] - 1 # subtract to make indexing pythonic
                post = a_line.split(':')[1]
                post = cdd_convert(post)
                post = [vert_index - 1 for vert_index in post] # subtract to make indexing pythonic
                incidences.append([pre,post])
    # add facial adjacencies, skipping linearities
    ieq_index = 0
    for index in range(adj_index+2,len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('end') != -1: break
        if a_line.find(':') != -1:
            ieq_index += 1
            if linearity_indices.count(ieq_index) == 0:
                pre = a_line.split(':')[0]
                pre = cdd_convert(pre)[0] - 1 # subtract to make indexing pythonic
                post = a_line.split(':')[1]
                post = cdd_convert(post)
                post = [face_index - 1 for face_index in post] # subtract to make indexing pythonic
                adjacencies.append([pre,post])
    poly_obj = Polyhedron(vertices = vertices, cdd_type = cdd_type)
    poly_obj._linearity = linearities
    poly_obj._inequalities = hieqs
    poly_obj._facial_incidences = incidences
    poly_obj._facial_adjacencies = adjacencies
    return poly_obj

def ieq_to_vert(in_list, linearities = [], cdd_type = 'rational', verbose = False):
    """
    Converts a list of inequalities (and optionally some equalities
    (aka linearities)) to a vertex/ray description.

    INPUT:
        in_list -- a list of inequalities.  Each inequality is
            represented by a list [b,a1,a2,...,an].  Suppose your
            variables are x1, x2,..., xn.  Then the inequality
            represented is b + a1*x1 + a2*x2 + ... + an*xn >= 0.
        linearities -- optional.  A list of equalities, in the same
            format as the inequality list described above (but with
            equality of course).
        cdd_type -- optional.  Defaults to rational, and the exact
            arithmetic of testcdd1_gmp.  Using real will invoke the
            floating-point program testcdd1 instead.  This is much
            faster but possibly will not give the correct results for
            highly degenerate polytopes in high dimensions.
        verbose -- optional.  Only interesting for debugging purposes,
            if True this may print out extra information.

    OUTPUT:
        A Polyhedron containing information about the vertices,
        faces, incidences, and adjancencies of the polyhedron.

    EXAMPLES:
        sage: i_list = [[1,0,0,-1],[1,0,-1,0],[1,-1,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
        sage: cube_info = ieq_to_vert(i_list)
        sage: cube_info._vertices
        [[1, 1, 0],
        [0, 1, 0],
        [0, 0, 0],
        [1, 0, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1]]
        sage: e_list = [[3,-2,-2,-2]]
        sage: cut_cube = ieq_to_vert(i_list, e_list)
        sage: cut_cube._vertices #vertices of cut cube, a hexagon.
        [[1/2, 1, 0], [0, 1, 1/2], [1, 1/2, 0], [1, 0, 1/2], [1/2, 0, 1], [0, 1/2, 1]]
    """
    # first we create the input for cddlib:
    in_length = len(in_list[0])
    in_str = 'H-representation\n'
    if linearities != []:
	in_str += 'linearity ' + str(len(linearities)) + ' '
        for index in range(len(linearities)):
	    in_str += str(index+1) + ' '
    in_str += '\n'
    in_str += 'begin\n'
    in_str += str(len(in_list)+len(linearities)) + ' ' + str(in_length) + ' ' + cdd_type + '\n'
    if linearities != []:
        for lin in linearities:
    	    for numb in lin:
    	        in_str += str(numb) + ' '
            in_str += '\n'
    for ieq in in_list:
        for numb in ieq:
            in_str += str(numb) + ' '
        in_str += '\n'
    in_str += 'end\n'
    in_filename = tmp_filename()
    in_file = file(in_filename,'w')
    in_file.write(in_str)
    in_file.close()
    if verbose: print in_str
    # now send the input to the appropriate cddlib function
    if cdd_type == 'real':
        cdd_procs = Popen('testcdd1',stdin = PIPE, stdout=PIPE, stderr=PIPE)
        cast = RR
    else:
        cdd_procs = Popen('testcdd1_gmp',stdin = PIPE, stdout=PIPE, stderr=PIPE)
        cast = Rational
    ans, err = cdd_procs.communicate(input = in_filename)
    # and now parse the output into Sage objects
    ans_lines = ans.splitlines()
    verts = []
    rays = []
    incidences = []
    adjacencies = []
    if verbose: print ans
    for index in range(len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('V-representation') != -1:
            vert_index = index+1
        if a_line.find('incidence') != -1:
            inc_index = index
        if a_line.find('adjacency') != -1:
            adj_index = index
    # read the vertices and rays:
    for index in range(vert_index,len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('end') != -1: break
        if a_line.find('begin') == -1 and a_line.find('r') == -1 and a_line.find('i') == -1:
            ray_or_vertex = cdd_convert(a_line)
            if ray_or_vertex[0] == 1:
                verts.append(ray_or_vertex[1:])
            elif ray_or_vertex[0] == 0:
                rays.append(ray_or_vertex[1:])
            else:
                raise RuntimeError
    # read the vertex incidences
    for index in range(inc_index+2,len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('end') != -1: break
        if a_line.find(':') != -1:
            pre = a_line.split(':')[0]
            pre = cdd_convert(pre)[0] - 1 # subtract to make indexing pythonic
            post = a_line.split(':')[1]
            post = cdd_convert(post)
            post = [face_index - 1 for face_index in post] # subtract to make indexing pythonic
            incidences.append([pre,post])
    # read the vertex-vertex adjacencies (edges)
    for index in range(adj_index+2,len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('end') != -1: break
        if a_line.find(':') != -1:
            pre = a_line.split(':')[0]
            pre = cdd_convert(pre)[0] - 1 # subtract to make indexing pythonic
            post = a_line.split(':')[1]
            post = cdd_convert(post)
            post = [vert_index - 1 for vert_index in post] # subtract to make indexing pythonic
            adjacencies.append([pre,post])
    poly_obj = Polyhedron()
    poly_obj._linearities = linearities
    poly_obj._inequalities = in_list
    poly_obj._vertex_incidences = incidences
    poly_obj._vertex_adjacencies = adjacencies
    poly_obj._vertices = verts
    poly_obj._rays = rays
    return poly_obj

def orthonormal_1(dim_n=5):
    """
    A matrix of rational approximations to orthonormal vectors to (1,...,1).

    INPUT:
        dim_n - the dimension of the vectors

    OUTPUT:
        a matrix over QQ whose rows are close to an orthonormal basis to the
        subspace normal to (1,...,1).

    EXAMPLES:
        sage: from sage.geometry.polyhedra import orthonormal_1
        sage: m = orthonormal_1(5)
        sage: m
        [ 70711/100000   -7071/10000             0             0             0]
        [    1633/4000     1633/4000 -81649/100000             0             0]
        [   7217/25000    7217/25000    7217/25000  -43301/50000             0]
        [ 22361/100000  22361/100000  22361/100000  22361/100000  -44721/50000]
    """
    pb = []
    for i in range(0,dim_n-1):
        pb.append([1.0/(i+1)]*(i+1) + [-1] + [0]*(dim_n-i-2))
    m = matrix(RDF,pb)
    new_m = []
    for i in range(0,dim_n-1):
        new_m.append([RDF(100000*q/norm(m[i])).ceil()/100000 for q in m[i]])
    return matrix(QQ,new_m)

def project_1(fpoint):
    """
    Takes a ndim-dimensional point and projects it onto the plane perpendicular
    to (1,1,...,1).

    INPUT:
        fpoint -- a list of ndim numbers

    EXAMPLES:
        sage: from sage.geometry.polyhedra import project_1
        sage: project_1([1,1,1,1,2])
        [1/100000, 1/100000, 1/50000, -559/625]
    """
    dim_n = len(fpoint)
    p_basis = [list(q) for q in orthonormal_1(dim_n)]
    out_v = []
    for v in p_basis:
        out_v.append(sum([fpoint[ind]*v[ind] for ind in range(dim_n)]))
    return out_v

def _pfunc(i,j,perm):
    """
    An internal utility function for constructing the Birkhoff polytopes.

    EXAMPLES:
        sage: from sage.geometry.polyhedra import _pfunc
        sage: _pfunc(1,2,permutations(3)[0])
        0
    """
    if perm[i-1] == j:
        return 1
    else:
        return 0

class Polytopes():
    """
    A class of constructors for commonly used, famous, or interesting polytopes.
    """

    def Birkhoff_polytope(self, n):
        """
        Returns the Birkhoff polytope with n! vertices.  Each vertex
        is a (flattened) n by n permutation matrix.

        INPUT:
            n -- a positive integer giving the size of the permutation matrices.

        EXAMPLES:
            sage: b3 = polytopes.Birkhoff_polytope(3)
            sage: b3.n_vertices()
            6
        """
        perms = permutations(range(1,n+1))
        verts = []
        for p in perms:
            verts = verts + [[_pfunc(i,j,p) for j in range(1,n+1) for i in range(1,n+1)]]
        return Polyhedron(vertices = verts)

    def n_simplex(self, dim_n):
        """
        Returns a rational approximation to a regular simplex in dimension dim_n.

        INPUT:
            dim_n -- The dimension of the cross-polytope, a positive integer.

        OUTPUT:
            A Polyhedron object of the dim_n-dimensional simplex.

        EXAMPLES:
            sage: s5 = polytopes.n_simplex(5)
            sage: s5.dim()
            5
        """
        verts = permutations([0 for i in range(dim_n)] + [1])
        verts = [project_1(x) for x in verts]
        return Polyhedron(vertices = verts)

    def icosahedron(self):
        """
        Returns an icosahedron with edge length 1.  The vertices are
        rational, so a rational approximation of the golden ratio
        is used.

        OUTPUT:
            A Polyhedron object of a rational approximation to the regular 3D icosahedrone.

        EXAMPLES:
            sage: ico = polytopes.icosahedron()
            sage: sum([len(x[1]) for x in ico.vertex_adjacencies()])/2
            30
        """
        g = QQ(1618033)/1000000 # Golden ratio approximation
        r12 = QQ(1)/2
        verts = [i([0,r12,g/2]) for i in AlternatingGroup(3)]
        verts = verts + [i([0,r12,-g/2]) for i in AlternatingGroup(3)]
        verts = verts + [i([0,-r12,g/2]) for i in AlternatingGroup(3)]
        verts = verts + [i([0,-r12,-g/2]) for i in AlternatingGroup(3)]
        return Polyhedron(vertices = verts)

    def small_rhombicuboctahedron(self):
        """
        An Archimedean solid with 24 vertices and 26 faces.
        """
        verts = [[-3/2, -1/2, -1/2], [-3/2, -1/2, 1/2], [-3/2, 1/2, -1/2], [-3/2, 1/2,
1/2], [-1/2, -3/2, -1/2], [-1/2, -3/2, 1/2], [-1/2, -1/2, -3/2], [-1/2,
-1/2, 3/2], [-1/2, 1/2, -3/2], [-1/2, 1/2, 3/2], [-1/2, 3/2, -1/2],
[-1/2, 3/2, 1/2], [1/2, -3/2, -1/2], [1/2, -3/2, 1/2], [1/2, -1/2,
-3/2], [1/2, -1/2, 3/2], [1/2, 1/2, -3/2], [1/2, 1/2, 3/2], [1/2, 3/2,
-1/2], [1/2, 3/2, 1/2], [3/2, -1/2, -1/2], [3/2, -1/2, 1/2], [3/2, 1/2,
-1/2], [3/2, 1/2, 1/2]]
        return Polyhedron(vertices = verts)

    def great_rhombicuboctahedron(self):
        """
        An Archimedean solid with 48 vertices and 26 faces.
        """
        sqr2 = floor(100000*n(sqrt(2)))/100000
        verts = list(Permutations([1,1+sqr2,1+2*sqr2]))
        verts = verts + [[x[0],x[1],-x[2]] for x in verts]
        verts = verts + [[x[0],-x[1],x[2]] for x in verts]
        verts = verts + [[-x[0],x[1],x[2]] for x in verts]
        return Polyhedron(vertices = verts)

    def twenty_four_cell(self):
        """
        Returns the standard 24-cell polytope.

        OUTPUT:
            A Polyhedron object of the 4-dimensional 24-cell, a regular polytope.
            The coordinates of this polytope are exact.

        EXAMPLES:
            sage: p24 = polytopes.twenty_four_cell()
            sage: len(p24.vertex_adjacencies()[0][1])
            8
        """
        verts = []
        q12 = QQ(1)/2
        base = [q12,q12,q12,q12]
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        verts.append([x for x in base])
                        base[3] = base[3]*(-1)
                    base[2] = base[2]*(-1)
                base[1] = base[1]*(-1)
            base[0] = base[0]*(-1)
        verts = verts + permutations([0,0,0,1])
        verts = verts + permutations([0,0,0,-1])
        return Polyhedron(vertices = verts)

    def six_hundred_cell(self):
        """
        Returns the standard 600-cell polytope.

        OUTPUT:
            A Polyhedron object of the 4-dimensional 600-cell, a regular polytope.
            In many ways this is an analogue of the icosahedron.
            The coordinates of this polytope are rational approximations of the true
            coordinates of the 600-cell, some of which involve the (irrational) golden
            ratio.

        EXAMPLES:
            sage: p600 = polytopes.six_hundred_cell() # long time
            sage: len(p600.vertex_adjacencies()) # long time
            120
        """
        verts = []
        q12 = QQ(1)/2
        base = [q12,q12,q12,q12]
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        verts.append([x for x in base])
                        base[3] = base[3]*(-1)
                    base[2] = base[2]*(-1)
                base[1] = base[1]*(-1)
            base[0] = base[0]*(-1)
        for x in permutations([0,0,0,1]):
            verts.append(x)
        for x in permutations([0,0,0,-1]):
            verts.append(x)
        g = QQ(1618033)/1000000 # Golden ratio approximation
        verts = verts + [i([q12,g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([q12,g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([q12,-g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([q12,-g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,-g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,-g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        return Polyhedron(vertices = verts)

    def cyclic_polytope(self, dim_n, points_n):
        """
        Returns a cyclic polytope in dimension=dim_n and with points_n number
        of points.

        INPUT:
            dim_n -- a positive integer, the dimension of the polytope.
            points_n -- the number of vertices.

        OUTPUT:
            A cyclic polytope of dim_n with points_n vertices on the moment
            curve (t,t^2,...,t^n), as Polyhedron object.

        EXAMPLES:
            sage: c = polytopes.cyclic_polytope(4,10)
            sage: len(c.ieqs())
            35
        """
        verts = [[t**i for i in range(1,dim_n+1)] for t in range(points_n)]
        return Polyhedron(vertices = verts)

    def hypersimplex(self, dim_n, k, project = True):
        """
        The hypersimplex in dimension dim_n with d choose k vertices,
        projected into (dim_n - 1) dimensions.

        INPUT:
            n -- the numbers (1,...,n) are permuted
            project -- optional, default = True.  If false the polyhedron is
        left in dimension n.

        OUTPUT:
            A Polyhedron object representing the hypersimplex.

        EXAMPLES:
            sage: h_4_2 = polytopes.hypersimplex(4,2) # combinatorially equivalent to octahedron
            sage: union([len(x[1]) for x in h_4_2.facial_adjacencies()])
            [3]
        """
        vert0 = [0]*(dim_n-k) + [1]*k
        verts = permutations(vert0)
        if project:
            verts = [project_1(x) for x in verts]
        return Polyhedron(vertices = verts)

    def permutahedron(self, n, project = True):
        """
        The standard permutahedron of (1,...,n) projected into n-1 dimensions.

        INPUT:
            n -- the numbers (1,...,n) are permuted
            project -- optional, default = True.  If false the polyhedron is
        left in dimension n.

        OUTPUT:
            A Polyhedron object representing the permutahedron.

        EXAMPLES:
            sage: perm4 = polytopes.permutahedron(4)
            sage: len(perm4.ieqs())
            14
            sage: union([len(x[1]) for x in perm4.facial_adjacencies()]) # two kinds of faces
            [4, 6]
        """
        verts = range(1,n+1)
        verts = permutations(verts)
        if project:
            verts = [project_1(x) for x in verts]
        p = Polyhedron(vertices = verts)
        return p

    def n_cube(self, dim_n):
        """
        Returns a cube in dimension dim_n.

        INPUT:
            dim_n -- The dimension of the cube, a positive integer.

        OUTPUT:
            A Polyhedron object of the dim_n-dimensional cube, with exact coordinates.

        EXAMPLES:
            sage: four_cube = polytopes.n_cube(4)
            sage: four_cube.is_simple()
            True
        """
        if dim_n == 1:
            verts = [[1],[-1]]
            ieqs = [[1,-1],[1,1]]
            return Polyhedron(vertices = verts, ieqs = ieqs)
        pre_cube = polytopes.n_cube(dim_n-1)
        pre_verts = pre_cube._vertices
        verts = [[1] + vert for vert in pre_verts] + [[-1] + vert for vert in pre_verts]
        ieqs = []
        for index in range(1,dim_n+1):
            temp = [1] + [0]*(dim_n-1)
            temp.insert(index,1)
            ieqs.append([q for q in temp])
            temp = [1] + [0]*(dim_n-1)
            temp.insert(index,-1)
            ieqs.append([q for q in temp])
        return Polyhedron(vertices = verts, ieqs = ieqs)

    def cross_polytope(self, dim_n):
        """
        Returns a cross-polytope in dimension dim_n.  These are
        the generalization of the octahedron.

        INPUT:
            dim_n -- The dimension of the cross-polytope, a positive integer.

        OUTPUT:
            A Polyhedron object of the dim_n-dimensional cross-polytope, with exact coordinates.

        EXAMPLES:
            sage: four_cross = polytopes.cross_polytope(4)
            sage: four_cross.is_simple()
            False
            sage: four_cross.n_vertices()
            8
        """
        verts = permutations([0 for i in range(dim_n-1)] + [1])
        verts = verts + permutations([0 for i in range(dim_n-1)] + [-1])
        return Polyhedron(vertices = verts)

polytopes = Polytopes()
