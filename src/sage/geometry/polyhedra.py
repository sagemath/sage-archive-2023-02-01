"""
Polyhedra

This module contains functions for computing with exact (rational)
or floating-point polyhedra in arbitrary dimensions.  The bulk of
this functionality is provided through the cddlib library of Komei
Fukuda.

There seems to be some inconsistency in the use of the word polyhedra.
In this module, a polyhedron is a convex (possibly unbounded) set
defined either as the convex hull of a set of vertices or rays or as
the intersection of a set of half-planes and hyperplanes.  The
half-planes are also referred to as inequalities, and abbreviated
ieq, and the hyperplanes are also referred to as linearities.

AUTHOR:
    -- Marshall Hampton: first version 2008
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.misc.all import SAGE_TMP, tmp_filename
from sage.rings.all import Integer, QQ
from sage.structure.sage_object import SageObject
from sage.rings.rational import Rational
from sage.rings.real_mpfr import RR
from sage.plot.plot import line
from random import randint
#from sage.plot.plot3d.platonic import * # needed in future for render_solid

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
        sage: len(poly_spam_and_eggs.vertices())
        12
    """
    answer = Polyhedron()
    temp_vertex_list = []
    for vertex1 in polyhedra_list[0].vertices():
        for vertex2 in polyhedra_list[1].vertices():
            temp_vertex_list.append([vertex1[i]+vertex2[i] for i in range(len(vertex1))])
    answer._vertices = temp_vertex_list
    answer.remove_redundant_vertices()
    if verbose:
	print 'Another pair of polytopes computed, sum has ' + str(len(answer._vertices)) + ' vertices'
    if len(polyhedra_list) == 2:
        return answer
    else:
        return mink_sum([answer] + polyhedra_list[2:], verbose = verbose)

class Polyhedron(SageObject):

    def __init__(self, vertices = [], rays = [], ieqs = [], linearities = [], cdd_type = 'rational'):
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
        """
	self._cdd_type = cdd_type
        self._vertices = vertices
        self._ieqs = ieqs
        self._linearities = linearities
        self._rays = rays
	if self._vertices != []:
	    self.remove_redundant_vertices(cdd_type = self._cdd_type)

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

    def is_simple(self):
        """
        Tests for simplicity of a polytope.

        EXAMPLES:
            sage: p = Polyhedron([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
            sage: p.remove_redundant_vertices()
            sage: p.is_simple()
            True
            sage: p = Polyhedron([[0,0,0],[4,4,0],[4,0,0],[0,4,0],[2,2,2]])
            sage: p.is_simple()
            False
        """
        self.remove_redundant_vertices()
        d = self.dim()
        for edge_data in self.vertex_adjacencies():
            if len(edge_data[1]) != d:
                return False
        return True

    def __add__(self, other):
        """
        Addition of two polyhedra is defined as their Minkowski sum.

        EXAMPLES:
            sage: four_cube = n_cube(4)
            sage: four_simplex = Polyhedron(vertices = [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: unholy_union = four_cube + four_simplex
            sage: unholy_union.dim()
            4
        """
        if not isinstance(other, Polyhedron):
            raise TypeError, "other (=%s) must be a Polyhedron"%other
        return mink_sum([self,other])

    def vertices(self, force_from_ieqs = False):
        """
        Returns the vertices of the polyhedron.

        EXAMPLES:
            sage: a_simplex = Polyhedron(ieqs = [[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]], linearities = [[1,-1,-1,-1,-1]])
            sage: a_simplex.vertices()
            [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]
        """
        if (self._vertices == [] and self._ieqs != []) or force_from_ieqs:
            temp_poly_object = ieq_to_vert(self._ieqs, linearities = self._linearities)
            self._vertices = temp_poly_object._vertices
            self._vertex_incidences = temp_poly_object._vertex_incidences
            self._vertex_adjacencies = temp_poly_object._vertex_adjacencies
            self._rays = temp_poly_object._rays
            self._linearities = temp_poly_object._linearities
            self._linearities_checked = True
            self._checked_rays = True
            self.n_vertices = len(self._vertices)
            return self._vertices
        return self._vertices

    def ieqs(self, force_from_vertices = False):
        """
        Return the inequalities (half-planes) defining the polyhedron.

        EXAMPLES:
            sage: permuto4 = Polyhedron(vertices = permutations([1,2,3,4,5]))
            sage: ieqs = permuto4.ieqs()
            sage: ieqs[0]
            [-1, 0, 0, 0, 1, 0]
            sage: ieqs[-1]
            [5, -1, 0, 0, 0, 0]
        """
        if (self._ieqs == [] and self._vertices != []) or force_from_vertices:
            temp_poly_object = vert_to_ieq(self._vertices, rays = self._rays, cdd_type = self._cdd_type)
            self._linearities = temp_poly_object._linearity
            self._ieqs = temp_poly_object._inequalities
            #temp_verts = temp_poly_object.
            self._facial_incidences = temp_poly_object._facial_incidences
            self._facial_adjacencies = temp_poly_object._facial_adjacencies
            return self._ieqs
        else:
            return self._ieqs

    def vertex_adjacencies(self):
        """
        Returns a list of vertex indices and their adjacent vertices.

        EXAMPLES:
            sage: permuto3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: permuto3.vertex_adjacencies()[0:3]
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
            sage: permuto3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: permuto3.facial_adjacencies()[0:3]
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

    def remove_redundant_vertices(self, cdd_type = 'rational'):
        """
        Removes vertices from the description of the polyhedron which
        do not affect the convex hull.  Now done automatically when a
        polyhedron is defined by vertices.

        EXAMPLES:
            sage: a_triangle = Polyhedron([[0,0,0],[4,0,0],[0,4,0],[1,1,0]])
            sage: a_triangle.vertices()
            [[0, 0, 0], [4, 0, 0], [0, 4, 0]]
            sage: a_triangle.remove_redundant_vertices()
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

    def triangulated_facial_incidences(self):
        """
        WARNING: EXPERIMENTAL
        Returns a list of the form [face_index, [v_i_0, v_i_1,,v_i_2]]
        where the face_index refers to the original defining inequality.
        For a given face, the collection of triangles formed by each
        triple of v_i should triangulate that face.
        This is computed by randomly lifting each face up a dimension; this does not always work!

        NOTE:
            The doctest below is very weak, since it does not have to create a triangulation.

        EXAMPLES:
            sage: p = Polyhedron(vertices = [[5,0,0],[0,5,0],[5,5,0],[2,2,5]])
            sage: len(p.triangulated_facial_incidences())
            4
        """
        try:
            return self._triangulated_facial_incidences
        except AttributeError:
            t_fac_incs = []
            for a_face in self.facial_incidences():
                lifted_verts = []
                vert_number = len(a_face[1])
                if vert_number == 3:
                    t_fac_incs.append(a_face)
                else:
                    for vert_index in a_face[1]:
                        lifted_verts.append(self.vertices()[vert_index] + [randint(0,vert_index + vert_number*2)])
                    temp_poly = Polyhedron(vertices = lifted_verts)
                    for t_face in temp_poly.facial_incidences():
                        if len(t_face[1]) != 3:
                            print a_face
                            print t_face
                            print temp_poly.ieqs()
                            print temp_poly.vertices()
                            print temp_poly.facial_incidences()
                            print temp_poly.linearities()
                            raise RuntimeError, "triangulation failed"
                        normal_fdir = temp_poly.ieqs()[t_face[0]][-1]
                        if normal_fdir >= 0:
                            t_fac_verts = [temp_poly.vertices()[i] for i in t_face[1]]
                            proj_verts = [q[0:3] for q in t_fac_verts]
                            t_fac_incs.append([a_face[0], [self.vertices().index(q) for q in proj_verts]])
        self._triangulated_facial_incidences = t_fac_incs
        return t_fac_incs

    #def face_lattice(self):
    #    return "Not implemented yet"

    def render_wireframe(self):
        """
        For polytopes in 2 or 3 dimensions, returns the edges
        as a list of lines.

        EXAMPLES:
            sage: p = Polyhedron([[1,2,],[1,1],[0,0]])
            sage: p_wireframe = p.render_wireframe()
            sage: p_wireframe._Graphics__objects
            [Line defined by 2 points, Line defined by 2 points, Line defined by 2 points]
        """
        if self.ieqs() == self.vertices() == []: return Graphics()
        edges = []
        verts = self.vertices()
        for adj in self.vertex_adjacencies():
            for vert in adj[1]:
                if vert > adj[0]:
                    edges.append([verts[adj[0]],verts[vert]])
        return sum([line(an_edge) for an_edge in edges])

    # Work in progress; needs better triangulation code.
    #def render_solid(self):
    #    """
    #    Returns solid 3d rendering of a 3d polytope.
    #    """
    #    tri_faces = self.triangulated_facial_incidences()
    #    try:
    #        return index_face_set([q[1] for q in tri_faces], self.vertices(), enclosed = True)
    #    except:
    #        print tri_faces
    #        print [q[1] for q in tri_faces], self.vertices()

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
                pre = cdd_convert(pre,a_type = cast)[0] - 1 # subtract to make indexing pythonic
                post = a_line.split(':')[1]
                post = cdd_convert(post,a_type = cast)
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
                pre = cdd_convert(pre,a_type = cast)[0] - 1 # subtract to make indexing pythonic
                post = a_line.split(':')[1]
                post = cdd_convert(post,a_type = cast)
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
            ray_or_vertex = cdd_convert(a_line, a_type = cast)
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
            pre = cdd_convert(pre,a_type = cast)[0] - 1 # subtract to make indexing pythonic
            post = a_line.split(':')[1]
            post = cdd_convert(post,a_type = cast)
            post = [face_index - 1 for face_index in post] # subtract to make indexing pythonic
            incidences.append([pre,post])
    # read the vertex-vertex adjacencies (edges)
    for index in range(adj_index+2,len(ans_lines)):
        a_line = ans_lines[index]
        if a_line.find('end') != -1: break
        if a_line.find(':') != -1:
            pre = a_line.split(':')[0]
            pre = cdd_convert(pre,a_type = cast)[0] - 1 # subtract to make indexing pythonic
            post = a_line.split(':')[1]
            post = cdd_convert(post,a_type = cast)
            post = [vert_index - 1 for vert_index in post] # subtract to make indexing pythonic
            adjacencies.append([pre,post])
    poly_obj = Polyhedron()
    poly_obj._linearity = linearities
    poly_obj._inequalities = in_list
    poly_obj._vertex_incidences = incidences
    poly_obj._vertex_adjacencies = adjacencies
    poly_obj._vertices = verts
    poly_obj._rays = rays
    return poly_obj

def n_cube(dim_n):
    """
    Returns a cube in dimension dim_n.

    EXAMPLES:
        sage: four_cube = n_cube(4)
        sage: four_cube.is_simple()
        True
    """
    if dim_n == 1:
        verts = [[1],[-1]]
        ieqs = [[1,-1],[1,1]]
        return Polyhedron(vertices = verts, ieqs = ieqs)
    pre_cube = n_cube(dim_n-1)
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

