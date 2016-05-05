"""
Indexed Face Sets

Graphics3D object that consists of a list of polygons, also used for
triangulations of other objects.

Usually these objects are not created directly by users.

AUTHORS:

- Robert Bradshaw (2007-08-26): initial version
- Robert Bradshaw (2007-08-28): significant optimizations

.. TODO::

    Smooth triangles using vertex normals

"""
#*****************************************************************************
#      Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
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

include "cysignals/signals.pxi"
include "cysignals/memory.pxi"

cdef extern from *:
    void memset(void *, int, Py_ssize_t)
    void memcpy(void * dest, void * src, Py_ssize_t n)
    int sprintf_3d "sprintf" (char*, char*, double, double, double)
    int sprintf_3i "sprintf" (char*, char*, int, int, int)
    int sprintf_4i "sprintf" (char*, char*, int, int, int, int)
    int sprintf_5i "sprintf" (char*, char*, int, int, int, int, int)
    int sprintf_6i "sprintf" (char*, char*, int, int, int, int, int, int)
    int sprintf_7i "sprintf" (char*, char*, int, int, int, int, int, int, int)
    int sprintf_9d "sprintf" (char*, char*, double, double, double, double, double, double, double, double, double)

# import the double infinity constant
cdef extern from "math.h":
    enum: INFINITY

from cpython.list cimport *
from cpython.string cimport *

include "point_c.pxi"


from math import sin, cos, sqrt
from random import randint

from sage.rings.real_double import RDF

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

from sage.plot.colors import Color, float_to_integer
from sage.plot.plot3d.base import Graphics3dGroup

from transform cimport Transformation


# --------------------------------------------------------------------
# Fast routines for generating string representations of the polygons.
# --------------------------------------------------------------------

cdef inline format_tachyon_texture(color_c rgb):
    cdef char rs[200]
    cdef Py_ssize_t cr = sprintf_3d(rs,
                                   "TEXTURE\n AMBIENT 0.3 DIFFUSE 0.7 SPECULAR 0 OPACITY 1.0\n COLOR %g %g %g \n TEXFUNC 0",
                                   rgb.r, rgb.g, rgb.b)
    return PyString_FromStringAndSize(rs, cr)


cdef inline format_tachyon_triangle(point_c P, point_c Q, point_c R):
    cdef char ss[250]
    # PyString_FromFormat doesn't do floats?
    cdef Py_ssize_t r = sprintf_9d(ss,
                                   "TRI V0 %g %g %g V1 %g %g %g V2 %g %g %g",
                                   P.x, P.y, P.z,
                                   Q.x, Q.y, Q.z,
                                   R.x, R.y, R.z )
    return PyString_FromStringAndSize(ss, r)


cdef inline format_json_vertex(point_c P):
    cdef char ss[100]
    cdef Py_ssize_t r = sprintf_3d(ss, "{x:%g,y:%g,z:%g}", P.x, P.y, P.z)
    return PyString_FromStringAndSize(ss, r)

cdef inline format_json_face(face_c face):
    return "[{}]".format(",".join([str(face.vertices[i])
                                   for i from 0 <= i < face.n]))

cdef inline format_obj_vertex(point_c P):
    cdef char ss[100]
    # PyString_FromFormat doesn't do floats?
    cdef Py_ssize_t r = sprintf_3d(ss, "v %g %g %g", P.x, P.y, P.z)
    return PyString_FromStringAndSize(ss, r)

cdef inline format_obj_face(face_c face, int off):
    cdef char ss[100]
    cdef Py_ssize_t r, i
    if face.n == 3:
        r = sprintf_3i(ss, "f %d %d %d", face.vertices[0] + off, face.vertices[1] + off, face.vertices[2] + off)
    elif face.n == 4:
        r = sprintf_4i(ss, "f %d %d %d %d", face.vertices[0] + off, face.vertices[1] + off, face.vertices[2] + off, face.vertices[3] + off)
    else:
        return "f " + " ".join([str(face.vertices[i] + off) for i from 0 <= i < face.n])
    # PyString_FromFormat is almost twice as slow
    return PyString_FromStringAndSize(ss, r)

cdef inline format_obj_face_back(face_c face, int off):
    cdef char ss[100]
    cdef Py_ssize_t r, i
    if face.n == 3:
        r = sprintf_3i(ss, "f %d %d %d", face.vertices[2] + off, face.vertices[1] + off, face.vertices[0] + off)
    elif face.n == 4:
        r = sprintf_4i(ss, "f %d %d %d %d", face.vertices[3] + off, face.vertices[2] + off, face.vertices[1] + off, face.vertices[0] + off)
    else:
        return "f " + " ".join([str(face.vertices[i] + off) for i from face.n > i >= 0])
    return PyString_FromStringAndSize(ss, r)

cdef inline format_pmesh_vertex(point_c P):
    cdef char ss[100]
    # PyString_FromFormat doesn't do floats?
    cdef Py_ssize_t r = sprintf_3d(ss, "%g %g %g", P.x, P.y, P.z)
    return PyString_FromStringAndSize(ss, r)

cdef inline format_pmesh_face(face_c face, int has_color):
    cdef char ss[100]
    cdef Py_ssize_t r, i
    cdef int color
    # if the face has an individual color, has_color is -1
    # otherwise it is 1
    if has_color == -1:
        color = float_to_integer(face.color.r,
                                 face.color.g,
                                 face.color.b)
        # it seems that Jmol does not like the 0 color at all
        if color == 0:
            color = 1

    if face.n == 3:
        if has_color == 1:
            r = sprintf_5i(ss, "%d\n%d\n%d\n%d\n%d", has_color * 4,
                           face.vertices[0],
                           face.vertices[1],
                           face.vertices[2],
                           face.vertices[0])
        else:
            r = sprintf_6i(ss, "%d\n%d\n%d\n%d\n%d\n%d", has_color * 4,
                           face.vertices[0],
                           face.vertices[1],
                           face.vertices[2],
                           face.vertices[0], color)
    elif face.n == 4:
        if has_color == 1:
            r = sprintf_6i(ss, "%d\n%d\n%d\n%d\n%d\n%d", has_color * 5,
                           face.vertices[0],
                           face.vertices[1],
                           face.vertices[2],
                           face.vertices[3],
                           face.vertices[0])
        else:
            r = sprintf_7i(ss, "%d\n%d\n%d\n%d\n%d\n%d\n%d", has_color * 5,
                           face.vertices[0],
                           face.vertices[1],
                           face.vertices[2],
                           face.vertices[3],
                           face.vertices[0], color)
    else:
        # Naive triangulation
        all = []
        if has_color == 1:
            for i from 1 <= i < face.n - 1:
                r = sprintf_5i(ss, "%d\n%d\n%d\n%d\n%d", has_color * 4,
                               face.vertices[0],
                               face.vertices[i],
                               face.vertices[i + 1],
                               face.vertices[0])
                PyList_Append(all, PyString_FromStringAndSize(ss, r))
        else:
            for i from 1 <= i < face.n - 1:
                r = sprintf_6i(ss, "%d\n%d\n%d\n%d\n%d\n%d", has_color * 4,
                               face.vertices[0],
                               face.vertices[i],
                               face.vertices[i + 1],
                               face.vertices[0], color)
                PyList_Append(all, PyString_FromStringAndSize(ss, r))
        return "\n".join(all)
    # PyString_FromFormat is almost twice as slow
    return PyString_FromStringAndSize(ss, r)


cdef class IndexFaceSet(PrimitiveObject):
    """
    Graphics3D object that consists of a list of polygons, also used for
    triangulations of other objects.

    Polygons (mostly triangles and quadrilaterals) are stored in the
    c struct ``face_c`` (see transform.pyx). Rather than storing
    the points directly for each polygon, each face consists a list
    of pointers into a common list of points which are basically triples
    of doubles in a ``point_c``.

    Moreover, each face has an attribute ``color`` which is used to
    store color information when faces are colored. The red/green/blue
    components are then available as floats between 0 and 1 using
    ``color.r,color.g,color.b``.

    Usually these objects are not created directly by users.

    EXAMPLES::

        sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
        sage: S = IndexFaceSet([[(1,0,0),(0,1,0),(0,0,1)],[(1,0,0),(0,1,0),(0,0,0)]])
        sage: S.face_list()
        [[(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)], [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)]]
        sage: S.vertex_list()
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 0.0)]

        sage: def make_face(n): return [(0,0,n),(0,1,n),(1,1,n),(1,0,n)]
        sage: S = IndexFaceSet([make_face(n) for n in range(10)])
        sage: S.show()

        sage: point_list = [(1,0,0),(0,1,0)] + [(0,0,n) for n in range(10)]
        sage: face_list = [[0,1,n] for n in range(2,10)]
        sage: S = IndexFaceSet(face_list, point_list, color='red')
        sage: S.face_list()
        [[(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 0.0)],
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 2.0)],
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 3.0)],
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 4.0)],
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 5.0)],
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 6.0)],
        [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 7.0)]]
        sage: S.show()

    A simple example of colored IndexFaceSet (:trac:`12212`)::

        sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
        sage: from sage.plot.plot3d.texture import Texture
        sage: point_list = [(2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0)]
        sage: face_list = [[0,4,5],[3,4,5],[2,3,4],[1,3,5]]
        sage: col = rainbow(10, 'rgbtuple')
        sage: t_list = [Texture(col[i]) for i in range(10)]
        sage: S = IndexFaceSet(face_list, point_list, texture_list=t_list)
        sage: S.show(viewer='tachyon')

    """
    def __cinit__(self):
        self.vs = <point_c *>NULL
        self.face_indices = <int *>NULL
        self._faces = <face_c *>NULL
    def __init__(self, faces, point_list=None,
                 enclosed=False, texture_list=None, **kwds):
        PrimitiveObject.__init__(self, **kwds)

        self.global_texture = (texture_list is None)

        self.enclosed = enclosed

        if point_list is None:
            face_list = faces
            faces = []
            point_list = []
            point_index = {}
            for face in face_list:
                iface = []
                for p in face:
                    try:
                        ix = point_index[p]
                    except KeyError:
                        ix = len(point_list)
                        point_index[p] = ix
                        point_list.append(p)
                    iface.append(ix)
                faces.append(iface)

        cdef Py_ssize_t i
        cdef Py_ssize_t index_len = 0
        for i from 0 <= i < len(faces):
            index_len += len(faces[i])

        self.vcount = len(point_list)
        self.fcount = len(faces)
        self.icount = index_len

        self.realloc(self.vcount, self.fcount, index_len)

        for i from 0 <= i < self.vcount:
            self.vs[i].x, self.vs[i].y, self.vs[i].z = point_list[i]

        cdef int cur_pt = 0
        for i from 0 <= i < self.fcount:
            self._faces[i].n = len(faces[i])
            self._faces[i].vertices = &self.face_indices[cur_pt]
            if self.global_texture:
                self._faces[i].color.r, self._faces[i].color.g, self._faces[i].color.b = self.texture.color
            else:
                self._faces[i].color.r, self._faces[i].color.g, self._faces[i].color.b = texture_list[i].color
            for ix in faces[i]:
                self.face_indices[cur_pt] = ix
                cur_pt += 1

    cdef realloc(self, vcount, fcount, icount):
        r"""
        Allocates memory for vertices, faces, and face indices.  Can
        only be called from Cython, so the doctests must be indirect.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: G = implicit_plot3d(x^2+y^2+z^2 - 1, (x, -2, 2), (y, -2, 2), (z, -2, 2), plot_points=6)
            sage: G.triangulate()  # indirect doctest
            sage: len(G.face_list())
            44
            sage: len(G.vertex_list())
            132
            sage: G = implicit_plot3d(x^2+y^2+z^2 - 100, (x, -2, 2), (y, -2, 2), (z, -2, 2), plot_points=6)
            sage: G.triangulate()  # indirect doctest
            sage: len(G.face_list())
            0
            sage: len(G.vertex_list())
            0
        """
        self.vs = <point_c *>check_reallocarray(self.vs, vcount, sizeof(point_c))
        self._faces = <face_c *>check_reallocarray(self._faces, fcount, sizeof(face_c))
        self.face_indices = <int *>check_reallocarray(self.face_indices, icount, sizeof(int))

    def _clean_point_list(self):
        # TODO: There is still wasted space where quadrilaterals were
        # converted to triangles...  but it's so little it's probably
        # not worth bothering with
        cdef int* point_map = <int *>check_calloc(self.vcount, sizeof(int))
        cdef Py_ssize_t i, j
        cdef face_c *face
        for i from 0 <= i < self.fcount:
            face = &self._faces[i]
            for j from 0 <= j < face.n:
                point_map[face.vertices[j]] += 1
        ix = 0
        for i from 0 <= i < self.vcount:
            if point_map[i] > 0:
                point_map[i] = ix
                self.vs[ix] = self.vs[i]
                ix += 1
        if ix != self.vcount:
            for i from 0 <= i < self.fcount:
                face = &self._faces[i]
                for j from 0 <= j < face.n:
                    face.vertices[j] = point_map[face.vertices[j]]
            self.realloc(ix, self.fcount, self.icount)
            self.vcount = ix
        sig_free(point_map)

    def _seperate_creases(self, threshold):
        """
        Some rendering engines Gouraud shading, which is great for smooth
        surfaces but looks bad if one actually has a polyhedron.

        INPUT:

        ``threshold`` -- the minimum cosine of the angle between adjacent
        faces a higher threshold separates more, all faces if >= 1, no
        faces if <= -1
        """
        cdef Py_ssize_t i, j, k
        cdef face_c *face
        cdef int v, count, total = 0
        cdef int* point_counts = <int *>check_calloc(self.vcount * 2 + 1, sizeof(int))
        # For each vertex, get number of faces
        cdef int* running_point_counts = &point_counts[self.vcount]
        for i from 0 <= i < self.fcount:
            face = &self._faces[i]
            total += face.n
            for j from 0 <= j < face.n:
                point_counts[face.vertices[j]] += 1
        # Running used as index into face list
        cdef int running = 0
        cdef int max = 0
        for i from 0 <= i < self.vcount:
            running_point_counts[i] = running
            running += point_counts[i]
            if point_counts[i] > max:
                max = point_counts[i]
        running_point_counts[self.vcount] = running
        # Create an array, indexed by running_point_counts[v], to the list of faces containing that vertex.
        cdef face_c** point_faces
        try:
            point_faces = <face_c **>check_allocarray(total, sizeof(face_c*))
        except MemoryError:
            sig_free(point_counts)
            raise
        sig_on()
        memset(point_counts, 0, sizeof(int) * self.vcount)
        for i from 0 <= i < self.fcount:
            face = &self._faces[i]
            for j from 0 <= j < face.n:
                v = face.vertices[j]
                point_faces[running_point_counts[v]+point_counts[v]] = face
                point_counts[v] += 1
        # Now, for each vertex, see if all faces are close enough,
        # or if it is a crease.
        cdef face_c** faces
        cdef int start = 0
        cdef bint any
        # We compare against face 0, and if it's not flat enough we push it to the end.
        # Then we come around again to compare everything that was put at the end, possibly
        # pushing stuff to the end again (until no further changes are needed).
        while start < self.vcount:
            ix = self.vcount
            # Find creases
            for i from 0 <= i < self.vcount - start:
                faces = &point_faces[running_point_counts[i]]
                any = 0
                for j from point_counts[i] > j >= 1:
                    if cos_face_angle(faces[0][0], faces[j][0], self.vs) < threshold:
                        any = 1
                        face = faces[j]
                        point_counts[i] -= 1
                        if j != point_counts[i]:
                            faces[j] = faces[point_counts[i]] # swap
                            faces[point_counts[i]] = face
                if any:
                    ix += 1
            # Reallocate room for vertices at end
            if ix > self.vcount:
                try:
                    self.vs = <point_c *>check_reallocarray(self.vs, ix, sizeof(point_c))
                except MemoryError:
                    sig_free(point_counts)
                    sig_free(point_faces)
                    self.vcount = self.fcount = self.icount = 0 # so we don't get segfaults on bad points
                    sig_off()
                    raise
                ix = self.vcount
                running = 0
                for i from 0 <= i < self.vcount - start:
                    if point_counts[i] != running_point_counts[i+1] - running_point_counts[i]:
                        # We have a new vertex
                        self.vs[ix] = self.vs[i+start]
                        # Update the point_counts and point_faces arrays for the next time around.
                        count = running_point_counts[i+1] - running_point_counts[i] - point_counts[i]
                        faces = &point_faces[running]
                        for j from 0 <= j < count:
                            faces[j] = point_faces[running_point_counts[i] + point_counts[i] + j]
                            face = faces[j]
                            for k from 0 <= k < face.n:
                                if face.vertices[k] == i + start:
                                    face.vertices[k] = ix
                        point_counts[ix-self.vcount] = count
                        running_point_counts[ix-self.vcount] = running
                        running += count
                        ix += 1
                running_point_counts[ix-self.vcount] = running
            start = self.vcount
            self.vcount = ix

        sig_free(point_counts)
        sig_free(point_faces)
        sig_off()

    def _mem_stats(self):
        return self.vcount, self.fcount, self.icount

    def __dealloc__(self):
        sig_free(self.vs)
        sig_free(self._faces)
        sig_free(self.face_indices)

    def is_enclosed(self):
        """
        Whether or not it is necessary to render the back sides of the polygons.

        One is assuming, of course, that they have the correct orientation.

        This is may be passed in on construction. It is also
        calculated in
        :class:`sage.plot.plot3d.parametric_surface.ParametricSurface`
        by verifying the opposite edges of the rendered domain either
        line up or are pinched together.

        EXAMPLES::

            sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
            sage: IndexFaceSet([[(0,0,1),(0,1,0),(1,0,0)]]).is_enclosed()
            False
        """
        return self.enclosed

    def index_faces(self):
        """
        Return the list over all faces of the indices of the vertices.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Box(1,2,3)
            sage: S.index_faces()
            [[0, 1, 2, 3],
             [0, 4, 5, 1],
             [0, 3, 6, 4],
             [5, 4, 6, 7],
             [6, 3, 2, 7],
             [2, 1, 5, 7]]
        """
        cdef Py_ssize_t i, j
        return [[self._faces[i].vertices[j]
                 for j from 0 <= j < self._faces[i].n]
                for i from 0 <= i < self.fcount]

    def has_local_colors(self):
        """
        Return ``True`` if and only if every face has an individual color.

        EXAMPLES::

            sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
            sage: from sage.plot.plot3d.texture import Texture
            sage: point_list = [(2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0)]
            sage: face_list = [[0,4,5],[3,4,5],[2,3,4],[1,3,5]]
            sage: col = rainbow(10, 'rgbtuple')
            sage: t_list=[Texture(col[i]) for i in range(10)]
            sage: S = IndexFaceSet(face_list, point_list, texture_list=t_list)
            sage: S.has_local_colors()
            True

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Box(1,2,3)
            sage: S.has_local_colors()
            False
        """
        return not(self.global_texture)
    
    def index_faces_with_colors(self):
        """
        Return the list over all faces of (indices of the vertices, color).

        This only works if every face has its own color.

        .. SEEALSO::

            :meth:`has_local_colors`

        EXAMPLES:

        A simple colored one::

            sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
            sage: from sage.plot.plot3d.texture import Texture
            sage: point_list = [(2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0)]
            sage: face_list = [[0,4,5],[3,4,5],[2,3,4],[1,3,5]]
            sage: col = rainbow(10, 'rgbtuple')
            sage: t_list=[Texture(col[i]) for i in range(10)]
            sage: S = IndexFaceSet(face_list, point_list, texture_list=t_list)
            sage: S.index_faces_with_colors()
            [([0, 4, 5], '#ff0000'),
            ([3, 4, 5], '#ff9900'),
            ([2, 3, 4], '#cbff00'),
            ([1, 3, 5], '#33ff00')]

        When the texture is global, an error is raised::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Box(1,2,3)
            sage: S.index_faces_with_colors()
            Traceback (most recent call last):
            ...
            ValueError: the texture is global
        """
        cdef Py_ssize_t i, j
        if self.global_texture:
            raise ValueError('the texture is global')
        return [([self._faces[i].vertices[j]
                  for j from 0 <= j < self._faces[i].n],
                 Color(self._faces[i].color.r,
                       self._faces[i].color.g,
                       self._faces[i].color.b).html_color())
                for i from 0 <= i < self.fcount]

    def faces(self):
        """
        An iterator over the faces.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Box(1,2,3)
            sage: list(S.faces()) == S.face_list()
            True
        """
        return FaceIter(self)

    def face_list(self):
        """
        Return the list of faces.

        Every face is given as a tuple of vertices.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Box(1,2,3)
            sage: S.face_list()[0]
            [(1.0, 2.0, 3.0), (-1.0, 2.0, 3.0), (-1.0, -2.0, 3.0), (1.0, -2.0, 3.0)]
        """
        points = self.vertex_list()
        cdef Py_ssize_t i, j
        return [[points[self._faces[i].vertices[j]]
                 for j from 0 <= j < self._faces[i].n]
                for i from 0 <= i < self.fcount]

    def edges(self):
        """
        An iterator over the edges.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Box(1,2,3)
            sage: list(S.edges())[0]
            ((1.0, -2.0, 3.0), (1.0, 2.0, 3.0))
        """
        return EdgeIter(self)

    def edge_list(self):
        """
        Return the list of edges.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Box(1,2,3)
            sage: S.edge_list()[0]
            ((1.0, -2.0, 3.0), (1.0, 2.0, 3.0))
        """
        return list(self.edges())

    def vertices(self):
        """
        An iterator over the vertices.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Cone(1,1)
            sage: list(S.vertices()) == S.vertex_list()
            True
        """
        return VertexIter(self)

    def vertex_list(self):
        """
        Return the list of vertices.

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = polygon([(0,0,1), (1,1,1), (2,0,1)])
            sage: S.vertex_list()[0]
            (0.0, 0.0, 1.0)
        """
        cdef Py_ssize_t i
        return [(self.vs[i].x, self.vs[i].y, self.vs[i].z) for i from 0 <= i < self.vcount]

    def x3d_geometry(self):
        """
        Return the x3d data.

        EXAMPLES:

        A basic test with a triangle::

            sage: G = polygon([(0,0,1), (1,1,1), (2,0,1)])
            sage: print G.x3d_geometry()
            <BLANKLINE>
            <IndexedFaceSet coordIndex='0,1,2,-1'>
              <Coordinate point='0.0 0.0 1.0,1.0 1.0 1.0,2.0 0.0 1.0'/>
            </IndexedFaceSet>
            <BLANKLINE>

        A simple colored one::

            sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
            sage: from sage.plot.plot3d.texture import Texture
            sage: point_list = [(2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0)]
            sage: face_list = [[0,4,5],[3,4,5],[2,3,4],[1,3,5]]
            sage: col = rainbow(10, 'rgbtuple')
            sage: t_list=[Texture(col[i]) for i in range(10)]
            sage: S = IndexFaceSet(face_list, point_list, texture_list=t_list)
            sage: print S.x3d_geometry()
            <BLANKLINE>
            <IndexedFaceSet solid='False' colorPerVertex='False' coordIndex='0,4,5,-1,3,4,5,-1,2,3,4,-1,1,3,5,-1'>
              <Coordinate point='2.0 0.0 0.0,0.0 2.0 0.0,0.0 0.0 2.0,0.0 1.0 1.0,1.0 0.0 1.0,1.0 1.0 0.0'/>
              <Color color='1.0 0.0 0.0,1.0 0.6 0.0,0.8 1.0 0.0,0.2 1.0 0.0' />
            </IndexedFaceSet>
            <BLANKLINE>
        """
        cdef Py_ssize_t i
        points = ",".join(["%s %s %s" % (self.vs[i].x,
                                         self.vs[i].y,
                                         self.vs[i].z)
                           for i from 0 <= i < self.vcount])
        coordIndex = ",-1,".join([",".join([str(self._faces[i].vertices[j])
                                            for j from 0 <= j < self._faces[i].n])
                                  for i from 0 <= i < self.fcount])
        if not self.global_texture:
            colorIndex = ",".join([str(self._faces[i].color.r) + " "
                                   + str(self._faces[i].color.g) + " "
                                   + str(self._faces[i].color.b)
                                   for i from 0 <= i < self.fcount])
            return """
<IndexedFaceSet solid='False' colorPerVertex='False' coordIndex='%s,-1'>
  <Coordinate point='%s'/>
  <Color color='%s' />
</IndexedFaceSet>
""" % (coordIndex, points, colorIndex)
        return """
<IndexedFaceSet coordIndex='%s,-1'>
  <Coordinate point='%s'/>
</IndexedFaceSet>
""" % (coordIndex, points)

    def bounding_box(self):
        r"""
        Calculate the bounding box for the vertices in this object
        (ignoring infinite or NaN coordinates).

        OUTPUT:

        a tuple ( (low_x, low_y, low_z), (high_x, high_y, high_z)),
        which gives the coordinates of opposite corners of the
        bounding box.

        EXAMPLE::

            sage: x,y = var('x,y')
            sage: p = plot3d(sqrt(sin(x)*sin(y)), (x,0,2*pi),(y,0,2*pi))
            sage: p.bounding_box()
            ((0.0, 0.0, -0.0), (6.283185307179586, 6.283185307179586, 0.9991889981715697))
        """
        if self.vcount == 0:
            return ((0,0,0),(0,0,0))

        cdef Py_ssize_t i
        cdef point_c low
        cdef point_c high

        low.x, low.y, low.z = INFINITY, INFINITY, INFINITY
        high.x, high.y, high.z = -INFINITY, -INFINITY, -INFINITY

        for i in range(0,self.vcount):
            point_c_update_finite_lower_bound(&low, self.vs[i])
            point_c_update_finite_upper_bound(&high, self.vs[i])
        return ((low.x, low.y, low.z), (high.x, high.y, high.z))

    def partition(self, f):
        r"""
        Partition the faces of ``self``.

        The partition is done according to the value of a map
        `f: \RR^3 \rightarrow \ZZ` applied to the center of each face.

        INPUT:

        - `f` -- a function from `\RR^3` to `\ZZ`

        EXAMPLES::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Box(1,2,3)
            sage: len(S.partition(lambda x,y,z : floor(x+y+z)))
            6
        """
        cdef Py_ssize_t i, j, ix, face_ix
        cdef int part
        cdef point_c P
        cdef face_c *face
        cdef face_c *new_face
        cdef IndexFaceSet face_set

        cdef int *partition = <int *>check_allocarray(self.fcount, sizeof(int))

        part_counts = {}
        for i from 0 <= i < self.fcount:
            face = &self._faces[i]
            P = self.vs[face.vertices[0]]
            for j from 1 <= j < face.n:
                point_c_add(&P, P, self.vs[face.vertices[j]])
            point_c_mul(&P, P, 1.0/face.n)
            partition[i] = part = f(P.x, P.y, P.z)
            try:
                count = part_counts[part]
            except KeyError:
                part_counts[part] = count = [0,0]
            count[0] += 1
            count[1] += face.n
        all = {}
        for part, count in part_counts.iteritems():
            face_set = IndexFaceSet([])
            face_set.realloc(self.vcount, count[0], count[1])
            face_set.vcount = self.vcount
            face_set.fcount = count[0]
            face_set.icount = count[1]
            memcpy(face_set.vs, self.vs, sizeof(point_c) * self.vcount)
            face_ix = 0
            ix = 0
            for i from 0 <= i < self.fcount:
                if partition[i] == part:
                    face = &self._faces[i]
                    new_face = &face_set._faces[face_ix]
                    new_face.n = face.n
                    new_face.vertices = &face_set.face_indices[ix]
                    for j from 0 <= j < face.n:
                        new_face.vertices[j] = face.vertices[j]
                    face_ix += 1
                    ix += face.n
            face_set._clean_point_list()
            all[part] = face_set
        sig_free(partition)
        return all

    def tachyon_repr(self, render_params):
        """
        Return a tachyon object for ``self``.

        EXAMPLES:

        A basic test with a triangle::

            sage: G = polygon([(0,0,1), (1,1,1), (2,0,1)])
            sage: s = G.tachyon_repr(G.default_render_params()); s
            ['TRI V0 0 0 1 V1 1 1 1 V2 2 0 1', ...]

        A simple colored one::

            sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
            sage: from sage.plot.plot3d.texture import Texture
            sage: point_list = [(2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0)]
            sage: face_list = [[0,4,5],[3,4,5],[2,3,4],[1,3,5]]
            sage: col = rainbow(10, 'rgbtuple')
            sage: t_list=[Texture(col[i]) for i in range(10)]
            sage: S = IndexFaceSet(face_list, point_list, texture_list=t_list)
            sage: S.tachyon_repr(S.default_render_params())
            ['TRI V0 2 0 0 V1 1 0 1 V2 1 1 0',
            'TEXTURE... AMBIENT 0.3 DIFFUSE 0.7 SPECULAR 0 OPACITY 1.0... COLOR 1 0 0 ... TEXFUNC 0',...]
        """
        cdef Transformation transform = render_params.transform
        lines = []
        cdef point_c P, Q, R
        cdef face_c face
        cdef Py_ssize_t i, k
        sig_on()
        for i from 0 <= i < self.fcount:
            face = self._faces[i]
            if transform is not None:
                transform.transform_point_c(&P, self.vs[face.vertices[0]])
                transform.transform_point_c(&Q, self.vs[face.vertices[1]])
                transform.transform_point_c(&R, self.vs[face.vertices[2]])
            else:
                P = self.vs[face.vertices[0]]
                Q = self.vs[face.vertices[1]]
                R = self.vs[face.vertices[2]]
            PyList_Append(lines, format_tachyon_triangle(P, Q, R))
            if self.global_texture:
                PyList_Append(lines, self.texture.id)
            else:
                PyList_Append(lines, format_tachyon_texture(face.color))
            if face.n > 3:
                for k from 3 <= k < face.n:
                    Q = R
                    if transform is not None:
                        transform.transform_point_c(&R, self.vs[face.vertices[k]])
                    else:
                        R = self.vs[face.vertices[k]]
                    PyList_Append(lines, format_tachyon_triangle(P, Q, R))
                    if self.global_texture:
                        PyList_Append(lines, self.texture.id)
                    else:
                        PyList_Append(lines, format_tachyon_texture(face.color))
        sig_off()

        return lines

    def json_repr(self, render_params):
        """
        Return a json representation for ``self``.

        TESTS:

        A basic test with a triangle::

            sage: G = polygon([(0,0,1), (1,1,1), (2,0,1)])
            sage: G.json_repr(G.default_render_params())
            ["{vertices:[{x:0,y:0,z:1},{x:1,y:1,z:1},{x:2,y:0,z:1}],faces:[[0,1,2]],color:'#0000ff'}"]

        A simple colored one::

            sage: from sage.plot.plot3d.index_face_set import IndexFaceSet
            sage: from sage.plot.plot3d.texture import Texture
            sage: point_list = [(2,0,0),(0,2,0),(0,0,2),(0,1,1),(1,0,1),(1,1,0)]
            sage: face_list = [[0,4,5],[3,4,5],[2,3,4],[1,3,5]]
            sage: col = rainbow(10, 'rgbtuple')
            sage: t_list=[Texture(col[i]) for i in range(10)]
            sage: S = IndexFaceSet(face_list, point_list, texture_list=t_list)
            sage: S.json_repr(S.default_render_params())
            ["{vertices:[{x:2,y:0,z:0},{x:0,y:2,z:0},{x:0,y:0,z:2},{x:0,y:1,z:1},{x:1,y:0,z:1},{x:1,y:1,z:0}],faces:[[0,4,5],[3,4,5],[2,3,4],[1,3,5]],face_colors:['#ff0000','#ff9900','#cbff00','#33ff00']}"]
        """
        cdef Transformation transform = render_params.transform
        cdef point_c res

        if transform is None:
            vertices_str = "[{}]".format(
                ",".join([format_json_vertex(self.vs[i])
                          for i from 0 <= i < self.vcount]))
        else:
            vertices_str = "["
            for i from 0 <= i < self.vcount:
                transform.transform_point_c(&res, self.vs[i])
                if i > 0:
                    vertices_str += ","
                vertices_str += format_json_vertex(res)
            vertices_str += "]"

        faces_str = "[{}]".format(",".join([format_json_face(self._faces[i])
                                            for i from 0 <= i < self.fcount]))
        if self.global_texture:
            color_str = "'#{}'".format(self.texture.hex_rgb())
            return ["{vertices:%s,faces:%s,color:%s}" %
                    (vertices_str, faces_str, color_str)]
        else:
            color_str = "[{}]".format(",".join(["'{}'".format(
                    Color(self._faces[i].color.r,
                          self._faces[i].color.g,
                          self._faces[i].color.b).html_color())
                                            for i from 0 <= i < self.fcount]))
            return ["{vertices:%s,faces:%s,face_colors:%s}" %
                    (vertices_str, faces_str, color_str)]

    def obj_repr(self, render_params):
        """
        Return an obj representation for ``self``.

        TESTS::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Cylinder(1,1)
            sage: s = S.obj_repr(S.default_render_params())
        """
        cdef Transformation transform = render_params.transform
        cdef int off = render_params.obj_vertex_offset
        cdef Py_ssize_t i
        cdef point_c res

        sig_on()
        if transform is None:
            points = [format_obj_vertex(self.vs[i]) for i from 0 <= i < self.vcount]
        else:
            points = []
            for i from 0 <= i < self.vcount:
                transform.transform_point_c(&res, self.vs[i])
                PyList_Append(points, format_obj_vertex(res))

        faces = [format_obj_face(self._faces[i], off) for i from 0 <= i < self.fcount]
        if not self.enclosed:
            back_faces = [format_obj_face_back(self._faces[i], off) for i from 0 <= i < self.fcount]
        else:
            back_faces = []

        render_params.obj_vertex_offset += self.vcount
        sig_off()

        return ["g " + render_params.unique_name('obj'),
                "usemtl " + self.texture.id,
                points,
                faces,
                back_faces]

    def jmol_repr(self, render_params):
        """
        Return a jmol representation for ``self``.

        TESTS::

            sage: from sage.plot.plot3d.shapes import *
            sage: S = Cylinder(1,1)
            sage: S.show(viewer='jmol')   # indirect doctest
        """
        cdef Transformation transform = render_params.transform
        cdef Py_ssize_t i
        cdef point_c res

        self._seperate_creases(render_params.crease_threshold)

        sig_on()
        if transform is None:
            points = [format_pmesh_vertex(self.vs[i])
                      for i from 0 <= i < self.vcount]
        else:
            points = []
            for i from 0 <= i < self.vcount:
                transform.transform_point_c(&res, self.vs[i])
                PyList_Append(points, format_pmesh_vertex(res))

        # activation of coloring in jmol
        if self.global_texture:
            faces = [format_pmesh_face(self._faces[i], 1)
                     for i from 0 <= i < self.fcount]
        else:
            faces = [format_pmesh_face(self._faces[i], -1)
                     for i from 0 <= i < self.fcount]

        # If a face has more than 4 vertices, it gets chopped up in
        # format_pmesh_face
        cdef Py_ssize_t extra_faces = 0
        for i from 0 <= i < self.fcount:
            if self._faces[i].n >= 5:
                extra_faces += self._faces[i].n-3

        sig_off()

        all = [str(self.vcount),
               points,
               str(self.fcount + extra_faces),
               faces]

        from base import flatten_list
        name = render_params.unique_name('obj')
        all = flatten_list(all)
        if render_params.output_archive:
            filename = "%s.pmesh" % (name)
            render_params.output_archive.writestr(filename, '\n'.join(all))
        else:
            filename = "%s-%s.pmesh" % (render_params.output_file, name)
            f = open(filename, 'w')
            for line in all:
                f.write(line)
                f.write('\n')
            f.close()

        if self.global_texture:
            s = 'pmesh {} "{}"\n{}'.format(name, filename,
                                           self.texture.jmol_str("pmesh"))
        else:
            s = 'pmesh {} "{}"'.format(name, filename)            

        # Turn on display of the mesh lines or dots?
        if render_params.mesh:
            s += '\npmesh %s mesh\n' % name
        if render_params.dots:
            s += '\npmesh %s dots\n' % name
        return [s]

    def dual(self, **kwds):
        """
        Return the dual.

        EXAMPLES::

            sage: S = cube()
            sage: T = S.dual()
            sage: len(T.vertex_list())
            6

        """
        cdef point_c P
        cdef face_c *face
        cdef Py_ssize_t i, j, ix, ff
        cdef IndexFaceSet dual = IndexFaceSet([], **kwds)
        cdef int incoming, outgoing

        dual.realloc(self.fcount, self.vcount, self.icount)

        sig_on()
        # is using dicts overly-heavy?
        dual_faces = [{} for i from 0 <= i < self.vcount]

        for i from 0 <= i < self.fcount:
            # Let the vertex be centered on the face according to a simple average
            face = &self._faces[i]
            dual.vs[i] = self.vs[face.vertices[0]]
            for j from 1 <= j < face.n:
                point_c_add(&dual.vs[i], dual.vs[i], self.vs[face.vertices[j]])
            point_c_mul(&dual.vs[i], dual.vs[i], 1.0/face.n)

            # Now compute the new face
            for j from 0 <= j < face.n:
                if j == 0:
                    incoming = face.vertices[face.n-1]
                else:
                    incoming = face.vertices[j-1]
                if j == face.n-1:
                    outgoing = face.vertices[0]
                else:
                    outgoing = face.vertices[j+1]
                dd = dual_faces[face.vertices[j]]
                dd[incoming] = i, outgoing

        i = 0
        ix = 0
        for dd in dual_faces:
            face = &dual._faces[i]
            face.n = len(dd)
            if face.n == 0: # skip unused vertices
                continue
            face.vertices = &dual.face_indices[ix]
            ff, next_ = next(dd.itervalues())
            face.vertices[0] = ff
            for j from 1 <= j < face.n:
                ff, next_ = dd[next_]
                face.vertices[j] = ff
            i += 1
            ix += face.n

        dual.vcount = self.fcount
        dual.fcount = i
        dual.icount = ix
        sig_off()

        return dual

    def stickers(self, colors, width, hover):
        """
        Return a group of IndexFaceSets.

        INPUT:

        - ``colors`` -- list of colors/textures to use (in cyclic order)

        - ``width`` -- offset perpendicular into the edge (to create a border)
          may also be negative

        - ``hover`` -- offset normal to the face (usually have to float above
          the original surface so it shows, typically this value is very
          small compared to the actual object

        OUTPUT:

        Graphics3dGroup of stickers

        EXAMPLE::

            sage: from sage.plot.plot3d.shapes import Box
            sage: B = Box(.5,.4,.3, color='black')
            sage: S = B.stickers(['red','yellow','blue'], 0.1, 0.05)
            sage: S.show()
            sage: (S+B).show()

        """
        all = []
        n = self.fcount
        ct = len(colors)
        for k in range(len(colors)):
            if colors[k]:
                all.append(self.sticker(range(k, n, ct), width, hover,
                                        texture=colors[k]))
        return Graphics3dGroup(all)

    def sticker(self, face_list, width, hover, **kwds):
        """
        Return a sticker on the chosen faces.
        """
        if not isinstance(face_list, (list, tuple)):
            face_list = (face_list,)
        faces = self.face_list()
        all = []
        for k in face_list:
            all.append(sticker(faces[k], width, hover))
        return IndexFaceSet(all, **kwds)


cdef class FaceIter:
    """
    A class for iteration over faces

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import *
        sage: S = Box(1,2,3)
        sage: len(list(S.faces())) == 6   # indirect doctest
        True
    """
    def __init__(self, face_set):
        """
        """
        self.set = face_set
        self.i = 0

    def __iter__(self):
        return self

    def __next__(self):
        cdef point_c P
        if self.i >= self.set.fcount:
            raise StopIteration
        else:
            face = []
            for j from 0 <= j < self.set._faces[self.i].n:
                P = self.set.vs[self.set._faces[self.i].vertices[j]]
                PyList_Append(face, (P.x, P.y, P.z))
            self.i += 1
            return face


cdef class EdgeIter:
    """
    A class for iteration over edges

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import *
        sage: S = Box(1,2,3)
        sage: len(list(S.edges())) == 12   # indirect doctest
        True
    """
    def __init__(self, face_set):
        self.set = face_set
        if not self.set.enclosed:
            raise TypeError("Must be closed to use the simple iterator.")
        self.i = 0
        self.j = 0
        self.seen = {}

    def __iter__(self):
        return self

    def __next__(self):
        cdef point_c P, Q
        cdef face_c face = self.set._faces[self.i]
        while self.i < self.set.fcount:
            if self.j == face.n:
                self.i += 1
                self.j = 0
                if self.i < self.set.fcount:
                    face = self.set._faces[self.i]
            else:
                if self.j == 0:
                    P = self.set.vs[face.vertices[face.n-1]]
                else:
                    P = self.set.vs[face.vertices[self.j-1]]
                Q = self.set.vs[face.vertices[self.j]]
                self.j += 1
                if self.set.enclosed:  # Every edge appears exactly twice, once in each orientation.
                    if point_c_cmp(P, Q) < 0:
                        return ((P.x, P.y, P.z), (Q.x, Q.y, Q.z))
                else:
                    if point_c_cmp(P, Q) > 0:
                        P,Q = Q,P
                    edge = ((P.x, P.y, P.z), (Q.x, Q.y, Q.z))
                    if not edge in self.seen:
                        self.seen[edge] = edge
                        return edge
        raise StopIteration


cdef class VertexIter:
    """
    A class for iteration over vertices

    EXAMPLES::

        sage: from sage.plot.plot3d.shapes import *
        sage: S = Box(1,2,3)
        sage: len(list(S.vertices())) == 8   # indirect doctest
        True
    """
    def __init__(self, face_set):
        self.set = face_set
        self.i = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.i >= self.set.vcount:
            raise StopIteration
        else:
            self.i += 1
            return (self.set.vs[self.i-1].x, self.set.vs[self.i-1].y, self.set.vs[self.i-1].z)


def len3d(v):
    """
    Return the norm of a vector in three dimensions.

    EXAMPLES::

        sage: from sage.plot.plot3d.index_face_set import len3d
        sage: len3d((1,2,3))
        3.7416573867739413
    """
    return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def sticker(face, width, hover):
    """
    Return a sticker over the given face.
    """
    n = len(face)
    edges = []
    for i from 0 <= i < n:
        edges.append(vector(RDF, [face[i-1][0] - face[i][0],
                                  face[i-1][1] - face[i][1],
                                  face[i-1][2] - face[i][2]]))
    sticker = []
    for i in range(n):
        v = -edges[i]
        w = edges[i - 1]
        N = v.cross_product(w)
        lenN = N.norm()
        dv = v * (width * w.norm() / lenN)
        dw = w * (width * v.norm() / lenN)
        sticker.append(tuple(vector(RDF, face[i-1]) + dv + dw + N*(hover/lenN)))
    return sticker
