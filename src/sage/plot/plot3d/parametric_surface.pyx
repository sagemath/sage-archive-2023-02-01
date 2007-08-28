include "../../ext/stdsage.pxi"
include "../../ext/interrupt.pxi"


cdef inline bint cmp_point_c(point_c P, point_c Q):
    return P.x == Q.x and P.y == Q.y and P.z == Q.z

cdef inline bint smash_edge(point_c* vs, face_c* f, int a, int b):
    if cmp_point_c(vs[f.vertices[a]], vs[f.vertices[b]]):
        f.vertices[b] = f.vertices[a]
        f.n = 3
        return 1
    else:
        return 0

cdef class ParametricSurface(IndexFaceSet):

    def __init__(self, f=None, **kwds):
        self.f = f
        IndexFaceSet.__init__(self, [], [], **kwds)

    def get_grid(self, ds):
        raise NotImplementedError, "You must override the get_grid method."

    def triangulate(self, ds):
        """
        Call self.eval() for all (u,v) in urange \times vrange
        to construct this surface.

        The most complicated part of this code is identifying shared
        vertices and shrinking trivial edges. This is not done so much
        to save memory, rather it is needed so normals of the triangles
        can be calculated correctly.
        """
        urange, vrange = self.get_grid(ds)
        if self.render_grid == (urange, vrange):
            return

        cdef Py_ssize_t i, j
        cdef Py_ssize_t n = len(urange) - 1
        cdef Py_ssize_t m = len(vrange) - 1
        cdef double u, v
        cdef Py_ssize_t ix = 0

        _sig_on
        try:
            self.realloc((m+1)*(n+1), m*n, 4*m*n)
            for u in urange:
                for v in vrange:
                    self.eval_c(&self.vs[ix], u, v)
                    ix += 1
        except:
            _sig_off
            self.fcount = self.vcount = 0
            self.render_grid = None
            raise

        # face_c.vertices:
        #
        #   0 - 1
        #   |   |
        #   3 - 2

        cdef face_c *face, *last_face

        for i from 0 <= i < n:
            for j from 0 <= j < m:
                ix = i*m+j
                face = &self._faces[ix]
                face.n = 4
                face.vertices = &self.face_indices[4*ix]

                # Connect to the i-1 row
                if i == 0:
                    if j == 0:
                        face.vertices[0] = 0
                    else:
                        face.vertices[0] = self._faces[ix-1].vertices[1]
                    face.vertices[1] = j+1
                    smash_edge(self.vs, face, 0, 1)
                else:
                    face.vertices[0] = self._faces[ix-m].vertices[3]
                    face.vertices[1] = self._faces[ix-m].vertices[2]

                # Connect to the j-1 col
                if j == 0:
                    face.vertices[3] = (i+1)*(m+1)
                    smash_edge(self.vs, face, 0, 3)
                else:
                    face.vertices[3] = self._faces[ix-1].vertices[2]

                # This is the newly-seen vertex, identify if its a triangle
                face.vertices[2] = (i+1)*(m+1)+j+1
                smash_edge(self.vs, face, 1, 2) or smash_edge(self.vs, face, 3, 2)

        # Now we see if it wraps around or is otherwise enclosed
        cdef bint enclosed = 1

        cdef face_c *first, *last
        for j from 0 <= j < m:
            first = &self._faces[j]
            last  = &self._faces[(n-1)*m+j]
            if cmp_point_c(self.vs[first.vertices[0]], self.vs[last.vertices[3]]):
                last.vertices[3] = first.vertices[0]
            elif first.vertices[0] != first.vertices[1] and last.vertices[3] != last.vertices[2]:
                enclosed = 0
            if cmp_point_c(self.vs[first.vertices[1]], self.vs[last.vertices[2]]):
                last.vertices[2] = first.vertices[1]
            elif first.vertices[0] != first.vertices[1] and last.vertices[3] != last.vertices[2]:
                enclosed = 0
                enclosed = 0

        for i from 0 <= i < n:
            first = &self._faces[i*m]
            last  = &self._faces[i*m + m-1]
            if cmp_point_c(self.vs[first.vertices[0]], self.vs[last.vertices[1]]):
                last.vertices[1] = first.vertices[0]
            elif first.vertices[0] != first.vertices[3] and last.vertices[1] != last.vertices[2]:
                enclosed = 0
            if cmp_point_c(self.vs[first.vertices[3]], self.vs[last.vertices[2]]):
                last.vertices[2] = first.vertices[3]
            elif first.vertices[0] != first.vertices[3] and last.vertices[1] != last.vertices[2]:
                enclosed = 0

        self.enclosed = enclosed

        # make sure we deleted the correct point from the triangles
        for ix from 0 <= ix < n*m:
            face = &self._faces[ix]
            if face.n == 3:
                if face.vertices[3] == face.vertices[2] or face.vertices[3] == face.vertices[0]:
                    pass
                else:
                    if face.vertices[0] == face.vertices[1]:
                        face.vertices[1] = face.vertices[2]
                    # face.vertices[1] == face.vertices[2]
                    face.vertices[2] = face.vertices[3]

        _sig_off

        self.vcount = (n+1)*(m+1)
        self.fcount = n*m
        self._clean_point_list()

        self.render_grid = urange, vrange


    cdef eval_c(self, point_c *res, double u, double v):
        p = self.eval(u, v)
        res.x, res.y, res.z = tuple(p)

    def eval(self, u, v):
        if self.f is None:
            raise NotImplementedError
        else:
            return self.f(u,v)