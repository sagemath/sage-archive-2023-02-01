r"""
Lattice and reflexive polytopes

This module provides tools for work with lattice and reflexive
polytopes. A \emph{convex polytope} is the convex hull of finitely
many points in $\RR^n$.  The dimension $n$ of a polytope is the
smallest $n$ such that the polytope can be embedded in $\RR^n$. A
\emph{lattice polytope} is a polytope whose vertices have integer
coordiantes.  A \emph{reflexive polytope} is

It uses Package for Analyzing Lattice Polytopes (PALP), which is a
program written in C by Maximilian Kreuzer and Harald Skarke, which is
freely available under the GNU licence terms at
\url{http://tph16.tuwien.ac.at/~kreuzer/CY/}.  Moreover, PALP is
included standard with SAGE.

PALP is described in the paper math.SC/0204356. Its distributive also
contains the application nef.x, which was created by Erwin Riegler and
computes nef partitions and Hodge data for toric complete
intersections.

ACKNOWLEDGMENT: polytope.py module written by William Stein was used
as an example of organizing an interface between an external program
and SAGE.  William Stein also helped Andrey Novoseltsev with debugging
and tuning of this module.

IMPORTANT NOTE: PALP requires some parameters to be determined during
compilation time, i.e., the maximum dimension of polytopes, the
maximum number of points, etc.  These limitations may lead to errors
during calls to different functions of these module. Currently, a
ValueError exception will be raised if the output of poly.x or nef.x
is empty or contains the exclamation mark. The error message will
contain the exact command that caused an error, the description and
vertices of the polytope, and the obtained output.

Data obtained from PALP and some other data is cached and most
returned values are immutable. In particular, you cannot change the
vertices of the polytope or their order after creation of the
polytope.

AUTHORS:
    -- Andrey Novoseltsev (2007-01-11): initial version of this module
    -- Maximilian Kreuzer and Harald Skarke: authors of PALP
    -- Erwin Riegler: the author of nef.x

"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#                     2007 Andrey Novoseltsev <novoselt@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.all import Integer, ZZ, QQ, gcd
from sage.matrix.all import matrix, is_Matrix
from sage.structure.all import Sequence
from sage.structure.sage_object import SageObject
from sage.sets.set import Set_generic
from sage.modules.all import vector
from sage.interfaces.all import maxima
from sage.plot.tachyon import Tachyon
from sage.plot.plot import hue
import random

import os
import StringIO

class SetOfAllLatticePolytopesClass(Set_generic):
    def _repr_(self):
        return "Set of all Lattice Polytopes"

    def __call__(self, x):
        if isinstance(x, LatticePolytopeClass):
            return x
        raise TypeError


SetOfAllLatticePolytopes = SetOfAllLatticePolytopesClass()


def LatticePolytope(data, desc=None, compute_vertices=False,
                    copy_vertices=True, n=0):
    r"""
    Construct a lattice polytope.

    LatticePolytope(data, [desc], [compute_vertices], [copy_vertices], [n])

    INPUT:
        data -- a matrix, whose columns are vertices of the polytope;
                a file with matrix data, open for reading;
                or a filename of such a file.
                See \code{read_palp_matrix} for the file format.
                Points of the given matrix must span the space.
        desc -- (default: "A lattice polytope") description of the polytope.
        compute_vertices -- (default: False) if True, the convex hull of
                the given points will be computed for determining vertices.
                Otherwise, the given points are vertices.
        copy_vertices -- (default: True) if False and \code{data} is a
                matrix of vertices, it will be made immutable.
        n -- (default: 0) if \code{data} is a name of a file, that contains
                data blocks for several polytopes, the n-th block will be used.
                \emph{NUMERATION STARTS WITH ZERO}.

    OUTPUT:
        a lattice polytope

    EXAMPLES:
    Here we construct a polytope from a matrix whose columns are
    vertices in 3-dimensional space. In the first case a copy of the
    given matrix is made during construction, in the second one the
    matrix is made immutable and used as a matrix of vertices.

        sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0],
        ...                   [0, 1, 0,  0, -1,  0],
        ...                   [0, 0, 1,  0,  0, -1]])
        ...
        sage: p = LatticePolytope(m)
        sage: p
        A lattice polytope: 3-dimensional, 6 vertices.
        sage: m.is_mutable()
        True
        sage: m is p.vertices()
        False
        sage: p = LatticePolytope(m, copy_vertices=False)
        sage: m.is_mutable()
        False
        sage: m is p.vertices()
        True

    We draw a pretty picture of the polytype in 3-dimensional space:
        sage: p.show()

    Now we add an extra point to the matrix...
        sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0, 0],
        ...                   [0, 1, 0,  0, -1,  0, 0],
        ...                   [0, 0, 1,  0,  0, -1, 0]])
        ...
        sage: p = LatticePolytope(m, "A lattice polytope with WRONG vertices")
        sage: p
        A lattice polytope with WRONG vertices: 3-dimensional, 7 vertices.

    The above construction is WRONG since the origin is an interior point of the
    octahedron. If you don't know in advance that your data are vertices of the
    polytope, use \code{compute_vertices} option as below:
        sage: p = LatticePolytope(m, "A lattice polytope constructed from 7 points",
        ...                         compute_vertices=True)
        ...
        sage: p
        A lattice polytope constructed from 7 points: 3-dimensional, 6 vertices.

    Points of the given matrix must always span the space:
        sage: m = matrix(ZZ, [[1, 0, -1,  0],
        ...                   [0, 1,  0, -1],
        ...                   [0, 0,  0,  0]])
        ...
        sage: p = LatticePolytope(m)
        Traceback (most recent call last):
        ...
        ValueError: Points must span the space!

    """
    if isinstance(data, LatticePolytopeClass):
        return data
    else:
        return LatticePolytopeClass(data, desc, compute_vertices, copy_vertices, n)


class LatticePolytopeClass(SageObject):

    r"""Class of lattice/reflexive polytopes.

    Use \code{LatticePolytope} for constructing a polytope.

    """

    def __init__(self, data, desc, compute_vertices, copy_vertices=True, n=0):
        r"""
        Construct a lattice polytope. See \code{LatticePolytope}.
        """
        if is_Matrix(data):
            if desc == None:
                self._desc = "A lattice polytope"
            else:
                self._desc = desc
            if data.rank() != data.nrows():
                raise ValueError, "Points must span the space!"
            if compute_vertices:
                self._vertices = data   # for using _poly_x
                self._vertices = read_palp_matrix(self.poly_x("v"))
            else:
                if copy_vertices:
                    self._vertices = data.copy()
                else:
                    self._vertices = data
            self._vertices.set_immutable()

        elif isinstance(data, file) or isinstance(data, StringIO.StringIO):
            m = read_palp_matrix(data)
            self.__init__(m, desc, compute_vertices, False)
        elif isinstance(data,str):
            f = open(data)
            skip_palp_matrix(f, n)
            self.__init__(f, desc, compute_vertices)
            f.close()
        else:
            raise TypeError, \
                "Cannot make a polytope from given data!\nData:\n%s" % data

    def _compute_faces(self):
        r"""Compute and cache faces of this polytope.

        If this polytope is reflexive and the polar polytope was already
        computed, computes faces of both in order to save time and preserve the
        one-to-one correspondence between the faces of this polytope of
        dimension d and the faces of the polar polytope of codimension d+1.

        """
        try:
            if self._constructed_as_polar:
                # "Polar of polar polytope" computed by poly.x may have the
                # order of vertices different from the original polytope. Thus,
                # in order to have consistent enumeration of vertices and faces
                # we must run poly.x on the original polytope.
                self._faces = Sequence([], cr=True)
                for polar_d_faces in reversed(self._polar.faces()):
                    self._faces.append([PolytopeFace(f._facets, f._vertices)
                                                    for f in polar_d_faces])
                self._faces.set_immutable()
                return
        except AttributeError:
            pass
        data = self.poly_x("i").splitlines()[2:]
        d = self.dim()
        v = _read_poly_x_incidences(data[:d])
        f = _read_poly_x_incidences(data[1+d:1+2*d])
        self._faces = Sequence([], cr=True)
        for i in range(len(v)):
            self._faces.append([PolytopeFace(vertices, facets)
                                    for vertices, facets in zip(v[i], f[i])])
        self._faces.set_immutable()

    def _palp(self, command):
        r"""Run \code{command} on vertices of this polytope.

        Returns the output of \code{command} as a string.

        """
        if self.dim() == 0:
            raise ValueError, ("Cannot run \"%s\" for the zero-dimensional "
                + "polytope!\nPolytope: %s") % (command, self)
        stdin, stdout, stderr = os.popen3(command)
        write_palp_matrix(self._vertices, stdin)
        stdin.close()
        err = stderr.read()
        if len(err) > 0:
            raise RuntimeError, ("Error executing \"%s\" for the given polytope!"
                + "\nPolytope: %s\nVertices:\n%s\nOutput:\n%s") % (command,
                self, self.vertices(), err)
        result = stdout.read()
        if len(result) == 0 or result.find("!") != -1:
            raise ValueError, ("Error executing \"%s\" for the given polytope!"
                + "\nPolytope: %s\nVertices:\n%s\nOutput:\n%s") % (command,
                self, self.vertices(), result)
        return result

    def _repr_(self):
        r"""Return a string representation of this polytope."""
        s = str(self._desc)
        s += ": %d-dimensional, %d vertices." % (self.dim(), self.nvertices())
        return s

    def dim(self):
        r"""Return the dimenstion of this polytope.

        EXAMPLES:
        We create a 3-dimensional octahedron and check its dimension:
            sage: o = octahedron(3)
            sage: o.dim()
            3

        """
        return self._vertices.nrows()

    def distances(self):
        r"""Return the matrix of distances for this polytope.

        This matrix \code{m} gives distances \code{m[i,j]} between the
        \code{i}-th facet (which is also the \code{i}-th vertex of the polar
        polytope in the reflexive case) and \code{j}-th point of this polytope.

        EXAMPLES:
        The matrix of distances for a 3-dimensional octahedron:
            sage: o = octahedron(3)
            sage: o.distances()
            [0 0 2 2 2 0 1]
            [2 0 2 0 2 0 1]
            [0 2 2 2 0 0 1]
            [2 2 2 0 0 0 1]
            [0 0 0 2 2 2 1]
            [2 0 0 0 2 2 1]
            [0 2 0 2 0 2 1]
            [2 2 0 0 0 2 1]

        """
        try:
            return self._distances
        except AttributeError:
            if self.is_reflexive():
                self._distances = (self._polar._vertices.transpose()
                                                        * self.points())
                for i in range(self._distances.nrows()):
                    for j in range(self._distances.ncols()):
                        self._distances[i,j] += 1
            else:
                self._distances = self._facet_normals * self.points()
                for i in range(self._distances.nrows()):
                    for j in range(self._distances.ncols()):
                        self._distances[i,j] += self._facet_constants[i]
            self._distances.set_immutable()
            return self._distances

    def face_boundary_points(self, face):
        r"""Return a sequence of indices of boundary lattice points of \code{face}.

        EXAMPLES:
        Boundary lattice points of one of the facets of a 3-dimensional cube:
            sage: o = octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: cube.face_boundary_points(face)
            [0, 1, 2, 3, 11, 15, 21, 25]

        """
        try:
            return face._boundary_points
        except AttributeError:
            face._boundary_points = [i for i in self.face_points(face)
                                    if not i in self.face_interior_points(face)]
            return face._boundary_points

    def face_facets(self, face):
        r"""Return a sequence of indices of facets containing \code{face}.

        EXAMPLES:
        Facets containing one of the edges of a 3-dimensional octahedron:
            sage: o = octahedron(3)
            sage: edge = o.faces(dim=1)[0]
            sage: o.face_facets(edge)
            [0, 1]

        Thus \code{edge} is the intersection of facets 0 and 1:
            sage: edge
            [1, 5]
            sage: o.facets()[0]
            [0, 1, 5]
            sage: o.facets()[1]
            [1, 3, 5]

        """
        return face._facets

    def face_interior_points(self, face):
        r"""Return a sequence of indices of interior lattice points of \code{face}.

        EXAMPLES:
        Interior lattice points of one of the facets of a 3-dimensional cube:
            sage: o = octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: cube.face_interior_points(face)
            [18]

        """
        try:
            return face._interior_points
        except AttributeError:
            non_vertices = self.face_points(face)[self.face_nvertices(face):]
            distances = self.distances()
            other_facets = [i for i in range(self.nfacets())
                                     if not i in face._facets]
            face._interior_points = Sequence([], int)
            for p in non_vertices:
                face._interior_points.append(p)
                for f in other_facets:
                    if distances[f, p] == 0:
                        face._interior_points.pop()
                        break
            face._interior_points.set_immutable()
            return face._interior_points

    def face_nboundary_points(self, face):
        r"""Return the number of boundary lattice points of \code{face}.

        EXAMPLES:
        The number of boundary lattice points of one of the facets of
        a 3-dimensional cube:
            sage: o = octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: cube.face_nboundary_points(face)
            8
        """
        return self.face_npoints(face) - self.face_ninterior_points(face)

    def face_nfacets(self, face):
        r"""Return the number of facets containing \code{face}.

        EXAMPLES:
        The number of facets containing one of the edges of
        a 3-dimensional octahedron:
            sage: o = octahedron(3)
            sage: edge = o.faces(dim=1)[0]
            sage: o.face_nfacets(edge)
            2

        """
        return len(face._facets)

    def face_ninterior_points(self, face):
        r"""Return the number of interior lattice points of \code{face}.

        EXAMPLES:
        The number of interior lattice points of one of the facets of
        a 3-dimensional cube:
            sage: o = octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: cube.face_ninterior_points(face)
            1

        """
        return len(self.face_interior_points(face))

    def face_npoints(self, face):
        r"""Return the number of lattice points of \code{face}.

        EXAMPLES:
        The number of lattice points of one of the facets of
        a 3-dimensional cube:
            sage: o = octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: cube.face_npoints(face)
            9

        """
        return len(self.face_points(face))

    def face_nvertices(self, face):
        r"""Return the number of vertices generating \code{face}.

        EXAMPLES:
        The number of vertices generating one of the facets of
        a 3-dimensional cube:
            sage: o = octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: cube.face_nvertices(face)
            4

        """
        return len(face._vertices)

    def face_points(self, face):
        r"""Return a sequence of indices of lattice points of \code{face}.

        EXAMPLES:
        The lattice points of one of the facets of a 3-dimensional cube:
            sage: o = octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: cube.face_points(face)
            [0, 1, 2, 3, 11, 15, 18, 21, 25]

        """
        try:
            return face._points
        except AttributeError:
            m = self.distances().matrix_from_rows(face._facets)
            cols = m.columns(copy=False)
            face._points = [i for i, col in enumerate(cols) if sum(col) == 0]
            return face._points

    def face_vertices(self, face):
        r"""Return a sequence of indices of vertices generating \code{face}.

        EXAMPLES:
        The vertices generating one of the facets of a 3-dimensional cube:
            sage: o = octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: cube.face_vertices(face)
            [0, 1, 2, 3]

        """
        return face._vertices

    def faces(self, dim=None, codim=None):
        r"""Return the sequence of faces of this polytope.

        If \code{dim} or \code{codim} are specified, returns a sequence of
        faces of the corresponding dimension or codimension. Otherwise returns
        the sequence of such sequences for all dimensions. Each face is
        an instance of class \code{PolytopeFace} and must be used as
        an argument of \code{face_xxx} functions of THE SAME polytope,
        otherwise results will be invalid.

        EXAMPLES:
        All faces of the 3-dimensional octahedron:
            sage: o = octahedron(3)
            sage: o.faces()
            [
            [[5], [1], [0], [3], [4], [2]],
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]],
            [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]
            ]

        Its faces of dimension one (i.e., edges):
            sage: o.faces(dim=1)
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]]

        Its faces of codimension two (also edges):
            sage: o.faces(codim=2)
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]]

        It is an error to specify both dimension and codimension at the same
        time, even if they do agree:
            sage: o.faces(dim=1, codim=2)
            Traceback (most recent call last):
            ...
            ValueError: Both dim and codim are given!

        """
        try:
            if dim == None and codim == None:
                return self._faces
            elif dim != None and codim == None:
                return self._faces[dim]
            elif dim == None and codim != None:
                return self._faces[self.dim()-codim]
            else:
                raise ValueError, "Both dim and codim are given!"
        except AttributeError:
            self._compute_faces()
            return self.faces(dim, codim)

    def facets(self):
        r"""Return the sequence of facets of this polytope.

        Each facet is an instance of class \code{PolytopeFace} and must be used
        as an argument of \code{face_xxx} functions of THE SAME polytope,
        otherwise results will be invalid.

        EXAMPLES:
        All facets of the 3-dimensional octahedron:
            sage: o = octahedron(3)
            sage: o.facets()
            [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]

        Facets are the same as faces of codimension one:
            sage: o.facets() is o.faces(codim=1)
            True

        """
        return self.faces(codim=1)

    def is_reflexive(self):
        r"""Return True if this polytope is reflexive.

        EXAMPLES:
        The 3-dimensional octahedron is reflexive (and 4318 other 3-polytopes):
            sage: o = octahedron(3)
            sage: o.is_reflexive()
            True

        But not all polytopes are reflexive:
            sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0],
            ...                   [0, 1, 0,  0, -1,  0],
            ...                   [0, 0, 0,  0,  0, -1]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.is_reflexive()
            False

        """
        try:
            return self._is_reflexive
        except AttributeError:
            # Determine if the polytope is reflexive by computing vertices
            # of the dual polytope and save all obtained information.
            result = self.poly_x("e")
            self._is_reflexive = (result.find("Vertices of P-dual") != -1)
            if self._is_reflexive:
                self._polar = LatticePolytope(read_palp_matrix(result),
                                "A polytope polar to " + str(self._desc),
                                copy_vertices=False)
                self._polar._is_reflexive = True
                self._polar._constructed_as_polar = True
                self._polar._polar = self
            else:
                normals = []
                constants = []
                d = self.dim()
                for line in result.splitlines()[1:]:
                    numbers = [int(number) for number in line.split()]
                    normals.append(numbers[:d])
                    constants.append(numbers[d])
                self._facet_normals = matrix(ZZ,normals)
                self._facet_constants = vector(ZZ,constants)
                self._facet_normals.set_immutable()
                # SAGE does not support the following as of 2006-12-26
                #self._facet_constants.set_immutable()
            return self._is_reflexive

    def mif(self, partition):
        r"""Return all vectors $m_{i,f}$, grouped into matrices.

        INPUT:
            partition -- NEF-partition (instance of class NEFPartition)

        OUTPUT:
            A sequence of matrices, one for each facet f of this polytope.
            Each row of each matrix corresponds to a part of the NEF-partition.

        EXAMPLES:
        We compute $m_{i,f}$ matrices for one of the nef-partitions of
        the 3-dimensional octahedron:
            sage: o = octahedron(3)
            sage: nefp = o.nef_partitions()[0]
            sage: o.mif(nefp)
            [
            [ 0  0  0]
            [-1 -1  1],
            [ 1  0  0]
            [ 0 -1  1],
            [ 0  1  0]
            [-1  0  1],
            [1 1 0]
            [0 0 1],
            [ 0  0 -1]
            [-1 -1  0],
            [ 1  0 -1]
            [ 0 -1  0],
            [ 0  1 -1]
            [-1  0  0],
            [ 1  1 -1]
            [ 0  0  0]
            ]

        """
        if not isinstance(partition, NEFPartition):
            raise TypeError, "Expect a NEF-partition, got %s!" % type(partition)
        result = Sequence([], cr=True)
        for f in self.facets():
            f_v = self.face_vertices(f)
            a = self.vertices().matrix_from_columns(f_v)
            b = matrix(ZZ,partition.nparts(), len(f_v))
            for j, v in enumerate(f_v):
                b[partition[v], j] = -1
            pivots = a.pivots()
            m = (b.matrix_from_columns(pivots)
                 * a.matrix_from_columns(pivots)**(-1))
            if m*a != b:
                raise RuntimeError, ("Not a nef-partition!\n"
                    + "Partition: %s \nFacet: %s" % (partition, f)
                    + "\na: \n%s \nb: \n%s \nm: \n%s" % (a, b, m))
            result.append(m)
        return result

    def nef_partitions(self):
        r"""Return the sequence of NEF-partitions for this polytope.

        EXAMPLES:
        NEF-partitions of the 3-dimensional octahedron:
            sage: o = octahedron(3)
            sage: o.nef_partitions()
            [
            [1, 1, 0, 0, 0, 1],
            [0, 1, 1, 0, 0, 1],
            [1, 1, 1, 0, 0, 1]
            ]

        NEF-partitions can be computed only for reflexive polytopes:
            sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0],
            ...                   [0, 1, 0,  0, -1,  0],
            ...                   [0, 0, 2,  0,  0, -1]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.nef_partitions()
            Traceback (most recent call last):
            ...
            ValueError: The given polytope is not reflexive!
            Polytope: A lattice polytope: 3-dimensional, 6 vertices.

        """
        if not self.is_reflexive():
            raise ValueError, ("The given polytope is not reflexive!\n"
                                + "Polytope: %s") % self
        try:
            return self._nef_partitions
        except AttributeError:
            result = self.nef_x("-N -Lv -p")
            partitions =_read_nef_x_partitions(result)
            # It is possible that nef.x changed the order of vertices
            nef_vertices = read_palp_matrix(result[result.find("\n")+1:])
            if self.vertices() != nef_vertices:
                trans = [self.vertices().columns().index(v)
                            for v in nef_vertices.columns()]
                for i, p in enumerate(partitions):
                    partitions[i] = [trans[v] for v in p]
            # Convert to the input format of nef_partition class
            for i, p in enumerate(partitions):
                new_p = [1]*self.nvertices()
                for v in p:
                    new_p[v] = 0
                partitions[i] = new_p
            self._nef_partitions = Sequence(partitions,NEFPartition,cr=True)
            self._nef_partitions.set_immutable()
            return self._nef_partitions

    def nef_x(self, keys):
        r"""Run nef.x wilth given \code{keys} on vertices of this polytope.

        INPUT:
            keys -- a string of options passed to nef.x. The key "-f" is added
                    automatically.

        OUTPUT:
            the output of nef.x as a string.

        EXAMPLES:
        This call is used internally for computing NEF-partitions:
            sage: o = octahedron(3)
            sage: print o.nef_x("-N -Lv -p")    # output contains random time
            M:27 8 N:7 6  #part=5
            3 6 Vertices in N-lattice:
                0    0    0    1   -1    0
                0    0    1    0    0   -1
               -1    1    0    0    0    0
            ------------------------------
                1    1    0    0    0    0  d=2  codim=2
                0    0    1    0    0    1  d=2  codim=2
                0    0    0    1    1    0  d=2  codim=2
             P:0 V:1 4 5   (1 1) (1 1) (1 1)     0sec  0cpu
             P:1 V:3 4 5   (0 2) (1 1) (2 0)     0sec  0cpu
             P:3 V:4 5   (0 2) (1 1) (1 1)     0sec  0cpu
            np=3 d:1 p:1    0sec     0cpu

        """
        return self._palp("nef.x -f " + keys)

    def nfacets(self):
        r"""Return the number of facets of this polytope.

        EXAMPLES:
        The number of facets of the 3-dimensional octahedron:
            sage: o = octahedron(3)
            sage: o.nfacets()
            8

        """
        return len(self.facets())

    def npoints(self):
        r"""Return the number of lattice points of this polytope.

        EXAMPLES:
        The number of lattice points of the 3-dimensional octahedron and its
        polar cube:
            sage: o = octahedron(3)
            sage: o.npoints()
            7
            sage: cube = o.polar()
            sage: cube.npoints()
            27

        """
        return self.points().ncols()

    def nvertices(self):
        r"""Return the number of vertices of this polytope.

        EXAMPLES:
        The number of vertices of the 3-dimensional octahedron and its
        polar cube:
            sage: o = octahedron(3)
            sage: o.nvertices()
            6
            sage: cube = o.polar()
            sage: cube.nvertices()
            8

        """
        return self._vertices.ncols()

    def parent(self):
        """
        Return the set of all lattice polytopes.

        EXAMPLES:
            sage:
        """
        return SetOfAllLatticePolytopes

    def plot(self, camera_center=None):
        """
        Draw a 3d picture of the polytope, which is assumed to be the
        convex hole of its vertices.

        INPUT:
            self -- polytope of dimension 3.
            camera_center -- (default: random) location of center of
            the camera (i.e., viewer)

        OUTPUT:
            -- a tachyon 3d raytracer plot of the polytope

        The face colors are random.

        AUTHORS:
            -- William Stein and Tom Boothby

        EXAMPLES:
            sage:
        """
        if self.dim() != 3:
            raise ValueError, "polytope must have dimension 3"
        m = self._vertices
        r = random.random
        if camera_center is None:
            x = max(m.list())
            camera_center = [(1+r())*x,(1+r())*x,(1+r())*x]
            for i in range(3):
                if r() > 0.5:
                    camera_center[i] *= -1
        t = Tachyon(camera_center=camera_center)
        t.light(tuple([2*v for v in camera_center]), .1, (1,1,1))

        cols = m.columns()
        n = len(cols)
        for i in range(n):
            for j in range(i,n):
                for k in range(j,n):
                    s = "T%010d%010d%010d"%(i,j,k)
                    t.texture(s, color=hue(r()))
                    t.triangle(cols[i],cols[j],cols[k], s)
        return t

    def show(self, camera_center=None):
        """
        Draw a 3d picture of the polytope, which is assumed to be the
        convex hole of its vertices.

        See self.plot? for more details.
        """

        self.plot(camera_center=camera_center).show()

    def points(self):
        r"""
        Return all lattice points of this polytope as columns of a matrix.

        EXAMPLES:
        The lattice points of the 3-dimensional octahedron and its polar cube:
            sage: o = octahedron(3)
            sage: o.points()
            [ 1  0  0 -1  0  0  0]
            [ 0  1  0  0 -1  0  0]
            [ 0  0  1  0  0 -1  0]
            sage: cube = o.polar()
            sage: cube.points()
            [-1  1 -1  1 -1  1 -1  1 -1 -1 -1 -1 -1  0  0  0  0  0  0  0  0  0  1  1  1  1  1]
            [-1 -1  1  1 -1 -1  1  1 -1  0  0  0  1 -1 -1 -1  0  0  0  1  1  1 -1  0  0  0  1]
            [ 1  1  1  1 -1 -1 -1 -1  0 -1  0  1  0 -1  0  1 -1  0  1 -1  0  1  0 -1  0  1  0]

        """
        try:
            return self._points
        except AttributeError:
            self._points = read_palp_matrix(self.poly_x("p"))
            return self._points

    def polar(self):
        r"""
        Return the polar polytope, if this polytope is reflexive.

        EXAMPLES:
        The polar polytope to the 3-dimensional octahedron:
            sage: o = octahedron(3)
            sage: cube = o.polar()
            sage: cube
            A polytope polar to An octahedron: 3-dimensional, 8 vertices.

        The polar polytope "remembers" the original one:
            sage: cube.polar()
            An octahedron: 3-dimensional, 6 vertices.
            sage: cube.polar().polar() is cube
            True

        Only reflexive polytopes have polars:
            sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0],
            ...                   [0, 1, 0,  0, -1,  0],
            ...                   [0, 0, 2,  0,  0, -1]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.polar()
            Traceback (most recent call last):
            ...
            ValueError: The given polytope is not reflexive!
            Polytope: A lattice polytope: 3-dimensional, 6 vertices.

        """
        if self.is_reflexive():
            return self._polar
        else:
            raise ValueError, ("The given polytope is not reflexive!\n"
                                + "Polytope: %s") % self

    def poly_x(self, keys):
        r"""
        Run poly.x with given \code{keys} on vertices of this polytope.

        INPUT:
            keys -- a string of options passed to poly.x. The key "f" is added
                    automatically.

        OUTPUT:
            the output of poly.x as a string.

        EXAMPLES:
        This call is used for determining if a polytope is reflexive or not:
            sage: o = octahedron(3)
            sage: print o.poly_x("e")
            8 3  Vertices of P-dual <-> Equations of P
              -1  -1   1
               1  -1   1
              -1   1   1
               1   1   1
              -1  -1  -1
               1  -1  -1
              -1   1  -1
               1   1  -1

        Since PALP has limits on different parameters determined during
        compilation, the following code is likely to fail, unless you change
        default settings of PALP:
            sage: BIGO = octahedron(7)
            sage: BIGO
            An octahedron: 7-dimensional, 14 vertices.
            sage: BIGO.poly_x("e")      # random output depending on your system
            Traceback (most recent call last):
            ...
            ValueError: Error executing "poly.x -fe" for the given polytope!
            Polytope: An octahedron: 7-dimensional, 14 vertices.
            Vertices:
            [ 1  0  0  0  0  0  0 -1  0  0  0  0  0  0]
            [ 0  1  0  0  0  0  0  0 -1  0  0  0  0  0]
            [ 0  0  1  0  0  0  0  0  0 -1  0  0  0  0]
            [ 0  0  0  1  0  0  0  0  0  0 -1  0  0  0]
            [ 0  0  0  0  1  0  0  0  0  0  0 -1  0  0]
            [ 0  0  0  0  0  1  0  0  0  0  0  0 -1  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0 -1]
            Output:
            increase POLY_Dmax!

        """
        return self._palp("poly.x -f" + keys)

    def vertices(self):
        r"""
        Return vertices of this polytope as columns of a matrix.

        EXAMPLES:
        The lattice points of the 3-dimensional octahedron and its polar cube:
            sage: o = octahedron(3)
            sage: o.vertices()
            [ 1  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: cube = o.polar()
            sage: cube.vertices()
            [-1  1 -1  1 -1  1 -1  1]
            [-1 -1  1  1 -1 -1  1  1]
            [ 1  1  1  1 -1 -1 -1 -1]

        """
        return self._vertices


class NEFPartition(Sequence):

    r"""
    Construct a NEF-partition from the given list.

    A NEF-partition with $k$ parts, $V = V_0 \cap V_1 \cap \dots \cap V_{k-1}$,
    is represented by a single list of lenght $n = #V$, in which the $i$-th
    entry is the part number of the $i$-th vertex of a polytope.

    NOTE THAT NUMERATON OF PARTS STARTS WITH ZERO.

    EXAMPLES:
    All elements of the list will be coerced to integers, so it is OK to use
    either a list of numbers or a list of strings:
        sage: NEFPartition([1, 1, 0, 0, 0, 1])
        [1, 1, 0, 0, 0, 1]
        sage: NEFPartition(['1', '1', '0', '0', '0', '1'])
        [1, 1, 0, 0, 0, 1]

    """

    def __init__(self, data):
        r"""Construct a NEF-partition."""
        Sequence.__init__(self, data, int)
        self._n = max(self) + 1
        self.set_immutable()

    def nparts(self):
        r"""
        Return the number of parts of this partitions.

        EXAMPLES:
            sage: nefp = NEFPartition([1, 1, 0, 0, 0, 1])
            sage: nefp.nparts()
            2

        """
        return self._n

    def part(self, i):
        r"""
        Return the \code{i}-th part of the partition.

        NUMERATON OF PARTS STARTS WITH ZERO.

        EXAMPLES:
            sage: nefp = NEFPartition([1, 1, 0, 0, 0, 1])
            sage: nefp.part(0)
            [2, 3, 4]
            sage: nefp.part(1)
            [0, 1, 5]

        """
        return [j for j, el in enumerate(self) if el == i]

    def part_of_vertex(self, i):
        r"""
        Return the index of the part containing the \code{i}-vertex.

        EXAMPLES:
        \code{nefp.part_of_vertex(i)} is equivalent to \code{nefp[i]}.
            sage: nefp = NEFPartition([1, 1, 0, 0, 0, 1])
            sage: nefp.part_of_vertex(3)
            0
            sage: nefp.part_of_vertex(5)
            1
            sage: nefp[3]
            0
            sage: nefp[5]
            1

        You cannot change a NEF-partition once it is constructed:
            sage: nefp[3] = 1
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.

        """
        return self[i]


class PolytopeFace(SageObject):

    r"""
    PolytopeFace(vertices, facets)

    Construct a polytope face.

    POLYTOPE FACES SHOULD NOT BE CONSTRUCTED OUTSIDE OF LATTICE POLYTOPES!

    INPUT:
        vertices -- a sequence of indices of generating vertices.
        facets -- a sequence of indices of facets containing this face.

    """

    def __init__(self, vertices, facets):
        r"""Construct a face."""
        self._vertices = vertices
        self._facets = facets

    def _repr_(self):
        r"""Return a string representation of this face."""
        return str(self._vertices)


def _create_octahedron(dim):
    r"""
    Create an octahedron of the given dimension."""
    m = matrix(ZZ, dim, 2*dim)
    for i in range(dim):
        m[i,i] = 1
        m[i,dim+i] = -1
    return LatticePolytope(m, "An octahedron", copy_vertices=False)


_octahedrons = dict()       # Dictionary for storing created octahedrons


def _read_nef_x_partitions(data):
    r"""
    Read all nef-partitions for one polytope from a string or an open file.

    \code{data} should be an output of nef.x.

    Returns the sequence of NEF-partitions. Each NEF-partition is given as
    a sequence of integers.

    If there are no NEF-partitions, returns the empty sequence.
    If the string is empty or EOF is reached, returns None.

    """
    if isinstance(data, str):
        f = StringIO.StringIO(data)
        partitions = _read_nef_x_partitions(f)
        f.close()
        return partitions
    # If data is not a string, try to treat it as a file.
    line = data.readline()
    if line == "":
        return None
    partitions = []
    while len(line) > 0 and line.find("np=") == -1:
        if line.find("V:") == -1:
            line = data.readline()
            continue
        start = line.find("V:") + 2
        end = line.find("  ", start)  # Find DOUBLE space
        partitions.append(Sequence(line[start:end].split(),int))
        line = data.readline()
    # Compare the number of found partitions with np in data.
    start = line.find("np=")
    if start != -1:
        start += 3
        end = line.find(" ", start)
        np = int(line[start:end])
        if np != len(partitions):
            raise ValueError, ("Found %d partitions, expected %d!" %
                                 (len(partitions), np))
    else:
        raise ValueError, "Wrong data format, cannot find \"np=\"!"
    partitions = Sequence(partitions)
    return partitions


def _read_poly_x_incidences(data):
    r"""Convert incidence data from binary numbers to sequences.

    INPUT:
        data -- list of lines. Each line contains incidence information for
            all faces of one dimension, the first word of each line is
            considered as a comment and is dropped.

    OUTPUT:
        a sequence F, such that F[d][i] is a sequence of vertices or facets
        corresponding to the i-th d-dimensional face.

    """
    n = len(data[0].split()[1])     # Number of vertices or facets
    result = Sequence([], cr=True)
    for d in data:
        subr = Sequence([])
        for el in d.split()[1:]:
            f = Sequence([j for j in range(n) if el[n-1-j] == '1'])
            f.set_immutable()
            subr.append(f)
        subr.set_immutable()
        result.append(subr)
    result.set_immutable()
    return result


def filter_polytopes(f, polytopes, subseq=None, print_numbers=False):
    r"""Use the function \code{f} to filter polytopes in a list.

    INPUT:
        f -- filtering function, it must take one argument, a lattice polytope,
            and return True or False.
        polytopes -- list of polytopes.
        subseq -- (default: None) list of integers. If it is specified, only
            polytopes with these numbers will be considered.
        print_numbers -- (default: False) if True, the number of the current
            polytope will be printed on the screen before calling \code{f}.

    OUTPUT:
        a list of integers --- numbers of polytopes in the given list, that
        satisfy the given condition (i.e. function \code{f} returns True)
        and are elements of subseq, if it is given.

    EXAMPLES:
    Consider a sequence of octahedrons:
        sage: polytopes = Sequence([octahedron(n) for n in range(2, 7)], cr=True)
        sage: polytopes
        [
        An octahedron: 2-dimensional, 4 vertices.,
        An octahedron: 3-dimensional, 6 vertices.,
        An octahedron: 4-dimensional, 8 vertices.,
        An octahedron: 5-dimensional, 10 vertices.,
        An octahedron: 6-dimensional, 12 vertices.
        ]

    This filters octahedrons of dimension at least 4:
        sage: filter_polytopes(lambda p: p.dim() >= 4, polytopes)
        [2, 3, 4]

    For long tests you can see the current progress:
        sage: filter_polytopes(lambda p: p.nvertices() >= 10, polytopes, \
        sage: print_numbers=True)
        0
        1
        2
        3
        4
        [3, 4]

    Here we consider only some of the polytopes:
        sage: filter_polytopes(lambda p: p.nvertices() >= 10, polytopes, \
        sage: [2, 3, 4], print_numbers=True)
        2
        3
        4
        [3, 4]

    """
    if subseq == []:
        return []
    elif subseq == None:
        subseq = range(len(polytopes))
    result = []
    for n in subseq:
        if print_numbers:
            print n
            os.sys.stdout.flush()
        if f(polytopes[n]):
            result.append(n)
    return result


def octahedron(dim):
    r"""Return an octahedron of the given dimension.

    EXAMPLES:
    Here are 3- and 4-dimensional octahedrons:
        sage: o = octahedron(3)
        sage: o
        An octahedron: 3-dimensional, 6 vertices.
        sage: o.vertices()
        [ 1  0  0 -1  0  0]
        [ 0  1  0  0 -1  0]
        [ 0  0  1  0  0 -1]
        sage: o = octahedron(4)
        sage: o
        An octahedron: 4-dimensional, 8 vertices.
        sage: o.vertices()
        [ 1  0  0  0 -1  0  0  0]
        [ 0  1  0  0  0 -1  0  0]
        [ 0  0  1  0  0  0 -1  0]
        [ 0  0  0  1  0  0  0 -1]

    There exists only one octahedron of each dimension:
        sage: o is octahedron(4)
        True

    """
    if _octahedrons.has_key(dim):
        return _octahedrons[dim]
    else:
        _octahedrons[dim] = _create_octahedron(dim)
        return _octahedrons[dim]


def positive_integer_relations(points):
    r"""Return relations between given points.

    INPUT:
        points -- lattice points given as columns of a matrix

    OUTPUT:
        matrix of relations between given points with non-negative integer
        coefficients

    EXAMPLES:
    This is a 3-dimensional reflexive polytope:
        sage: m = matrix(ZZ,[[1, 0, -1, 0, -1],
        ...                  [0, 1, -1, 0,  0],
        ...                  [0, 0,  0, 1, -1]])
        ...
        sage: p = LatticePolytope(m)
        sage: p.points()
        [ 1  0 -1  0 -1  0]
        [ 0  1 -1  0  0  0]
        [ 0  0  0  1 -1  0]

    We can compute linear relations between its points in the following way:
        sage: p.points().transpose().kernel().echelonized_basis_matrix()
        [ 1  0  0  1  1  0]
        [ 0  1  1 -1 -1  0]
        [ 0  0  0  0  0  1]

    However, the above relations may contain negative and rational numbers. This
    function transforms them in such a way, that all coefficients are
    non-negative integers:
        sage: positive_integer_relations(p.points())
        [1 0 0 1 1 0]
        [1 1 1 0 0 0]
        [0 0 0 0 0 1]

    """
    points = points.transpose().base_extend(QQ)
    relations = points.kernel().echelonized_basis_matrix()
    nonpivots = relations.nonpivots()
    nonpivot_relations = relations.matrix_from_columns(nonpivots)
    n_nonpivots = len(nonpivots)
    n = nonpivot_relations.nrows()
    a = matrix(QQ,n_nonpivots,n_nonpivots)
    for i in range(n_nonpivots):
        a[i, i] = -1
    a = nonpivot_relations.stack(a).transpose()
    a = sage_matrix_to_maxima(a)
    maxima.load("simplex")

    new_relations = []
    for i in range(n_nonpivots):
        # Find a non-negative linear combination of relations,
        # such that all components are non-negative and the i-th one is 1
        b = [0]*i + [1] + [0]*(n_nonpivots - i - 1)
        c = [0]*(n+i) + [1] + [0]*(n_nonpivots - i - 1)
        x = maxima.linear_program(a, b, c)
        if x.str() == r'?Problem\not\feasible\!':
            raise ValueError, "cannot find required relations"
        x = x.sage()[0][:n]
        v = relations.linear_combination_of_rows(x)
        new_relations.append(v)

    relations = relations.stack(matrix(QQ, new_relations))
    # Use the new relation to remove negative entries in non-pivot columns
    for i in range(n_nonpivots):
        for j in range(n):
            coef = relations[j,nonpivots[i]]
            if coef < 0:
                relations.add_multiple_of_row(j, n+i, -coef)
    # Get a new basis
    relations = relations.matrix_from_rows(relations.transpose().pivots())
    # Switch to integers
    for i in range(n):
        relations.rescale_row(i, 1/gcd(relations[i].list()))
    return relations.change_ring(ZZ)


def projective_space(dim):
    r"""Return a simplex of the given dimension, corresponding to $P_{dim}$.

    EXAMPLES:
    We construct 3- and 4-dimensional simplexes:
        sage: p = projective_space(3)
        sage: p
        A simplex: 3-dimensional, 4 vertices.
        sage: p.vertices()
        [ 1  0  0 -1]
        [ 0  1  0 -1]
        [ 0  0  1 -1]
        sage: p = projective_space(4)
        sage: p
        A simplex: 4-dimensional, 5 vertices.
        sage: p.vertices()
        [ 1  0  0  0 -1]
        [ 0  1  0  0 -1]
        [ 0  0  1  0 -1]
        [ 0  0  0  1 -1]

    """
    m = matrix(ZZ, dim, dim+1)
    for i in range(dim):
        m[i,i] = 1
        m[i,dim] = -1
    return LatticePolytope(m, "A simplex", copy_vertices=False)


def read_all_polytopes(file_name, desc=None):
    r"""Read all polytopes from the given file.

    INPUT:
        file_name -- the name of a file with vertices of polytopes
        desc -- a string, that will be used for creating polytope descriptions.
            By default it will be set to 'A lattice polytope #%d from "filename"'
            and will be used as \code{desc % n} where \code{n} is the number of
            the polytope in the file (STARTING WITH ZERO).

    OUTPUT:
        a sequence of polytopes

    """
    if desc == None:
        desc = r'A lattice polytope #%d from "'+file_name+'"'
    polytopes = Sequence([], LatticePolytope, cr=True)
    f = open(file_name)
    n = 0
    p = LatticePolytope(f, desc % n)
    while p.dim() != 0:
        polytopes.append(p)
        n += 1
        p = LatticePolytope(f, desc % n)
    f.close()
    return polytopes


def read_palp_matrix(data):
    r"""
    Read and return an integer matrix from a string or an opened file.

    First input line must start with two integers m and n, the number of rows
    and columns of the matrix. The rest of the first line is ignored. The next
    m lines must contain n numbers each.

    If m>n, returns the transposed matrix.
    If the string is empty or EOF is reached, returns the empty matrix,
    constructed by \code{matrix()}.

    EXAMPLES:
        sage: read_palp_matrix("2 3 comment \n 1 2 3 \n 4 5 6")
        [1 2 3]
        [4 5 6]
        sage: read_palp_matrix("3 2 Will be transposed \n 1 2 \n 3 4 \n 5 6")
        [1 3 5]
        [2 4 6]

    """
    if isinstance(data,str):
        f = StringIO.StringIO(data)
        mat = read_palp_matrix(f)
        f.close()
        return mat
    # If data is not a string, try to treat it as a file.
    line = data.readline()
    if line == "":
        return matrix()
    line = line.split()
    m = int(line[0])
    n = int(line[1])
    seq = []
    for i in range(m):
        seq.extend(int(el) for el in data.readline().split())
    mat = matrix(ZZ,m,n,seq)
    if m <= n:
        return mat
    else:
        return mat.transpose()


def sage_matrix_to_maxima(m):
    r"""Convert a SAGE matrix to the string representation of Maxima."""
    return maxima("matrix("+",".join(str(v.list()) for v in m.rows())+")")


def skip_palp_matrix(data, n=1):
    r"""
    Skip matrix data in a file.

    INPUT:
        data -- opened file with blocks of matrix data in the following format:
                A block consisting of m+1 lines has the number m as the first
                element of its first line.
        n -- (default: 1) integer, specifies how many blocks should be skipped

    If EOF is reached during the process, raises ValueError exception.

    """
    for i in range(n):
        line = data.readline()
        if line == "":
            raise ValueError, "There are not enough data to skip!"
        for j in range(int(line.split()[0])):
            if data.readline() == "":
                raise ValueError, "There are not enough data to skip!"


def write_palp_matrix(m, ofile=os.sys.stdout, comment="", format=None):
    r"""Write a matrix into a file.

    INPUT:
        m -- a matrix over integers.
        ofile -- a file opened for writing (default: stdout)
        comment -- a string (default: empty) see output description
        format -- a format string used to print matrix entries. By default,
                  "%nd" will be used, where n is the maximum entry length.

    OUTPUT:
        First line: number_of_rows number_of_columns comment
        Next number_of_rows lines: rows of the matrix.

    EXAMPLES:
    This functions is used for writing polytope vertices in PALP format:
        sage: o = octahedron(3)
        sage: write_palp_matrix(o.vertices(), comment="3D Octahedron")
        3 6 3D Octahedron
         1  0  0 -1  0  0
         0  1  0  0 -1  0
         0  0  1  0  0 -1
        sage: write_palp_matrix(o.vertices(), format="%4d")
        3 6
           1    0    0   -1    0    0
           0    1    0    0   -1    0
           0    0    1    0    0   -1

    """
    if format == None:
        n = max(len(str(m[i,j]))
                for i in range(m.nrows()) for j in range(m.ncols()))
        format = "%" + str(n) + "d"
    ofile.write("%d %d %s\n" % (m.nrows(),m.ncols(),comment))
    for i in range(m.nrows()):
        ofile.write(" ".join(format % m[i,j] for j in range(m.ncols()))+"\n")
