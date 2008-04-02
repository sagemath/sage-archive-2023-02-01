r"""
Lattice and reflexive polytopes

This module provides tools for work with lattice and reflexive
polytopes. A \emph{convex polytope} is the convex hull of finitely
many points in $\RR^n$.  The dimension $n$ of a polytope is the
smallest $n$ such that the polytope can be embedded in $\RR^n$.

A \emph{lattice polytope} is a polytope whose vertices all have
integer coordinates.

If $L$ is a lattice polytope, the dual polytope of $L$ is
$$
  \{y \in \Z^n :   x\cdot y \geq -1 \text{ all } x \in L\}
$$
A \emph{reflexive polytope} is a lattice polytope, such that its polar
is also a lattice polytope, i.e. has vertices with integer coordinates.

This SAGE module uses Package for Analyzing Lattice Polytopes (PALP),
which is a program written in C by Maximilian Kreuzer and Harald
Skarke, which is freely available under the GNU licence terms at
\url{http://tph16.tuwien.ac.at/~kreuzer/CY/}.  Moreover, PALP is
included standard with SAGE.

PALP is described in the paper math.SC/0204356. Its distribution also
contains the application nef.x, which was created by Erwin Riegler and
computes nef partitions and Hodge data for toric complete
intersections.

ACKNOWLEDGMENT: polytope.py module written by William Stein was used
as an example of organizing an interface between an external program
and SAGE.  William Stein also helped Andrey Novoseltsev with debugging
and tuning of this module.

Robert Bradshaw helped Andrey Novoseltsev to realize plot3d function.

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

If you are going to  work with large sets of data, take a look at all_*
functions in this module. They precompute different data for sequences
of polynomials with a few runs of external programs. This can significantly
affect the time of future computations. You can also use dump/load, but not
all data will be stored (currently only faces and the number of their internal
and boundary points are stored, in addition to polytope vertices and its polar).

AUTHORS:
    -- Andrey Novoseltsev (2007-01-11): initial version of this module
    -- Andrey Novoseltsev (2007-01-15): all_* functions
    -- Andrey Novoseltsev (2008-04-01): second version, including:
            - dual NEF-partitions and necessary convex_hull and minkowski_sum
            - built-in sequences of 2- and 3-dimensional reflexive polytopes
            - plot3d, skeleton_show

    -- Maximilian Kreuzer and Harald Skarke: authors of PALP (which was also
            used to obtain the list of 3-dimensional reflexive polytopes)
    -- Erwin Riegler: the author of nef.x
"""

#*****************************************************************************
#       Copyright (C) 2007-2008 William Stein <wstein@gmail.com>
#                     2007-2008 Andrey Novoseltsev <novoselt@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.databases.db import DB_HOME
from sage.graphs.all import Graph
from sage.interfaces.all import maxima
from sage.matrix.all import matrix, is_Matrix
from sage.misc.all import tmp_filename
from sage.modules.all import vector
from sage.plot.plot import hue
from sage.plot.plot3d.index_face_set import IndexFaceSet
from sage.plot.plot3d.all import line3d, point3d
from sage.plot.plot3d.shapes2 import text3d
from sage.plot.tachyon import Tachyon
from sage.rings.all import Integer, ZZ, QQ, gcd
from sage.sets.set import Set_generic
from sage.structure.all import Sequence
from sage.structure.sage_object import SageObject, load, loads

import copy_reg
import os
import pickle
import sage.misc.prandom as random
import StringIO
import zlib


data_location = DB_HOME + '/reflexive_polytopes/'


class SetOfAllLatticePolytopesClass(Set_generic):
    def _repr_(self):
        r"""
        Return a string representation.

        TESTS:
            sage: lattice_polytope.SetOfAllLatticePolytopesClass()._repr_()
            'Set of all Lattice Polytopes'
        """
        return "Set of all Lattice Polytopes"

    def __call__(self, x):
        r"""
        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: lattice_polytope.SetOfAllLatticePolytopesClass().__call__(o)
            An octahedron: 3-dimensional, 6 vertices.
        """
        if isinstance(x, LatticePolytopeClass):
            return x
        raise TypeError


SetOfAllLatticePolytopes = SetOfAllLatticePolytopesClass()


def LatticePolytope(data, desc=None, compute_vertices=True,
                    copy_vertices=True, n=0):
    r"""
    Construct a lattice polytope.

    LatticePolytope(data, [desc], [compute_vertices], [copy_vertices], [n])

    INPUT:
        data -- matrix whose columns are vertices of the polytope
                (unless \code{compute_vertices} is True);
                a file with matrix data, open for reading;
                or a filename of such a file.
                See \code{read_palp_matrix} for the file format.
                Points of the given matrix must span the space.
        desc -- (default: "A lattice polytope") description of the polytope.
        compute_vertices -- (default: True) if True, the convex hull of
                the given points will be computed for determining vertices.
                Otherwise, the given points are vertices.
        copy_vertices -- (default: True) if False, and compute_vertices is False,
                and \code{data} is a matrix of vertices, it will be made immutable.
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
        sage: p = LatticePolytope(m, compute_vertices=False, copy_vertices=False)
        sage: m.is_mutable()
        False
        sage: m is p.vertices()
        True

    We draw a pretty picture of the polytype in 3-dimensional space:
        sage: p.plot().show()

    Now we add an extra point, which is in the interiour of the polytope...
        sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0, 0],
        ...                   [0, 1, 0,  0, -1,  0, 0],
        ...                   [0, 0, 1,  0,  0, -1, 0]])
        ...
        sage: p = LatticePolytope(m, "A lattice polytope constructed from 7 points")
        sage: p
        A lattice polytope constructed from 7 points: 3-dimensional, 6 vertices.

    You can suppress vertex computation for speed but this can lead to mistakes:
        sage: p = LatticePolytope(m, "A lattice polytope with WRONG vertices",
        ...                         compute_vertices=False)
        ...
        sage: p
        A lattice polytope with WRONG vertices: 3-dimensional, 7 vertices.

    Points of the given matrix must always span the space, this conditions will
    be checked unless you specify \code{compute_vertices=False} option:
        sage: m = matrix(ZZ, [[1, 0, -1,  0],
        ...                   [0, 1,  0, -1],
        ...                   [0, 0,  0,  0]])
        ...
        sage: p = LatticePolytope(m)
        Traceback (most recent call last):
        ...
        ValueError: Points must span the space!
        Given:
        [ 1  0 -1  0]
        [ 0  1  0 -1]
        [ 0  0  0  0]
    """
    if isinstance(data, LatticePolytopeClass):
        return data
    else:
        return LatticePolytopeClass(data, desc, compute_vertices, copy_vertices, n)

copy_reg.constructor(LatticePolytope)   # "safe for unpickling"

def ReflexivePolytope(dim, n):
    r"""
    Return n-th reflexive polytope from the database of 2- or 3-dimensional
    reflexive polytopes.

    NOTES:
    1. Numeration starts with zero: 0<=n<=15 for dim=2 and 0<=n<=4318 for dim=3.
    2. During the first call, all reflexive polytopes of requested dimension are
    loaded and cached for future use, so the first call for 3-dimensional
    polytopes can take several seconds, but all consecutive calls are fast.
    3. Equivalent to ReflexivePolytopes(dim)[n] but checks bounds first.

    EXAMPLES:
    The 3rd 2-dimensional polytope is "the diamond:"
      sage: ReflexivePolytope(2,3)
      Reflexive polygon  3: 2-dimensional, 4 vertices.
      sage: lattice_polytope.ReflexivePolytope(2,3).vertices()
      [ 1  0  0 -1]
      [ 0  1 -1  0]

    There are 16 reflexive polygons and numeration starts with 0:
      sage: ReflexivePolytope(2,16)
      Traceback (most recent call last):
      ...
      ValueError: there are only 16 reflexive polygons!

    It is not possible to load a 4-dimensional polytope in this way:
      sage: ReflexivePolytope(4,16)
      Traceback (most recent call last):
      ...
      NotImplementedError: only 2- and 3-dimensional reflexive polytopes are available!

    """

    if dim == 2:
        if n > 15:
            raise ValueError, "there are only 16 reflexive polygons!"
        return ReflexivePolytopes(2)[n]
    elif dim == 3:
        if n > 4318:
            raise ValueError, "there are only 4319 reflexive 3-polytopes!"
        return ReflexivePolytopes(3)[n]
    else:
        raise NotImplementedError, "only 2- and 3-dimensional reflexive polytopes are available!"

# Sequences of reflexive polytopes
_rp2d = None
_rp3d = None

def ReflexivePolytopes(dim):
    r"""
    Return the sequence of all 2- or 3-dimensional reflexive polytopes.

    NOTE: During the first call the database is loaded and cached for future use,
    so repetitive calls will return the same object in memory.

    TODO: load a database with precomputed faces. Requires better pickling,
    otherwise takes forever!

    EXAMPLES:
    There are 16 reflexive polygons:
      sage: len(ReflexivePolytopes(2))
      16

    It is not possible to load 4-dimensional polytopes in this way:
      sage: ReflexivePolytopes(4)
      Traceback (most recent call last):
      ...
      NotImplementedError: only 2- and 3-dimensional reflexive polytopes are available!

    """
    global _rp2d, _rp3d
    if dim == 2:
        if _rp2d == None:
            _rp2d = load(data_location+"reflexive_polytopes_2d")
        return _rp2d
    elif dim == 3:
        if _rp3d == None:
            _rp3d = load(data_location+"reflexive_polytopes_3d")
        return _rp3d
    else:
        raise NotImplementedError, "only 2- and 3-dimensional reflexive polytopes are available!"


class LatticePolytopeClass(SageObject):
    r"""
    Class of lattice/reflexive polytopes.

    Use \code{LatticePolytope} for constructing a polytope.
    """

    def __init__(self, data, desc, compute_vertices, copy_vertices=True, n=0):
        r"""
        Construct a lattice polytope. See \code{LatticePolytope}.

        Most tests/examples are also done in LatticePolytope.

        TESTS:
            sage: LatticePolytope(matrix(ZZ,[[1,2,3],[4,5,6]]))
            Traceback (most recent call last):
            ...
            ValueError: Points must span the space!
            Given:
            [1 2 3]
            [4 5 6]
        """
        if is_Matrix(data):
            if desc == None:
                self._desc = "A lattice polytope"
            else:
                self._desc = desc
            if compute_vertices:
                shifted = matrix(ZZ,[data.column(i, from_list=True) - data.column(0, from_list=True)
                        for i in range(1, data.ncols())])
                if shifted.rank() != data.nrows():
                    raise ValueError, "Points must span the space!\nGiven:\n%s" % data
                self._vertices = data   # for using poly_x
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

    def __reduce__(self):
        r"""
        Reduction function. Does not store data that can be relatively fast
        recomputed.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.vertices() == loads(o.dumps()).vertices()
            True
        """
        state = self.__dict__.copy()
        state.pop('_vertices')
        state.pop('_desc')
        state.pop('_distances', None)
        state.pop('_skeleton', None)
        if state.has_key('_points'):
            state['_npoints'] = state.pop('_points').ncols()
        return (LatticePolytope, (self._vertices, self._desc, False, False), state)

    def __setstate__(self, state):
        r"""
        Restores the state of pickled polytope.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.vertices() == loads(o.dumps()).vertices()
            True
        """
        self.__dict__.update(state)
        if state.has_key('_faces'):     # Faces do not remember polytopes
            for d_faces in self._faces:
                for face in d_faces:
                    face._polytope = self
        if state.has_key('_nef_partitions'): # NEF-partitions do not remember polytopes
            for partition in self._nef_partitions:
                partition._polytope = self

    def _compute_faces(self):
        r"""
        Compute and cache faces of this polytope.

        If this polytope is reflexive and the polar polytope was already
        computed, computes faces of both in order to save time and preserve the
        one-to-one correspondence between the faces of this polytope of
        dimension d and the faces of the polar polytope of codimension d+1.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: v = o.__dict__.pop("_faces", None) # faces may be cached already
            sage: o.__dict__.has_key("_faces")
            False
            sage: o._compute_faces()
            sage: o.__dict__.has_key("_faces")
            True
        """
        try:
            if self._constructed_as_polar:
                # "Polar of polar polytope" computed by poly.x may have the
                # order of vertices different from the original polytope. Thus,
                # in order to have consistent enumeration of vertices and faces
                # we must run poly.x on the original polytope.
                self._faces = Sequence([], cr=True)
                for polar_d_faces in reversed(self._polar.faces()):
                    self._faces.append([_PolytopeFace(self, f._facets,
                                        f._vertices) for f in polar_d_faces])
                self._faces.set_immutable()
                return
        except AttributeError:
            self._read_faces(self.poly_x("i"))

    def _compute_hodge_numbers(self):
        r"""
        Compute Hodge numbers for the current nef_partitions.

        This function (currently) always raises an exception directing to use
        another way for computing Hodge numbers.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: o._compute_hodge_numbers()
            Traceback (most recent call last):
            ...
            NotImplementedError: use nef_partitions(hodge_numbers=True)!
        """
        raise NotImplementedError, "use nef_partitions(hodge_numbers=True)!"

    def _face_compute_points(self, face):
        r"""
        Compute and cache lattice points of the given \code{face} of this polytope.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: e = o.faces(dim=1)[0]
            sage: v = e.__dict__.pop("_points", None) # points may be cached already
            sage: e.__dict__.has_key("_points")
            False
            sage: o._face_compute_points(e)
            sage: e.__dict__.has_key("_points")
            True
        """
        m = self.distances().matrix_from_rows(face._facets)
        cols = m.columns(copy=False)
        points = [i for i, col in enumerate(cols) if sum(col) == 0]
        face._points = Sequence(points, int, check=False)
        face._points.set_immutable()

    def _face_split_points(self, face):
        r"""
        Compute and cache boundary and interior lattice points of \code{face}.

        TESTS:
            sage: c = lattice_polytope.octahedron(3).polar()
            sage: f = c.facets()[0]
            sage: v = f.__dict__.pop("_interior_points", None)
            sage: f.__dict__.has_key("_interior_points")
            False
            sage: v = f.__dict__.pop("_boundary_points", None)
            sage: f.__dict__.has_key("_boundary_points")
            False
            sage: c._face_split_points(f)
            sage: f._interior_points
            [18]
            sage: f._boundary_points
            [0, 1, 2, 3, 11, 15, 21, 25]
            sage: f.points()
            [0, 1, 2, 3, 11, 15, 18, 21, 25]
        """
        face._interior_points = Sequence([], int, check=False)
        face._boundary_points = Sequence(face.points()[:face.nvertices()], int,
                                                                    check=False)
        non_vertices = face.points()[face.nvertices():]
        distances = self.distances()

        other_facets = [i for i in range(self.nfacets())
                                     if not i in face._facets]
        for p in non_vertices:
            face._interior_points.append(p)
            for f in other_facets:
                if distances[f, p] == 0:
                    face._interior_points.pop()
                    face._boundary_points.append(p)
                    break
        face._interior_points.set_immutable()
        face._boundary_points.set_immutable()

    def _palp(self, command):
        r"""
        Run \code{command} on vertices of this polytope.

        Returns the output of \code{command} as a string.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: o._palp("poly.x -f")
            'M:7 6 N:27 8 Pic:17 Cor:0\n'
            sage: print o._palp("nef.x -f -N -p") # random time information
            M:27 8 N:7 6  codim=2 #part=5
            H:[0] P:0 V:2 4 5       0sec  0cpu
            H:[0] P:2 V:3 4 5       0sec  0cpu
            H:[0] P:3 V:4 5       0sec  0cpu
            np=3 d:1 p:1    0sec     0cpu
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

    def _read_equations(self, data):
        r"""
        Read equations of facets/vertices of polar polytope from string or file.

        TESTS:
        For a reflexive polytope construct the polar polytope:
            sage: p = LatticePolytope(matrix(ZZ,2,3,[1,0,-1,0,1,-1]))
            sage: p.vertices()
            [ 1  0 -1]
            [ 0  1 -1]
            sage: s = p.poly_x("e")
            sage: print s
            3 2  Vertices of P-dual <-> Equations of P
               2  -1
              -1   2
              -1  -1
            sage: p.__dict__.has_key("_polar")
            False
            sage: p._read_equations(s)
            sage: p._polar._vertices
            [ 2 -1 -1]
            [-1  2 -1]

        For a non-reflexive polytope cache facet equations:
            sage: p = LatticePolytope(matrix(ZZ,2,3,[1,0,-1,0,2,-3]))
            sage: p.vertices()
            [ 1  0 -1]
            [ 0  2 -3]
            sage: p.__dict__.has_key("_facet_normals")
            False
            sage: p.__dict__.has_key("_facet_constants")
            False
            sage: s = p.poly_x("e")
            sage: print s
            3 2  Equations of P
               5  -1     2
              -2  -1     2
              -3   2     3
            sage: p._read_equations(s)
            sage: p._facet_normals
            [ 5 -1]
            [-2 -1]
            [-3  2]
            sage: p._facet_constants
            (2, 2, 3)
        """
        if isinstance(data, str):
            f = StringIO.StringIO(data)
            self._read_equations(f)
            f.close()
            return
        try:
            self._is_reflexive
            # If the above line didn't cause an exception, there is no need to
            # read equations of this polytope. Moreover, doing so can corrupt
            # data if this polytope was constructed as polar. Just skip data:
            skip_palp_matrix(data)
        except AttributeError:
            pos = data.tell()
            line = data.readline()
            self._is_reflexive = line.find("Vertices of P-dual") != -1
            if self._is_reflexive:
                data.seek(pos)
                self._polar = LatticePolytope(read_palp_matrix(data),
                                "A polytope polar to " + str(self._desc),
                                copy_vertices=False)
                self._polar._is_reflexive = True
                self._polar._constructed_as_polar = True
                self._polar._polar = self
            else:
                normals = []
                constants = []
                d = self.dim()
                for i in range(int(line.split()[0])):
                    line = data.readline()
                    numbers = [int(number) for number in line.split()]
                    constants.append(numbers.pop())
                    normals.append(numbers)
                self._facet_normals = matrix(ZZ, normals)
                self._facet_constants = vector(ZZ, constants)
                self._facet_normals.set_immutable()

    def _read_faces(self, data):
        r"""
        Read faces informations from string or file.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: s = o.poly_x("i")
            sage: print s
            Incidences as binary numbers [F-vector=(6 12 8)]:
            v[d][i]: sum_j Incidence(i'th dim-d-face, j-th vertex) x 2^j
            v[0]: 100000 000010 000001 001000 010000 000100
            v[1]: 100010 100001 000011 101000 001010 110000 010001 011000 000110 000101 001100 010100
            v[2]: 100011 101010 110001 111000 000111 001110 010101 011100
            f[d][i]: sum_j Incidence(i'th dim-d-face, j-th facet) x 2^j
            f[0]: 00001111 00110011 01010101 10101010 11001100 11110000
            f[1]: 00000011 00000101 00010001 00001010 00100010 00001100 01000100 10001000 00110000 01010000 10100000 11000000
            f[2]: 00000001 00000010 00000100 00001000 00010000 00100000 01000000 10000000
            sage: v = o.__dict__.pop("_faces", None)
            sage: o.__dict__.has_key("_faces")
            False
            sage: o._read_faces(s)
            sage: o._faces
            [
            [[5], [1], [0], [3], [4], [2]],
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]],
            [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]
            ]

        Cannot be used for "polar polytopes," their faces are constructed from faces of the original one to preserve facial duality.
            sage: c = o.polar()
            sage: s = c.poly_x("i")
            sage: print s
            Incidences as binary numbers [F-vector=(8 12 6)]:
            v[d][i]: sum_j Incidence(i'th dim-d-face, j-th vertex) x 2^j
            v[0]: 00010000 00000001 01000000 00000100 00100000 00000010 10000000 00001000
            v[1]: 00010001 01010000 00000101 01000100 00110000 00000011 00100010 11000000 10100000 00001100 00001010 10001000
            v[2]: 01010101 00110011 11110000 00001111 11001100 10101010
            f[d][i]: sum_j Incidence(i'th dim-d-face, j-th facet) x 2^j
            f[0]: 000111 001011 010101 011001 100110 101010 110100 111000
            f[1]: 000011 000101 001001 010001 000110 001010 100010 010100 100100 011000 101000 110000
            f[2]: 000001 000010 000100 001000 010000 100000
            sage: c._read_faces(s)
            Traceback (most recent call last):
            ...
            ValueError: Cannot read face structure for a polytope constructed as polar, use _compute_faces!
        """
        if isinstance(data, str):
            f = StringIO.StringIO(data)
            self._read_faces(f)
            f.close()
            return
        try:
            if self._constructed_as_polar:
                raise ValueError, ("Cannot read face structure for a polytope "
                    + "constructed as polar, use _compute_faces!")
        except AttributeError:
            pass
        data.readline()
        v = _read_poly_x_incidences(data, self.dim())
        f = _read_poly_x_incidences(data, self.dim())
        self._faces = Sequence([], cr=True)
        for i in range(len(v)):
            self._faces.append([_PolytopeFace(self, vertices, facets)
                                    for vertices, facets in zip(v[i], f[i])])
        self._faces.set_immutable()

    def _read_nef_partitions(self, data):
        r"""
        Read NEF-partitions from string or file.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: s = o.nef_x("-p -N -Lv")
            sage: print s # random time values
            M:27 8 N:7 6  codim=2 #part=5
            3 6 Vertices in N-lattice:
                1    0    0   -1    0    0
                0    1    0    0   -1    0
                0    0    1    0    0   -1
            ------------------------------
                1    0    0    1    0    0  d=2  codim=2
                0    1    0    0    1    0  d=2  codim=2
                0    0    1    0    0    1  d=2  codim=2
             P:0 V:2 4 5   (0 2) (1 1) (2 0)     0sec  0cpu
             P:2 V:3 4 5   (1 1) (1 1) (1 1)     0sec  0cpu
             P:3 V:4 5   (0 2) (1 1) (1 1)     0sec  0cpu
            np=3 d:1 p:1    0sec     0cpu
            sage: v = o.__dict__.pop("_nef_partitions", None)
            sage: o.__dict__.has_key("_nef_partitions")
            False
            sage: o._read_nef_partitions(s)
            sage: o._nef_partitions
            [
            [1, 1, 0, 1, 0, 0],
            [1, 1, 1, 0, 0, 0],
            [1, 1, 1, 1, 0, 0]
            ]
        """
        if isinstance(data, str):
            f = StringIO.StringIO(data)
            self._read_nef_partitions(f)
            f.close()
            return
        nvertices = self.nvertices()
        data.readline() # Skip M/N information
        nef_vertices = read_palp_matrix(data)
        if self.vertices() != nef_vertices:
            # It seems that we SHOULD worry about this...
            # raise RunTimeError, "nef.x changed the order of vertices!"
            trans = [self.vertices().columns().index(v)
                        for v in nef_vertices.columns()]
            for i, p in enumerate(partitions):
                partitions[i] = [trans[v] for v in p]
        line = data.readline()
        if line == "":
            raise ValueError, "more data expected!"
        partitions = Sequence([], cr=True)
        while len(line) > 0 and line.find("np=") == -1:
            if line.find("V:") == -1:
                line = data.readline()
                continue
            start = line.find("V:") + 2
            end = line.find("  ", start)  # Find DOUBLE space
            pvertices = Sequence(line[start:end].split(),int)
            partition = [1] * nvertices
            for v in pvertices:
                partition[v] = 0
            partition = NEFPartition(partition, self)
            partition._is_product = line.find(" D ") != -1
            partition._is_projection = line.find(" DP ") != -1
            # Set the stuff
            start = line.find("H:")
            if start != -1:
                start += 2
                end = line.find("[", start)
                partition._hodge_numbers = Sequence(line[start:end].split(),
                        int, cr=False)
                partition._hodge_numbers.set_immutable()
            partitions.append(partition)
            line = data.readline()
        start = line.find("np=")
        if start == -1:
            raise ValueError, """Wrong data format, cannot find "np="!"""
#         The following block seems to be unnecessary (and requires taking into
#         account projections/products)
#         # Compare the number of found partitions with statistic.
#         start += 3
#         end = line.find(" ", start)
#         np = int(line[start:end])
#         if False and np != len(partitions):
#             raise ValueError, ("Found %d partitions, expected %d!" %
#                                  (len(partitions), np))
        partitions.set_immutable()
        self._nef_partitions = partitions

    def _repr_(self):
        r"""
        Return a string representation of this polytope.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: o._repr_()
            'An octahedron: 3-dimensional, 6 vertices.'
        """
        s = str(self._desc)
        s += ": %d-dimensional, %d vertices." % (self.dim(), self.nvertices())
        return s

    def dim(self):
        r"""
        Return the dimension of this polytope.

        EXAMPLES:
        We create a 3-dimensional octahedron and check its dimension:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.dim()
            3
        """
        return self._vertices.nrows()

    def distances(self, point=None):
        r"""
        Return the matrix of distances for this polytope or distances for the
        given point.

        The matrix of distances m gives distances m[i,j] between the i-th facet
        (which is also the i-th vertex of the polar polytope in the reflexive
        case) and j-th point of this polytope.

        If point is specified, integral distances from the point to all facets
        of this polytope will be computed.

        EXAMPLES:
        The matrix of distances for a 3-dimensional octahedron:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.distances()
            [0 0 2 2 2 0 1]
            [2 0 2 0 2 0 1]
            [0 2 2 2 0 0 1]
            [2 2 2 0 0 0 1]
            [0 0 0 2 2 2 1]
            [2 0 0 0 2 2 1]
            [0 2 0 2 0 2 1]
            [2 2 0 0 0 2 1]

        Distances from facets to the point (1,2,3):
            sage: o.distances([1,2,3])
            (1, 3, 5, 7, -5, -3, -1, 1)

        It is OK to use RATIONAL coordinates:
            sage: o.distances([1,2,3/2])
            (-1/2, 3/2, 7/2, 11/2, -7/2, -3/2, 1/2, 5/2)
            sage: o.distances([1,2,sqrt(2)])
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sqrt(2) to a rational

        """
        if point != None:
            if self.is_reflexive():
                return (self._polar._vertices.transpose() * vector(QQ, point)
                        + vector(QQ, [1]*self.nfacets()))
            else:
                return self._facet_normals*vector(QQ, point)+self._facet_constants
        try:
            return self._distances
        except AttributeError:
            if self.is_reflexive():
                V = self._polar._vertices.transpose()
                P = self.points()
                self._distances = V * P
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

    def edges(self):
        r"""
        Return the sequence of edges of this polytope (i.e. faces of dimension 1).

        EXAMPLES:
        The octahedron has 12 edges:
            sage: o = lattice_polytope.octahedron(3)
            sage: len(o.edges())
            12
            sage: o.edges()
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]]

        """
        return self.faces(dim=1)

    def faces(self, dim=None, codim=None):
        r"""
        Return the sequence of faces of this polytope.

        If \code{dim} or \code{codim} are specified, returns a sequence of
        faces of the corresponding dimension or codimension. Otherwise returns
        the sequence of such sequences for all dimensions.

        EXAMPLES:
        All faces of the 3-dimensional octahedron:
            sage: o = lattice_polytope.octahedron(3)
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
        r"""
        Return the sequence of facets of this polytope (i.e. faces of codimension 1).

        EXAMPLES:
        All facets of the 3-dimensional octahedron:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.facets()
            [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]

        Facets are the same as faces of codimension one:
            sage: o.facets() is o.faces(codim=1)
            True

        """
        return self.faces(codim=1)

    # Dictionaries of normal forms
    _rp2d_dict = None
    _rp3d_dict = None

    def index(self):
        r"""
        Return the index of this polytope in the internal database of 2- or
        3-dimensional reflexive polytopes. Databases are stored in the directory
        of the package.

        NOTE: The first call to this function for each dimension can take
        a second or so while the dictionary of all polytopes is loaded, but
        after that it is cached and fast.

        EXAMPLES:
        We check what is the index of the "diamond" in the database:
            sage: o = lattice_polytope.octahedron(2)
            sage: o.index()
            3

        Note that polytopes with the same index are not necessarily the same:
            sage: o.vertices()
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: lattice_polytope.ReflexivePolytope(2,3).vertices()
            [ 1  0  0 -1]
            [ 0  1 -1  0]

        But they are in the same GL(Z^n) orbit and have the same normal form:
            sage: o.normal_form()
            [ 1  0  0 -1]
            [ 0  1 -1  0]
            sage: lattice_polytope.ReflexivePolytope(2,3).normal_form()
            [ 1  0  0 -1]
            [ 0  1 -1  0]

        """
        try:
            return self._index
        except AttributeError:
            if not self.is_reflexive():
                raise NotImplementedError, "only reflexive polytopes can be indexed!"
            if self.dim() == 2:
                if LatticePolytopeClass._rp2d_dict == None:
                    LatticePolytopeClass._rp2d_dict = load(data_location+"reflexive_polytopes_2d_dict")
#                     f = open(data_location+"reflexive_polytopes_2d_dict", "rb")
#                     LatticePolytopeClass._rp2d_dict = loads(f.read())
#                     f.close()
                self._index = LatticePolytopeClass._rp2d_dict[str(self.normal_form())]
            elif self.dim() == 3:
                if LatticePolytopeClass._rp3d_dict == None:
                    LatticePolytopeClass._rp3d_dict = load(data_location+"reflexive_polytopes_3d_dict")
#                     f = open(data_location+"reflexive_polytopes_3d_dict", "rb")
#                     LatticePolytopeClass._rp3d_dict = loads(f.read())
#                     f.close()
                self._index = LatticePolytopeClass._rp3d_dict[str(self.normal_form())]
            else:
                raise NotImplementedError, "only 2- and 3-dimensional polytopes can be indexed!"
            return self._index

    def is_reflexive(self):
        r"""
        Return True if this polytope is reflexive.

        EXAMPLES:
        The 3-dimensional octahedron is reflexive (and 4318 other 3-polytopes):
            sage: o = lattice_polytope.octahedron(3)
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
            self._read_equations(self.poly_x("e"))
            return self._is_reflexive

    def mif(self, partition):
        r"""
        Return all vectors $m_{i,f}$, grouped into matrices.

        INPUT:
            partition -- NEF-partition (instance of class NEFPartition)

        OUTPUT:
            A sequence of matrices, one for each facet f of this polytope.
            Each row of each matrix corresponds to a part of the NEF-partition.

        EXAMPLES:
        We compute $m_{i,f}$ matrices for one of the nef-partitions of
        the 3-dimensional octahedron:
            sage: o = lattice_polytope.octahedron(3)
            sage: nefp = o.nef_partitions()[0]
            sage: o.mif(nefp)
            [
            [ 0  0  1]
            [-1 -1  0],
            [ 0  0  1]
            [ 1 -1  0],
            [ 0  1  1]
            [-1  0  0],
            [0 1 1]
            [1 0 0],
            [ 0  0 -1]
            [-1 -1  0],
            [ 0  0 -1]
            [ 1 -1  0],
            [ 0  1 -1]
            [-1  0  0],
            [ 0  1 -1]
            [ 1  0  0]
            ]
        """
        if not isinstance(partition, NEFPartition):
            raise TypeError, "Expect a NEF-partition, got %s!" % type(partition)
        result = Sequence([], cr=True)
        for f in self.facets():
            a = self.vertices().matrix_from_columns(f.vertices())
            b = matrix(ZZ,partition.nparts(), f.nvertices())
            for j, v in enumerate(f.vertices()):
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

    def nef_partitions(self, keep_symmetric=False, keep_products=True,
        keep_projections=True, hodge_numbers=False):
        r"""
        Return the sequence of NEF-partitions for this polytope.

        INPUT:
            keep_symmetric -- (default: False) if True, "-s" option will be
                    passed to nef.x in order to keep symmetric partitions;
            keep_products -- (default: True) if True, "-D" option will be
                    passed to nef.x in order to keep product partitions;
            keep_projections -- (default: True) if True, "-P" option will be
                    passed to nef.x in order to keep projection partitions;
            hodge_numbers -- (default: False) if False, "-p" option will be
                    passed to nef.x in order to skip Hodge numbers computation,
                    which takes a lot of time.

        EXAMPLES:
        NEF-partitions of the 4-dimensional octahedron:
            sage: o = lattice_polytope.octahedron(4)
            sage: o.nef_partitions()
            [
            [1, 1, 0, 0, 1, 1, 0, 0] (direct product),
            [1, 1, 1, 0, 1, 0, 0, 0],
            [1, 1, 1, 0, 1, 1, 0, 0],
            [1, 1, 1, 0, 1, 1, 1, 0] (direct product),
            [1, 1, 1, 1, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 0] (projection)
            ]

        Now we omit projections:
            sage: o.nef_partitions(keep_projections=False)
            [
            [1, 1, 0, 0, 1, 1, 0, 0] (direct product),
            [1, 1, 1, 0, 1, 0, 0, 0],
            [1, 1, 1, 0, 1, 1, 0, 0],
            [1, 1, 1, 0, 1, 1, 1, 0] (direct product),
            [1, 1, 1, 1, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 0, 0]
            ]

        At present, Hodge numbers cannot be computed for a given NEF-partition:
            sage: o.nef_partitions()[1].hodge_numbers()
            Traceback (most recent call last):
            ...
            NotImplementedError: use nef_partitions(hodge_numbers=True)!

        But they can be obtained from nef.x for all partitions at once.
        Partitions will be exactly the same:
            sage: o.nef_partitions(hodge_numbers=True)
            [
            [1, 1, 0, 0, 1, 1, 0, 0] (direct product),
            [1, 1, 1, 0, 1, 0, 0, 0],
            [1, 1, 1, 0, 1, 1, 0, 0],
            [1, 1, 1, 0, 1, 1, 1, 0] (direct product),
            [1, 1, 1, 1, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 0, 0, 0],
            [1, 1, 1, 1, 1, 1, 0, 0],
            [1, 1, 1, 1, 1, 1, 1, 0] (projection)
            ]

        Now it is possible to get Hodge numbers:
            sage: o.nef_partitions(hodge_numbers=True)[1].hodge_numbers()
            [20]

        Since NEF-partitions are cached, Hodge numbers are accessible after the
        first request, even if you do not specify hodge_numbers anymore:
            sage: o.nef_partitions()[1].hodge_numbers()
            [20]

        We illustrate removal of symmetric partitions on a diamond:
            sage: o = lattice_polytope.octahedron(2)
            sage: o.nef_partitions()
            [
            [1, 0, 1, 0] (direct product),
            [1, 1, 0, 0],
            [1, 1, 1, 0] (projection)
            ]
            sage: o.nef_partitions(keep_symmetric=True)
            [
            [1, 1, 0, 1] (projection),
            [1, 0, 1, 1] (projection),
            [1, 0, 0, 1],
            [0, 1, 1, 1] (projection),
            [0, 1, 0, 1] (direct product),
            [0, 0, 1, 1],
            [1, 1, 1, 0] (projection)
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
        keys = "-N -V"
        if keep_symmetric:
            keys += " -s"
        if keep_products:
            keys += " -D"
        if keep_projections:
            keys += " -P"
        if not hodge_numbers:
            keys += " -p"
        try:
            oldkeys = self._npkeys
            if oldkeys == keys:
                return self._nef_partitions
            if not (hodge_numbers and oldkeys.find("-p") != -1
                or keep_symmetric and oldkeys.find("-s") == -1
                or not keep_symmetric and oldkeys.find("-s") != -1
                or keep_projections and oldkeys.find("-P") == -1
                or keep_products and oldkeys.find("-D") == -1):
                # Select only necessary partitions
                return Sequence([p for p in self._nef_partitions if
                    (keep_projections or not p._is_projection)
                    and (keep_products or not p._is_product)], cr=True)
        except AttributeError:
            pass
        self._read_nef_partitions(self.nef_x(keys))
        self._npkeys = keys
        return self._nef_partitions

    def nef_x(self, keys):
        r"""
        Run nef.x wilth given \code{keys} on vertices of this polytope.

        INPUT:
            keys -- a string of options passed to nef.x. The key "-f" is added
                    automatically.

        OUTPUT:
            the output of nef.x as a string.

        EXAMPLES:
        This call is used internally for computing NEF-partitions:
            sage: o = lattice_polytope.octahedron(3)
            sage: s = o.nef_x("-N -V -p")
            sage: s                      # output contains random time
            M:27 8 N:7 6  codim=2 #part=5
            3 6  Vertices of P:
                1    0    0   -1    0    0
                0    1    0    0   -1    0
                0    0    1    0    0   -1
             P:0 V:2 4 5       0sec  0cpu
             P:2 V:3 4 5       0sec  0cpu
             P:3 V:4 5       0sec  0cpu
            np=3 d:1 p:1    0sec     0cpu
        """
        return self._palp("nef.x -f " + keys)

    def nfacets(self):
        r"""
        Return the number of facets of this polytope.

        EXAMPLES:
        The number of facets of the 3-dimensional octahedron:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.nfacets()
            8
        """
        try:
            return self._nfacets
        except AttributeError:
            if self.is_reflexive():
                self._nfacets = self.polar().nvertices()
            else:
                self._nfacets = len(self._facet_normals)
            return self._nfacets

    def normal_form(self):
        r"""Return the normal form of vertices of the polytope.

        Two lattice polytopes are in the same GL(Z)-orbit if and only if their
        normal forms are the same.

        EXAMPLES:
        We compute the normal form of the "diamond":
            sage: o = lattice_polytope.octahedron(2)
            sage: o.vertices()
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: o.normal_form()
            [ 1  0  0 -1]
            [ 0  1 -1  0]

        The diamond is the 3rd polytope in the internal database...
            sage: o.index()
            3
            sage: lattice_polytope.ReflexivePolytope(2,3).vertices()
            [ 1  0  0 -1]
            [ 0  1 -1  0]
        """
        try:
            return self._normal_form
        except AttributeError:
            self._normal_form = read_palp_matrix(self.poly_x("N"))
            self._normal_form.set_immutable()
            return self._normal_form

    def npoints(self):
        r"""
        Return the number of lattice points of this polytope.

        EXAMPLES:
        The number of lattice points of the 3-dimensional octahedron and its
        polar cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.npoints()
            7
            sage: cube = o.polar()
            sage: cube.npoints()
            27
        """
        try:
            return self._npoints
        except AttributeError:
            return self.points().ncols()

    def nvertices(self):
        r"""
        Return the number of vertices of this polytope.

        EXAMPLES:
        The number of vertices of the 3-dimensional octahedron and its
        polar cube:
            sage: o = lattice_polytope.octahedron(3)
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
            sage: o = lattice_polytope.octahedron(3)
            sage: o.parent()
            Set of all Lattice Polytopes
        """
        return SetOfAllLatticePolytopes

    def plot(self, camera_center=None):
        """
        Draw a 3d picture of the polytope.

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
            sage: o = lattice_polytope.octahedron(3)
            sage: o.plot().show()
        """
        if self.dim() != 3:
            raise ValueError, "Polytope must have dimension 3!"
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

    def plot3d(self,
            show_facets=True, facet_opacity=0.5, facet_color=(0,1,0),
            show_edges=True, edge_thickness=3, edge_color=(0.5,0.5,0.5),
            show_vertices=True, vertex_size=10, vertex_color=(1,0,0),
            show_points=True, point_size=10, point_color=(0,0,1),
            show_vindices=None, vindex_color=(0,0,0),
            show_pindices=None, pindex_color=(0,0,0),
            index_shift=1.1):
        r"""
        Return a 3d-plot of a 3-dimensional polytope.

        By default, everything is shown with more or less pretty combination
        of size and color parameters.

        INPUT:
        Most of the parameters are self-explanatory:
            show_facets     -- (default:True)
            facet_opacity   -- (default:0.5)
            facet_color     -- (default:(0,1,0))
            show_edges      -- (default:True) whether to draw edges as lines
            edge_thickness  -- (default:3)
            edge_color      -- (default:(0.5,0.5,0.5))
            show_vertices   -- (default:True) whether to draw vertices as balls
            vertex_size     -- (default:10)
            vertex_color    -- (default:(1,0,0))
            show_points     -- (default:True) whether to draw other poits as balls
            point_size      -- (default:10)
            point_color     -- (default:(0,0,1))
            show_vindices   -- (default:same as show_vertices) whether to show indices of vertices
            vindex_color    -- (default:(0,0,0)) color for vertex labels
            show_pindices   -- (default:same as show_points) whether to show indices of other points
            pindex_color    -- (default:(0,0,0)) color for point labels
            index_shift     -- (default:1.1)) if 1, labels are placed exactly at
                                the corresponding points. Otherwise the label position
                                is computed as a multiple of the point position vector.

        EXAMPLES:
        The default plot of a cube:
            sage: c = lattice_polytope.octahedron(3).polar()
            sage: c.plot3d()

        Plot without facets and points, shown without the frame:
            sage: c.plot3d(show_facets=false,show_points=false).show(frame=False)
        """
        if self.dim() != 3:
            raise ValueError, "Only 3-dimensional polytopes can be plotted in 3D!"
        pplot = 0
        if show_facets:
            pplot += IndexFaceSet([f.traverse_boundary() for f in self.facets()],
                    self.vertices().columns(copy=False), opacity=facet_opacity, rgbcolor=facet_color)
        if show_edges:
            for e in self.edges():
                pplot += line3d([self.vertex(e.vertices()[0]), self.vertex(e.vertices()[1])],
                        thickness=edge_thickness, rgbcolor=edge_color)
        if show_vertices:
            pplot += point3d(self.vertices().columns(copy=False), size=vertex_size, rgbcolor=vertex_color)
        if show_vindices == None:
            show_vindices = show_vertices
        if show_vindices:
            for i,v in enumerate(self.vertices().columns(copy=False)):
                pplot += text3d(i, index_shift*v, rgbcolor=vindex_color)
        if show_points:
            pplot += point3d(self.points().columns(copy=False)[self.nvertices():],
                        size=point_size, rgbcolor=point_color)
        if show_pindices == None:
            show_pindices = show_points
        if show_pindices:
            for i in range(self.nvertices(),self.npoints()):
                pplot += text3d(i, index_shift*self.point(i), rgbcolor=pindex_color)
        return pplot

    def show(self, camera_center=None):
        """
        Draw a 3d picture of the polytope.

        See self.plot? for more details.

        EXAMPLES:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.show()
        """
        self.plot(camera_center=camera_center).show()

    def point(self, i):
        r"""
        Return the i-th point of this polytope, i.e. the i-th column of the
        matrix returned by points().

        EXAMPLES:
        First few points are actually vertices:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.vertices()
            [ 1  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: o.point(1)
            (0, 1, 0)

        The only other point in the octahedron is the origin:
            sage: o.point(6)
            (0, 0, 0)
            sage: o.points()
            [ 1  0  0 -1  0  0  0]
            [ 0  1  0  0 -1  0  0]
            [ 0  0  1  0  0 -1  0]
        """
        # Extra checks are made to compensate for a bug in SAGE - column accepts any number.
        if i < 0:
            raise ValueError, "polytopes don't have negative points!"
        elif i < self.nvertices():
            return self.vertex(i)
        elif i < self.npoints():
            return self.points().column(i, from_list=True)
        else:
            raise ValueError, "the polytope does not have %d points!" % (i+1)

    def points(self):
        r"""
        Return all lattice points of this polytope as columns of a matrix.

        EXAMPLES:
        The lattice points of the 3-dimensional octahedron and its polar cube:
            sage: o = lattice_polytope.octahedron(3)
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
            sage: o = lattice_polytope.octahedron(3)
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
            sage: o = lattice_polytope.octahedron(3)
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
            sage: BIGO = lattice_polytope.octahedron(7)
            sage: BIGO
            An octahedron: 7-dimensional, 14 vertices.
            sage: BIGO.poly_x("e")      # possibly different output depending on your system
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

    def skeleton(self):
        r"""
        Return the graph of the one-skeleton of this polytope.

        EXAMPLES:
        We construct the one-skeleton graph for the "diamond":
            sage: o = lattice_polytope.octahedron(2)
            sage: g = o.skeleton()
            sage: g
            Graph on 4 vertices
            sage: g.edges()
            [(0, 1, None), (0, 3, None), (1, 2, None), (2, 3, None)]
        """
        try:
            return self._skeleton
        except AttributeError:
            skeleton = Graph()
            skeleton.add_vertices(self.skeleton_points(1))
            for e in self.faces(dim=1):
                points = e.ordered_points()
                for i in range(len(points)-1):
                    skeleton.add_edge(points[i], points[i+1])
            self._skeleton = skeleton
            return skeleton

    def skeleton_points(self, k=1):
        r"""
        Return the increasing list of indices of lattice points in k-skeleton
        of the polytope (k is 1 by default).

        EXAMPLES:
        We compute all skeleton points for the cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: c = o.polar()
            sage: c.skeleton_points()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 19, 21, 22, 23, 25, 26]

        The default was 1-skeleton:
            sage: c.skeleton_points(k=1)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 19, 21, 22, 23, 25, 26]

        0-skeleton just lists all vertices:
            sage: c.skeleton_points(k=0)
            [0, 1, 2, 3, 4, 5, 6, 7]

        2-skeleton lists all points except for the origin (point #17):
            sage: c.skeleton_points(k=2)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26]

        3-skeleton includes all points:
            sage: c.skeleton_points(k=3)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]

        It is OK to compute higher dimensional skeletons - you will get the list of all points:
            sage: c.skeleton_points(k=100)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]
        """
        if k >= self.dim():
            return range(self.npoints())
        skeleton = set([])
        for face in self.faces(dim=k):
            skeleton.update(face.points())
        skeleton = list(skeleton)
        skeleton.sort()
        return skeleton

    def skeleton_show(self, normal=None):
        r"""Show the graph of one-skeleton of this polytope.

        Works only for 3-dimensional polytopes.

        INPUT:
            normal -- a 3-dimensional vector (can be given as a list), which
                should be perpendicular to the screen. If not given, will be
                selected randomly (new each time and it may be far from "nice").

        EXAMPLES:
        Show a pretty picture of the octahedron:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.skeleton_show([1,2,4])
        """
        if normal == None:
            normal = [ZZ.random_element(20),ZZ.random_element(20),ZZ.random_element(20)]
        normal = matrix(QQ,3,1,list(normal))
        projectionm = normal.kernel().basis_matrix()
        positions = dict(enumerate([list(c) for c in (projectionm*self.points()).columns(copy=False)]))
        self.skeleton().show(pos=positions)

    def vertex(self, i):
        r"""
        Return the i-th vertex of this polytope, i.e. the i-th column of the
        matrix returned by vertices().

        EXAMPLES:
        Note that numeration starts with zero:
            sage: o = lattice_polytope.octahedron(3)
            sage: o.vertices()
            [ 1  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: o.vertex(3)
            (-1, 0, 0)
        """
        # The check is added to compensate for a bug in SAGE - column works for any numbers
        if i < 0:
            raise ValueError, "polytopes don't have negative vertices!"
        elif i > self.nvertices():
            raise ValueError, "the polytope does not have %d vertices!" % (i+1)
        else:
            return self._vertices.column(i, from_list=True)

    def vertices(self):
        r"""
        Return vertices of this polytope as columns of a matrix.

        EXAMPLES:
        The lattice points of the 3-dimensional octahedron and its polar cube:
            sage: o = lattice_polytope.octahedron(3)
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
    Construct a NEF-partition from the given list. If the second argument is given,
    it must be the polytope of the NEF-partition.

    A NEF-partition with $k$ parts, $V = V_0 \cap V_1 \cap \dots \cap V_{k-1}$,
    is represented by a single list of lenght $n = \#V$, in which the $i$-th
    entry is the part number of the $i$-th vertex of a polytope.

    \emph{NOTE THAT NUMERATON OF PARTS STARTS WITH ZERO.}

    EXAMPLES:
    All elements of the list will be coerced to integers, so it is OK to use
    either a list of numbers or a list of strings:
        sage: lattice_polytope.NEFPartition([1, 1, 0, 0, 0, 1])
        [1, 1, 0, 0, 0, 1]
        sage: lattice_polytope.NEFPartition(['1', '1', '0', '0', '0', '1'])
        [1, 1, 0, 0, 0, 1]
    """

    def __init__(self, data, polytope=None):
        r"""
        Construct a NEF-partition.

        TESTS:
            sage: lattice_polytope.NEFPartition([1,0,1,0])
            [1, 0, 1, 0]
        """
        Sequence.__init__(self, data, int, cr=False)
        self._n = max(self) + 1
        self._polytope = polytope
        self.set_immutable()

    def __reduce__(self):
        r"""
        Reduction function. Does not store data that can be relatively fast
        recomputed.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np == loads(np.dumps())
            True
        """
        state = self.__dict__.copy()
        # Reference to the polytope is not pickles - the polytope will restore it
        state.pop('_polytope')
        state.pop('_n')
        return (NEFPartition, (Sequence(self), None), state)

    def _repr_(self):
        r"""
        Return a string representation of this NEF-partition.

        TESTS:
            sage: np = lattice_polytope.NEFPartition([1,1,0,0])
            sage: np._repr_()
            '[1, 1, 0, 0]'
        """
        s = Sequence._repr_(self)
        try:
            # We may or may not know the type of the partition
            if self._is_projection:
                s += " (projection)"
            if self._is_product:
                s += " (direct product)"
        except AttributeError:
            pass
        return s

    def dual(self):
        r"""
        Compute dual nef-partition.

        Computation follows the article "Combinatorial aspects of mirror symmetry"
        by V.Batyrev and B.Nill, initial polytope corresponds to Nabla*, output to Delta*.

        EXAMPLES:
        We compute the dual nef-partition for the first nef-partition of the octahedron:
            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            [1, 1, 0, 1, 0, 0]
            sage: npd = np.dual()
            sage: npd
            [1, 0, 1, 0, 0, 1, 0, 1]
            sage: npd.polytope().vertices()
            [ 1  0  1  0  0 -1  0 -1]
            [-1  1  0  0  0 -1  1  0]
            [ 0  1  0  1 -1  0 -1  0]
        """
        Nabla_p = self._polytope
        r = self.nparts()
        delta=[]
        for i in range(r):
            delta.append([[Integer(0)]*Nabla_p.dim()])
        for p in range(Nabla_p.nvertices()):
            delta[self[p]].append(Nabla_p.point(p))
        sum = delta[0]
        for i in range(1, r):
            sum = minkowski_sum(sum, delta[i])
        Delta_p = LatticePolytope(matrix(ZZ, sum).transpose(), compute_vertices=False).polar()
        np_dual = [-1]*Delta_p.nvertices()
        for i in range(r):
            A = Nabla_p.vertices().matrix_from_columns(self.part(i)).transpose() * Delta_p.vertices()
            for j, cj in enumerate(A.columns(copy=False)):
                if min(cj) == -1:
                    if np_dual[j] != -1:
                        raise RunTimeError, "one vertex in many parts!"
                    np_dual[j] = i
        return NEFPartition(np_dual, Delta_p)

    def hodge_numbers(self):
        r"""
        Return Hodge numbers corresponding to the given nef-partition.

        NOTE: In order to use this function you must use
        nef_partitions(hodge_numbers=True) for computing nef-partitions.

        EXAMPLE:
        Currently, you need to request hodge_numbers when you compute nef-partitions.
        (The following weird construction ensures that there are no cached numbers.)
            sage: o = LatticePolytope(lattice_polytope.octahedron(4).vertices())
            sage: np = o.nef_partitions()[0]
            sage: np.hodge_numbers()
            Traceback (most recent call last):
            ...
            NotImplementedError: use nef_partitions(hodge_numbers=True)!
            sage: np = o.nef_partitions(hodge_numbers=True)[0]
            sage: np.hodge_numbers()
            [4]
        """
        try:
            return self._hodge_numbers
        except AttributeError:
            self._polytope._compute_hodge_numbers()
            return self._hodge_numbers

    def nparts(self):
        r"""
        Return the number of parts of this partitions.

        EXAMPLES:
            sage: nefp = lattice_polytope.NEFPartition([1, 1, 0, 0, 0, 1])
            sage: nefp.nparts()
            2
        """
        return self._n

    def part(self, i):
        r"""
        Return the \code{i}-th part of the partition.

        \emph{NUMERATON OF PARTS STARTS WITH ZERO.}

        EXAMPLES:
            sage: nefp = lattice_polytope.NEFPartition([1, 1, 0, 0, 0, 1])
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
            sage: nefp = lattice_polytope.NEFPartition([1, 1, 0, 0, 0, 1])
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

    def polytope(self):
        r"""
        Return the polytope of this NEF-partition.

        EXAMPLE:
            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np.polytope() == o
            True
        """
        return self._polytope


class _PolytopeFace(SageObject):
    r"""
    _PolytopeFace(polytope, vertices, facets)

    Construct a polytope face.

    POLYTOPE FACES SHOULD NOT BE CONSTRUCTED OUTSIDE OF LATTICE POLYTOPES!

    INPUT:
        polytope -- a polytope whose face is being constructed.
        vertices -- a sequence of indices of generating vertices.
        facets -- a sequence of indices of facets containing this face.
    """
    def __init__(self, polytope, vertices, facets):
        r"""
        Construct a face.

        TESTS:
            sage: o = lattice_polytope.octahedron(2)
            sage: o.faces()
            [
            [[3], [0], [2], [1]],
            [[0, 3], [2, 3], [0, 1], [1, 2]]
            ]
        """
        self._polytope = polytope
        self._vertices = vertices
        self._facets = facets

    def __reduce__(self):
        r"""
        Reduction function. Does not store data that can be relatively fast
        recomputed.

        TESTS:
            sage: o = lattice_polytope.octahedron(2)
            sage: f = o.facets()[0]
            sage: fl = loads(f.dumps())
            sage: f.vertices() == fl.vertices()
            True
            sage: f.facets() == fl.facets()
            True
        """
        state = self.__dict__.copy()
        state.pop('_polytope')
        state.pop('_vertices')
        state.pop('_facets')
        if state.has_key('_points'):
            state['_npoints'] = len(state.pop('_points'))
        if state.has_key('_interior_points'):
            state['_ninterior_points'] = len(state.pop('_interior_points'))
            state.pop('_boundary_points')
        # Reference to the polytope is not pickled - the polytope will restore it
        return (_PolytopeFace, (None, self._vertices, self._facets), state)

    def _repr_(self):
        r"""
        Return a string representation of this face.

        TESTS:
            sage: o = lattice_polytope.octahedron(3)
            sage: f = o.facets()[0]
            sage: f._repr_()
            '[0, 1, 5]'
        """
        return repr(self._vertices)

    def boundary_points(self):
        r"""
        Return a sequence of indices of boundary lattice points of this face.

        EXAMPLES:
        Boundary lattice points of one of the facets of the 3-dimensional cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.boundary_points()
            [0, 1, 2, 3, 11, 15, 21, 25]
        """
        try:
            return self._boundary_points
        except AttributeError:
            self._polytope._face_split_points(self)
            return self._boundary_points

    def facets(self):
        r"""
        Return a sequence of indices of facets containing this face.

        EXAMPLES:
        Facets containing one of the edges of the 3-dimensional octahedron:
            sage: o = lattice_polytope.octahedron(3)
            sage: edge = o.faces(dim=1)[0]
            sage: edge.facets()
            [0, 1]

        Thus \code{edge} is the intersection of facets 0 and 1:
            sage: edge
            [1, 5]
            sage: o.facets()[0]
            [0, 1, 5]
            sage: o.facets()[1]
            [1, 3, 5]
        """
        return self._facets

    def interior_points(self):
        r"""
        Return a sequence of indices of interior lattice points of this face.

        EXAMPLES:
        Interior lattice points of one of the facets of the 3-dimensional cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.interior_points()
            [18]
        """
        try:
            return self._interior_points
        except AttributeError:
            self._polytope._face_split_points(self)
            return self._interior_points

    def nboundary_points(self):
        r"""
        Return the number of boundary lattice points of this face.

        EXAMPLES:
        The number of boundary lattice points of one of the facets of
        the 3-dimensional cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.nboundary_points()
            8
        """
        return self.npoints() - self.ninterior_points()

    def nfacets(self):
        r"""
        Return the number of facets containing this face.

        EXAMPLES:
        The number of facets containing one of the edges of
        the 3-dimensional octahedron:
            sage: o = lattice_polytope.octahedron(3)
            sage: edge = o.faces(dim=1)[0]
            sage: edge.nfacets()
            2
        """
        return len(self._facets)

    def ninterior_points(self):
        r"""
        Return the number of interior lattice points of this face.

        EXAMPLES:
        The number of interior lattice points of one of the facets of
        the 3-dimensional cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.ninterior_points()
            1
        """
        try:
            return self._ninterior_points
        except AttributeError:
            return len(self.interior_points())

    def npoints(self):
        r"""
        Return the number of lattice points of this face.

        EXAMPLES:
        The number of lattice points of one of the facets of
        the 3-dimensional cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.npoints()
            9
        """
        try:
            return self._npoints
        except AttributeError:
            return len(self.points())

    def nvertices(self):
        r"""
        Return the number of vertices generating this face.

        EXAMPLES:
        The number of vertices generating one of the facets of
        the 3-dimensional cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.nvertices()
            4
        """
        return len(self._vertices)

    def ordered_points(self):
        r"""
        Return the list of indices of lattice points on the edge in their
        geometric order, from one vertex to other.

        Works only for edges, i.e. faces generated by exactly two vertices.

        EXAMPLE:
        We find all points along an edge of the cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: c = o.polar()
            sage: e = c.edges()[0]
            sage: e.vertices()
            [0, 1]
            sage: e.ordered_points()
            [0, 15, 1]
        """
        if len(self.vertices()) != 2:
            raise ValueError, "Order of points is defined for edges only!"
        pcol = self._polytope.points().columns(copy=False)
        start = pcol[self.vertices()[0]]
        end = pcol[self.vertices()[1]]
        primitive = end - start
        primitive = primitive * (1/gcd(primitive.list()))
        result = [self.vertices()[0]]
        start = start + primitive
        while start != end:
            for i in self.points():
                if start == pcol[i]:
                    result.append(i)
                    break
            start = start + primitive
        result.append(self.vertices()[1])
        return result

    def points(self):
        r"""
        Return a sequence of indices of lattice points of this face.

        EXAMPLES:
        The lattice points of one of the facets of the 3-dimensional cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.points()
            [0, 1, 2, 3, 11, 15, 18, 21, 25]
        """
        try:
            return self._points
        except AttributeError:
            self._polytope._face_compute_points(self)
            return self._points

    def traverse_boundary(self):
        r"""
        Return a list of indices of vertices of a 2-face in their boundary order.

        Needed for plot3d function of polytopes.

        EXAMPLES:
            sage: c = lattice_polytope.octahedron(3).polar()
            sage: f = c.facets()[0]
            sage: f.vertices()
            [0, 1, 2, 3]
            sage: f.traverse_boundary()
            [0, 2, 3, 1]
        """
        if self not in self._polytope.faces(dim=2):
            raise ValueError, "Boundary can be traversed only for 2-faces!"
        edges = [e for e in self._polytope.edges() if e.vertices()[0] in self.vertices() and
                e.vertices()[1] in self.vertices()]
        start = self.vertices()[0]
        l = [start]
        for e in edges:
            if start in e.vertices():
                next = e.vertices()[0] if e.vertices()[0] != start else e.vertices()[1]
        l.append(next)
        prev = start
        while len(l) < self.nvertices():
            for e in edges:
                if next in e.vertices() and prev not in e.vertices():
                    prev = next
                    next = e.vertices()[0] if e.vertices()[0] != next else e.vertices()[1]
                    l.append(next)
                    break
        return l

    def vertices(self):
        r"""
        Return a sequence of indices of vertices generating this face.

        EXAMPLES:
        The vertices generating one of the facets of the 3-dimensional cube:
            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.vertices()
            [0, 1, 2, 3]
        """
        return self._vertices


def _create_octahedron(dim):
    r"""
    Create an octahedron of the given dimension.

    Since we know that we are passing vertices, we suppress their computation.

    TESTS:
        sage: o = lattice_polytope._create_octahedron(3)
        sage: o.vertices()
        [ 1  0  0 -1  0  0]
        [ 0  1  0  0 -1  0]
        [ 0  0  1  0  0 -1]
    """
    m = matrix(ZZ, dim, 2*dim)
    for i in range(dim):
        m[i,i] = 1
        m[i,dim+i] = -1
    return LatticePolytope(m, "An octahedron", compute_vertices=False, copy_vertices=False)

_octahedrons = dict()       # Dictionary for storing created octahedrons

def _palp(command, polytopes):
    r"""
    Run \code{command} on vertices of given \code{polytopes}.

    Returns the name of the file conatining the output of \code{command}. You
    should delete it after using.

    TESTS:
        sage: o = lattice_polytope.octahedron(3)
        sage: result_name = lattice_polytope._palp("poly.x -f", [o])
        sage: f = open(result_name)
        sage: f.readlines()
        ['M:7 6 N:27 8 Pic:17 Cor:0\n']
        sage: f.close()
        sage: os.remove(result_name)
    """
    input_file_name = tmp_filename()
    input_file = open(input_file_name, "w")
    for p in polytopes:
        write_palp_matrix(p._vertices, input_file)
    input_file.close()
    output_file_name = tmp_filename()
    stdin, stdout, stderr = os.popen3("%s <%s >%s" % (command, input_file_name,
                                                             output_file_name))
    err = stderr.read()
    if len(err) > 0:
        raise RuntimeError, ("Error executing \"%s\" for a polytope sequence!"
            + "\nOutput:\n%s") % (command, err)
    os.remove(input_file_name)
    return output_file_name

def _read_nef_x_partitions(data):
    r"""
    Read all nef-partitions for one polytope from a string or an open file.

    \code{data} should be an output of nef.x.

    Returns the sequence of NEF-partitions. Each NEF-partition is given as
    a sequence of integers.

    If there are no NEF-partitions, returns the empty sequence.
    If the string is empty or EOF is reached, raises ValueError.

    TESTS:
        sage: o = lattice_polytope.octahedron(3)
        sage: s = o.nef_x("-N -p")
        sage: print s
        M:27 8 N:7 6  codim=2 #part=5
         P:0 V:2 4 5       0sec  0cpu
         P:2 V:3 4 5       0sec  0cpu
         P:3 V:4 5       0sec  0cpu
        np=3 d:1 p:1    0sec     0cpu
        sage: lattice_polytope._read_nef_x_partitions(s)
        [[2, 4, 5], [3, 4, 5], [4, 5]]
    """
    if isinstance(data, str):
        f = StringIO.StringIO(data)
        partitions = _read_nef_x_partitions(f)
        f.close()
        return partitions
    line = data.readline()
    if line == "":
        raise ValueError, "Empty file!"
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
        if False and np != len(partitions):
            raise ValueError, ("Found %d partitions, expected %d!" %
                                 (len(partitions), np))
    else:
        raise ValueError, "Wrong data format, cannot find \"np=\"!"
    return partitions

def _read_poly_x_incidences(data, dim):
    r"""
    Convert incidence data from binary numbers to sequences.

    INPUT:
        data -- an opened file with incidence information. The first line will
            be skipped, each consecutive line contains incidence information for
            all faces of one dimension, the first word of each line is a comment
            and is dropped.
        dim -- dimension of the polytope.

    OUTPUT:
        a sequence F, such that F[d][i] is a sequence of vertices or facets
        corresponding to the i-th d-dimensional face.

    TESTS:
        sage: o = lattice_polytope.octahedron(2)
        sage: result_name = lattice_polytope._palp("poly.x -fi", [o])
        sage: f = open(result_name)
        sage: f.readlines()
        ['Incidences as binary numbers [F-vector=(4 4)]:\n', "v[d][i]: sum_j Incidence(i'th dim-d-face, j-th vertex) x 2^j\n", 'v[0]: 1000 0001 0100 0010 \n', 'v[1]: 1001 1100 0011 0110 \n', "f[d][i]: sum_j Incidence(i'th dim-d-face, j-th facet) x 2^j\n", 'f[0]: 0011 0101 1010 1100 \n', 'f[1]: 0001 0010 0100 1000 \n']
        sage: f.close()
        sage: f = open(result_name)
        sage: l = f.readline()
        sage: lattice_polytope._read_poly_x_incidences(f, 2)
        [[[3], [0], [2], [1]], [[0, 3], [2, 3], [0, 1], [1, 2]]]
        sage: f.close()
        sage: os.remove(result_name)
    """
    data.readline()
    lines = [data.readline().split() for i in range(dim)]
    if len(lines) != dim:
        raise ValueError, "Not enough data!"
    n = len(lines[0][1])     # Number of vertices or facets
    result = []
    for line in lines:
        line.pop(0)
        subr = []
        for e in line:
            f = Sequence([j for j in range(n) if e[n-1-j] == '1'], int, check=False)
            f.set_immutable()
            subr.append(f)
        result.append(subr)
    return result

def all_cached_data(polytopes):
    r"""
    Compute all cached data for all given \code{polytopes} and their polars.

    This functions does it MUCH faster than member functions of
    code{LatticePolytope} during the first run. So it is recommended to use
    this functions if you work with big sets of data. None of the polytopes in
    the given sequence should be constructed as the polar polytope to another one.

    INPUT:
        a sequence of lattice polytopes.

    EXAMPLES:
    This function has no output, it is just a fast way to work with long
    sequences of polytopes. Of course, you can use short sequences as well:
        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.all_cached_data([o])
        sage: o.faces()
        [
        [[5], [1], [0], [3], [4], [2]],
        [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]],
        [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]
        ]

    However, you cannot use it for polytopes that are constructed as polar
    polytopes of others:
        sage: lattice_polytope.all_cached_data([o.polar()])
        Traceback (most recent call last):
        ...
        ValueError: Cannot read face structure for a polytope constructed as polar, use _compute_faces!
    """
    all_polars(polytopes)
    all_points(polytopes)
    all_faces(polytopes)
    reflexive = [p for p in polytopes if p.is_reflexive()]
    all_nef_partitions(reflexive)
    polar = [p.polar() for p in reflexive]
    all_points(polar)
    all_nef_partitions(polar)
    for p in polytopes + polar:
        for d_faces in p.faces():
            for face in d_faces:
                face.boundary_points()

def all_faces(polytopes):
    r"""
    Compute faces for all given \code{polytopes}.

    This functions does it MUCH faster than member functions of
    code{LatticePolytope} during the first run. So it is recommended to use
    this functions if you work with big sets of data.

    INPUT:
        a sequence of lattice polytopes.

    EXAMPLES:
    This function has no output, it is just a fast way to work with long
    sequences of polytopes. Of course, you can use short sequences as well:
        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.all_faces([o])
        sage: o.faces()
        [
        [[5], [1], [0], [3], [4], [2]],
        [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]],
        [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]
        ]

    However, you cannot use it for polytopes that are constructed as polar
    polytopes of others:
        sage: lattice_polytope.all_faces([o.polar()])
        Traceback (most recent call last):
        ...
        ValueError: Cannot read face structure for a polytope constructed as polar, use _compute_faces!
    """
    result_name = _palp("poly.x -fi", polytopes)
    result = open(result_name)
    for p in polytopes:
        p._read_faces(result)
    result.close()
    os.remove(result_name)

def all_nef_partitions(polytopes, keep_symmetric=False):
    r"""
    Compute NEF-partitions for all given \code{polytopes}.

    This functions does it MUCH faster than member functions of
    code{LatticePolytope} during the first run. So it is recommended to use
    this functions if you work with big sets of data.

    Note: member function \code{is_reflexive} will be called separately for each
    polytope. It is strictly recommended to call \code{all_polars} on the
    sequence of \code{polytopes} before using this function.

    INPUT:
        a sequence of lattice polytopes.

    EXAMPLES:
    This function has no output, it is just a fast way to work with long
    sequences of polytopes. Of course, you can use short sequences as well:
        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.all_nef_partitions([o])
        sage: o.nef_partitions()
        [
        [1, 1, 0, 1, 0, 0],
        [1, 1, 1, 0, 0, 0],
        [1, 1, 1, 1, 0, 0]
        ]

    You cannot use this function for non-reflexive polytopes:
        sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0],
        ...                   [0, 1, 0,  0, -1,  0],
        ...                   [0, 0, 2,  0,  0, -1]])
        ...
        sage: p = LatticePolytope(m)
        sage: lattice_polytope.all_nef_partitions([o, p])
        Traceback (most recent call last):
        ...
        ValueError: The given polytope is not reflexive!
        Polytope: A lattice polytope: 3-dimensional, 6 vertices.
    """
    keys = "-N -V -p"
    if keep_symmetric:
        keys += " -s"
    result_name = _palp("nef.x -f " + keys, polytopes)
    result = open(result_name)
    for p in polytopes:
        if not p.is_reflexive():
            raise ValueError, "The given polytope is not reflexive!\nPolytope: %s" % p
        p._read_nef_partitions(result)
        p._nef_partitions_s = keep_symmetric
    result.close()
    os.remove(result_name)

def all_points(polytopes):
    r"""
    Compute lattice points for all given \code{polytopes}.

    This functions does it MUCH faster than member functions of
    code{LatticePolytope} during the first run. So it is recommended to use
    this functions if you work with big sets of data.

    INPUT:
        a sequence of lattice polytopes.

    EXAMPLES:
    This function has no output, it is just a fast way to work with long
    sequences of polytopes. Of course, you can use short sequences as well:
        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.all_points([o])
        sage: o.points()
        [ 1  0  0 -1  0  0  0]
        [ 0  1  0  0 -1  0  0]
        [ 0  0  1  0  0 -1  0]
    """
    result_name = _palp("poly.x -fp", polytopes)
    result = open(result_name)
    for p in polytopes:
        p._points = read_palp_matrix(result)
        if p._points.nrows() == 0:
            raise RuntimeError, ("Cannot read points of a polytope!"
                                                        +"\nPolytope: %s" % p)
    result.close()
    os.remove(result_name)

def all_polars(polytopes):
    r"""
    Compute polar polytopes for all reflexive and equations of facets for all
    non-reflexive \code{polytopes}.

    This functions does it MUCH faster than member functions of
    code{LatticePolytope} during the first run. So it is recommended to use
    this functions if you work with big sets of data.

    INPUT:
        a sequence of lattice polytopes.

    EXAMPLES:
    This function has no output, it is just a fast way to work with long
    sequences of polytopes. Of course, you can use short sequences as well:
        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.all_polars([o])
        sage: o.polar()
        A polytope polar to An octahedron: 3-dimensional, 8 vertices.
    """
    result_name = _palp("poly.x -fe", polytopes)
    result = open(result_name)
    for p in polytopes:
        p._read_equations(result)
    result.close()
    os.remove(result_name)

def convex_hull(points):
    r"""
    Compute the convex hull of the given points.

    NOTE: points might not span the space.
    NOTE: fails for large numbers of vertices in dimensions 4 or greater

    INPUT:
        points -- a list that can be converted into vectors of the same dimension over ZZ.

    OUTPUT:
        list of vertices of the convex hull of the given points (as vectors).

    EXAMPLES:
    Let's compute the convex hull of several points on a line in the plane:
        sage: lattice_polytope.convex_hull([[1,2],[3,4],[5,6],[7,8]])
        [(1, 2), (7, 8)]
    """
    if len(points) == 0:
        return []
    vpoints = []
    for p in points:
        v = vector(ZZ,p)
        if not v in vpoints:
            vpoints.append(v)
    p0 = vpoints[0]
    vpoints = [p-p0 for p in vpoints]
    N = ZZ**p0.degree()
    H = N.submodule(vpoints)
    if H.rank() == 0:
        return [p0]
    elif H.rank() == N.rank():
        vpoints = LatticePolytope(matrix(ZZ, vpoints).transpose(), compute_vertices=True).vertices().columns(copy=False)
    else:
        H_points = [H.coordinates(p) for p in vpoints]
        H_polytope = LatticePolytope(matrix(ZZ, H_points).transpose(), compute_vertices=True)
        vpoints = (H.basis_matrix().transpose() * H_polytope.vertices()).columns(copy=False)
    vpoints = [p+p0 for p in vpoints]
    return vpoints

def filter_polytopes(f, polytopes, subseq=None, print_numbers=False):
    r"""
    Use the function \code{f} to filter polytopes in a list.

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
        sage: polytopes = Sequence([lattice_polytope.octahedron(n) for n in range(2, 7)], cr=True)
        sage: polytopes
        [
        An octahedron: 2-dimensional, 4 vertices.,
        An octahedron: 3-dimensional, 6 vertices.,
        An octahedron: 4-dimensional, 8 vertices.,
        An octahedron: 5-dimensional, 10 vertices.,
        An octahedron: 6-dimensional, 12 vertices.
        ]

    This filters octahedrons of dimension at least 4:
        sage: lattice_polytope.filter_polytopes(lambda p: p.dim() >= 4, polytopes)
        [2, 3, 4]

    For long tests you can see the current progress:
        sage: lattice_polytope.filter_polytopes(lambda p: p.nvertices() >= 10, polytopes, print_numbers=True)
        0
        1
        2
        3
        4
        [3, 4]

    Here we consider only some of the polytopes:
        sage: lattice_polytope.filter_polytopes(lambda p: p.nvertices() >= 10, polytopes, [2, 3, 4], print_numbers=True)
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


def minkowski_sum(points1, points2):
    r"""
    Compute the Minkowski sum of two convex polytopes.

    NOTE: polytopes might not be of maximal dimension.

    INPUT:
        points1, points2 -- lists of objects that can be converted into vectors
            of the same dimension, treated as vertices of two polytopes.

    OUTPUT:
        list of vertices of the Minkowski sum, given as vectors.

    EXAMPLES:
    Let's compute the Minkovski sum of two line segments:
        sage: lattice_polytope.minkowski_sum([[1,0],[-1,0]],[[0,1],[0,-1]])
        [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    """
    points1 = [vector(p) for p in points1]
    points2 = [vector(p) for p in points2]
    points = []
    for p1 in points1:
        for p2 in points2:
            points.append(p1+p2)
    return convex_hull(points)


def octahedron(dim):
    r"""
    Return an octahedron of the given dimension.

    EXAMPLES:
    Here are 3- and 4-dimensional octahedrons:
        sage: o = lattice_polytope.octahedron(3)
        sage: o
        An octahedron: 3-dimensional, 6 vertices.
        sage: o.vertices()
        [ 1  0  0 -1  0  0]
        [ 0  1  0  0 -1  0]
        [ 0  0  1  0  0 -1]
        sage: o = lattice_polytope.octahedron(4)
        sage: o
        An octahedron: 4-dimensional, 8 vertices.
        sage: o.vertices()
        [ 1  0  0  0 -1  0  0  0]
        [ 0  1  0  0  0 -1  0  0]
        [ 0  0  1  0  0  0 -1  0]
        [ 0  0  0  1  0  0  0 -1]

    There exists only one octahedron of each dimension:
        sage: o is lattice_polytope.octahedron(4)
        True
    """
    if _octahedrons.has_key(dim):
        return _octahedrons[dim]
    else:
        _octahedrons[dim] = _create_octahedron(dim)
        return _octahedrons[dim]


def positive_integer_relations(points):
    r"""
    Return relations between given points.

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
        sage: lattice_polytope.positive_integer_relations(p.points())
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
        relations.rescale_row(i, 1/gcd(relations.row(i).list()))
    return relations.change_ring(ZZ)


def projective_space(dim):
    r"""
    Return a simplex of the given dimension, corresponding to $P_{dim}$.

    EXAMPLES:
    We construct 3- and 4-dimensional simplexes:
        sage: p = lattice_polytope.projective_space(3)
        sage: p
        A simplex: 3-dimensional, 4 vertices.
        sage: p.vertices()
        [ 1  0  0 -1]
        [ 0  1  0 -1]
        [ 0  0  1 -1]
        sage: p = lattice_polytope.projective_space(4)
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
    r"""
    Read all polytopes from the given file.

    INPUT:
        file_name -- the name of a file with VERTICES of polytopes
        desc -- a string, that will be used for creating polytope descriptions.
            By default it will be set to 'A lattice polytope \#\%d from "filename"'
+           and will be used as \code{desc \% n} where \code{n} is the number of
            the polytope in the file (\emph{STARTING WITH ZERO}).

    OUTPUT:
        a sequence of polytopes

    EXAMPLES:
    We use poly.x to compute polar polytopes of 2- and 3-octahedrons and read them:
        sage: o2 = lattice_polytope.octahedron(2)
        sage: o3 = lattice_polytope.octahedron(3)
        sage: result_name = lattice_polytope._palp("poly.x -fe", [o2, o3])
        sage: f = open(result_name)
        sage: f.readlines()
        ['4 2  Vertices of P-dual <-> Equations of P\n', '  -1   1\n', '   1   1\n', '  -1  -1\n', '   1  -1\n', '8 3  Vertices of P-dual <-> Equations of P\n', '  -1  -1   1\n', '   1  -1   1\n', '  -1   1   1\n', '   1   1   1\n', '  -1  -1  -1\n', '   1  -1  -1\n', '  -1   1  -1\n', '   1   1  -1\n']
        sage: f.close()
        sage: lattice_polytope.read_all_polytopes(result_name, desc="Polytope from file %d")
        [
        Polytope from file 0: 2-dimensional, 4 vertices.,
        Polytope from file 1: 3-dimensional, 8 vertices.
        ]
        sage: os.remove(result_name)
    """
    if desc == None:
        desc = r'A lattice polytope #%d from "'+file_name+'"'
    polytopes = Sequence([], LatticePolytope, cr=True)
    f = open(file_name)
    n = 0
    p = LatticePolytope(f, desc=(desc % n), compute_vertices=False)
    while p.dim() != 0:
        polytopes.append(p)
        n += 1
        p = LatticePolytope(f, desc=(desc % n), compute_vertices=False)
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
        sage: lattice_polytope.read_palp_matrix("2 3 comment \n 1 2 3 \n 4 5 6")
        [1 2 3]
        [4 5 6]
        sage: lattice_polytope.read_palp_matrix("3 2 Will be transposed \n 1 2 \n 3 4 \n 5 6")
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
    r"""
    Convert a SAGE matrix to the string representation of Maxima.

    EXAMPLE:
        sage: m = matrix(ZZ,2)
        sage: lattice_polytope.sage_matrix_to_maxima(m)
        matrix([0,0],[0,0])
    """
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

    EXAMPLE:
    We create a file with vertices of the square and the cube, but read only
    the second set:
        sage: o2 = lattice_polytope.octahedron(2)
        sage: o3 = lattice_polytope.octahedron(3)
        sage: result_name = lattice_polytope._palp("poly.x -fe", [o2, o3])
        sage: f = open(result_name)
        sage: f.readlines()
        ['4 2  Vertices of P-dual <-> Equations of P\n', '  -1   1\n', '   1   1\n', '  -1  -1\n', '   1  -1\n', '8 3  Vertices of P-dual <-> Equations of P\n', '  -1  -1   1\n', '   1  -1   1\n', '  -1   1   1\n', '   1   1   1\n', '  -1  -1  -1\n', '   1  -1  -1\n', '  -1   1  -1\n', '   1   1  -1\n']
        sage: f.close()
        sage: f = open(result_name)
        sage: lattice_polytope.skip_palp_matrix(f)
        sage: lattice_polytope.read_palp_matrix(f)
        [-1  1 -1  1 -1  1 -1  1]
        [-1 -1  1  1 -1 -1  1  1]
        [ 1  1  1  1 -1 -1 -1 -1]
        sage: f.close()
        sage: os.remove(result_name)
    """
    for i in range(n):
        line = data.readline()
        if line == "":
            raise ValueError, "There are not enough data to skip!"
        for j in range(int(line.split()[0])):
            if data.readline() == "":
                raise ValueError, "There are not enough data to skip!"


def write_palp_matrix(m, ofile=None, comment="", format=None):
    r"""
    Write a matrix into a file.

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
        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.write_palp_matrix(o.vertices(), comment="3D Octahedron")
        3 6 3D Octahedron
         1  0  0 -1  0  0
         0  1  0  0 -1  0
         0  0  1  0  0 -1
        sage: lattice_polytope.write_palp_matrix(o.vertices(), format="%4d")
        3 6
           1    0    0   -1    0    0
           0    1    0    0   -1    0
           0    0    1    0    0   -1
    """
    if format == None:
        n = max(len(str(m[i,j]))
                for i in range(m.nrows()) for j in range(m.ncols()))
        format = "%" + str(n) + "d"
    s = "%d %d %s\n" % (m.nrows(),m.ncols(),comment)
    if ofile is None:
        print s
    else:
        ofile.write(s)
    for i in range(m.nrows()):
        s = " ".join(format % m[i,j] for j in range(m.ncols()))+"\n"
        if ofile is None:
            print s
        else:
            ofile.write(s)
