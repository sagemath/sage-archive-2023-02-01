r"""
Lattice and reflexive polytopes

This module provides tools for work with lattice and reflexive
polytopes. A *convex polytope* is the convex hull of finitely many
points in `\RR^n`. The dimension `n` of a
polytope is the smallest `n` such that the polytope can be
embedded in `\RR^n`.

A *lattice polytope* is a polytope whose vertices all have integer
coordinates.

If `L` is a lattice polytope, the dual polytope of
`L` is

.. math::

       \{y \in \ZZ^n :   x\cdot y \geq -1 \text{ all } x \in L\}


A *reflexive polytope* is a lattice polytope, such that its polar
is also a lattice polytope, i.e. it is bounded and has vertices with
integer coordinates.

This Sage module uses Package for Analyzing Lattice Polytopes
(PALP), which is a program written in C by Maximilian Kreuzer and
Harald Skarke, which is freely available under the GNU license
terms at http://hep.itp.tuwien.ac.at/~kreuzer/CY/. Moreover, PALP is
included standard with Sage.

PALP is described in the paper :arxiv:`math.SC/0204356`. Its distribution
also contains the application nef.x, which was created by Erwin
Riegler and computes nef-partitions and Hodge data for toric
complete intersections.

ACKNOWLEDGMENT: polytope.py module written by William Stein was
used as an example of organizing an interface between an external
program and Sage. William Stein also helped Andrey Novoseltsev with
debugging and tuning of this module.

Robert Bradshaw helped Andrey Novoseltsev to realize plot3d
function.

.. note::

   IMPORTANT: PALP requires some parameters to be determined during
   compilation time, i.e., the maximum dimension of polytopes, the
   maximum number of points, etc. These limitations may lead to errors
   during calls to different functions of these module.  Currently, a
   ValueError exception will be raised if the output of poly.x or
   nef.x is empty or contains the exclamation mark. The error message
   will contain the exact command that caused an error, the
   description and vertices of the polytope, and the obtained output.

Data obtained from PALP and some other data is cached and most
returned values are immutable. In particular, you cannot change the
vertices of the polytope or their order after creation of the
polytope.

If you are going to work with large sets of data, take a look at
``all_*`` functions in this module. They precompute different data
for sequences of polynomials with a few runs of external programs.
This can significantly affect the time of future computations. You
can also use dump/load, but not all data will be stored (currently
only faces and the number of their internal and boundary points are
stored, in addition to polytope vertices and its polar).

AUTHORS:

- Andrey Novoseltsev (2007-01-11): initial version

- Andrey Novoseltsev (2007-01-15): ``all_*`` functions

- Andrey Novoseltsev (2008-04-01): second version, including:

    - dual nef-partitions and necessary convex_hull and minkowski_sum

    - built-in sequences of 2- and 3-dimensional reflexive polytopes

    - plot3d, skeleton_show

- Andrey Novoseltsev (2009-08-26): dropped maximal dimension requirement

- Andrey Novoseltsev (2010-12-15): new version of nef-partitions

- Maximilian Kreuzer and Harald Skarke: authors of PALP (which was
  also used to obtain the list of 3-dimensional reflexive polytopes)

- Erwin Riegler: the author of nef.x

"""

#*****************************************************************************
#       Copyright (C) 2007-2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2007-2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.graphs.graph import Graph
from sage.interfaces.all import maxima
from sage.matrix.all import matrix
from sage.matrix.matrix import is_Matrix
from sage.misc.all import tmp_filename
from sage.misc.misc import SAGE_SHARE
from sage.modules.all import vector
from sage.plot.plot3d.index_face_set import IndexFaceSet
from sage.plot.plot3d.all import line3d, point3d
from sage.plot.plot3d.shapes2 import text3d
from sage.rings.all import Integer, ZZ, QQ, gcd, lcm
from sage.sets.set import Set_generic
from sage.structure.all import Sequence
from sage.structure.sequence import Sequence_generic
from sage.structure.sage_object import SageObject

import collections
import copy_reg
import os
import subprocess
import StringIO


data_location = os.path.join(SAGE_SHARE,'reflexive_polytopes')


class SetOfAllLatticePolytopesClass(Set_generic):
    def _repr_(self):
        r"""
        Return a string representation.

        TESTS::

            sage: lattice_polytope.SetOfAllLatticePolytopesClass()._repr_()
            'Set of all Lattice Polytopes'
        """
        return "Set of all Lattice Polytopes"

    def __call__(self, x):
        r"""
        TESTS::

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

    LatticePolytope(data, [desc], [compute_vertices],
    [copy_vertices], [n])

    INPUT:

    - ``data`` -- The points (which must be vertices if
      ``compute_vertices`` is False) spanning the lattice polytope,
      specified as one of:

        * a matrix whose columns are vertices of the polytope.

        * an iterable of iterables (for example, a list of vectors)
          defining the point coordinates.

        * a file with matrix data, open for reading, or

        * a filename of such a file. See ``read_palp_matrix`` for the
          file format.

    -  ``desc`` -- (default: "A lattice polytope")
       description of the polytope.

    - ``compute_vertices`` -- boolean (default: True) if True, the
       convex hull of the given points will be computed for
       determining vertices. Otherwise, the given points must be
       vertices.

    - ``copy_vertices`` -- boolean (default: True) if False, and
       compute_vertices is False, and ``data`` is a matrix of
       vertices, it will be made immutable.

    - ``n`` -- integer (default: 0) if ``data`` is a name of a file,
       that contains data blocks for several polytopes, the n-th block
       will be used. *ENUMERATION STARTS WITH ZERO*.

    OUTPUT: a lattice polytope

    EXAMPLES: Here we construct a polytope from a matrix whose columns
    are vertices in 3-dimensional space. In the first case a copy of
    the given matrix is made during construction, in the second one the
    matrix is made immutable and used as a matrix of vertices.

    ::

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

    We draw a pretty picture of the polytype in 3-dimensional space::

        sage: p.plot3d().show()

    Now we add an extra point, which is in the interior of the
    polytope...

    ::

        sage: p = LatticePolytope(m.columns() + [(0,0,0)], "A lattice polytope constructed from 7 points")
        sage: p
        A lattice polytope constructed from 7 points: 3-dimensional, 6 vertices.

    You can suppress vertex computation for speed but this can lead to
    mistakes::

        sage: p = LatticePolytope(m.columns() + [(0,0,0)], "A lattice polytope with WRONG vertices",
        ...                         compute_vertices=False)
        ...
        sage: p
        A lattice polytope with WRONG vertices: 3-dimensional, 7 vertices.

    Given points must be in the lattice::

        sage: LatticePolytope(matrix([1/2, 3/2]))
        Traceback (most recent call last):
        ...
        ValueError: Points must be integral!
        Given:
        [1/2 3/2]

    But it is OK to create polytopes of non-maximal dimension::

        sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0, 0],
        ...                   [0, 1, 0,  0, -1,  0, 0],
        ...                   [0, 0, 0,  0,  0,  0, 0]])
        ...
        sage: p = LatticePolytope(m)
        sage: p
        A lattice polytope: 2-dimensional, 4 vertices.
        sage: p.vertices()
        [ 1  0 -1  0]
        [ 0  1  0 -1]
        [ 0  0  0  0]

    An empty lattice polytope can be specified by a matrix with zero columns:

        sage: p = LatticePolytope(matrix(ZZ,3,0)); p
        A lattice polytope: -1-dimensional, 0 vertices.
        sage: p.ambient_dim()
        3
        sage: p.npoints()
        0
        sage: p.nfacets()
        0
        sage: p.points()
        []
        sage: p.faces()
        []
    """
    if isinstance(data, LatticePolytopeClass):
        return data

    if not is_Matrix(data) and not isinstance(data,(file,basestring)):
        try:
            data = matrix(ZZ,data).transpose()
        except ValueError:
            pass

    return LatticePolytopeClass(data, desc, compute_vertices, copy_vertices, n)

copy_reg.constructor(LatticePolytope)   # "safe for unpickling"

def ReflexivePolytope(dim, n):
    r"""
    Return n-th reflexive polytope from the database of 2- or
    3-dimensional reflexive polytopes.

    .. note::

       #. Numeration starts with zero: `0 \leq n \leq 15` for `{\rm dim} = 2`
          and `0 \leq n \leq 4318` for `{\rm dim} = 3`.

       #. During the first call, all reflexive polytopes of requested
          dimension are loaded and cached for future use, so the first
          call for 3-dimensional polytopes can take several seconds,
          but all consecutive calls are fast.

       #. Equivalent to ``ReflexivePolytopes(dim)[n]`` but checks bounds
          first.

    EXAMPLES: The 3rd 2-dimensional polytope is "the diamond:"

    ::

        sage: ReflexivePolytope(2,3)
        Reflexive polytope    3: 2-dimensional, 4 vertices.
        sage: lattice_polytope.ReflexivePolytope(2,3).vertices()
        [ 1  0  0 -1]
        [ 0  1 -1  0]

    There are 16 reflexive polygons and numeration starts with 0::

        sage: ReflexivePolytope(2,16)
        Traceback (most recent call last):
        ...
        ValueError: there are only 16 reflexive polygons!

    It is not possible to load a 4-dimensional polytope in this way::

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
_rp = [None]*4

def ReflexivePolytopes(dim):
    r"""
    Return the sequence of all 2- or 3-dimensional reflexive polytopes.

    .. note::

       During the first call the database is loaded and cached for
       future use, so repetitive calls will return the same object in
       memory.

    :param dim: dimension of required reflexive polytopes
    :type dim: 2 or 3
    :rtype: list of lattice polytopes

    EXAMPLES: There are 16 reflexive polygons::

        sage: len(ReflexivePolytopes(2))
        16

    It is not possible to load 4-dimensional polytopes in this way::


        sage: ReflexivePolytopes(4)
        Traceback (most recent call last):
        ...
        NotImplementedError: only 2- and 3-dimensional reflexive polytopes are available!
    """
    global _rp
    if dim not in [2, 3]:
        raise NotImplementedError, "only 2- and 3-dimensional reflexive polytopes are available!"
    if _rp[dim] == None:
        rp = read_all_polytopes(
            os.path.join(data_location,"reflexive_polytopes_%dd"%dim),
            desc="Reflexive polytope %4d")
        for n, p in enumerate(rp):
            # Data files have normal form of reflexive polytopes
            p._normal_form = p._vertices
            # Prevents dimension computation later
            p._dim = dim
        # Compute "fast" data in one call to PALP
        all_polars(rp)
        all_points(rp + [p._polar for p in rp])
        # TODO: improve faces representation so that we can uncomment
        # all_faces(rp)
        # It adds ~10s for dim=3, which is a bit annoying to wait for.
        _rp[dim] = rp
    return _rp[dim]


def is_LatticePolytope(x):
    r"""
    Check if ``x`` is a lattice polytope.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a :class:`lattice polytope <LatticePolytopeClass>`,
      ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.lattice_polytope import is_LatticePolytope
        sage: is_LatticePolytope(1)
        False
        sage: p = LatticePolytope([(1,0), (0,1), (-1,-1)])
        sage: p
        A lattice polytope: 2-dimensional, 3 vertices.
        sage: is_LatticePolytope(p)
        True
    """
    return isinstance(x, LatticePolytopeClass)


class LatticePolytopeClass(SageObject, collections.Hashable):
    r"""
    Class of lattice/reflexive polytopes.

    Use ``LatticePolytope`` for constructing a polytope.
    """

    def __init__(self, data, desc, compute_vertices, copy_vertices=True, n=0):
        r"""
        Construct a lattice polytope. See ``LatticePolytope``.

        Most tests/examples are also done in LatticePolytope.

        TESTS::

            sage: LatticePolytope(matrix(ZZ,[[1,2,3],[4,5,6]]))
            A lattice polytope: 1-dimensional, 2 vertices.
        """
        if is_Matrix(data):
            if data.base_ring() != ZZ:
                try:
                    data = matrix(ZZ, data)
                except TypeError:
                    raise ValueError("Points must be integral!\nGiven:\n%s" %data)
            if desc == None:
                self._desc = "A lattice polytope"
            else:
                self._desc = desc
            if compute_vertices:
                self._vertices = data       # Temporary assignment
                self._compute_dim(compute_vertices=True)
            else:
                if copy_vertices:
                    self._vertices = data.__copy__()
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

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        INPUT:

        - ``other`` -- anything.

        OUTPUT:

        - ``True`` if ``other`` is a :class:`lattice polytope
          <LatticePolytopeClass>` equal to ``self``, ``False`` otherwise.

        .. NOTE::

            Two lattice polytopes are if they have the same vertices listed in
            the same order.

        TESTS::

            sage: p1 = LatticePolytope([(1,0), (0,1), (-1,-1)])
            sage: p2 = LatticePolytope([(1,0), (0,1), (-1,-1)])
            sage: p3 = LatticePolytope([(0,1), (1,0), (-1,-1)])
            sage: p1 == p1
            True
            sage: p1 == p2
            True
            sage: p1 is p2
            False
            sage: p1 == p3
            False
            sage: p1 == 0
            False
        """
        return (is_LatticePolytope(other)
                and self._vertices == other._vertices)

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        OUTPUT:

        - an integer.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: hash(o) == hash(o)
            True
        """
        try:
            return self._hash
        except AttributeError:
            self._hash = hash(self._vertices)
            return self._hash

    def __ne__(self, other):
        r"""
        Compare ``self`` with ``other``.

        INPUT:

        - ``other`` -- anything.

        OUTPUT:

        - ``False`` if ``other`` is a :class:`lattice polytope
          <LatticePolytopeClass>` equal to ``self``, ``True`` otherwise.

        .. NOTE::

            Two lattice polytopes are if they have the same vertices listed in
            the same order.

        TESTS::

            sage: p1 = LatticePolytope([(1,0), (0,1), (-1,-1)])
            sage: p2 = LatticePolytope([(1,0), (0,1), (-1,-1)])
            sage: p3 = LatticePolytope([(0,1), (1,0), (-1,-1)])
            sage: p1 != p1
            False
            sage: p1 != p2
            False
            sage: p1 is p2
            False
            sage: p1 != p3
            True
            sage: p1 != 0
            True
        """
        return not (self == other)

    def __reduce__(self):
        r"""
        Reduction function. Does not store data that can be relatively fast
        recomputed.

        TESTS::

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

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.vertices() == loads(o.dumps()).vertices()
            True
        """
        self.__dict__.update(state)
        if state.has_key('_faces'):     # Faces do not remember polytopes
            for d_faces in self._faces:
                for face in d_faces:
                    face._polytope = self

    def _compute_dim(self, compute_vertices):
        r"""
        Compute the dimension of this polytope and its vertices, if necessary.

        If ``compute_vertices`` is ``True``, then ``self._vertices`` should
        contain points whose convex hull will be computed and placed back into
        ``self._vertices``.

        If the dimension of this polytope is not equal to its ambient dimension,
        auxiliary polytope will be created and stored for using PALP commands.

        TESTS::

            sage: p = LatticePolytope(matrix([1, 2, 3]), compute_vertices=False)
            sage: p.vertices() # wrong, since these were not vertices
            [1 2 3]
            sage: hasattr(p, "_dim")
            False
            sage: p._compute_dim(compute_vertices=True)
            sage: p.vertices()
            [1 3]
            sage: p._dim
            1
        """
        if hasattr(self, "_dim"):
            return
        if self._vertices.ncols()==0:  # the empty lattice polytope
            self._dim = -1
            return
        if compute_vertices:
            points = []
            for point in self._vertices.columns(copy=False):
                if not point in points:
                    points.append(point)
        else:
            points = self._vertices.columns(copy=False)
        p0 = points[0]
        if p0 != 0:
            points = [point - p0 for point in points]
        N = ZZ**self.ambient_dim()
        H = N.submodule(points)
        self._dim = H.rank()
        if self._dim == 0:
            if compute_vertices:
                self._vertices = matrix(p0).transpose()
        elif self._dim == self.ambient_dim():
            if compute_vertices:
                if len(points) < self._vertices.ncols():  # Repeated vertices
                    if p0 != 0:
                        points = [point + p0 for point in points]
                    self._vertices = matrix(ZZ, points).transpose()
                self._vertices = read_palp_matrix(self.poly_x("v"))
        else:
            # Setup auxiliary polytope and maps
            H = H.saturation()
            H_points = [H.coordinates(point) for point in points]
            H_polytope = LatticePolytope(matrix(ZZ, H_points).transpose(),
                compute_vertices=True)
            self._sublattice = H
            self._sublattice_polytope = H_polytope
            self._embedding_matrix = H.basis_matrix().transpose()
            self._shift_vector = p0
            if compute_vertices:
                self._vertices = self._embed(H_polytope.vertices())
            # In order to use facet normals obtained from subpolytopes, we
            # need the following (see Trac #9188).
            M = self._embedding_matrix
            # Basis for the ambient space with spanned subspace in front
            basis = M.columns() + M.integer_kernel().basis()
            # Let's represent it as columns of a matrix
            basis = matrix(basis).transpose()
            # Absolute value helps to keep normals "inner"
            self._dual_embedding_scale = abs(basis.det())
            dualbasis = matrix(ZZ, self._dual_embedding_scale * basis.inverse())
            self._dual_embedding_matrix = dualbasis.submatrix(0,0,M.ncols())

    def _compute_faces(self):
        r"""
        Compute and cache faces of this polytope.

        If this polytope is reflexive and the polar polytope was already
        computed, computes faces of both in order to save time and preserve
        the one-to-one correspondence between the faces of this polytope of
        dimension d and the faces of the polar polytope of codimension
        d+1.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: v = o.__dict__.pop("_faces", None) # faces may be cached already
            sage: o.__dict__.has_key("_faces")
            False
            sage: o._compute_faces()
            sage: o.__dict__.has_key("_faces")
            True

        Check that Trac 8934 is fixed::

            sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0, 0],
            ...                   [0, 1, 0,  0, -1,  0, 0],
            ...                   [0, 0, 0,  0,  0,  0, 0]])
            ...
            sage: p = LatticePolytope(m)
            sage: p._compute_faces()
            sage: p.facets()
            [[0, 3], [2, 3], [0, 1], [1, 2]]
        """
        if hasattr(self, "_constructed_as_polar"):
                # "Polar of polar polytope" computed by poly.x may have the
                # order of vertices different from the original polytope. Thus,
                # in order to have consistent enumeration of vertices and faces
                # we must run poly.x on the original polytope.
                self._copy_faces(self._polar, reverse=True)
        elif hasattr(self, "_constructed_as_affine_transform"):
                self._copy_faces(self._original)
        elif self.dim() <= 0:
            self._faces = []
        else:
            self._read_faces(self.poly_x("i", reduce_dimension=True))

    def _compute_hodge_numbers(self):
        r"""
        Compute Hodge numbers for the current nef_partitions.

        This function (currently) always raises an exception directing to
        use another way for computing Hodge numbers.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: o._compute_hodge_numbers()
            Traceback (most recent call last):
            ...
            NotImplementedError: use nef_partitions(hodge_numbers=True)!
        """
        raise NotImplementedError, "use nef_partitions(hodge_numbers=True)!"

    def _copy_faces(self, other, reverse=False):
        r"""
        Copy facial structure of another polytope.

        This may be necessary for preserving natural correspondence of faces,
        e.g. between this polytope and its multiple or translation. In case of
        reflexive polytopes, faces of this polytope and its polar are in
        inclusion reversing bijection.

        .. note::

            This function does not perform any checks that this operation makes
            sense.

        INPUT:

        - ``other`` - another LatticePolytope, whose facial structure will be
          copied

        -  ``reverse`` - (default: False) if True, the facial structure of the
          other polytope will be reversed, i.e. vertices will correspond to
          facets etc.

        TESTS::

            sage: o = lattice_polytope.octahedron(2)
            sage: p = LatticePolytope(o.vertices())
            sage: p._copy_faces(o)
            sage: str(o.faces()) == str(p.faces())
            True
            sage: c = o.polar()
            sage: p._copy_faces(c, reverse=True)
            sage: str(o.faces()) == str(p.faces())
            True
        """
        self._faces = Sequence([], cr=True)
        if reverse:
            for d_faces in reversed(other.faces()):
                self._faces.append([_PolytopeFace(self, f._facets, f._vertices)
                                     for f in d_faces])
        else:
            for d_faces in other.faces():
                self._faces.append([_PolytopeFace(self, f._vertices, f._facets)
                                     for f in d_faces])
        self._faces.set_immutable()

    def _embed(self, data):
        r"""
        Embed given point(s) into the ambient space of this polytope.

        INPUT:

        - ``data`` - point or matrix of points (as columns) in the affine
          subspace spanned by this polytope

        OUTPUT: The same point(s) in the coordinates of the ambient space of
        this polytope.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: o._embed(o.vertices()) == o.vertices()
            True
            sage: m = matrix(ZZ, 3)
            sage: m[0, 0] = 1
            sage: m[1, 1] = 1
            sage: p = o.affine_transform(m)
            sage: p._embed((0,0))
            (1, 0, 0)
        """
        if self.ambient_dim() == self.dim():
            return data
        elif is_Matrix(data):
            r = self._embedding_matrix * data
            for i, col in enumerate(r.columns(copy=False)):
                r.set_column(i, col + self._shift_vector)
            return r
        else:
            return self._embedding_matrix * vector(QQ, data) + self._shift_vector

    def _face_compute_points(self, face):
        r"""
        Compute and cache lattice points of the given ``face``
        of this polytope.

        TESTS::

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
        Compute and cache boundary and interior lattice points of
        ``face``.

        TESTS::

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
            [10]
            sage: f._boundary_points
            [0, 2, 4, 6, 8, 9, 11, 12]
            sage: f.points()
            [0, 2, 4, 6, 8, 9, 10, 11, 12]

        Vertices don't have boundary::

            sage: f = c.faces(dim=0)[0]
            sage: c._face_split_points(f)
            sage: len(f._interior_points)
            1
            sage: len(f._boundary_points)
            0
        """
        if face.npoints() == 1: # Vertex
            face._interior_points = face.points()
            face._boundary_points = Sequence([], int, check=False)
        else:
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

    def _latex_(self):
        r"""
        Return the latex representation of self.

        OUTPUT:

        - string

        EXAMPLES:

        Arbitrary lattice polytopes are printed as `\Delta^d`, where `d` is
        the (actual) dimension of the polytope::

            sage: LatticePolytope(matrix(2, [1, 0, 1, 0]))._latex_()
            '\\Delta^{1}'

        For reflexive polytopes the output is the same... ::

            sage: o = lattice_polytope.octahedron(2)
            sage: o._latex_()
            '\\Delta^{2}'

        ... unless they are written in the normal from in which case the index
        in the internal database is printed as a subscript::

            sage: o.vertices()
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: o = LatticePolytope(o.normal_form())
            sage: o.vertices()
            [ 1  0  0 -1]
            [ 0  1 -1  0]
            sage: o._latex_()
            '\\Delta^{2}_{3}'
        """
        if (self.dim() in (2, 3)
            and self.is_reflexive()
            and self.normal_form() == self.vertices()):
            return r"\Delta^{%d}_{%d}" % (self.dim(), self.index())
        else:
            return r"\Delta^{%d}" % self.dim()

    def _palp(self, command, reduce_dimension=False):
        r"""
        Run ``command`` on vertices of this polytope.

        Returns the output of ``command`` as a string.

        .. note::

          PALP cannot be called for polytopes that do not span the ambient space.
          If you specify ``reduce_dimension=True`` argument, PALP will be
          called for vertices of this polytope in some basis of the affine space
          it spans.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: o._palp("poly.x -f")
            'M:7 6 N:27 8 Pic:17 Cor:0\n'
            sage: print o._palp("nef.x -f -N -p") # random time information
            M:27 8 N:7 6  codim=2 #part=5
            H:[0] P:0 V:2 4 5       0sec  0cpu
            H:[0] P:2 V:3 4 5       0sec  0cpu
            H:[0] P:3 V:4 5       0sec  0cpu
            np=3 d:1 p:1    0sec     0cpu

            sage: p = LatticePolytope(matrix([1]))
            sage: p._palp("poly.x -f")
            Traceback (most recent call last):
            ...
            ValueError: Cannot run "poly.x -f" for the zero-dimensional polytope!
            Polytope: A lattice polytope: 0-dimensional, 1 vertices.

            sage: m = matrix(ZZ, [[1, 0, -1,  0],
            ...                   [0, 1,  0, -1],
            ...                   [0, 0,  0,  0]])
            ...
            sage: p = LatticePolytope(m)
            sage: p._palp("poly.x -f")
            Traceback (most recent call last):
            ...
            ValueError: Cannot run PALP for a 2-dimensional polytope in a 3-dimensional space!
            sage: p._palp("poly.x -f", reduce_dimension=True)
            'M:5 4 F:4\n'
        """
        if self.dim() <= 0:
            raise ValueError, ("Cannot run \"%s\" for the zero-dimensional "
                + "polytope!\nPolytope: %s") % (command, self)
        if self.dim() < self.ambient_dim() and not reduce_dimension:
            raise ValueError(("Cannot run PALP for a %d-dimensional polytope " +
            "in a %d-dimensional space!") % (self.dim(), self.ambient_dim()))
        if _always_use_files:
            fn = _palp(command, [self], reduce_dimension)
            f = open(fn)
            result = f.read()
            f.close()
            os.remove(fn)
        else:
            if _palp_dimension is not None:
                dot = command.find(".")
                command = (command[:dot] + "-%dd" % _palp_dimension +
                           command[dot:])
            p = subprocess.Popen(command, shell=True, bufsize=2048,
                                 stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 close_fds=True)
            stdin, stdout, stderr = (p.stdin, p.stdout, p.stderr)
            write_palp_matrix(self._pullback(self._vertices), stdin)
            stdin.close()
            err = stderr.read()
            if len(err) > 0:
                raise RuntimeError, ("Error executing \"%s\" for the given polytope!"
                    + "\nPolytope: %s\nVertices:\n%s\nOutput:\n%s") % (command,
                    self, self.vertices(), err)
            result = stdout.read()
            # We program around an issue with subprocess (this so far seems to
            # only be an issue on Cygwin).
            try:
                p.terminate()
            except OSError:
                pass
        if (not result or
            "!" in result or
            "failed." in result or
            "increase" in result or
            "Unable" in result):
            raise ValueError, ("Error executing \"%s\" for the given polytope!"
                + "\nPolytope: %s\nVertices:\n%s\nOutput:\n%s") % (command,
                self, self.vertices(), result)
        return result

    def _pullback(self, data):
        r"""
        Pull back given point(s) to the affine subspace spanned by this polytope.

        INPUT:

        - ``data`` - rational point or matrix of points (as columns) in the
        ambient space

        OUTPUT: The same point(s) in the coordinates of the affine subspace
        space spanned by this polytope.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: o._pullback(o.vertices()) == o.vertices()
            True
            sage: m = matrix(ZZ, 3)
            sage: m[0, 0] = 1
            sage: m[1, 1] = 1
            sage: p = o.affine_transform(m)
            sage: p._pullback((0, 0, 0))
            [-1, 0]
        """
        if self.ambient_dim() == self.dim():
            return data
        elif data is self._vertices:
            return self._sublattice_polytope._vertices
        elif is_Matrix(data):
            r = matrix([self._pullback(col)
                    for col in data.columns(copy=False)]).transpose()
            return r
        else:
            data = vector(QQ, data)
            return self._sublattice.coordinates(data - self._shift_vector)

    def _read_equations(self, data):
        r"""
        Read equations of facets/vertices of polar polytope from string or
        file.

        TESTS: For a reflexive polytope construct the polar polytope::

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

        For a non-reflexive polytope cache facet equations::

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
        if hasattr(self, "_is_reflexive"):
            # If it is already known that this polytope is reflexive, its
            # polar (whose vertices are equations of facets of this one)
            # is already computed and there is no need to read equations
            # of facets of this polytope. Moreover, doing so can corrupt
            # data if this polytope was constructed as polar. Skip input.
            skip_palp_matrix(data)
            return
        pos = data.tell()
        line = data.readline()
        self._is_reflexive = line.find("Vertices of P-dual") != -1
        if self._is_reflexive:
            data.seek(pos)
            self._polar = LatticePolytope(read_palp_matrix(data),
                            "A polytope polar to " + str(self._desc),
                            copy_vertices=False, compute_vertices=False)
            self._polar._dim = self._dim
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
        Read faces information from string or file.

        TESTS::

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
            [[0], [1], [2], [3], [4], [5]],
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]],
            [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]
            ]

        Cannot be used for "polar polytopes," their faces are constructed
        from faces of the original one to preserve facial duality.

        ::

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
        # Zero-dimensional faces (i.e. vertices) from poly.x can be in "random"
        # order, so that the lists of corresponding facets are in increasing
        # order.
        # While this may be convenient for something, it is quite confusing to
        # have p.faces(dim=0)[0].vertices() == [5], which means "the 5th vertex
        # spans the 0th 0-dimensional face" and, on the polar side, "the 0th
        # facet is described by the 5th equation."
        # The next line sorts 0-dimensional faces to make these enumerations
        # more transparent.
        self._faces[0].sort(cmp = lambda x,y: cmp(x._vertices[0], y._vertices[0]))
        self._faces.set_immutable()

    def _read_nef_partitions(self, data):
        r"""
        Read nef-partitions of ``self`` from ``data``.

        INPUT:

        - ``data`` -- a string or a file.

        OUTPUT:

        - none.

        TESTS::

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

        We make a copy of the octahedron since direct use of this function may
        destroy cache integrity and lead so strange effects in other doctests::

            sage: o_copy = LatticePolytope(o.vertices())
            sage: o_copy.__dict__.has_key("_nef_partitions")
            False
            sage: o_copy._read_nef_partitions(s)
            sage: o_copy._nef_partitions
            [
            Nef-partition {0, 1, 3} U {2, 4, 5},
            Nef-partition {0, 1, 2} U {3, 4, 5},
            Nef-partition {0, 1, 2, 3} U {4, 5}
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
            partition = [0] * nvertices
            for v in pvertices:
                partition[v] = 1
            partition = NefPartition(partition, self)
            partition._is_product = line.find(" D ") != -1
            partition._is_projection = line.find(" DP ") != -1
            # Set the stuff
            start = line.find("H:")
            if start != -1:
                start += 2
                end = line.find("[", start)
                partition._hodge_numbers = tuple(int(h)
                                            for h in line[start:end].split())
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

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: o._repr_()
            'An octahedron: 3-dimensional, 6 vertices.'
        """
        s = str(self._desc)
        s += ": %d-dimensional, %d vertices." % (self.dim(), self.nvertices())
        return s

    def affine_transform(self, a=1, b=0):
        r"""
        Return a*P+b, where P is this lattice polytope.

        .. note::

          #. While ``a`` and ``b`` may be rational, the final result must be a
          lattice polytope, i.e. all vertices must be integral.

          #. If the transform (restricted to this polytope) is bijective, facial
          structure will be preserved, e.g. the first facet of the image will be
          spanned by the images of vertices which span the first facet of the
          original polytope.

        INPUT:

        - ``a`` - (default: 1) rational scalar or matrix

        - ``b`` - (default: 0) rational scalar or vector, scalars are
          interpreted as vectors with the same components

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(2)
            sage: o.vertices()
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: o.affine_transform(2).vertices()
            [ 2  0 -2  0]
            [ 0  2  0 -2]
            sage: o.affine_transform(1,1).vertices()
            [2 1 0 1]
            [1 2 1 0]
            sage: o.affine_transform(b=1).vertices()
            [2 1 0 1]
            [1 2 1 0]
            sage: o.affine_transform(b=(1, 0)).vertices()
            [ 2  1  0  1]
            [ 0  1  0 -1]
            sage: a = matrix(QQ, 2, [1/2, 0, 0, 3/2])
            sage: o.polar().vertices()
            [-1  1 -1  1]
            [ 1  1 -1 -1]
            sage: o.polar().affine_transform(a, (1/2, -1/2)).vertices()
            [ 0  1  0  1]
            [ 1  1 -2 -2]

        While you can use rational transformation, the result must be integer::

            sage: o.affine_transform(a)
            Traceback (most recent call last):
            ...
            ValueError: Points must be integral!
            Given:
            [ 1/2    0 -1/2    0]
            [   0  3/2    0 -3/2]
        """
        new_vertices = a * self.vertices()
        if b in QQ:
            b = vector(QQ, [b]*new_vertices.nrows())
        else:
            b = vector(QQ, b)
        for i, col in enumerate(new_vertices.columns(copy=False)):
            new_vertices.set_column(i, col + b)
        r = LatticePolytope(new_vertices, compute_vertices=True)
        if (a in QQ and a != 0) or r.dim() == self.dim():
            r._constructed_as_affine_transform = True
            if hasattr(self, "_constructed_as_affine_transform"):
                # Prevent long chains of "original-transform"
                r._original = self._original
            else:
                r._original = self
        return r

    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space of this polytope.

        EXAMPLES: We create a 3-dimensional octahedron and check its
        ambient dimension::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.ambient_dim()
            3
        """
        return self._vertices.nrows()

    def dim(self):
        r"""
        Return the dimension of this polytope.

        EXAMPLES: We create a 3-dimensional octahedron and check its
        dimension::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.dim()
            3

        Now we create a 2-dimensional diamond in a 3-dimensional space::

            sage: m = matrix(ZZ, [[1, 0, -1,  0],
            ...                   [0, 1,  0, -1],
            ...                   [0, 0,  0,  0]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.dim()
            2
            sage: p.ambient_dim()
            3
        """
        if not hasattr(self, "_dim"):
            self._compute_dim(compute_vertices=False)
        return self._dim

    def distances(self, point=None):
        r"""
        Return the matrix of distances for this polytope or distances for
        the given point.

        The matrix of distances m gives distances m[i,j] between the i-th
        facet (which is also the i-th vertex of the polar polytope in the
        reflexive case) and j-th point of this polytope.

        If point is specified, integral distances from the point to all
        facets of this polytope will be computed.

        This function CAN be used for polytopes whose dimension is smaller
        than the dimension of the ambient space. In this case distances are
        computed in the affine subspace spanned by the polytope and if the
        point is given, it must be in this subspace.

        EXAMPLES: The matrix of distances for a 3-dimensional octahedron::

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

        Distances from facets to the point (1,2,3)::

            sage: o.distances([1,2,3])
            (1, 3, 5, 7, -5, -3, -1, 1)

        It is OK to use RATIONAL coordinates::

            sage: o.distances([1,2,3/2])
            (-1/2, 3/2, 7/2, 11/2, -7/2, -3/2, 1/2, 5/2)
            sage: o.distances([1,2,sqrt(2)])
            Traceback (most recent call last):
            ...
            TypeError: unable to convert sqrt(2) to a rational

        Now we create a non-spanning polytope::

            sage: m = matrix(ZZ, [[1, 0, -1,  0],
            ...                   [0, 1,  0, -1],
            ...                   [0, 0,  0,  0]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.distances()
            [0 2 2 0 1]
            [2 2 0 0 1]
            [0 0 2 2 1]
            [2 0 0 2 1]
            sage: p.distances((1/2, 3, 0))
            (7/2, 9/2, -5/2, -3/2)
            sage: p.distances((1, 1, 1))
            Traceback (most recent call last):
            ...
            ArithmeticError: vector is not in free module
        """
        if point != None:
            if self.dim() < self.ambient_dim():
                return self._sublattice_polytope.distances(self._pullback(point))
            elif self.is_reflexive():
                return (self._polar._vertices.transpose() * vector(QQ, point)
                        + vector(QQ, [1]*self.nfacets()))
            else:
                return self._facet_normals*vector(QQ, point)+self._facet_constants
        if not hasattr(self, "_distances"):
            if self.dim() < self.ambient_dim():
                return self._sublattice_polytope.distances()
            elif self.is_reflexive():
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
        Return the sequence of edges of this polytope (i.e. faces of
        dimension 1).

        EXAMPLES: The octahedron has 12 edges::

            sage: o = lattice_polytope.octahedron(3)
            sage: len(o.edges())
            12
            sage: o.edges()
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]]
        """
        return self.faces(dim=1)

    def faces(self, dim=None, codim=None):
        r"""
        Return the sequence of proper faces of this polytope.

        If ``dim`` or ``codim`` are specified,
        returns a sequence of faces of the corresponding dimension or
        codimension. Otherwise returns the sequence of such sequences for
        all dimensions.

        EXAMPLES: All faces of the 3-dimensional octahedron::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.faces()
            [
            [[0], [1], [2], [3], [4], [5]],
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]],
            [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]
            ]

        Its faces of dimension one (i.e., edges)::

            sage: o.faces(dim=1)
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]]

        Its faces of codimension two (also edges)::

            sage: o.faces(codim=2)
            [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]]

        It is an error to specify both dimension and codimension at the
        same time, even if they do agree::

            sage: o.faces(dim=1, codim=2)
            Traceback (most recent call last):
            ...
            ValueError: Both dim and codim are given!

        The only faces of a zero-dimensional polytope are the empty set and
        the polytope itself, i.e. it has no proper faces at all::

            sage: p = LatticePolytope(matrix([[1]]))
            sage: p.vertices()
            [1]
            sage: p.faces()
            []

        In particular, you an exception will be raised if you try to access
        faces of the given dimension or codimension, including edges and
        facets::

            sage: p.facets()
            Traceback (most recent call last):
            ...
            IndexError: list index out of range
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

    def facet_constant(self, i):
        r"""
        Return the constant in the ``i``-th facet inequality of this polytope.

        The i-th facet inequality is given by
        self.facet_normal(i) * X + self.facet_constant(i) >= 0.

        INPUT:

        - ``i`` - integer, the index of the facet

        OUTPUT:

        - integer -- the constant in the ``i``-th facet inequality.

        EXAMPLES:

        Let's take a look at facets of the octahedron and some polytopes
        inside it::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.vertices()
            [ 1  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: o.facet_normal(0)
            (-1, -1, 1)
            sage: o.facet_constant(0)
            1
            sage: m = copy(o.vertices())
            sage: m[0,0] = 0
            sage: m
            [ 0  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: p = LatticePolytope(m)
            sage: p.facet_normal(0)
            (-1, 0, 0)
            sage: p.facet_constant(0)
            0
            sage: m[0,3] = 0
            sage: m
            [ 0  0  0  0  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: p = LatticePolytope(m)
            sage: p.facet_normal(0)
            (0, -1, 1)
            sage: p.facet_constant(0)
            1

        This is a 2-dimensional lattice polytope in a 4-dimensional space::

            sage: m = matrix([(1,-1,1,3), (-1,-1,1,3), (0,0,0,0)])
            sage: p = LatticePolytope(m.transpose())
            sage: p
            A lattice polytope: 2-dimensional, 3 vertices.
            sage: p.vertices()
            [ 1 -1  0]
            [-1 -1  0]
            [ 1  1  0]
            [ 3  3  0]
            sage: fns = [p.facet_normal(i) for i in range(p.nfacets())]
            sage: fns
            [(11, -1, 1, 3), (-11, -1, 1, 3), (0, 1, -1, -3)]
            sage: fcs = [p.facet_constant(i) for i in range(p.nfacets())]
            sage: fcs
            [0, 0, 11]

        Now we manually compute the distance matrix of this polytope. Since it
        is a triangle, each line (corresponding to a facet) should have two
        zeros (vertices of the corresponding facet) and one positive number
        (since our normals are inner)::

            sage: matrix([[fns[i] * p.vertex(j) + fcs[i]
            ...            for j in range(p.nvertices())]
            ...           for i in range(p.nfacets())])
            [22  0  0]
            [ 0 22  0]
            [ 0  0 11]
        """
        if self.is_reflexive():
            return 1
        elif self.ambient_dim() == self.dim():
            return self._facet_constants[i]
        else:
            return (self._sublattice_polytope.facet_constant(i)
                    * self._dual_embedding_scale
                    - self.facet_normal(i) * self._shift_vector)

    def facet_normal(self, i):
        r"""
        Return the inner normal to the ``i``-th facet of this polytope.

        If this polytope is not full-dimensional, facet normals will be
        orthogonal to the integer kernel of the affine subspace spanned by
        this polytope.

        INPUT:

        - ``i`` -- integer, the index of the facet

        OUTPUT:

        - vectors -- the inner normal of the ``i``-th facet

        EXAMPLES:

        Let's take a look at facets of the octahedron and some polytopes
        inside it::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.vertices()
            [ 1  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: o.facet_normal(0)
            (-1, -1, 1)
            sage: o.facet_constant(0)
            1
            sage: m = copy(o.vertices())
            sage: m[0,0] = 0
            sage: m
            [ 0  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: p = LatticePolytope(m)
            sage: p.facet_normal(0)
            (-1, 0, 0)
            sage: p.facet_constant(0)
            0
            sage: m[0,3] = 0
            sage: m
            [ 0  0  0  0  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: p = LatticePolytope(m)
            sage: p.facet_normal(0)
            (0, -1, 1)
            sage: p.facet_constant(0)
            1

        Here is an example of a 3-dimensional polytope in a 4-dimensional
        space::

            sage: m = matrix([(0,0,0,0), (1,1,1,3), (1,-1,1,3), (-1,-1,1,3)])
            sage: p = LatticePolytope(m.transpose())
            sage: p.vertices()
            [ 0  1  1 -1]
            [ 0  1 -1 -1]
            [ 0  1  1  1]
            [ 0  3  3  3]
            sage: ker = p.vertices().integer_kernel().matrix()
            sage: ker
            [ 0  0  3 -1]
            sage: [ker * p.facet_normal(i) for i in range(p.nfacets())]
            [(0), (0), (0), (0)]

        Now we manually compute the distance matrix of this polytope. Since it
        is a simplex, each line (corresponding to a facet) should consist of
        zeros (indicating generating vertices of the corresponding facet) and
        a single positive number (since our normals are inner)::

            sage: matrix([[p.facet_normal(i) * p.vertex(j)
            ...            + p.facet_constant(i)
            ...            for j in range(p.nvertices())]
            ...           for i in range(p.nfacets())])
            [ 0  0  0 20]
            [ 0  0 20  0]
            [ 0 20  0  0]
            [10  0  0  0]
        """
        if self.is_reflexive():
            return self.polar().vertex(i)
        elif self.ambient_dim() == self.dim():
            return self._facet_normals[i]
        else:
            return ( self._sublattice_polytope.facet_normal(i) *
                     self._dual_embedding_matrix )

    def facets(self):
        r"""
        Return the sequence of facets of this polytope (i.e. faces of
        codimension 1).

        EXAMPLES: All facets of the 3-dimensional octahedron::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.facets()
            [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]

        Facets are the same as faces of codimension one::

            sage: o.facets() is o.faces(codim=1)
            True
        """
        return self.faces(codim=1)

    # Dictionaries of normal forms
    _rp_dict = [None]*4

    def index(self):
        r"""
        Return the index of this polytope in the internal database of 2- or
        3-dimensional reflexive polytopes. Databases are stored in the
        directory of the package.

        .. note::

           The first call to this function for each dimension can take
           a few seconds while the dictionary of all polytopes is
           constructed, but after that it is cached and fast.

        :rtype: integer

        EXAMPLES: We check what is the index of the "diamond" in the
        database::

            sage: o = lattice_polytope.octahedron(2)
            sage: o.index()
            3

        Note that polytopes with the same index are not necessarily the
        same::

            sage: o.vertices()
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: lattice_polytope.ReflexivePolytope(2,3).vertices()
            [ 1  0  0 -1]
            [ 0  1 -1  0]

        But they are in the same `GL(Z^n)` orbit and have the same
        normal form::

            sage: o.normal_form()
            [ 1  0  0 -1]
            [ 0  1 -1  0]
            sage: lattice_polytope.ReflexivePolytope(2,3).normal_form()
            [ 1  0  0 -1]
            [ 0  1 -1  0]
        """
        if not hasattr(self, "_index"):
            if not self.is_reflexive():
                raise NotImplementedError, "only reflexive polytopes can be indexed!"
            dim = self.dim()
            if dim not in [2, 3]:
                raise NotImplementedError, "only 2- and 3-dimensional polytopes can be indexed!"
            if LatticePolytopeClass._rp_dict[dim] == None:
                rp_dict = dict()
                for n, p in enumerate(ReflexivePolytopes(dim)):
                    rp_dict[str(p.normal_form())] = n
                LatticePolytopeClass._rp_dict[dim] = rp_dict
            self._index = LatticePolytopeClass._rp_dict[dim][str(self.normal_form())]
        return self._index

    def is_reflexive(self):
        r"""
        Return True if this polytope is reflexive.

        EXAMPLES: The 3-dimensional octahedron is reflexive (and 4318 other
        3-polytopes)::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.is_reflexive()
            True

        But not all polytopes are reflexive::

            sage: m = matrix(ZZ, [[1, 0, 0, -1,  0,  0],
            ...                   [0, 1, 0,  0, -1,  0],
            ...                   [0, 0, 0,  0,  0, -1]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.is_reflexive()
            False

        Only full-dimensional polytopes can be reflexive (otherwise the polar
        set is not a polytope at all, since it is unbounded)::

            sage: m = matrix(ZZ, [[0, 1, 0, -1,  0],
            ...                   [0, 0, 1,  0, -1],
            ...                   [0, 0, 0,  0,  0]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.is_reflexive()
            False
        """
        # In the last doctest above the origin is added intentionall - it
        # causes the sublattice polytope to be reflexive, yet the original
        # still should not be one.
        if not hasattr(self, "_is_reflexive"):
            if self.dim() != self.ambient_dim():
                self._is_reflexive = False
            else:
                # Determine if the polytope is reflexive by computing vertices
                # of the dual polytope and save all obtained information.
                self._read_equations(self.poly_x("e"))
        return self._is_reflexive

    def nef_partitions(self, keep_symmetric=False, keep_products=True,
        keep_projections=True, hodge_numbers=False):
        r"""
        Return 2-part nef-partitions of ``self``.

        INPUT:

        - ``keep_symmetric`` -- (default: ``False``) if ``True``, "-s" option
          will be passed to ``nef.x`` in order to keep symmetric partitions,
          i.e. partitions related by lattice automorphisms preserving ``self``;

        - ``keep_products`` -- (default: ``True``) if ``True``, "-D" option
          will be passed to ``nef.x`` in order to keep product partitions,
          with corresponding complete intersections being direct products;

        - ``keep_projections`` -- (default: ``True``) if ``True``, "-P" option
          will be passed to ``nef.x`` in order to keep projection partitions,
          i.e. partitions with one of the parts consisting of a single vertex;

        - ``hodge_numbers`` -- (default: ``False``) if ``False``, "-p" option
          will be passed to ``nef.x`` in order to skip Hodge numbers
          computation, which takes a lot of time.

        OUTPUT:

        - a sequence of :class:`nef-partitions <NefPartition>`.

        Type ``NefPartition?`` for definitions and notation.

        EXAMPLES:

        Nef-partitions of the 4-dimensional octahedron::

            sage: o = lattice_polytope.octahedron(4)
            sage: o.nef_partitions()
            [
            Nef-partition {0, 1, 4, 5} U {2, 3, 6, 7} (direct product),
            Nef-partition {0, 1, 2, 4} U {3, 5, 6, 7},
            Nef-partition {0, 1, 2, 4, 5} U {3, 6, 7},
            Nef-partition {0, 1, 2, 4, 5, 6} U {3, 7} (direct product),
            Nef-partition {0, 1, 2, 3} U {4, 5, 6, 7},
            Nef-partition {0, 1, 2, 3, 4} U {5, 6, 7},
            Nef-partition {0, 1, 2, 3, 4, 5} U {6, 7},
            Nef-partition {0, 1, 2, 3, 4, 5, 6} U {7} (projection)
            ]

        Now we omit projections::

            sage: o.nef_partitions(keep_projections=False)
            [
            Nef-partition {0, 1, 4, 5} U {2, 3, 6, 7} (direct product),
            Nef-partition {0, 1, 2, 4} U {3, 5, 6, 7},
            Nef-partition {0, 1, 2, 4, 5} U {3, 6, 7},
            Nef-partition {0, 1, 2, 4, 5, 6} U {3, 7} (direct product),
            Nef-partition {0, 1, 2, 3} U {4, 5, 6, 7},
            Nef-partition {0, 1, 2, 3, 4} U {5, 6, 7},
            Nef-partition {0, 1, 2, 3, 4, 5} U {6, 7}
            ]

        Currently Hodge numbers cannot be computed for a given nef-partition::

            sage: o.nef_partitions()[1].hodge_numbers()
            Traceback (most recent call last):
            ...
            NotImplementedError: use nef_partitions(hodge_numbers=True)!

        But they can be obtained from ``nef.x`` for all nef-partitions at once.
        Partitions will be exactly the same::

            sage: o.nef_partitions(hodge_numbers=True)  # long time (2s on sage.math, 2011)
            [
            Nef-partition {0, 1, 4, 5} U {2, 3, 6, 7} (direct product),
            Nef-partition {0, 1, 2, 4} U {3, 5, 6, 7},
            Nef-partition {0, 1, 2, 4, 5} U {3, 6, 7},
            Nef-partition {0, 1, 2, 4, 5, 6} U {3, 7} (direct product),
            Nef-partition {0, 1, 2, 3} U {4, 5, 6, 7},
            Nef-partition {0, 1, 2, 3, 4} U {5, 6, 7},
            Nef-partition {0, 1, 2, 3, 4, 5} U {6, 7},
            Nef-partition {0, 1, 2, 3, 4, 5, 6} U {7} (projection)
            ]

        Now it is possible to get Hodge numbers::

            sage: o.nef_partitions(hodge_numbers=True)[1].hodge_numbers()
            (20,)

        Since nef-partitions are cached, their Hodge numbers are accessible
        after the first request, even if you do not specify
        ``hodge_numbers=True`` anymore::

            sage: o.nef_partitions()[1].hodge_numbers()
            (20,)

        We illustrate removal of symmetric partitions on a diamond::

            sage: o = lattice_polytope.octahedron(2)
            sage: o.nef_partitions()
            [
            Nef-partition {0, 2} U {1, 3} (direct product),
            Nef-partition {0, 1} U {2, 3},
            Nef-partition {0, 1, 2} U {3} (projection)
            ]
            sage: o.nef_partitions(keep_symmetric=True)
            [
            Nef-partition {0, 1, 3} U {2} (projection),
            Nef-partition {0, 2, 3} U {1} (projection),
            Nef-partition {0, 3} U {1, 2},
            Nef-partition {1, 2, 3} U {0} (projection),
            Nef-partition {1, 3} U {0, 2} (direct product),
            Nef-partition {2, 3} U {0, 1},
            Nef-partition {0, 1, 2} U {3} (projection)
            ]

        Nef-partitions can be computed only for reflexive polytopes::

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
        if hasattr(self, "_npkeys"):
            oldkeys = self._npkeys
            if oldkeys == keys:
                return self._nef_partitions
            if not (hodge_numbers and oldkeys.find("-p") != -1
                or keep_symmetric and oldkeys.find("-s") == -1
                or not keep_symmetric and oldkeys.find("-s") != -1
                or keep_projections and oldkeys.find("-P") == -1
                or keep_products and oldkeys.find("-D") == -1):
                # Select only necessary partitions
                return Sequence([p for p in self._nef_partitions
                                 if (keep_projections or not p._is_projection)
                                     and (keep_products or not p._is_product)],
                                cr=True, check=False)
        self._read_nef_partitions(self.nef_x(keys))
        self._npkeys = keys
        return self._nef_partitions

    def nef_x(self, keys):
        r"""
        Run nef.x with given ``keys`` on vertices of this
        polytope.

        INPUT:


        -  ``keys`` - a string of options passed to nef.x. The
           key "-f" is added automatically.


        OUTPUT: the output of nef.x as a string.

        EXAMPLES: This call is used internally for computing
        nef-partitions::

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

        EXAMPLES: The number of facets of the 3-dimensional octahedron::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.nfacets()
            8

        The number of facets of an interval is 2::

            sage: LatticePolytope(matrix([1,2])).nfacets()
            2

        Now consider a 2-dimensional diamond in a 3-dimensional space::

            sage: m = matrix(ZZ, [[1, 0, -1,  0],
            ...                   [0, 1,  0, -1],
            ...                   [0, 0,  0,  0]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.nfacets()
            4
        """
        if not hasattr(self, "_nfacets"):
            if self.is_reflexive():
                self._nfacets = self.polar().nvertices()
            elif self.dim() == self.ambient_dim():
                self._nfacets = self._facet_normals.nrows()
            elif self.dim() <= 0:
                self._nfacets = 0
            else:
                self._nfacets = self._sublattice_polytope.nfacets()
        return self._nfacets

    def normal_form(self):
        r"""
        Return the normal form of vertices of the polytope.

        Two full-dimensional lattice polytopes are in the same GL(Z)-orbit if
        and only if their normal forms are the same. Normal form is not defined
        and thus cannot be used for polytopes whose dimension is smaller than
        the dimension of the ambient space.

        EXAMPLES: We compute the normal form of the "diamond"::

            sage: o = lattice_polytope.octahedron(2)
            sage: o.vertices()
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: o.normal_form()
            [ 1  0  0 -1]
            [ 0  1 -1  0]

        The diamond is the 3rd polytope in the internal database...

        ::

            sage: o.index()
            3
            sage: lattice_polytope.ReflexivePolytope(2,3).vertices()
            [ 1  0  0 -1]
            [ 0  1 -1  0]

        It is not possible to compute normal forms for polytopes which do not
        span the space::

            sage: m = matrix(ZZ, [[1, 0, -1,  0],
            ...                   [0, 1,  0, -1],
            ...                   [0, 0,  0,  0]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.normal_form()
            Traceback (most recent call last):
            ...
            ValueError: Normal form is not defined for a 2-dimensional polytope in a 3-dimensional space!
        """
        if not hasattr(self, "_normal_form"):
            if self.dim() < self.ambient_dim():
                raise ValueError(
                ("Normal form is not defined for a %d-dimensional polytope " +
                "in a %d-dimensional space!") % (self.dim(), self.ambient_dim()))
            self._normal_form = read_palp_matrix(self.poly_x("N"))
            self._normal_form.set_immutable()
        return self._normal_form

    def npoints(self):
        r"""
        Return the number of lattice points of this polytope.

        EXAMPLES: The number of lattice points of the 3-dimensional
        octahedron and its polar cube::

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

        EXAMPLES: The number of vertices of the 3-dimensional octahedron
        and its polar cube::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.nvertices()
            6
            sage: cube = o.polar()
            sage: cube.nvertices()
            8
        """
        return self._vertices.ncols()

    def origin(self):
        r"""
        Return the index of the origin in the list of points of self.

        OUTPUT:

        - integer if the origin belongs to this polytope, ``None`` otherwise.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(2)
            sage: o.origin()
            4
            sage: o.point(o.origin())
            (0, 0)

            sage: p = LatticePolytope(matrix([[1,2]]))
            sage: p.points()
            [1 2]
            sage: print p.origin()
            None

        Now we make sure that the origin of non-full-dimensional polytopes can
        be identified correctly (Trac #10661)::

            sage: LatticePolytope([(1,0,0), (-1,0,0)]).origin()
            2
        """
        if "_origin" not in self.__dict__:
            origin = vector(ZZ, self.ambient_dim())
            points = self.points().columns(copy=False)
            self._origin = points.index(origin) if origin in points else None
        return self._origin

    def parent(self):
        """
        Return the set of all lattice polytopes.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.parent()
            Set of all Lattice Polytopes
        """
        return SetOfAllLatticePolytopes

    def plot3d(self,
            show_facets=True, facet_opacity=0.5, facet_color=(0,1,0),
            facet_colors=None,
            show_edges=True, edge_thickness=3, edge_color=(0.5,0.5,0.5),
            show_vertices=True, vertex_size=10, vertex_color=(1,0,0),
            show_points=True, point_size=10, point_color=(0,0,1),
            show_vindices=None, vindex_color=(0,0,0),
            vlabels=None,
            show_pindices=None, pindex_color=(0,0,0),
            index_shift=1.1):
        r"""
        Return a 3d-plot of this polytope.

        Polytopes with ambient dimension 1 and 2 will be plotted along x-axis
        or in xy-plane respectively. Polytopes of dimension 3 and less with
        ambient dimension 4 and greater will be plotted in some basis of the
        spanned space.

        By default, everything is shown with more or less pretty
        combination of size and color parameters.

        INPUT: Most of the parameters are self-explanatory:


        -  ``show_facets`` - (default:True)

        -  ``facet_opacity`` - (default:0.5)

        -  ``facet_color`` - (default:(0,1,0))

        -  ``facet_colors`` - (default:None) if specified, must be a list of
           colors for each facet separately, used instead of ``facet_color``

        -  ``show_edges`` - (default:True) whether to draw
           edges as lines

        -  ``edge_thickness`` - (default:3)

        -  ``edge_color`` - (default:(0.5,0.5,0.5))

        -  ``show_vertices`` - (default:True) whether to draw
           vertices as balls

        -  ``vertex_size`` - (default:10)

        -  ``vertex_color`` - (default:(1,0,0))

        -  ``show_points`` - (default:True) whether to draw
           other poits as balls

        -  ``point_size`` - (default:10)

        -  ``point_color`` - (default:(0,0,1))

        -  ``show_vindices`` - (default:same as
           show_vertices) whether to show indices of vertices

        -  ``vindex_color`` - (default:(0,0,0)) color for
           vertex labels

        -  ``vlabels`` - (default:None) if specified, must be a list of labels
           for each vertex, default labels are vertex indicies

        -  ``show_pindices`` - (default:same as show_points)
           whether to show indices of other points

        -  ``pindex_color`` - (default:(0,0,0)) color for
           point labels

        -  ``index_shift`` - (default:1.1)) if 1, labels are
           placed exactly at the corresponding points. Otherwise the label
           position is computed as a multiple of the point position vector.


        EXAMPLES: The default plot of a cube::

            sage: c = lattice_polytope.octahedron(3).polar()
            sage: c.plot3d()

        Plot without facets and points, shown without the frame::

            sage: c.plot3d(show_facets=false,show_points=false).show(frame=False)

        Plot with facets of different colors::

            sage: c.plot3d(facet_colors=rainbow(c.nfacets(), 'rgbtuple'))

        It is also possible to plot lower dimensional polytops in 3D (let's
        also change labels of vertices)::

            sage: lattice_polytope.octahedron(2).plot3d(vlabels=["A", "B", "C", "D"])

        TESTS::

            sage: m = matrix([[0,0,0],[0,1,1],[1,0,1],[1,1,0]]).transpose()
            sage: p = LatticePolytope(m, compute_vertices=True)
            sage: p.plot3d()
        """
        dim = self.dim()
        amb_dim = self.ambient_dim()
        if dim > 3:
            raise ValueError, "%d-dimensional polytopes can not be plotted in 3D!" % self.dim()
        elif amb_dim > 3:
            return self._sublattice_polytope.plot3d(
                show_facets, facet_opacity, facet_color,
                facet_colors,
                show_edges, edge_thickness, edge_color,
                show_vertices, vertex_size, vertex_color,
                show_points, point_size, point_color,
                show_vindices, vindex_color,
                vlabels,
                show_pindices, pindex_color,
                index_shift)
        elif dim == 3:
            vertices = self.vertices().columns(copy=False)
            if show_points or show_pindices:
                points = self.points().columns(copy=False)[self.nvertices():]
        else:
            vertices = [vector(ZZ, list(self.vertex(i))+[0]*(3-amb_dim))
                        for i in range(self.nvertices())]
            if show_points or show_pindices:
                points = [vector(ZZ, list(self.point(i))+[0]*(3-amb_dim))
                        for i in range(self.nvertices(), self.npoints())]
        pplot = 0
        if show_facets:
            if dim == 2:
                pplot += IndexFaceSet([self.traverse_boundary()],
                        vertices, opacity=facet_opacity, rgbcolor=facet_color)
            elif dim == 3:
                if facet_colors != None:
                    for i, f in enumerate(self.facets()):
                        pplot += IndexFaceSet([f.traverse_boundary()],
                            vertices, opacity=facet_opacity, rgbcolor=facet_colors[i])
                else:
                    pplot += IndexFaceSet([f.traverse_boundary() for f in self.facets()],
                        vertices, opacity=facet_opacity, rgbcolor=facet_color)
        if show_edges:
            if dim == 1:
                pplot += line3d(vertices, thickness=edge_thickness, rgbcolor=edge_color)
            else:
                for e in self.edges():
                    pplot += line3d([vertices[e.vertices()[0]], vertices[e.vertices()[1]]],
                            thickness=edge_thickness, rgbcolor=edge_color)
        if show_vertices:
            pplot += point3d(vertices, size=vertex_size, rgbcolor=vertex_color)
        if show_vindices == None:
            show_vindices = show_vertices
        if show_pindices == None:
            show_pindices = show_points
        if show_vindices or show_pindices:
            # Compute the barycenter and shift text of labels away from it
            bc = 1/Integer(len(vertices)) * sum(vertices)
        if show_vindices:
            if vlabels == None:
                vlabels = range(len(vertices))
            for i,v in enumerate(vertices):
                pplot += text3d(vlabels[i], bc+index_shift*(v-bc), rgbcolor=vindex_color)
        if show_points and len(points) > 0:
            pplot += point3d(points, size=point_size, rgbcolor=point_color)
        if show_pindices:
            for i, p in enumerate(points):
                pplot += text3d(i+self.nvertices(), bc+index_shift*(p-bc), rgbcolor=pindex_color)
        return pplot

    def show3d(self):
        """
        Show a 3d picture of the polytope with default settings and without axes or frame.

        See self.plot3d? for more details.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.show3d()
        """
        self.plot3d().show(axis=False, frame=False)

    def point(self, i):
        r"""
        Return the i-th point of this polytope, i.e. the i-th column of the
        matrix returned by points().

        EXAMPLES: First few points are actually vertices::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.vertices()
            [ 1  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: o.point(1)
            (0, 1, 0)

        The only other point in the octahedron is the origin::

            sage: o.point(6)
            (0, 0, 0)
            sage: o.points()
            [ 1  0  0 -1  0  0  0]
            [ 0  1  0  0 -1  0  0]
            [ 0  0  1  0  0 -1  0]
        """
        # Extra checks are made to compensate for a bug in Sage - column accepts any number.
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

        EXAMPLES: The lattice points of the 3-dimensional octahedron and
        its polar cube::

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

        Lattice points of a 2-dimensional diamond in a 3-dimensional space::

            sage: m = matrix(ZZ, [[1, 0, -1,  0],
            ...                   [0, 1,  0, -1],
            ...                   [0, 0,  0,  0]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.points()
            [ 1  0 -1  0  0]
            [ 0  1  0 -1  0]
            [ 0  0  0  0  0]

        We check that points of a zero-dimensional polytope can be computed::

            sage: p = LatticePolytope(matrix([[1]]))
            sage: p.vertices()
            [1]
            sage: p.points()
            [1]
        """
        if not hasattr(self, "_points"):
            if self.dim() <= 0:
                self._points = self._vertices
            else:
                self._points = self._embed(read_palp_matrix(
                                self.poly_x("p", reduce_dimension=True)))
                self._points.set_immutable()
        return self._points

    def polar(self):
        r"""
        Return the polar polytope, if this polytope is reflexive.

        EXAMPLES: The polar polytope to the 3-dimensional octahedron::

            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: cube
            A polytope polar to An octahedron: 3-dimensional, 8 vertices.

        The polar polytope "remembers" the original one::

            sage: cube.polar()
            An octahedron: 3-dimensional, 6 vertices.
            sage: cube.polar().polar() is cube
            True

        Only reflexive polytopes have polars::

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

    def poly_x(self, keys, reduce_dimension=False):
        r"""
        Run poly.x with given ``keys`` on vertices of this
        polytope.

        INPUT:


        -  ``keys`` - a string of options passed to poly.x. The
           key "f" is added automatically.

        -  ``reduce_dimension`` - (default: False) if ``True`` and this
           polytope is not full-dimensional, poly.x will be called for the
           vertices of this polytope in some basis of the spanned affine space.


        OUTPUT: the output of poly.x as a string.

        EXAMPLES: This call is used for determining if a polytope is
        reflexive or not::

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
        compilation, the following code is likely to fail, unless you
        change default settings of PALP::

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
            Please increase POLY_Dmax to at least 7

        You cannot call poly.x for polytopes that don't span the space (if you
        could, it would crush anyway)::

            sage: m = matrix(ZZ, [[1, 0, -1,  0],
            ...                   [0, 1,  0, -1],
            ...                   [0, 0,  0,  0]])
            ...
            sage: p = LatticePolytope(m)
            sage: p.poly_x("e")
            Traceback (most recent call last):
            ...
            ValueError: Cannot run PALP for a 2-dimensional polytope in a 3-dimensional space!

        But if you know what you are doing, you can call it for the polytope in
        some basis of the spanned space::

            sage: print p.poly_x("e", reduce_dimension=True)
            4 2  Equations of P
              -1   1     0
               1   1     2
              -1  -1     0
               1  -1     2
        """
        return self._palp("poly.x -f" + keys, reduce_dimension)

    def skeleton(self):
        r"""
        Return the graph of the one-skeleton of this polytope.

        EXAMPLES: We construct the one-skeleton graph for the "diamond"::

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
        Return the increasing list of indices of lattice points in
        k-skeleton of the polytope (k is 1 by default).

        EXAMPLES: We compute all skeleton points for the cube::

            sage: o = lattice_polytope.octahedron(3)
            sage: c = o.polar()
            sage: c.skeleton_points()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 19, 21, 22, 23, 25, 26]

        The default was 1-skeleton::

            sage: c.skeleton_points(k=1)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 19, 21, 22, 23, 25, 26]

        0-skeleton just lists all vertices::

            sage: c.skeleton_points(k=0)
            [0, 1, 2, 3, 4, 5, 6, 7]

        2-skeleton lists all points except for the origin (point #17)::

            sage: c.skeleton_points(k=2)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26]

        3-skeleton includes all points::

            sage: c.skeleton_points(k=3)
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]

        It is OK to compute higher dimensional skeletons - you will get the
        list of all points::

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
        Works only for polytopes in a 3-dimensional space.

        INPUT:


        -  ``normal`` - a 3-dimensional vector (can be given as
           a list), which should be perpendicular to the screen. If not given,
           will be selected randomly (new each time and it may be far from
           "nice").


        EXAMPLES: Show a pretty picture of the octahedron::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.skeleton_show([1,2,4])

        Does not work for a diamond at the moment::

            sage: o = lattice_polytope.octahedron(2)
            sage: o.skeleton_show()
            Traceback (most recent call last):
            ...
            ValueError: need a polytope in a 3-dimensional space! Got:
            An octahedron: 2-dimensional, 4 vertices.
        """
        if self.ambient_dim() != 3:
            raise ValueError("need a polytope in a 3-dimensional space! Got:\n%s" % self)
        if normal == None:
            normal = [ZZ.random_element(20),ZZ.random_element(20),ZZ.random_element(20)]
        normal = matrix(QQ,3,1,list(normal))
        projectionm = normal.kernel().basis_matrix()
        positions = dict(enumerate([list(c) for c in (projectionm*self.points()).columns(copy=False)]))
        self.skeleton().show(pos=positions)

    def traverse_boundary(self):
        r"""
        Return a list of indices of vertices of a 2-dimensional polytope in their boundary order.

        Needed for plot3d function of polytopes.

        EXAMPLES:

            sage: c = lattice_polytope.octahedron(2).polar()
            sage: c.traverse_boundary()
            [0, 1, 3, 2]
        """
        if self.dim() != 2:
            raise ValueError, "Boundary can be traversed only for 2-polytopes!"
        edges = self.edges()
        l = [0]
        for e in edges:
            if 0 in e.vertices():
                next = e.vertices()[0] if e.vertices()[0] != 0 else e.vertices()[1]
        l.append(next)
        prev = 0
        while len(l) < self.nvertices():
            for e in edges:
                if next in e.vertices() and prev not in e.vertices():
                    prev = next
                    next = e.vertices()[0] if e.vertices()[0] != next else e.vertices()[1]
                    l.append(next)
                    break
        return l

    def vertex(self, i):
        r"""
        Return the i-th vertex of this polytope, i.e. the i-th column of
        the matrix returned by vertices().

        EXAMPLES: Note that numeration starts with zero::

            sage: o = lattice_polytope.octahedron(3)
            sage: o.vertices()
            [ 1  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: o.vertex(3)
            (-1, 0, 0)
        """
        # The check is added to compensate for a bug in Sage - column works for any numbers
        if i < 0:
            raise ValueError, "polytopes don't have negative vertices!"
        elif i > self.nvertices():
            raise ValueError, "the polytope does not have %d vertices!" % (i+1)
        else:
            return self._vertices.column(i, from_list=True)

    def vertices(self):
        r"""
        Return vertices of this polytope as columns of a matrix.

        EXAMPLES: The lattice points of the 3-dimensional octahedron and
        its polar cube::

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


def is_NefPartition(x):
    r"""
    Check if ``x`` is a nef-partition.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a :class:`nef-partition <NefPartition>` and
      ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.lattice_polytope import is_NefPartition
        sage: is_NefPartition(1)
        False
        sage: o = lattice_polytope.octahedron(3)
        sage: np = o.nef_partitions()[0]
        sage: np
        Nef-partition {0, 1, 3} U {2, 4, 5}
        sage: is_NefPartition(np)
        True
    """
    return isinstance(x, NefPartition)


class NefPartition(SageObject,
                   collections.Hashable):
    r"""
    Create a nef-partition.

    INPUT:

    - ``data`` -- a list of integers, the $i$-th element of this list must be
      the part of the $i$-th vertex of ``Delta_polar`` in this nef-partition;

    - ``Delta_polar`` -- a :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`;

    - ``check`` -- by default the input will be checked for correctness, i.e.
      that ``data`` indeed specify a nef-partition. If you are sure that the
      input is correct, you can speed up construction via ``check=False``
      option.

    OUTPUT:

    - a nef-partition of ``Delta_polar``.

    Let $M$ and $N$ be dual lattices. Let $\Delta \subset M_\RR$ be a reflexive
    polytope with polar $\Delta^\circ \subset N_\RR$. Let $X_\Delta$ be the
    toric variety associated to the normal fan of $\Delta$. A **nef-partition**
    is a decomposition of the vertex set $V$ of $\Delta^\circ$ into a disjoint
    union $V = V_0 \sqcup V_1 \sqcup \dots \sqcup V_{k-1}$ such that divisors
    $E_i = \sum_{v\in V_i} D_v$ are Cartier (here $D_v$ are prime
    torus-invariant Weil divisors corresponding to vertices of $\Delta^\circ$).
    Equivalently, let $\nabla_i \subset N_\RR$ be the convex hull of vertices
    from $V_i$ and the origin. These polytopes form a nef-partition if their
    Minkowski sum $\nabla \subset N_\RR$ is a reflexive polytope.

    The **dual nef-partition** is formed by polytopes $\Delta_i \subset M_\RR$
    of $E_i$, which give a decomposition of the vertex set of $\nabla^\circ
    \subset M_\RR$ and their Minkowski sum is $\Delta$, i.e. the polar duality
    of reflexive polytopes switches convex hull and Minkowski sum for dual
    nef-partitions:

    .. MATH::

        \Delta^\circ
        &=
        \mathrm{Conv} \left(\nabla_0, \nabla_1, \dots, \nabla_{k-1}\right), \\
        \nabla^{\phantom{\circ}}
        &=
        \nabla_0 + \nabla_1 + \dots + \nabla_{k-1}, \\
        &
        \\
        \Delta^{\phantom{\circ}}
        &=
        \Delta_0 + \Delta_1 + \dots + \Delta_{k-1}, \\
        \nabla^\circ
        &=
        \mathrm{Conv} \left(\Delta_0, \Delta_1, \dots, \Delta_{k-1}\right).

    See Section 4.3.1 in [CK99]_ and references therein for further details, or
    [BN08]_ for a purely combinatorial approach.

    REFERENCES:

    ..  [BN08]
        Victor V. Batyrev and Benjamin Nill.
        Combinatorial aspects of mirror symmetry.
        In *Integer points in polyhedra --- geometry, number theory,
        representation theory, algebra, optimization, statistics*,
        volume 452 of *Contemp. Math.*, pages 35--66.
        Amer. Math. Soc., Providence, RI, 2008.
        arXiv:math/0703456v2 [math.CO].

    ..  [CK99]
        David A. Cox and Sheldon Katz.
        *Mirror symmetry and algebraic geometry*,
        volume 68 of *Mathematical Surveys and Monographs*.
        American Mathematical Society, Providence, RI, 1999.

    EXAMPLES:

    It is very easy to create a nef-partition for the octahedron, since for
    this polytope any decomposition of vertices is a nef-partition. We create a
    3-part nef-partition with the 0-th and 1-st vertices belonging to the 0-th
    part (recall that numeration in Sage starts with 0), the 2-nd and 5-th
    vertices belonging to the 1-st part, and 3-rd and 4-th vertices belonging
    to the 2-nd part::

        sage: o = lattice_polytope.octahedron(3)
        sage: np = NefPartition([0,0,1,2,2,1], o)
        sage: np
        Nef-partition {0, 1} U {2, 5} U {3, 4}

    The octahedron plays the role of `\Delta^\circ` in the above description::

        sage: np.Delta_polar() is o
        True

    The dual nef-partition (corresponding to the "mirror complete
    intersection") gives decomposition of the vertex set of `\nabla^\circ`::

        sage: np.dual()
        Nef-partition {4, 5, 6} U {1, 3} U {0, 2, 7}
        sage: np.nabla_polar().vertices()
        [ 1  0  0  0 -1  0 -1  1]
        [ 1  0  1  0 -1 -1  0  0]
        [ 0  1  0 -1  0  0  0  0]

    Of course, `\nabla^\circ` is `\Delta^\circ` from the point of view of the
    dual nef-partition::

        sage: np.dual().Delta_polar() is np.nabla_polar()
        True
        sage: np.Delta(1).vertices()
        [ 0  0]
        [ 0  0]
        [ 1 -1]
        sage: np.dual().nabla(1).vertices()
        [ 0  0]
        [ 0  0]
        [ 1 -1]

    Instead of constructing nef-partitions directly, you can request all 2-part
    nef-partitions of a given reflexive polytope (they will be computed using
    ``nef.x`` program from PALP)::

        sage: o.nef_partitions()
        [
        Nef-partition {0, 1, 3} U {2, 4, 5},
        Nef-partition {0, 1, 3, 4} U {2, 5} (direct product),
        Nef-partition {0, 1, 2} U {3, 4, 5},
        Nef-partition {0, 1, 2, 3} U {4, 5},
        Nef-partition {0, 1, 2, 3, 4} U {5} (projection)
        ]
    """

    def __init__(self, data, Delta_polar, check=True):
        r"""
        See :class:`NefPartition` for documentation.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: TestSuite(np).run()
        """
        if check and not Delta_polar.is_reflexive():
            raise ValueError("nef-partitions can be constructed for reflexive "
                             "polytopes ony!")
        self._vertex_to_part = tuple(int(el) for el in data)
        self._nparts = max(self._vertex_to_part) + 1
        self._Delta_polar = Delta_polar
        if check and not self.nabla().is_reflexive():
            raise ValueError("%s do not form a nef-partition!" % str(data))

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        INPUT:

        - ``other`` -- anything.

        OUTPUT:

        - ``True`` if ``other`` is a :class:`nef-partition <NefPartition>`
          equal to ``self``, ``False`` otherwise.

        .. NOTE::

            Two nef-partitions are equal if they correspond to equal polytopes
            and their parts are the same, including their order.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np == np
            True
            sage: np == o.nef_partitions()[1]
            False
            sage: np2 = NefPartition(np._vertex_to_part, o)
            sage: np2 is np
            False
            sage: np2 == np
            True
            sage: np == 0
            False
        """
        return (is_NefPartition(other)
                and self._Delta_polar == other._Delta_polar
                and self._vertex_to_part == other._vertex_to_part)

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        OUTPUT:

        - an integer.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: hash(np) == hash(np)
            True
        """
        try:
            return self._hash
        except AttributeError:
            self._hash = hash(self._vertex_to_part) + hash(self._Delta_polar)
            return self._hash

    def __ne__(self, other):
        r"""
        Compare ``self`` with ``other``.

        INPUT:

        - ``other`` -- anything.

        OUTPUT:

        - ``False`` if ``other`` is a :class:`nef-partition <NefPartition>`
          equal to ``self``, ``True`` otherwise.

        .. NOTE::

            Two nef-partitions are equal if they correspond to equal polytopes
            and their parts are the same, including their order.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np != np
            False
            sage: np != o.nef_partitions()[1]
            True
            sage: np2 = NefPartition(np._vertex_to_part, o)
            sage: np2 is np
            False
            sage: np2 != np
            False
            sage: np != 0
            True
        """
        return not (self == other)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: latex(np) # indirect doctest
            \text{Nef-partition } \{0, 1, 3\} \sqcup \{2, 4, 5\}
        """
        result = r"\text{Nef-partition } "
        for i, part in enumerate(self.parts()):
            if i != 0:
                result += " \sqcup "
            result += r"\{" + ", ".join("%d" % v for v in part) + r"\}"
        try:
            # We may or may not know the type of the partition
            if self._is_product:
                result += r" \text{ (direct product)}"
            if self._is_projection:
                result += r" \text{ (projection)}"
        except AttributeError:
            pass
        return result

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: repr(np)  # indirect doctest
            'Nef-partition {0, 1, 3} U {2, 4, 5}'
        """
        result = "Nef-partition "
        for i, part in enumerate(self.parts()):
            if i != 0:
                result += " U "
            result += "{" + ", ".join("%d" % v for v in part) + "}"
        try:
            # We may or may not know the type of the partition
            if self._is_product:
                result += " (direct product)"
            if self._is_projection:
                result += " (projection)"
        except AttributeError:
            pass
        return result

    def Delta(self, i=None):
        r"""
        Return the polytope $\Delta$ or $\Delta_i$ corresponding to ``self``.

        INPUT:

        - ``i`` -- an integer. If not given, $\Delta$ will be returned.

        OUTPUT:

        - a :class:`lattice polytope <LatticePolytopeClass>`.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.Delta().polar() is o
            True
            sage: np.Delta().vertices()
            [-1  1 -1  1 -1  1 -1  1]
            [-1 -1  1  1 -1 -1  1  1]
            [ 1  1  1  1 -1 -1 -1 -1]
            sage: np.Delta(0).vertices()
            [ 1  1 -1 -1]
            [-1  0 -1  0]
            [ 0  0  0  0]
        """
        if i is None:
            return self._Delta_polar.polar()
        else:
            return self.dual().nabla(i)

    def Delta_polar(self):
        r"""
        Return the polytope $\Delta^\circ$ corresponding to ``self``.

        OUTPUT:

        - a :class:`lattice polytope <LatticePolytopeClass>`.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLE::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.Delta_polar() is o
            True
        """
        return self._Delta_polar

    def Deltas(self):
        r"""
        Return the polytopes $\Delta_i$ corresponding to ``self``.

        OUTPUT:

        - a tuple of :class:`lattice polytopes <LatticePolytopeClass>`.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.Delta().vertices()
            [-1  1 -1  1 -1  1 -1  1]
            [-1 -1  1  1 -1 -1  1  1]
            [ 1  1  1  1 -1 -1 -1 -1]
            sage: [Delta_i.vertices() for Delta_i in np.Deltas()]
            [
            [ 1  1 -1 -1]  [ 0  0  0  0]
            [-1  0 -1  0]  [ 1  0  0  1]
            [ 0  0  0  0], [ 1  1 -1 -1]
            ]
            sage: np.nabla_polar().vertices()
            [ 1  0  1  0  0 -1  0 -1]
            [-1  1  0  0  0 -1  1  0]
            [ 0  1  0  1 -1  0 -1  0]
        """
        return self.dual().nablas()

    def dual(self):
        r"""
        Return the dual nef-partition.

        OUTPUT:

        - a :class:`nef-partition <NefPartition>`.

        See the class documentation for the definition.

        ALGORITHM:

        See Proposition 3.19 in [BN08]_.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.dual()
            Nef-partition {0, 2, 5, 7} U {1, 3, 4, 6}
            sage: np.dual().Delta() is np.nabla()
            True
            sage: np.dual().nabla(0) is np.Delta(0)
            True
        """
        try:
            return self._dual
        except AttributeError:
            # Delta and nabla are interchanged compared to [BN08]_.
            nabla_polar = self.nabla_polar()
            n = nabla_polar.nvertices()
            vertex_to_part = [-1] * n
            for i in range(self._nparts):
                A = nabla_polar.vertices().transpose()*self.nabla(i).vertices()
                for j in range(n):
                    if min(A[j]) == -1:
                        vertex_to_part[j] = i
            self._dual = NefPartition(vertex_to_part, nabla_polar)
            self._dual._dual = self
            self._dual._nabla = self.Delta() # For vertex order consistency
            return self._dual

    def hodge_numbers(self):
        r"""
        Return Hodge numbers corresponding to ``self``.

        OUTPUT:

        - a tuple of integers (produced by ``nef.x`` program from PALP).

        EXAMPLES:

        Currently, you need to request Hodge numbers when you compute
        nef-partitions::

            sage: o = lattice_polytope.octahedron(5)
            sage: np = o.nef_partitions()[0]  # long time (4s on sage.math, 2011)
            sage: np.hodge_numbers()  # long time
            Traceback (most recent call last):
            ...
            NotImplementedError: use nef_partitions(hodge_numbers=True)!
            sage: np = o.nef_partitions(hodge_numbers=True)[0]  # long time (13s on sage.math, 2011)
            sage: np.hodge_numbers()  # long time
            (19, 19)
        """
        try:
            return self._hodge_numbers
        except AttributeError:
            self._Delta_polar._compute_hodge_numbers()
            return self._hodge_numbers

    def nabla(self, i=None):
        r"""
        Return the polytope $\nabla$ or $\nabla_i$ corresponding to ``self``.

        INPUT:

        - ``i`` -- an integer. If not given, $\nabla$ will be returned.

        OUTPUT:

        - a :class:`lattice polytope <LatticePolytopeClass>`.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.Delta_polar().vertices()
            [ 1  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: np.nabla(0).vertices()
            [ 1  0 -1]
            [ 0  1  0]
            [ 0  0  0]
            sage: np.nabla().vertices()
            [ 1  1  1  0  0 -1 -1 -1]
            [ 0 -1  0  1  1  0 -1  0]
            [ 1  0 -1  1 -1  1  0 -1]
        """
        if i is None:
            try:
                return self._nabla
            except AttributeError:
                self._nabla = LatticePolytope(reduce(minkowski_sum,
                    (nabla.vertices().columns() for nabla in self.nablas())))
                return self._nabla
        else:
            return self.nablas()[i]

    def nabla_polar(self):
        r"""
        Return the polytope $\nabla^\circ$ corresponding to ``self``.

        OUTPUT:

        - a :class:`lattice polytope <LatticePolytopeClass>`.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.nabla_polar().vertices()
            [ 1  0  1  0  0 -1  0 -1]
            [-1  1  0  0  0 -1  1  0]
            [ 0  1  0  1 -1  0 -1  0]
            sage: np.nabla_polar() is np.dual().Delta_polar()
            True
        """
        return self.nabla().polar()

    def nablas(self):
        r"""
        Return the polytopes $\nabla_i$ corresponding to ``self``.

        OUTPUT:

        - a tuple of :class:`lattice polytopes <LatticePolytopeClass>`.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.Delta_polar().vertices()
            [ 1  0  0 -1  0  0]
            [ 0  1  0  0 -1  0]
            [ 0  0  1  0  0 -1]
            sage: [nabla_i.vertices() for nabla_i in np.nablas()]
            [
            [ 1  0 -1]  [ 0  0  0]
            [ 0  1  0]  [ 0 -1  0]
            [ 0  0  0], [ 1  0 -1]
            ]
        """
        try:
            return self._nablas
        except AttributeError:
            Delta_polar = self._Delta_polar
            origin = [[0] * Delta_polar.dim()]
            self._nablas = tuple(LatticePolytope(
                                [Delta_polar.vertex(j) for j in part] + origin)
                                for part in self.parts())
            return self._nablas

    def nparts(self):
        r"""
        Return the number of parts in ``self``.

        OUTPUT:

        - an integer.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.nparts()
            2
        """
        return self._nparts

    def part(self, i):
        r"""
        Return the ``i``-th part of ``self``.

        INPUT:

        - ``i`` -- an integer.

        OUTPUT:

        - a tuple of integers, indices of vertices of $\Delta^\circ$ belonging
          to $V_i$.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.part(0)
            (0, 1, 3)
        """
        return self.parts()[i]

    def parts(self):
        r"""
        Return all parts of ``self``.

        OUTPUT:

        - a tuple of tuples of integers. The $i$-th tuple contains indices of
          vertices of $\Delta^\circ$ belonging to $V_i$.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.parts()
            ((0, 1, 3), (2, 4, 5))
        """
        try:
            return self._parts
        except AttributeError:
            parts = []
            for part in range(self._nparts):
                parts.append([])
            for vertex, part in enumerate(self._vertex_to_part):
                parts[part].append(vertex)
            self._parts = tuple(tuple(part) for part in parts)
            return self._parts

    def part_of(self, i):
        r"""
        Return the index of the part containing the ``i``-th vertex.

        INPUT:

        - ``i`` -- an integer.

        OUTPUT:

        - an integer $j$ such that the ``i``-th vertex of $\Delta^\circ$
          belongs to $V_j$.

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES::

            sage: o = lattice_polytope.octahedron(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} U {2, 4, 5}
            sage: np.part_of(3)
            0
            sage: np.part_of(2)
            1
        """
        return self._vertex_to_part[i]

    def part_of_point(self, i):
        r"""
        Return the index of the part containing the ``i``-th point.

        INPUT:

        - ``i`` -- an integer.

        OUTPUT:

        - an integer `j` such that the ``i``-th point of `\Delta^\circ`
          belongs to `\nabla_j`.

        .. NOTE::

            Since a nef-partition induces a partition on the set of boundary
            lattice points of `\Delta^\circ`, the value of `j` is well-defined
            for all `i` but the one that corresponds to the origin, in which
            case this method will raise a ``ValueError`` exception. (The origin
            always belongs to all `\nabla_j`.)

        See :class:`nef-partition <NefPartition>` class documentation for
        definitions and notation.

        EXAMPLES:

        We consider a relatively complicated reflexive polytope #2252 (easily
        accessible in Sage as ``ReflexivePolytope(3, 2252)``, we create it here
        explicitly to avoid loading the whole database)::

            sage: p = LatticePolytope([(1,0,0), (0,1,0), (0,0,1), (0,1,-1),
            ...           (0,-1,1), (-1,1,0), (0,-1,-1), (-1,-1,0), (-1,-1,2)])
            sage: np = p.nef_partitions()[0]
            sage: np
            Nef-partition {1, 2, 5, 7, 8} U {0, 3, 4, 6}
            sage: p.nvertices()
            9
            sage: p.npoints()
            15

        We see that the polytope has 6 more points in addition to vertices. One
        of them is the origin::

            sage: p.origin()
            14
            sage: np.part_of_point(14)
            Traceback (most recent call last):
            ...
            ValueError: the origin belongs to all parts!

        But the remaining 5 are partitioned by ``np``::

            sage: [n for n in range(p.npoints())
            ...      if p.origin() != n and np.part_of_point(n) == 0]
            [1, 2, 5, 7, 8, 9, 11, 13]
            sage: [n for n in range(p.npoints())
            ...      if p.origin() != n and np.part_of_point(n) == 1]
            [0, 3, 4, 6, 10, 12]
        """
        try:
            ptp = self._point_to_part
        except AttributeError:
            ptp = [-1] * self._Delta_polar.npoints()
            for v, part in enumerate(self._vertex_to_part):
                ptp[v] = part
            self._point_to_part = ptp
        if ptp[i] > 0:
            return ptp[i]
        if i == self._Delta_polar.origin():
            raise ValueError("the origin belongs to all parts!")
        point = self._Delta_polar.point(i)
        for part, nabla in enumerate(self.nablas()):
            if min(nabla.distances(point)) >= 0:
                ptp[i] = part
                break
        return ptp[i]


class _PolytopeFace(SageObject):
    r"""
    _PolytopeFace(polytope, vertices, facets)

    Construct a polytope face.

    POLYTOPE FACES SHOULD NOT BE CONSTRUCTED OUTSIDE OF LATTICE
    POLYTOPES!

    INPUT:


    -  ``polytope`` - a polytope whose face is being
       constructed.

    -  ``vertices`` - a sequence of indices of generating
       vertices.

    -  ``facets`` - a sequence of indices of facets
       containing this face.
    """
    def __init__(self, polytope, vertices, facets):
        r"""
        Construct a face.

        TESTS::

            sage: o = lattice_polytope.octahedron(2)
            sage: o.faces()
            [
            [[0], [1], [2], [3]],
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

        TESTS::

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

        TESTS::

            sage: o = lattice_polytope.octahedron(3)
            sage: f = o.facets()[0]
            sage: f._repr_()
            '[0, 1, 5]'
        """
        return repr(self._vertices)

    def boundary_points(self):
        r"""
        Return a sequence of indices of boundary lattice points of this
        face.

        EXAMPLES: Boundary lattice points of one of the facets of the
        3-dimensional cube::

            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.boundary_points()
            [0, 2, 4, 6, 8, 9, 11, 12]
        """
        try:
            return self._boundary_points
        except AttributeError:
            self._polytope._face_split_points(self)
            return self._boundary_points

    def facets(self):
        r"""
        Return a sequence of indices of facets containing this face.

        EXAMPLES: Facets containing one of the edges of the 3-dimensional
        octahedron::

            sage: o = lattice_polytope.octahedron(3)
            sage: edge = o.faces(dim=1)[0]
            sage: edge.facets()
            [0, 1]

        Thus ``edge`` is the intersection of facets 0 and 1::

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
        Return a sequence of indices of interior lattice points of this
        face.

        EXAMPLES: Interior lattice points of one of the facets of the
        3-dimensional cube::

            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.interior_points()
            [10]
        """
        try:
            return self._interior_points
        except AttributeError:
            self._polytope._face_split_points(self)
            return self._interior_points

    def nboundary_points(self):
        r"""
        Return the number of boundary lattice points of this face.

        EXAMPLES: The number of boundary lattice points of one of the
        facets of the 3-dimensional cube::

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

        EXAMPLES: The number of facets containing one of the edges of the
        3-dimensional octahedron::

            sage: o = lattice_polytope.octahedron(3)
            sage: edge = o.faces(dim=1)[0]
            sage: edge.nfacets()
            2
        """
        return len(self._facets)

    def ninterior_points(self):
        r"""
        Return the number of interior lattice points of this face.

        EXAMPLES: The number of interior lattice points of one of the
        facets of the 3-dimensional cube::

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

        EXAMPLES: The number of lattice points of one of the facets of the
        3-dimensional cube::

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

        EXAMPLES: The number of vertices generating one of the facets of
        the 3-dimensional cube::

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

        Works only for edges, i.e. faces generated by exactly two
        vertices.

        EXAMPLE: We find all points along an edge of the cube::

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
        primitive = primitive * (1/integral_length(primitive.list()))
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

        EXAMPLES: The lattice points of one of the facets of the
        3-dimensional cube::

            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.points()
            [0, 2, 4, 6, 8, 9, 10, 11, 12]
        """
        try:
            return self._points
        except AttributeError:
            self._polytope._face_compute_points(self)
            return self._points

    def traverse_boundary(self):
        r"""
        Return a list of indices of vertices of a 2-face in their boundary
        order.

        Needed for plot3d function of polytopes.

        EXAMPLES::

            sage: c = lattice_polytope.octahedron(3).polar()
            sage: f = c.facets()[0]
            sage: f.vertices()
            [0, 2, 4, 6]
            sage: f.traverse_boundary()
            [0, 4, 6, 2]
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

        EXAMPLES: The vertices generating one of the facets of the
        3-dimensional cube::

            sage: o = lattice_polytope.octahedron(3)
            sage: cube = o.polar()
            sage: face = cube.facets()[0]
            sage: face.vertices()
            [0, 2, 4, 6]
        """
        return self._vertices


def _create_octahedron(dim):
    r"""
    Create an octahedron of the given dimension.

    Since we know that we are passing vertices, we suppress their
    computation.

    TESTS::

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

_palp_dimension = None

def _palp(command, polytopes, reduce_dimension=False):
    r"""
    Run ``command`` on vertices of given
    ``polytopes``.

    Returns the name of the file containing the output of
    ``command``. You should delete it after using.

    .. note::

      PALP cannot be called for polytopes that do not span the ambient space.
      If you specify ``reduce_dimension=True`` argument, PALP will be
      called for vertices of this polytope in some basis of the affine space
      it spans.

    TESTS::

        sage: o = lattice_polytope.octahedron(3)
        sage: result_name = lattice_polytope._palp("poly.x -f", [o])
        sage: f = open(result_name)
        sage: f.readlines()
        ['M:7 6 N:27 8 Pic:17 Cor:0\n']
        sage: f.close()
        sage: os.remove(result_name)

        sage: m = matrix(ZZ, [[1, 0, -1,  0],
        ...                   [0, 1,  0, -1],
        ...                   [0, 0,  0,  0]])
        ...
        sage: p = LatticePolytope(m)
        sage: lattice_polytope._palp("poly.x -f", [p])
        Traceback (most recent call last):
        ValueError: Cannot run PALP for a 2-dimensional polytope in a 3-dimensional space!

        sage: result_name = lattice_polytope._palp("poly.x -f", [p], reduce_dimension=True)
        sage: f = open(result_name)
        sage: f.readlines()
        ['M:5 4 F:4\n']
        sage: f.close()
        sage: os.remove(result_name)
    """
    if _palp_dimension is not None:
        dot = command.find(".")
        command = command[:dot] + "-%dd" % _palp_dimension + command[dot:]
    input_file_name = tmp_filename()
    input_file = open(input_file_name, "w")
    for p in polytopes:
        if p.dim() == 0:
            raise ValueError(("Cannot run \"%s\" for the zero-dimensional "
                + "polytope!\nPolytope: %s") % (command, p))
        if p.dim() < p.ambient_dim() and not reduce_dimension:
            raise ValueError(("Cannot run PALP for a %d-dimensional polytope " +
            "in a %d-dimensional space!") % (p.dim(), p.ambient_dim()))
        write_palp_matrix(p._pullback(p._vertices), input_file)
    input_file.close()
    output_file_name = tmp_filename()
    c = "%s <%s >%s" % (command, input_file_name, output_file_name)
    p = subprocess.Popen(c, shell=True, bufsize=2048,
                     stdin=subprocess.PIPE,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     close_fds=True)
    stdin, stdout, stderr = (p.stdin, p.stdout, p.stderr)
    err = stderr.read()
    if len(err) > 0:
        raise RuntimeError, ("Error executing \"%s\" for a polytope sequence!"
            + "\nOutput:\n%s") % (command, err)
    os.remove(input_file_name)
    try:
        p.terminate()
    except OSError:
        pass
    return output_file_name

def _read_nef_x_partitions(data):
    r"""
    Read all nef-partitions for one polytope from a string or an open
    file.

    ``data`` should be an output of nef.x.

    Returns the sequence of nef-partitions. Each nef-partition is given
    as a sequence of integers.

    If there are no nef-partitions, returns the empty sequence. If the
    string is empty or EOF is reached, raises ValueError.

    TESTS::

        sage: o = lattice_polytope.octahedron(3)
        sage: s = o.nef_x("-N -p")
        sage: print s # random
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


    -  ``data`` - an opened file with incidence
       information. The first line will be skipped, each consecutive line
       contains incidence information for all faces of one dimension, the
       first word of each line is a comment and is dropped.

    -  ``dim`` - dimension of the polytope.


    OUTPUT: a sequence F, such that F[d][i] is a sequence of vertices
    or facets corresponding to the i-th d-dimensional face.

    TESTS::

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
    Compute all cached data for all given ``polytopes`` and
    their polars.

    This functions does it MUCH faster than member functions of
    ``LatticePolytope`` during the first run. So it is recommended to
    use this functions if you work with big sets of data. None of the
    polytopes in the given sequence should be constructed as the polar
    polytope to another one.

    INPUT: a sequence of lattice polytopes.

    EXAMPLES: This function has no output, it is just a fast way to
    work with long sequences of polytopes. Of course, you can use short
    sequences as well::

        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.all_cached_data([o])
        sage: o.faces()
        [
        [[0], [1], [2], [3], [4], [5]],
        [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]],
        [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]
        ]

    However, you cannot use it for polytopes that are constructed as
    polar polytopes of others::

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
    Compute faces for all given ``polytopes``.

    This functions does it MUCH faster than member functions of
    ``LatticePolytope`` during the first run. So it is recommended to
    use this functions if you work with big sets of data.

    INPUT: a sequence of lattice polytopes.

    EXAMPLES: This function has no output, it is just a fast way to
    work with long sequences of polytopes. Of course, you can use short
    sequences as well::

        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.all_faces([o])
        sage: o.faces()
        [
        [[0], [1], [2], [3], [4], [5]],
        [[1, 5], [0, 5], [0, 1], [3, 5], [1, 3], [4, 5], [0, 4], [3, 4], [1, 2], [0, 2], [2, 3], [2, 4]],
        [[0, 1, 5], [1, 3, 5], [0, 4, 5], [3, 4, 5], [0, 1, 2], [1, 2, 3], [0, 2, 4], [2, 3, 4]]
        ]

    However, you cannot use it for polytopes that are constructed as
    polar polytopes of others::

        sage: lattice_polytope.all_faces([o.polar()])
        Traceback (most recent call last):
        ...
        ValueError: Cannot read face structure for a polytope constructed as polar, use _compute_faces!
    """
    result_name = _palp("poly.x -fi", polytopes, reduce_dimension=True)
    result = open(result_name)
    for p in polytopes:
        p._read_faces(result)
    result.close()
    os.remove(result_name)

def all_nef_partitions(polytopes, keep_symmetric=False):
    r"""
    Compute nef-partitions for all given ``polytopes``.

    This functions does it MUCH faster than member functions of
    ``LatticePolytope`` during the first run. So it is recommended to
    use this functions if you work with big sets of data.

    Note: member function ``is_reflexive`` will be called
    separately for each polytope. It is strictly recommended to call
    ``all_polars`` on the sequence of
    ``polytopes`` before using this function.

    INPUT: a sequence of lattice polytopes.

    EXAMPLES: This function has no output, it is just a fast way to
    work with long sequences of polytopes. Of course, you can use short
    sequences as well::

        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.all_nef_partitions([o])
        sage: o.nef_partitions()
        [
        Nef-partition {0, 1, 3} U {2, 4, 5},
        Nef-partition {0, 1, 3, 4} U {2, 5} (direct product),
        Nef-partition {0, 1, 2} U {3, 4, 5},
        Nef-partition {0, 1, 2, 3} U {4, 5},
        Nef-partition {0, 1, 2, 3, 4} U {5} (projection)
        ]

    You cannot use this function for non-reflexive polytopes::

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
    keys = "-N -V -D -P -p"
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
    Compute lattice points for all given ``polytopes``.

    This functions does it MUCH faster than member functions of
    ``LatticePolytope`` during the first run. So it is recommended to
    use this functions if you work with big sets of data.

    INPUT: a sequence of lattice polytopes.

    EXAMPLES: This function has no output, it is just a fast way to
    work with long sequences of polytopes. Of course, you can use short
    sequences as well::

        sage: o = lattice_polytope.octahedron(3)
        sage: lattice_polytope.all_points([o])
        sage: o.points()
        [ 1  0  0 -1  0  0  0]
        [ 0  1  0  0 -1  0  0]
        [ 0  0  1  0  0 -1  0]
    """
    result_name = _palp("poly.x -fp", polytopes, reduce_dimension=True)
    result = open(result_name)
    for p in polytopes:
        p._points = p._embed(read_palp_matrix(result))
        if p._points.nrows() == 0:
            raise RuntimeError, ("Cannot read points of a polytope!"
                                                        +"\nPolytope: %s" % p)
    result.close()
    os.remove(result_name)

def all_polars(polytopes):
    r"""
    Compute polar polytopes for all reflexive and equations of facets
    for all non-reflexive ``polytopes``.

    ``all_facet_equations`` and ``all_polars`` are synonyms.

    This functions does it MUCH faster than member functions of
    ``LatticePolytope`` during the first run. So it is recommended to
    use this functions if you work with big sets of data.

    INPUT: a sequence of lattice polytopes.

    EXAMPLES: This function has no output, it is just a fast way to
    work with long sequences of polytopes. Of course, you can use short
    sequences as well::

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

# Synonym for the above function
all_facet_equations = all_polars


_always_use_files = True


def always_use_files(new_state=None):
    r"""
    Set or get the way of using PALP for lattice polytopes.

    INPUT:

    - ``new_state`` - (default:None) if specified, must be ``True`` or ``False``.

    OUTPUT: The current state of using PALP. If ``True``, files are used
    for all calls to PALP, otherwise pipes are used for single polytopes.
    While the latter may have some advantage in speed, the first method
    is more reliable when working with large outputs. The initial state
    is ``True``.

    EXAMPLES::

        sage: lattice_polytope.always_use_files()
        True
        sage: p = LatticePolytope(matrix([1, 20]))
        sage: p.npoints()
        20

    Now let's use pipes instead of files::

        sage: lattice_polytope.always_use_files(False)
        False
        sage: p = LatticePolytope(matrix([1, 20]))
        sage: p.npoints()
        20
    """
    global _always_use_files
    if new_state != None:
        _always_use_files = new_state
    return _always_use_files


def convex_hull(points):
    r"""
    Compute the convex hull of the given points.

    .. note::

       ``points`` might not span the space. Also, it fails for large
       numbers of vertices in dimensions 4 or greater

    INPUT:


    -  ``points`` - a list that can be converted into
       vectors of the same dimension over ZZ.


    OUTPUT: list of vertices of the convex hull of the given points (as
    vectors).

    EXAMPLES: Let's compute the convex hull of several points on a line
    in the plane::

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
    Use the function ``f`` to filter polytopes in a list.

    INPUT:


    -  ``f`` - filtering function, it must take one
       argument, a lattice polytope, and return True or False.

    -  ``polytopes`` - list of polytopes.

    -  ``subseq`` - (default: None) list of integers. If it
       is specified, only polytopes with these numbers will be
       considered.

    -  ``print_numbers`` - (default: False) if True, the
       number of the current polytope will be printed on the screen before
       calling ``f``.


    OUTPUT: a list of integers -- numbers of polytopes in the given
    list, that satisfy the given condition (i.e. function
    ``f`` returns True) and are elements of subseq, if it
    is given.

    EXAMPLES: Consider a sequence of octahedrons::

        sage: polytopes = Sequence([lattice_polytope.octahedron(n) for n in range(2, 7)], cr=True)
        sage: polytopes
        [
        An octahedron: 2-dimensional, 4 vertices.,
        An octahedron: 3-dimensional, 6 vertices.,
        An octahedron: 4-dimensional, 8 vertices.,
        An octahedron: 5-dimensional, 10 vertices.,
        An octahedron: 6-dimensional, 12 vertices.
        ]

    This filters octahedrons of dimension at least 4::

        sage: lattice_polytope.filter_polytopes(lambda p: p.dim() >= 4, polytopes)
        [2, 3, 4]

    For long tests you can see the current progress::

        sage: lattice_polytope.filter_polytopes(lambda p: p.nvertices() >= 10, polytopes, print_numbers=True)
        0
        1
        2
        3
        4
        [3, 4]

    Here we consider only some of the polytopes::

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


def integral_length(v):
    """
    Compute the integral length of a given rational vector.

    INPUT:

    -  ``v`` - any object which can be converted to a list of rationals

    OUTPUT: Rational number ``r`` such that ``v = r u``, where ``u`` is the
    primitive integral vector in the direction of ``v``.

    EXAMPLES::

        sage: lattice_polytope.integral_length([1, 2, 4])
        1
        sage: lattice_polytope.integral_length([2, 2, 4])
        2
        sage: lattice_polytope.integral_length([2/3, 2, 4])
        2/3
    """
    data = [QQ(e) for e in list(v)]
    ns = [e.numerator() for e in data]
    ds = [e.denominator() for e in data]
    return gcd(ns)/lcm(ds)


def minkowski_sum(points1, points2):
    r"""
    Compute the Minkowski sum of two convex polytopes.

    .. note::

       Polytopes might not be of maximal dimension.

    INPUT:


    -  ``points1, points2`` - lists of objects that can be
       converted into vectors of the same dimension, treated as vertices
       of two polytopes.


    OUTPUT: list of vertices of the Minkowski sum, given as vectors.

    EXAMPLES: Let's compute the Minkowski sum of two line segments::

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

    EXAMPLES: Here are 3- and 4-dimensional octahedrons::

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

    There exists only one octahedron of each dimension::

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


    -  ``points`` - lattice points given as columns of a
       matrix


    OUTPUT: matrix of relations between given points with non-negative
    integer coefficients

    EXAMPLES: This is a 3-dimensional reflexive polytope::

        sage: m = matrix(ZZ,[[1, 0, -1, 0, -1],
        ...                  [0, 1, -1, 0,  0],
        ...                  [0, 0,  0, 1, -1]])
        ...
        sage: p = LatticePolytope(m)
        sage: p.points()
        [ 1  0 -1  0 -1  0]
        [ 0  1 -1  0  0  0]
        [ 0  0  0  1 -1  0]

    We can compute linear relations between its points in the following
    way::

        sage: p.points().transpose().kernel().echelonized_basis_matrix()
        [ 1  0  0  1  1  0]
        [ 0  1  1 -1 -1  0]
        [ 0  0  0  0  0  1]

    However, the above relations may contain negative and rational
    numbers. This function transforms them in such a way, that all
    coefficients are non-negative integers::

        sage: lattice_polytope.positive_integer_relations(p.points())
        [1 0 0 1 1 0]
        [1 1 1 0 0 0]
        [0 0 0 0 0 1]
        sage: lattice_polytope.positive_integer_relations(ReflexivePolytope(2,1).vertices())
        [2 1 1]
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
        relations.rescale_row(i, 1/integral_length(relations[i]))
    return relations.change_ring(ZZ)


def projective_space(dim):
    r"""
    Return a simplex of the given dimension, corresponding to
    `P_{dim}`.

    EXAMPLES: We construct 3- and 4-dimensional simplexes::

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


    -  ``file_name`` - the name of a file with VERTICES of
       polytopes

    -  ``desc`` - a string, that will be used for creating
       polytope descriptions. By default it will be set to 'A lattice
       polytope #%d from "filename"' + and will be used as ``desc %
       n`` where ``n`` is the number of the polytope in
       the file (*STARTING WITH ZERO*).


    OUTPUT: a sequence of polytopes

    EXAMPLES: We use poly.x to compute polar polytopes of 2- and
    3-octahedrons and read them::

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
    m = read_palp_matrix(f)
    while m.nrows() != 0:
        polytopes.append(LatticePolytope(m,
            desc=(desc % n), compute_vertices=False, copy_vertices=False))
        n += 1
        m = read_palp_matrix(f)
    f.close()
    return polytopes


def read_palp_matrix(data):
    r"""
    Read and return an integer matrix from a string or an opened file.

    First input line must start with two integers m and n, the number
    of rows and columns of the matrix. The rest of the first line is
    ignored. The next m lines must contain n numbers each.

    If m>n, returns the transposed matrix. If the string is empty or EOF
    is reached, returns the empty matrix, constructed by
    ``matrix()``.

    EXAMPLES::

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
    Convert a Sage matrix to the string representation of Maxima.

    EXAMPLE::

        sage: m = matrix(ZZ,2)
        sage: lattice_polytope.sage_matrix_to_maxima(m)
        matrix([0,0],[0,0])
    """
    return maxima("matrix("+",".join(str(v.list()) for v in m.rows())+")")


def set_palp_dimension(d):
    r"""
    Set the dimension for PALP calls to ``d``.

    INPUT:

    - ``d`` -- an integer from the list [4,5,6,11] or ``None``.

    OUTPUT:

    - none.

    PALP has many hard-coded limits, which must be specified before
    compilation, one of them is dimension. Sage includes several versions with
    different dimension settings (which may also affect other limits and enable
    certain features of PALP). You can change the version which will be used by
    calling this function. Such a change is not done automatically for each
    polytope based on its dimension, since depending on what you are doing it
    may be necessary to use dimensions higher than that of the input polytope.

    EXAMPLES:

    By default, it is not possible to create the 7-dimensional simplex with
    vertices at the basis of the 8-dimensional space::

        sage: LatticePolytope(identity_matrix(8))
        Traceback (most recent call last):
        ...
        ValueError: Error executing "poly.x -fv" for the given polytope!
        Polytope: A lattice polytope: 7-dimensional, 8 vertices.
        Vertices:
        [ 0 -1 -1 -1 -1 -1 -1 -1]
        [ 0  1  0  0  0  0  0  0]
        [ 0  0  1  0  0  0  0  0]
        [ 0  0  0  1  0  0  0  0]
        [ 0  0  0  0  1  0  0  0]
        [ 0  0  0  0  0  1  0  0]
        [ 0  0  0  0  0  0  1  0]
        Output:
        Please increase POLY_Dmax to at least 7

    However, we can work with this polytope by changing PALP dimension to 11::

        sage: lattice_polytope.set_palp_dimension(11)
        sage: LatticePolytope(identity_matrix(8))
        A lattice polytope: 7-dimensional, 8 vertices.

    Let's go back to default settings::

        sage: lattice_polytope.set_palp_dimension(None)
    """
    global _palp_dimension
    _palp_dimension = d


def skip_palp_matrix(data, n=1):
    r"""
    Skip matrix data in a file.

    INPUT:


    -  ``data`` - opened file with blocks of matrix data in
       the following format: A block consisting of m+1 lines has the
       number m as the first element of its first line.

    -  ``n`` - (default: 1) integer, specifies how many
       blocks should be skipped


    If EOF is reached during the process, raises ValueError exception.

    EXAMPLE: We create a file with vertices of the square and the cube,
    but read only the second set::

        sage: o2 = lattice_polytope.octahedron(2)
        sage: o3 = lattice_polytope.octahedron(3)
        sage: result_name = lattice_polytope._palp("poly.x -fe", [o2, o3])
        sage: f = open(result_name)
        sage: f.readlines()
        ['4 2  Vertices of P-dual <-> Equations of P\n',
         '  -1   1\n',
         '   1   1\n',
         '  -1  -1\n',
         '   1  -1\n',
         '8 3  Vertices of P-dual <-> Equations of P\n',
         '  -1  -1   1\n',
         '   1  -1   1\n',
         '  -1   1   1\n',
         '   1   1   1\n',
         '  -1  -1  -1\n',
         '   1  -1  -1\n',
         '  -1   1  -1\n',
         '   1   1  -1\n']
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


    -  ``m`` - a matrix over integers.

    -  ``ofile`` - a file opened for writing (default:
       stdout)

    -  ``comment`` - a string (default: empty) see output
       description

    -  ``format`` - a format string used to print matrix entries.


    OUTPUT: First line: number_of_rows number_of_columns comment
    Next number_of_rows lines: rows of the matrix.

    EXAMPLES: This functions is used for writing polytope vertices in
    PALP format::

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
        print s,
    else:
        ofile.write(s)
    for i in range(m.nrows()):
        s = " ".join(format % m[i,j] for j in range(m.ncols()))+"\n"
        if ofile is None:
            print s,
        else:
            ofile.write(s)
