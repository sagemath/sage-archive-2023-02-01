r"""
Point collections

This module was designed as a part of framework for toric varieties
(:mod:`~sage.schemes.generic.toric_variety`,
:mod:`~sage.schemes.generic.fano_toric_variety`).

AUTHORS:

- Andrey Novoseltsev (2011-04-25): initial version, based on cone module.

EXAMPLES:

The idea behind :class:`point collections <PointCollection>` is to have a
container for points of the same space that

* behaves like a tuple *without significant performance penalty*::

    sage: from sage.geometry.point_collection import PointCollection
    sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
    sage: c = PointCollection(c.rays(), c.lattice())
    sage: c[1]
    N(1, 0, 1)
    sage: for point in c: point
    N(0, 0, 1)
    N(1, 0, 1)
    N(0, 1, 1)
    N(1, 1, 1)

* prints like a matrix with points represented by columns and with some
  indication of the ambient space::

    sage: c
    [0 1 0 1]
    [0 0 1 1]
    [1 1 1 1]
    in 3-d lattice N

* allows (cached) access to alternative representations::

    sage: c.set()
    frozenset([N(0, 1, 1), N(1, 1, 1), N(0, 0, 1), N(1, 0, 1)])

* allows introduction of additional methods::

    sage: c.basis()
    [0 1 0]
    [0 0 1]
    [1 1 1]
    in 3-d lattice N

Examples of natural point collections include ray and line generators of cones,
vertices and points of polytopes, normals to facets, their subcollections, etc.

Using this class for all of the above cases allows for unified interface *and*
cache sharing. Suppose that `\Delta` is a reflexive polytope. Then the same
point collection can be linked as

#. vertices of `\Delta`;
#. facet normals of its polar `\Delta^\circ`;
#. ray generators of the face fan of `\Delta`;
#. ray generators of the normal fan of `\Delta`.

If all these objects are in use and, say, a matrix representation was computed
for one of them, it becomes available to all others as well, eliminating the
need to spend time and memory four times.
"""

#*****************************************************************************
#       Copyright (C) 2011 Andrey Novoseltsev <novoselt@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object cimport SageObject

from sage.matrix.all import column_matrix
from sage.misc.all import latex


def is_PointCollection(x):
    r"""
    Check if ``x`` is a :class:`point collection <PointCollection>`.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a point collection and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.point_collection import is_PointCollection
        sage: is_PointCollection(1)
        False
        sage: from sage.geometry.point_collection import PointCollection
        sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
        sage: is_PointCollection(c.rays()) # TODO: make it True
        False
        sage: c = PointCollection(c.rays(), c.lattice())
        sage: is_PointCollection(c)
        True
    """
    return isinstance(x, PointCollection)


cdef class PointCollection(SageObject):
    r"""
    Create a point collection.

    .. WARNING::

        No correctness check or normalization is performed on the input data.
        This class is designed for internal operations and you probably should
        not use it directly.

    Point collections are immutable, but cache most of the returned values.

    INPUT:

    - ``points`` -- an iterable structure of immutable elements of ``module``,
      if ``points`` are already accessible to you as a :class:`tuple`, it is
      preferable to use it for speed and memory consumption reasons;

    - ``module`` -- an ambient module for ``points``. If ``None``, it will be
      determined as :func:`parent` of the first point. Of course, this cannot
      be done if there are no points, so in this case you must give an
      appropriate ``module`` directly. Note that ``None`` is *not* the default
      value - you always *must* give this argument explicitly, even if it is
      ``None``.

    OUTPUT:

    - a point collection.
    """

    cdef:
        tuple _points
        object _module
        # cache attributes
        PointCollection _basis
        object _matrix
        frozenset _set

    def __init__(self, points, module=None):
        r"""
        See :class:`PointCollection` for documentation.

        TESTS::

            sage: from sage.geometry.point_collection import PointCollection
            sage: v = vector([1,0])
            sage: v.set_immutable()
            sage: c = PointCollection([v], ZZ^2)
            sage: c.module()
            Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: c = PointCollection([v], None)
            sage: c.module()  # Determined automatically
            Ambient free module of rank 2
            over the principal ideal domain Integer Ring
            sage: TestSuite(c).run()
        """
        super(PointCollection, self).__init__()
        self._points = tuple(points)
        self._module = self._points[0].parent() if module is None else module

    def __call__(self, *args):
        r"""
        Return a subcollection of ``self``.

        INPUT:

        - a list of integers (as a single or many arguments).

        OUTPUT:

        - a :class:`point collection <PointCollection>`.

        TESTS::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c()
            []
            in 3-d lattice N
            sage: c(2,1)
            [0 1]
            [1 0]
            [1 1]
            in 3-d lattice N
            sage: c(range(4)) == c
            True
        """
        if len(args) == 1:
            try:
                args = tuple(args[0])
            except TypeError:
                pass
        # Avoid creating a copy of self
        if len(args) == len(self) and args == tuple(range(len(self))):
            return self
        else:
            return PointCollection([self[i] for i in args], self._module)

    def __cmp__(self, right):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        - 0 if ``right`` is of the same type as ``self`` (i.e. it is another
          :class:`point collection <PointCollection>`), they have the same
          :meth:`module`, and their points are the same and listed in the same
          order. 1 or -1 otherwise.

        TESTS::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: cmp(c, c)
            0
            sage: cmp(c, 1) * cmp(1, c)
            -1
        """
        c = cmp(type(self), type(right))
        if c != 0:
            return c
        cdef PointCollection pc = right
        c = cmp(self._module, pc._module)
        if c != 0:
            return c
        return cmp(self._points, pc._points)

    def __getitem__(self, n):
        r"""
        Return the ``n``-th point of ``self``.

        INPUT:

        - ``n`` -- an integer.

        OUTPUT:

        - a point, an element of the ambient :meth:`module` of ``self``.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c[0]
            N(0, 0, 1)
        """
        return self._points[n]

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        OUTPUT:

        - an integer.

        TESTS::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: hash(c) == hash(c)
            True
        """
        return hash(self._points)

    def __iter__(self):
        r"""
        Return an iterator over points of ``self``.

        OUTPUT:

        - an iterator.

        TESTS::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: for point in c: print point
            N(0, 0, 1)
            N(1, 0, 1)
            N(0, 1, 1)
            N(1, 1, 1)
        """
        return iter(self._points)

    def __len__(self):
        r"""
        Return the number of points in ``self``.

        OUTPUT:

        - an integer.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: len(c)
            4
        """
        return len(self._points)

    def __reduce__(self):
        r"""
        Prepare ``self`` for pickling.

        OUTPUT:

        - a tuple, currently the class name and a tuple consisting of points
          and the ambient module.

        TESTS::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: loads(dumps(c))
            [0 1 0 1]
            [0 0 1 1]
            [1 1 1 1]
            in 3-d lattice N
            sage: loads(dumps(c)) == c
            True
        """
        return (PointCollection, (self._points, self._module))

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: print c._latex_()
            \left(\begin{array}{rrrr}
            0 & 1 & 0 & 1 \\
            0 & 0 & 1 & 1 \\
            1 & 1 & 1 & 1
            \end{array}\right)_{N}
        """
        return r"%s_{%s}" % (latex(self.matrix()), latex(self.module()))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: print c._repr_()
            [0 1 0 1]
            [0 0 1 1]
            [1 1 1 1]
            in 3-d lattice N
        """
        return "%s\nin %s" % (self.matrix(), self.module())

    def basis(self):
        r"""
        Return a linearly independent subset of points of ``self``.

        OUTPUT:

        - a :class:`point collection <PointCollection>` giving a random (but
          fixed) choice of an `\RR`-basis for the vector space spanned by the
          points of ``self``.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c.basis()
            [0 1 0]
            [0 0 1]
            [1 1 1]
            in 3-d lattice N

        Calling this method twice will always return *exactly the same* point
        collection::

            sage: c.basis().basis() is c.basis()
            True
        """
        if self._basis is None:
            self._basis = self(self.matrix().pivots())
        return self._basis

    def cardinality(self):
        r"""
        Return the number of points in ``self``.

        OUTPUT:

        - an integer.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c.cardinality()
            4
        """
        return len(self._points)

    def cartesian_product(self, other, module=None):
        r"""
        Return the Cartesian product of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a :class:`point collection <PointCollection>`;

        - ``module`` -- (optional) the ambient module for the result. By
          default, the direct sum of the ambient modules of ``self`` and
          ``other`` is constructed.

        OUTPUT:

        - a :class:`point collection <PointCollection>`.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c.cartesian_product(c)
            [0 1 0 1]
            [0 1 0 1]
            [1 1 1 1]
            [0 0 1 1]
            [0 0 1 1]
            [1 1 1 1]
            in 6-d lattice N+N
        """
        assert is_PointCollection(other)
        if module is None:
            module = self._module.direct_sum(other.module())
        P = [list(p) for p in self]
        Q = [list(q) for q in other]
        PQ = [module(p + q) for q in Q for p in P]
        for pq in PQ:
            pq.set_immutable()
        return PointCollection(PQ, module)

    def dimension(self):
        r"""
        Return the dimension of the space spanned by points of ``self``.

        .. note:: You can use either :meth:`dim` or :meth:`dimension`.

        OUTPUT:

        - an integer.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c.dimension()
            2
            sage: c.dim()
            2
        """
        return self.matrix().rank()

    dim = dimension

    def dual_module(self):
        r"""
        Return the dual of the ambient module of ``self``.

        OUTPUT:

        - a :class:`module <FreeModule_generic>`. If possible (that is, if the
          ambient :meth:`module` `M` of ``self`` has a ``dual()`` method), the
          dual module is returned. Otherwise, `R^n` is returned, where `n` is
          the dimension of `M` and `R` is its base ring.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c.dual_module()
            3-d lattice M
        """
        M = self._module
        try:
            return M.dual()
        except AttributeError:
            # TODO: add support for torsion modules as well?
            return M.base_ring() ** M.dimension()

    def matrix(self):
        r"""
        Return a matrix whose columns are points of ``self``.

        OUTPUT:

        - a :class:`matrix <Matrix>`.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c.matrix()
            [0 1 0 1]
            [0 0 1 1]
            [1 1 1 1]
        """
        if self._matrix is None:
            M = column_matrix(self._module.base_ring(), len(self._points),
                              self._module.degree(), self._points)
            M.set_immutable()
            self._matrix = M
        return self._matrix

    def module(self):
        r"""
        Return the ambient module of ``self``.

        OUTPUT:

        - a :class:`module <FreeModule_generic>`.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c.module()
            3-d lattice N
        """
        return self._module

    def set(self):
        r"""
        Return points of ``self`` as a :class:`frozenset`.

        OUTPUT:

        - a :class:`frozenset`.

        EXAMPLES::

            sage: from sage.geometry.point_collection import PointCollection
            sage: c = Cone([(0,0,1), (1,0,1), (0,1,1), (1,1,1)])
            sage: c = PointCollection(c.rays(), c.lattice())
            sage: c.set()
            frozenset([N(0, 1, 1), N(1, 1, 1), N(0, 0, 1), N(1, 0, 1)])
        """
        if self._set is None:
            self._set = frozenset(self._points)
        return self._set
