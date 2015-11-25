r"""
Abstract objects and sets

These are base classes for the manifolds.

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2015-11-25): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.parent import Parent
from sage.misc.fast_methods import Singleton, WithEqualityById
from sage.categories.fields import Fields
from sage.categories.manifolds import Manifolds
from sage.categories.sets_cat import Sets
from sage.rings.integer import Integer
from sage.manifolds.point import ManifoldPoint

class AbstractObject(object):
    """
    An abstract object.

    An abstract object is a variable name and latex name.
    """
    def __init__(self, name, latex_name=None):
        """
        Initialize ``self``.
        """
        if not isinstance(name, str):
            raise TypeError("{} is not a string".format(name))
        self._name = name

        if latex_name is None:
            self._latex_name = self._name
        else:
            if not isinstance(latex_name, str):
                raise TypeError("{} is not a string".format(latex_name))
            self._latex_name = latex_name

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return self._name

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: A._latex_()
            'A'
            sage: B = A.subset('B', latex_name=r'\mathcal{B}')
            sage: B._latex_()
            '\\mathcal{B}'
            sage: latex(B)  # indirect doctest
            \mathcal{B}

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M._latex_()
            'M'
            sage: latex(M)
            M
            sage: M = Manifold(3, 'M', latex_name=r'\mathcal{M}',
            ....:              structure='topological')
            sage: M._latex_()
            '\\mathcal{M}'
            sage: latex(M)
            \mathcal{M}
        """
        return self._latex_name

class AbstractSet(AbstractObject, WithEqualityById, Parent):
    """
    An abstract set.

    An abstract set is an :class:`AbstractObject` along with its known
    subsets, supersets, intersections, unions, and open covers.
    """
    def __init__(self, name, latex_name, base=None, category=None):
        r"""
        Initialize ``self``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.subset('A'); A
            Subset A of the 2-dimensional topological manifold M
        """
        AbstractObject.__init__(self, name, latex_name)

        category = Sets().or_subcategory(category)
        Parent.__init__(self, base=base, category=category)

        self._supersets = set([self]) # subsets containing self
        self._subsets = set([self]) # subsets of self
        self._top_subsets = set([self]) # subsets contained in self but not
                                        # in another strict subset of self
        self._intersections = {} # dict. of intersections with other subsets
                                 # (key: subset name)
        self._unions = {} # dict. of unions with other subsets (key: subset
                          # name)
        self._open_covers = set() # set of open covers of self

    #### Methods required for any Parent in the category of sets:

    def _element_constructor_(self, coords=None, chart=None, name=None,
                              latex_name=None, check_coords=True):
        r"""
        Construct a point in the subset from its coordinates in some chart.

        INPUT:

        - ``coords`` -- (default: ``None``) either (i) the point coordinates
          (as a tuple or a list) in the chart ``chart`` or (ii) another point
          in the subset
        - ``chart`` -- (default: ``None``) chart in which the coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the subset's default chart
        - ``name`` -- (default: ``None``) name given to the point
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          point; if none is provided, the LaTeX symbol is set to ``name``
        - ``check_coords`` -- (default: ``True``) determines whether ``coords``
          are valid coordinates for the chart ``chart``; for symbolic
          coordinates, it is recommended to set ``check_coords`` to ``False``.

        OUTPUT:

        - an instance of :class:`~sage.manifolds.point.ManifoldPoint`
          representing a point in the current subset.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M((-2,3)); p  # coord in the default chart
            Point on the 2-dimensional topological manifold M
            sage: X(p)
            (-2, 3)

        A generic subset has no default chart, so the chart must be explicited::

            sage: A = M.subset('A')
            sage: p = A((-2,3), chart=X); p
            Point on the 2-dimensional topological manifold M
            sage: X(p)
            (-2, 3)
            sage: p.containing_set()
            Subset A of the 2-dimensional topological manifold M
            sage: p in A
            True

        Coordinates in a chart with some coordinate restrictions::

            sage: Y.<u,v> = M.chart('u:(-1,1) v:(-1,1)')
            sage: p = A((0,1/2), chart=Y); p
            Point on the 2-dimensional topological manifold M
            sage: Y(p)
            (0, 1/2)
            sage: p = A((0,1/2), chart=Y, check_coords=False); p
            Point on the 2-dimensional topological manifold M
            sage: Y(p)
            (0, 1/2)
            sage: p = A((3,1/2), chart=Y)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates (3, 1/2) are not valid on the Chart (M, (u, v))

        Specifying the name of the point::

            sage: p = A((-2,3), chart=X, name='p'); p
            Point p on the 2-dimensional topological manifold M

        A point as entry::

            sage: q = A(p); q
            Point p on the 2-dimensional topological manifold M
            sage: X(q)
            (-2, 3)

        """
        if isinstance(coords, ManifoldPoint):
            point = coords # for readability
            # This should actually never happen by the coercion framework...
            if point.parent() is self:
                return point
            if point in self:
                resu = self.element_class(self, name=point._name,
                                          latex_name=point._latex_name)
                for chart, coords in point._coordinates.iteritems():
                    resu._coordinates[chart] = coords
                return resu
            else:
                raise ValueError("the {}".format(point) +
                                 " is not in {}".format(self))
        return self.element_class(self, coords=coords, chart=chart, name=name,
                                  latex_name=latex_name, check_coords=check_coords)

    def _an_element_(self):
        r"""
        Construct some point in the subset.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.subset('A')
            sage: p = A._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p in A
            True

        """
        #!# should be improved...
        return self.element_class(self)

    #### End of methods required for any Parent in the category of sets

    def open_covers(self):
        r"""
        Return the list of open covers of the current subset.

        If the current subset, `A` say, is a subset of the manifold `M`, an
        *open cover* of `A` is list (indexed set) `(U_i)_{i\in I}` of
        open subsets of `M` such that

        .. MATH::

            A \subset \bigcup_{i \in I} U_i

        If `A` is open, we ask that the above inclusion is actually an
        identity:

        .. MATH::

            A = \bigcup_{i \in I} U_i

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.open_covers()
            [[2-dimensional topological manifold M]]
            sage: U = M.open_subset('U')
            sage: U.open_covers()
            [[Open subset U of the 2-dimensional topological manifold M]]
            sage: A = U.open_subset('A')
            sage: B = U.open_subset('B')
            sage: U.declare_union(A,B)
            sage: U.open_covers()
            [[Open subset U of the 2-dimensional topological manifold M],
             [Open subset A of the 2-dimensional topological manifold M,
              Open subset B of the 2-dimensional topological manifold M]]
            sage: V = M.open_subset('V')
            sage: M.declare_union(U,V)
            sage: M.open_covers()
            [[2-dimensional topological manifold M],
             [Open subset U of the 2-dimensional topological manifold M,
              Open subset V of the 2-dimensional topological manifold M],
             [Open subset A of the 2-dimensional topological manifold M,
              Open subset B of the 2-dimensional topological manifold M,
              Open subset V of the 2-dimensional topological manifold M]]

        """
        return sorted(sorted(oc, key=lambda x: x._name) for oc in self._open_covers)

    def subsets(self):
        r"""
        Return the set of subsets that have been defined on the current subset.

        OUTPUT:

        - A Python set containing all the subsets that have been defined on
          the current subset.

        .. NOTE::

            To get the subsets as a list, used the method
            :meth:`list_of_subsets` instead.

        EXAMPLE:

        Subsets of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: M.subsets()  # random (set output)
            {Subset V of the 2-dimensional topological manifold M,
             2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M}
            sage: type(M.subsets())
            <type 'frozenset'>
            sage: U in M.subsets()
            True

        The method :meth:`list_of_subsets` returns a list (sorted
        alphabetically by the subset names) instead of a set::

            sage: M.list_of_subsets()
            [2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M,
             Subset V of the 2-dimensional topological manifold M]

        """
        return frozenset(self._subsets)

    # Do we need this?
    def list_of_subsets(self):
        r"""
        Return the list of subsets that have been defined on the current
        subset.

        The list is sorted by the alphabetical names of the subsets.

        OUTPUT:

        - A list containing all the subsets that have been defined on
          the current subset.

        .. NOTE::

            To get the subsets as a Python set, used the method
            :meth:`subsets` instead.

        EXAMPLE:

        Subsets of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: M.list_of_subsets()
            [2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M,
             Subset V of the 2-dimensional topological manifold M]

        The method :meth:`subsets` returns a set instead of a list::

            sage: M.subsets()  # random (set output)
            {Subset V of the 2-dimensional topological manifold M,
             2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M}

        """
        return sorted(self._subsets, key=lambda x: x._name)

    def subset(self, name, latex_name=None, is_open=False):
        r"""
        Create a subset of ``self``.

        INPUT:

        - ``name`` -- name given to the subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``is_open`` -- (default: ``False``) if ``True``, the created subset
          is assumed to be open with respect to the manifold's topology

        OUTPUT:

        - the subset, as an instance of :class:`ManifoldSubset`, or
          of the derived class
          :class:`~sage.manifolds.manifold.Manifold` if ``is_open``
          is ``True``.

        EXAMPLES:

        Creating a subset of a manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A'); a
            Subset A of the 2-dimensional topological manifold M

        Creating a subset of A::

            sage: b = a.subset('B', latex_name=r'\mathcal{B}'); b
            Subset B of the 2-dimensional topological manifold M
            sage: latex(b)
            \mathcal{B}

        We have then::

            sage: b.is_subset(a)
            True
            sage: b in a.subsets()
            True
        """
        if is_open:
            return self.open_subset(name, latex_name=latex_name)
        from sage.manifolds.subset import ManifoldSubset
        res = ManifoldSubset(self.manifold(), name,
                                        latex_name=latex_name)
        res._supersets.update(self._supersets)
        for sd in self._supersets:
            sd._subsets.add(res)
        self._top_subsets.add(res)
        return res

    def get_subset(self, name):
        r"""
        Get a subset by its name.

        The subset must have been previously created by the method
        :meth:`subset` (or
        :meth:`~sage.manifolds.manifold.Manifold.open_subset`)

        INPUT:

        - ``name`` -- (string) name of the subset

        OUTPUT:

        - instance of :class:`ManifoldSubset` (or
          of the derived class
          :class:`~sage.manifolds.manifold.Manifold` for an open
          subset) representing the subset whose name is ``name``.

        EXAMPLES::

            sage: M = Manifold(4, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = A.subset('B')
            sage: U = M.open_subset('U')
            sage: M.list_of_subsets()
            [Subset A of the 4-dimensional topological manifold M,
             Subset B of the 4-dimensional topological manifold M,
             4-dimensional topological manifold M,
             Open subset U of the 4-dimensional topological manifold M]
            sage: M.get_subset('A')
            Subset A of the 4-dimensional topological manifold M
            sage: M.get_subset('A') is A
            True
            sage: M.get_subset('B') is B
            True
            sage: A.get_subset('B') is B
            True
            sage: M.get_subset('U')
            Open subset U of the 4-dimensional topological manifold M
            sage: M.get_subset('U') is U
            True

        """
        for ss in self._subsets:
            if ss._name == name:
                return ss
        raise ValueError("no subset of name '{}' found".format(name))

    def is_subset(self, other):
        r"""
        Return ``True`` iff ``self`` is included in ``other``.

        EXAMPLES:

        Subsets on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A')
            sage: b = a.subset('B')
            sage: c = M.subset('C')
            sage: a.is_subset(M)
            True
            sage: b.is_subset(a)
            True
            sage: b.is_subset(M)
            True
            sage: a.is_subset(b)
            False
            sage: c.is_subset(a)
            False
        """
        return self in other._subsets

    def __contains__(self, point):
        r"""
        Check whether ``point`` is contained in ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.subset('A')
            sage: p = A((-2,3), chart=X); p
            Point on the 2-dimensional topological manifold M
            sage: A.__contains__(p)
            True
            sage: p in A  # indirect doctest
            True
            sage: A.__contains__(A.an_element())
            True
            sage: q = M((0,0), chart=X); q
            Point on the 2-dimensional topological manifold M
            sage: A.__contains__(q)
            False
        """
        # for efficiency, a quick test first:
        if point.parent() is self:
            return True
        if point.parent().is_subset(self):
            return True
        #!# should be improved once coordinate definition have been introduced
        # in ManifoldSubset
        return False

    def declare_union(self, dom1, dom2):
        r"""
        Declare that the current subset is the union of two subsets,
        i.e. that

        .. MATH::

            U = U_1 \cup U_2,

        where `U` is the current subset, `U_1 \subset U` and `U_2 \subset U`.

        INPUT:

        - ``dom1`` -- subset `U_1`
        - ``dom2`` -- subset `U_2`

        EXAMPLE::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: M.declare_union(A, B)
            sage: A.union(B)
            2-dimensional topological manifold M

        """
        if dom1 == dom2:
            if dom1 != self:
                raise ValueError("the union of two identical sets must be " +
                                 "this set")
            return
        if not dom1.is_subset(self):
            raise TypeError("the {} is not a subset of ".format(dom1) +
                            "the {}".format(self))
        if not dom2.is_subset(self):
            raise TypeError("the {} is not a subset of ".format(dom2) +
                            "the {}".format(self))
        dom1._unions[dom2._name] = self
        dom2._unions[dom1._name] = self
        for oc1 in dom1._open_covers:
            for oc2 in dom2._open_covers:
                oc = set(oc1) # Make a (shallow) copy
                for s in oc2:
                    if s not in oc:
                        oc.add(s)
            self._open_covers.add(frozenset(oc))

    def is_open(self):
        """
        Return if ``self`` is an open set.
        """
        return True

    def point(self, coords=None, chart=None, name=None, latex_name=None):
        r"""
        Define a point in ``self``.

        See :class:`~sage.manifolds.point.ManifoldPoint` for a
        complete documentation.

        INPUT:

        - ``coords`` -- the point coordinates (as a tuple or a list) in the
          chart specified by ``chart``
        - ``chart`` -- (default: ``None``) chart in which the point coordinates
          are given; if ``None``, the coordinates are assumed to refer to
          the default chart of the current subset
        - ``name`` -- (default: ``None``) name given to the point
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          point; if ``None``, the LaTeX symbol is set to ``name``

        OUTPUT:

        - the declared point, as an instance of
          :class:`~sage.manifolds.point.ManifoldPoint`.

        EXAMPLES:

        Points on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: p = M.point((1,2), name='p'); p
            Point p on the 2-dimensional topological manifold M
            sage: p in M
            True
            sage: a = M.open_subset('A')
            sage: c_uv.<u,v> = a.chart()
            sage: q = a.point((-1,0), name='q'); q
            Point q on the 2-dimensional topological manifold M
            sage: q in a
            True
            sage: p._coordinates
            {Chart (M, (x, y)): (1, 2)}
            sage: q._coordinates
            {Chart (A, (u, v)): (-1, 0)}
        """
        return self.element_class(self, coords=coords, chart=chart,
                                  name=name, latex_name=latex_name)

    Element = ManifoldPoint

