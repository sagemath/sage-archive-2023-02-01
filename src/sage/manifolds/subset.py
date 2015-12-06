r"""
Subsets of topological manifolds

The class :class:`TopologicalManifoldSubset` implements generic subsets of a
topological manifold. Open subsets are implemented by the class
:class:`~sage.manifolds.manifold.TopologicalManifold` (since an open subset of
a manifold is a manifold by itself), which inherits from
:class:`TopologicalManifoldSubset`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version

REFERENCES:

- [Lee11]_ J.M. Lee : *Introduction to Topological Manifolds*, 2nd ed.,
  Springer (New York) (2011)

EXAMPLES:

Two subsets on a manifold::

    sage: M = Manifold(2, 'M', structure='topological')
    sage: a = M.subset('A'); a
    Subset A of the 2-dimensional topological manifold M
    sage: b = M.subset('B'); b
    Subset B of the 2-dimensional topological manifold M
    sage: M.list_of_subsets()
    [Subset A of the 2-dimensional topological manifold M,
     Subset B of the 2-dimensional topological manifold M,
     2-dimensional topological manifold M]

The intersection of the two subsets::

    sage: c = a.intersection(b); c
    Subset A_inter_B of the 2-dimensional topological manifold M

Their union::

    sage: d = a.union(b); d
    Subset A_union_B of the 2-dimensional topological manifold M

Lists of subsets after the above operations::

    sage: M.list_of_subsets()
    [Subset A of the 2-dimensional topological manifold M,
     Subset A_inter_B of the 2-dimensional topological manifold M,
     Subset A_union_B of the 2-dimensional topological manifold M,
     Subset B of the 2-dimensional topological manifold M,
     2-dimensional topological manifold M]
    sage: a.list_of_subsets()
    [Subset A of the 2-dimensional topological manifold M,
     Subset A_inter_B of the 2-dimensional topological manifold M]
    sage: c.list_of_subsets()
    [Subset A_inter_B of the 2-dimensional topological manifold M]
    sage: d.list_of_subsets()
    [Subset A of the 2-dimensional topological manifold M,
     Subset A_inter_B of the 2-dimensional topological manifold M,
     Subset A_union_B of the 2-dimensional topological manifold M,
     Subset B of the 2-dimensional topological manifold M]

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
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
from sage.manifolds.point import TopologicalManifoldPoint

class TopologicalManifoldSubset(UniqueRepresentation, Parent):
    r"""
    Subset of a topological manifold.

    The class :class:`TopologicalManifoldSubset` inherits from the generic Sage
    class :class:`~sage.structure.parent.Parent` and is declared to belong to
    the category of facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`).
    The corresponding element class is
    :class:`~sage.manifolds.point.TopologicalManifoldPoint`. A subset acts
    as a facade for the true parent of its points, which is the whole manifold
    (see example below).

    Note that open subsets are not implemented directly by this class, but
    by the derived class :class:`~sage.manifolds.manifold.TopologicalManifold`
    (an open subset of a topological manifold being itself a topological
    manifold).

    INPUT:

    - ``manifold`` -- topological manifold on which the subset is defined
    - ``name`` -- string; name (symbol) given to the subset
    - ``latex_name`` --  (default: ``None``) string: LaTeX symbol to denote the
      subset; if none is provided, it is set to ``name``
    - ``category`` -- (default: ``None``) to specify the categeory; if ``None``,
      the category of sets (:class:`~sage.categories.sets_cat.Sets`) is used

    EXAMPLES:

    A subset of a manifold::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: from sage.manifolds.subset import TopologicalManifoldSubset
        sage: A = TopologicalManifoldSubset(M, 'A', latex_name=r'\mathcal{A}')
        sage: A
        Subset A of the 2-dimensional topological manifold M
        sage: latex(A)
        \mathcal{A}
        sage: A.is_subset(M)
        True

    Instead of importing :class:`TopologicalManifoldSubset` in the global
    namespace, it is recommended to use the method :meth:`subset` to create a
    new subset::

        sage: B = M.subset('B', latex_name=r'\mathcal{B}'); B
        Subset B of the 2-dimensional topological manifold M
        sage: M.list_of_subsets()
        [Subset A of the 2-dimensional topological manifold M,
         Subset B of the 2-dimensional topological manifold M,
         2-dimensional topological manifold M]

    The manifold is itself a subset::

        sage: isinstance(M, TopologicalManifoldSubset)
        True

    Instances of :class:`TopologicalManifoldSubset` are Sage's facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`):
    their elements are manifold points
    (class :class:`~sage.manifolds.point.TopologicalManifoldPoint`),
    which have the manifold (and not the subset) as parent::

        sage: isinstance(A, Parent)
        True
        sage: A.category()
        Category of facade sets
        sage: A.facade_for()
        (2-dimensional topological manifold M,)
        sage: p = A.an_element(); p
        Point on the 2-dimensional topological manifold M
        sage: p.parent()
        2-dimensional topological manifold M
        sage: p in A
        True
        sage: p in M
        True

    """

    Element = TopologicalManifoldPoint

    def __init__(self, manifold, name, latex_name=None, category=None):
        r"""
        Construct a manifold subset.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.subset('A'); A
            Subset A of the 2-dimensional topological manifold M

        """
        if category is None:
            base = None
            category = Sets()
        else:
            base = manifold._field
        # Except for the manifold itself, the subsets are facade sets:
        if self is manifold:
            Parent.__init__(self, base=base, category=category)
        else:
            Parent.__init__(self, base=base, category=category,
                            facade=manifold)
        self._manifold = manifold
        if not isinstance(name, str):
            raise TypeError("{} is not a string".format(name))
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            if not isinstance(latex_name, str):
                raise TypeError("{} is not a string".format(latex_name))
            self._latex_name = latex_name
        if self is not self._manifold:
            for dom in manifold._subsets:
                if name == dom._name:
                    raise ValueError("the name '" + name +
                                     "' is already used for another " +
                                     "subset of the {}".format(manifold))
            manifold._subsets.add(self)
        self._supersets = set([manifold, self]) # subsets containing self
        self._subsets = set([self]) # subsets of self
        self._top_subsets = set([self]) # subsets contained in self but not
                                        # in another strict subset of self
        self._intersections = {} # dict. of intersections with other subsets
                                 # (key: subset name)
        self._unions = {} # dict. of unions with other subsets (key: subset
                          # name)
        self._open_covers = [] # list of open covers of self
        self._is_open = False  # a priori (may be redifined by subclasses)

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

        - an instance of :class:`~sage.manifolds.point.TopologicalManifoldPoint`
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
        if isinstance(coords, TopologicalManifoldPoint):
            point = coords # for readability
            if point._subset is self:
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

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: A._repr_()
            'Subset A of the 2-dimensional topological manifold M'
            sage: repr(A)  # indirect doctest
            'Subset A of the 2-dimensional topological manifold M'

        """
        return "Subset {} of the {}".format(self._name, self._manifold)

    def _latex_(self):
        r"""
        LaTeX representation of the object.

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

        """
        return self._latex_name

    def manifold(self):
        r"""
        Return the manifold of which the current object is a subset.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: A.manifold()
            2-dimensional topological manifold M
            sage: B = A.subset('B')
            sage: B.manifold()
            2-dimensional topological manifold M
            sage: M.manifold() is M
            True

        """
        return self._manifold

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
        return self._open_covers

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
            <type 'set'>
            sage: U in M.subsets()
            True

        The method :meth:`list_of_subsets` returns a list (sorted
        alphabetically by the subset names) instead of a set::

            sage: M.list_of_subsets()
            [2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M,
             Subset V of the 2-dimensional topological manifold M]

        """
        return self._subsets

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
        return sorted(self._subsets, key = lambda x: x._name)

    def subset(self, name, latex_name=None, is_open=False):
        r"""
        Create a subset of the current subset.

        INPUT:

        - ``name`` -- name given to the subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``is_open`` -- (default: ``False``) if ``True``, the created subset
          is assumed to be open with respect to the manifold's topology

        OUTPUT:

        - the subset, as an instance of :class:`TopologicalManifoldSubset`, or
          of the derived class
          :class:`~sage.manifolds.manifold.TopologicalManifold` if ``is_open``
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
        res = TopologicalManifoldSubset(self._manifold, name,
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
        :meth:`~sage.manifolds.manifold.TopologicalManifold.open_subset`)

        INPUT:

        - ``name`` -- (string) name of the subset

        OUTPUT:

        - instance of :class:`TopologicalManifoldSubset` (or
          of the derived class
          :class:`~sage.manifolds.manifold.TopologicalManifold` for an open
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

    def superset(self, name, latex_name=None, is_open=False):
        r"""
        Create a superset of the current subset.

        A *superset* is a manifold subset in which the current subset is
        included.

        INPUT:

        - ``name`` -- name given to the superset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          superset; if none is provided, it is set to ``name``
        - ``is_open`` -- (default: ``False``) if ``True``, the created subset
          is assumed to be open with respect to the manifold's topology

        OUTPUT:

        - the superset, as an instance of :class:`TopologicalManifoldSubset` or
          of the derived class
          :class:`~sage.manifolds.manifold.TopologicalManifold` if ``is_open``
          is ``True``.

        EXAMPLES:

        Creating some superset of a given subset::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A')
            sage: b = a.superset('B'); b
            Subset B of the 2-dimensional topological manifold M
            sage: b.list_of_subsets()
            [Subset A of the 2-dimensional topological manifold M,
             Subset B of the 2-dimensional topological manifold M]
            sage: a._supersets # random (set output)
            {Subset B of the 2-dimensional topological manifold M,
             Subset A of the 2-dimensional topological manifold M,
             2-dimensional topological manifold M}

        The superset of the whole manifold is itself::

            sage: M.superset('SM') is M
            True

        Two supersets of a given subset are a priori different::

            sage: c = a.superset('C')
            sage: c == b
            False

        """
        if self is self._manifold:
            return self
        if is_open:
            res = self._manifold.open_subset(name, latex_name=latex_name)
        else:
            res = TopologicalManifoldSubset(self._manifold, name,
                                    latex_name=latex_name)
        res._subsets.update(self._subsets)
        for sd in self._subsets:
            sd._supersets.add(res)
        if is_open and self._is_open:
            res._atlas = list(self._atlas)
            res._top_charts = list(self._top_charts)
            res._coord_changes = dict(self._coord_changes)
            res._def_chart = self._def_chart
        return res

    def intersection(self, other, name=None, latex_name=None):
        r"""
        Return the intersection of the current subset with another subset.

        INPUT:

        - ``other`` -- another subset of the same manifold
        - ``name`` -- (default: ``None``) name given to the intersection in the
          case the latter has to be created; the default is
          ``self._name`` inter ``other._name``
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          intersection in the case the latter has to be created; the default
          is built upon the symbol `\cap`

        OUTPUT:

        - instance of :class:`TopologicalManifoldSubset` representing the
          subset that is the intersection of the current subset with ``other``

        EXAMPLES:

        Intersection of two subsets::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A')
            sage: b = M.subset('B')
            sage: c = a.intersection(b); c
            Subset A_inter_B of the 2-dimensional topological manifold M
            sage: a.list_of_subsets()
            [Subset A of the 2-dimensional topological manifold M,
             Subset A_inter_B of the 2-dimensional topological manifold M]
            sage: b.list_of_subsets()
            [Subset A_inter_B of the 2-dimensional topological manifold M,
             Subset B of the 2-dimensional topological manifold M]
            sage: c._supersets  # random (set output)
            {Subset B of the 2-dimensional topological manifold M,
             Subset A_inter_B of the 2-dimensional topological manifold M,
             Subset A of the 2-dimensional topological manifold M,
             2-dimensional topological manifold M}

        Some checks::

            sage: (a.intersection(b)).is_subset(a)
            True
            sage: (a.intersection(b)).is_subset(a)
            True
            sage: a.intersection(b) is b.intersection(a)
            True
            sage: a.intersection(a.intersection(b)) is a.intersection(b)
            True
            sage: (a.intersection(b)).intersection(a) is a.intersection(b)
            True
            sage: M.intersection(a) is a
            True
            sage: a.intersection(M) is a
            True

        """
        if other._manifold != self._manifold:
            raise ValueError(
                          "the two subsets do not belong to the same manifold")
        # Particular cases:
        if self is self._manifold:
            return other
        if other is self._manifold:
            return self
        if self in other._subsets:
            return self
        if other in self._subsets:
            return other
        # Generic case:
        if other._name in self._intersections:
            # the intersection has already been created:
            return self._intersections[other._name]
        else:
            # the intersection must be created:
            if latex_name is None:
                if name is None:
                    latex_name = self._latex_name + r'\cap ' + other._latex_name
                else:
                    latex_name = name
            if name is None:
                name = self._name + "_inter_" + other._name
            if self._is_open and other._is_open:
                res = self.open_subset(name, latex_name=latex_name)
            else:
                res = self.subset(name, latex_name=latex_name)
            res._supersets.update(other._supersets)
            for sd in other._supersets:
                sd._subsets.add(res)
            other._top_subsets.add(res)
            self._intersections[other._name] = res
            other._intersections[self._name] = res
            return res

    def union(self, other, name=None, latex_name=None):
        r"""
        Return the union of the current subset with another subset.

        INPUT:

        - ``other`` -- another subset of the same manifold
        - ``name`` -- (default: ``None``) name given to the union in the
          case the latter has to be created; the default is
          ``self._name`` union ``other._name``
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          union in the case the latter has to be created; the default
          is built upon the symbol `\cup`

        OUTPUT:

        - instance of :class:`TopologicalManifoldSubset` representing the
          subset that is the union of the current subset with ``other``

        EXAMPLES:

        Union of two subsets::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A')
            sage: b = M.subset('B')
            sage: c = a.union(b); c
            Subset A_union_B of the 2-dimensional topological manifold M
            sage: a._supersets  # random (set output)
            set([subset 'A_union_B' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])
            sage: b._supersets  # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A_union_B' of the 2-dimensional manifold 'M'])
            sage: c._subsets  # random (set output)
            set([subset 'A_union_B' of the 2-dimensional manifold 'M',
                subset 'A' of the 2-dimensional manifold 'M',
                subset 'B' of the 2-dimensional manifold 'M'])

        Some checks::

            sage: a.is_subset(a.union(b))
            True
            sage: b.is_subset(a.union(b))
            True
            sage: a.union(b) is b.union(a)
            True
            sage: a.union(a.union(b)) is a.union(b)
            True
            sage: (a.union(b)).union(a) is a.union(b)
            True
            sage: a.union(M) is M
            True
            sage: M.union(a) is M
            True

        """
        if other._manifold != self._manifold:
            raise ValueError(
                          "the two subsets do not belong to the same manifold")
        # Particular cases:
        if (self is self._manifold) or (other is self._manifold):
            return self._manifold
        if self in other._subsets:
            return other
        if other in self._subsets:
            return self
        # Generic case:
        if other._name in self._unions:
            # the union has already been created:
            return self._unions[other._name]
        else:
            # the union must be created:
            if latex_name is None:
                if name is None:
                    latex_name = self._latex_name + r'\cup ' + other._latex_name
                else:
                    latex_name = name
            if name is None:
                name = self._name + "_union_" + other._name
            res_open = self._is_open and other._is_open
            res = self.superset(name, latex_name, is_open=res_open)
            res._subsets.update(other._subsets)
            res._top_subsets.add(self)
            res._top_subsets.add(other)
            for sd in other._subsets:
                sd._supersets.add(res)
            if res._is_open:
                for chart in other._atlas:
                    if chart not in res._atlas:
                        res._atlas.append(chart)
                for chart in other._top_charts:
                    if chart not in res._top_charts:
                        res._top_charts.append(chart)
                res._coord_changes.update(other._coord_changes)
            self._unions[other._name] = res
            other._unions[self._name] = res
            # Open covers of the union:
            for oc1 in self._open_covers:
                for oc2 in other._open_covers:
                    oc = oc1[:]
                    for s in oc2:
                        if s not in oc:
                            oc.append(s)
                res._open_covers.append(oc)
            return res

    def declare_union(self, dom1, dom2):
        r"""
        Declare that the current subset is the union of two subsets, i.e.
        that

        .. MATH::

            U = U_1 \cup U_2

        where `U` is the current subset,  `U_1\subset U` and `U_2\subset U`.

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
                oc = oc1[:]
                for s in oc2:
                    if s not in oc:
                        oc.append(s)
            self._open_covers.append(oc)

    def is_subset(self, other):
        r"""
        Return ``True`` iff the current subset is included in ``other``.

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
        Check whether a point is contained in the current subset.

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
        if point._subset is self:
            return True
        if point._subset.is_subset(self):
            return True
        #!# should be improved once coordinate definition have been introduced
        # in TopologicalManifoldSubset
        return False

    def point(self, coords=None, chart=None, name=None, latex_name=None):
        r"""
        Define a point in the subset.

        See :class:`~sage.manifolds.point.TopologicalManifoldPoint` for a
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
          :class:`~sage.manifolds.point.TopologicalManifoldPoint`.

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
        return self.element_class(self, coords=coords, chart=chart, name=name,
                                  latex_name=latex_name)
