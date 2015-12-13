r"""
Subsets of topological manifolds

The class :class:`ManifoldSubset` implements generic subsets of a
topological manifold. Open subsets are implemented by the class
:class:`~sage.manifolds.manifold.Manifold` (since an open subset of
a manifold is a manifold by itself), which inherits from
:class:`ManifoldSubset`.

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

from sage.categories.sets_cat import Sets
from sage.manifolds.abstract import AbstractSet
from sage.manifolds.manifold import Manifold

class ManifoldSubset(AbstractSet):
    r"""
    Subset of a topological manifold.

    The class :class:`ManifoldSubset` inherits from the generic Sage
    class :class:`~sage.structure.parent.Parent` and is declared to belong to
    the category of facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`).
    The corresponding element class is
    :class:`~sage.manifolds.point.ManifoldPoint`. A subset acts
    as a facade for the true parent of its points, which is the whole manifold
    (see example below).

    Note that open subsets are not implemented directly by this class, but
    by the derived class :class:`~sage.manifolds.manifold.Manifold`
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
        sage: from sage.manifolds.subset import ManifoldSubset
        sage: A = ManifoldSubset(M, 'A', latex_name=r'\mathcal{A}')
        sage: A
        Subset A of the 2-dimensional topological manifold M
        sage: latex(A)
        \mathcal{A}
        sage: A.is_subset(M)
        True

    Instead of importing :class:`ManifoldSubset` in the global
    namespace, it is recommended to use the method :meth:`subset` to create a
    new subset::

        sage: B = M.subset('B', latex_name=r'\mathcal{B}'); B
        Subset B of the 2-dimensional topological manifold M
        sage: M.list_of_subsets()
        [Subset A of the 2-dimensional topological manifold M,
         Subset B of the 2-dimensional topological manifold M,
         2-dimensional topological manifold M]

    The manifold is itself a subset::

        sage: isinstance(M, ManifoldSubset)
        True

    Instances of :class:`ManifoldSubset` are Sage's facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`):
    their elements are manifold points
    (class :class:`~sage.manifolds.point.ManifoldPoint`),
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

    def __init__(self, manifold, name, latex_name=None, category=None):
        r"""
        Construct a manifold subset.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.subset('A'); A
            Subset A of the 2-dimensional topological manifold M
        """
        category = Sets().Subobjects()
        AbstractSet.__init__(self, name=name, latex_name=latex_name,
                                  category=category)

        self._manifold = manifold
        for dom in manifold._subsets:
            if name == dom._name:
                raise ValueError("the name '" + name +
                                 "' is already used for another " +
                                 "subset of the {}".format(manifold))
        manifold._subsets.add(self)

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

    ambient = manifold

    def is_open(self):
        """
        Return if ``self`` is an open set.
        """
        return False

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

        - the superset, as an instance of :class:`ManifoldSubset` or
          of the derived class
          :class:`~sage.manifolds.manifold.Manifold` if ``is_open``
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
        res = self._manifold.subset(name, latex_name, is_open)
        res._subsets.update(self._subsets)
        for sd in self._subsets:
            sd._supersets.add(res)
        return res

    def intersection(self, other, name=None, latex_name=None):
        r"""
        Return the intersection of the current subset with another subset.

        INPUT:

        - ``other`` -- another subset of the same manifold
        - ``name`` -- (default: ``None``) name given to the intersection in the
          case the latter has to be created; the default is
          ``self._name`` inter ``other._name``
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          intersection in the case the latter has to be created; the default
          is built upon the symbol `\cap`

        OUTPUT:

        - instance of :class:`ManifoldSubset` representing the
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
        # Particular cases:
        if other is self._manifold:
            return self
        if self in other._subsets:
            return self
        if other in self._subsets:
            return other
        if other._manifold != self._manifold:
            raise ValueError("the two subsets do not belong to the same manifold")

        # Generic case:
        if other._name in self._intersections:
            # the intersection has already been created:
            return self._intersections[other._name]

        # TODO: Check to see if we've already created a union of ``self`` and ``other``

        # the intersection must be created:
        if latex_name is None:
            if name is None:
                latex_name = self._latex_name + r'\cap ' + other._latex_name
            else:
                latex_name = name
        if name is None:
            name = self._name + "_inter_" + other._name
        if self.is_open() and other.is_open():
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

        - instance of :class:`ManifoldSubset` representing the
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
        # Particular cases:
        if other is self._manifold:
            return self._manifold
        if self in other._subsets:
            return other
        if other in self._subsets:
            return self
        if other.manifold() != self._manifold:
            raise ValueError("the two subsets do not belong to the same manifold")

        # Generic case:
        if other._name in self._unions:
            # the union has already been created:
            return self._unions[other._name]

        # TODO: Check to see if we've already created a union of ``self`` and ``other``

        # the union must be created:
        if latex_name is None:
            if name is None:
                latex_name = self._latex_name + r'\cup ' + other._latex_name
            else:
                latex_name = name
        if name is None:
            name = self._name + "_union_" + other._name
        res_open = self.is_open() and other.is_open()
        res = self.superset(name, latex_name, is_open=res_open)
        res._subsets.update(other._subsets)
        res._top_subsets.add(self)
        res._top_subsets.add(other)
        for sd in other._subsets:
            sd._supersets.add(res)
        if res_open:
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

class TopologicalSubmanifold(ManifoldSubset, Manifold):
    """
    A submanifold of a manifold, which is any open subset of a manifold.
    """
    def __init__(self, ambient, name, latex_name=None, category=None):
        """
        Initialize ``self``.
        """
        # TODO: This test probably needs a better place, like open_subset
        if not isinstance(ambient, Manifold):
            raise TypeError("the argument 'manifold' must be " +
                            "a topological manifold")

        # This is copied from ManifoldSubset to avoid twice
        #   initializing AbstractSet
        self._manifold = ambient
        for dom in ambient._subsets:
            if name == dom._name:
                raise ValueError("the name '" + name +
                                 "' is already used for another " +
                                 "subset of the {}".format(ambient))
        ambient._subsets.add(self)

        category = ambient.category().Subobjects()
        Manifold.__init__(self, ambient.dim(), name, latex_name,
                          ambient._field, ambient._structure,
                          ambient._sindex, category)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Open subset {} of the {}".format(self._name, self._manifold)


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

        - the superset, as an instance of :class:`ManifoldSubset` or
          of the derived class
          :class:`~sage.manifolds.manifold.Manifold` if ``is_open``
          is ``True``.

        EXAMPLES:
        """
        res = ManifoldSubset.superset(self, name, latex_name, is_open)
        if is_open:
            res._atlas = list(self._atlas)
            res._top_charts = list(self._top_charts)
            res._coord_changes = dict(self._coord_changes)
            res._def_chart = self._def_chart
        return res

    def is_open(self):
        """
        Return if ``self`` is an open set.
        """
        return True
