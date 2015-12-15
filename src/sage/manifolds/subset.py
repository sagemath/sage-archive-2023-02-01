r"""
Subsets of topological manifolds

The class :class:`ManifoldSubset` implements generic subsets of a
topological manifold. Open subsets are implemented by the class
:class:`~sage.manifolds.subset.OpenTopologicalSubmanifold`, which inherits from
both :class:`ManifoldSubset` and
:class:`~sage.manifolds.manifold.TopologicalManifold`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version
- Travis Scrimshaw (2015): inheritance from
  :class:`~sage.manifolds.abstract.AbstractSet`

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
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.sets_cat import Sets
from sage.manifolds.abstract import AbstractSet
from sage.manifolds.manifold import TopologicalManifold

class ManifoldSubset(AbstractSet):
    r"""
    Subset of a topological manifold.

    The class :class:`ManifoldSubset` inherits (via
    :class:`~sage.manifolds.abstract.AbstractSet`) from the generic
    class :class:`~sage.structure.parent.Parent`.
    The corresponding element class is
    :class:`~sage.manifolds.point.ManifoldPoint`.

    Note that open subsets are not implemented directly by this class, but
    by the derived class
    :class:`~sage.manifolds.subset.OpenTopologicalSubmanifold`.

    INPUT:

    - ``manifold`` -- topological manifold on which the subset is defined
    - ``name`` -- string; name (symbol) given to the subset
    - ``latex_name`` --  (default: ``None``) string: LaTeX symbol to denote the
      subset; if none is provided, it is set to ``name``

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
    namespace, it is recommended to use the method
    :meth:`~sage.manifolds.abstract.AbstractSet.subset` to create a new subset::

        sage: B = M.subset('B', latex_name=r'\mathcal{B}'); B
        Subset B of the 2-dimensional topological manifold M
        sage: M.list_of_subsets()
        [Subset A of the 2-dimensional topological manifold M,
         Subset B of the 2-dimensional topological manifold M,
         2-dimensional topological manifold M]

    Instances of :class:`ManifoldSubset` are parents::

        sage: isinstance(A, Parent)
        True
        sage: A.category()
        Category of subobjects of sets
        sage: p = A.an_element(); p
        Point on the 2-dimensional topological manifold M
        sage: p.parent()
        Subset A of the 2-dimensional topological manifold M
        sage: p in A
        True
        sage: p in M
        True
    """

    def __init__(self, manifold, name, latex_name=None):
        r"""
        Construct a manifold subset.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A'); A
            Subset A of the 2-dimensional topological manifold M
            sage: type(A)
            <class 'sage.manifolds.subset.ManifoldSubset_with_category'>
            sage: A.category()
            Category of subobjects of sets
            sage: TestSuite(A).run(skip=['_test_elements', \
                                         '_test_not_implemented_methods'])

        .. TODO::

            implement method ``lift`` so that ``_test_not_implemented_methods``
            is passed.
            NB: ``_test_elements`` cannot be passed without a proper
            coordinate definition of the subset.

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
            sage: A.manifold() is M
            True
            sage: B = A.subset('B')
            sage: B.manifold() is M
            True

        An alias is ``ambient``::

            sage: A.ambient() is A.manifold()
            True

        """
        return self._manifold

    ambient = manifold

    def is_open(self):
        """
        Return if ``self`` is an open set.

        This method always returns ``False``, since open subsets must be
        constructed as instances of the subclass
        :class:`~sage.manifolds.subset.OpenTopologicalSubmanifold`
        (which redefines ``is_open``)

        EXAMPLE::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: A.is_open()
            False

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
          :class:`~sage.manifolds.subset.OpenTopologicalSubmanifold` if ``is_open``
          is ``True``.

        EXAMPLES:

        Creating some superset of a given subset::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = A.superset('B'); B
            Subset B of the 2-dimensional topological manifold M
            sage: B.list_of_subsets()
            [Subset A of the 2-dimensional topological manifold M,
             Subset B of the 2-dimensional topological manifold M]
            sage: A._supersets # random (set output)
            {Subset B of the 2-dimensional topological manifold M,
             Subset A of the 2-dimensional topological manifold M,
             2-dimensional topological manifold M}

        Two supersets of a given subset are a priori different::

            sage: C = A.superset('C')
            sage: C == B
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

#*****************************************************************************

class OpenTopologicalSubmanifold(ManifoldSubset, TopologicalManifold):
    """
    An open submanifold of a topological manifold, which is any open subset
    of the manifold.

    EXAMPLE:

    The unit ball of the Euclidean 2-plane::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart()
        sage: B = M.open_subset('B', coord_def={X: x^2+y^2<1})
        sage: B
        Open subset B of the 2-dimensional topological manifold M
        sage: type(B)
        <class 'sage.manifolds.subset.OpenTopologicalSubmanifold_with_category'>
        sage: B.category()
        Join of Category of subobjects of sets and Category of manifolds over
         Real Field with 53 bits of precision

    An open subset of a topological manifold being a topological manifold
    by itself, ``B`` inherits of all the methods of class
    :class:`~sage.manifolds.manifold.TopologicalManifold`::

        sage: isinstance(B, sage.manifolds.manifold.TopologicalManifold)
        True
        sage: dim(B)
        2
        sage: B.base_field()
        Real Field with 53 bits of precision

    Given the definition of ``B``, it is automatically endowed with a chart,
    which is the restriction of ``X`` to ``B``::

        sage: B.atlas()
        [Chart (B, (x, y))]
        sage: B.default_chart()
        Chart (B, (x, y))
        sage: B.default_chart() is X.restrict(B)
        True

    An point in ``B``::

        sage: p = B.an_element(); p
        Point on the 2-dimensional topological manifold M
        sage: X(p)  # the coordinates (x,y) of p
        (0, 0)
        sage: p in B
        True

    Checking whether various points, defined by their coordinates w.r.t.
    chart ``X``,  are in ``B``::

        sage: M((0,1/2)) in B
        True
        sage: M((0,1)) in B
        False
        sage: M((1/2,1)) in B
        False
        sage: M((-1/2,1/3)) in B
        True

    """
    def __init__(self, ambient, name, latex_name=None, category=None):
        """
        Initialize ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.open_subset('A', coord_def={X: x^2+y^2<1}); A
            Open subset A of the 2-dimensional topological manifold M
            sage: type(A)
            <class 'sage.manifolds.subset.OpenTopologicalSubmanifold_with_category'>
            sage: A.category()
            Join of Category of subobjects of sets and Category of manifolds
             over Real Field with 53 bits of precision
            sage: TestSuite(A).run(skip='_test_not_implemented_methods')

        .. TODO::

            implement method ``lift`` so that ``_test_not_implemented_methods``
            is passed.

        """
        # TODO: This test probably needs a better place, like open_subset
        if not isinstance(ambient, TopologicalManifold):
            raise TypeError("the argument 'ambient' must be " +
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
        TopologicalManifold.__init__(self, ambient.dim(), name, latex_name,
                                     ambient._field, ambient._structure,
                                     ambient._sindex, category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TEST::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: A = M.open_subset('A')
            sage: A._repr_()
            'Open subset A of the 3-dimensional topological manifold M'

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
          :class:`~sage.manifolds.subset.OpenTopologicalSubmanifold` if ``is_open``
          is ``True``.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: A = M.open_subset('A')
            sage: B = A.superset('B'); B
            Subset B of the 3-dimensional topological manifold M
            sage: B.list_of_subsets()
            [Open subset A of the 3-dimensional topological manifold M,
             Subset B of the 3-dimensional topological manifold M]
            sage: A._supersets  # random (set output)
            {Subset B of the 3-dimensional topological manifold M,
             Open subset A of the 3-dimensional topological manifold M,
             3-dimensional topological manifold M}

        Creating an open superset::

            sage: C = A.superset('C', is_open=True); C
            Open subset C of the 3-dimensional topological manifold M
            sage: C.list_of_subsets()
            [Open subset A of the 3-dimensional topological manifold M,
             Open subset C of the 3-dimensional topological manifold M]
            sage: A._supersets  # random (set output)
            {Subset B of the 3-dimensional topological manifold M,
             Open subset A of the 3-dimensional topological manifold M,
             Open subset C of the 3-dimensional topological manifold M,
             3-dimensional topological manifold M}

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

        In the present case (open submanifold), always return ``True``.

        EXAMPLE::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: A = M.open_subset('A')
            sage: A.is_open()
            True

        """
        return True
