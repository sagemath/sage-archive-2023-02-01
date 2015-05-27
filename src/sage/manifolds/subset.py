r"""
Subsets of topological manifolds

The class :class:`TopManifoldSubset` implements generic subsets of a topological
manifold. The subclass :class:`TopManifoldOpenSubset` is devoted
to open subsets, with respect to the manifold topology.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version

REFERENCES:

- J.M. Lee : *Introduction to Topological Manifolds*, 2nd ed., Springer (New
  York) (2011)

EXAMPLES:

Two subsets on a manifold::

    sage: M = Manifold(2, 'M')
    sage: a = M.subset('A') ; a
    subset 'A' of the 2-dimensional manifold 'M'
    sage: b = M.subset('B') ; b
    subset 'B' of the 2-dimensional manifold 'M'
    sage: M.subsets()  # random (set output)
    {subset 'A' of the 2-dimensional manifold 'M',
     subset 'B' of the 2-dimensional manifold 'M',
     2-dimensional manifold 'M'}

The intersection of the two subsets::

    sage: c = a.intersection(b) ; c
    subset 'A_inter_B' of the 2-dimensional manifold 'M'

Their union::

    sage: d = a.union(b) ; d
    subset 'A_union_B' of the 2-dimensional manifold 'M'

State of various data members after the above operations::

    sage: M.subsets()  # random (set output)
    {subset 'A' of the 2-dimensional manifold 'M',
     subset 'B' of the 2-dimensional manifold 'M',
     subset 'A_inter_B' of the 2-dimensional manifold 'M',
     2-dimensional manifold 'M'}
    sage: a.subsets()  # random (set output)
    set([subset 'A' of the 2-dimensional manifold 'M',
         subset 'A_inter_B' of the 2-dimensional manifold 'M'])
    sage: a._supersets  # random (set output)
    set([subset 'A_union_B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M',
         subset 'A' of the 2-dimensional manifold 'M'])
    sage: c._supersets  # random (set output)
    set([subset 'B' of the 2-dimensional manifold 'M',
         subset 'A_union_B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M',
         subset 'A' of the 2-dimensional manifold 'M',
         subset 'A_inter_B' of the 2-dimensional manifold 'M'])
    sage: c.subsets()  # random (set output)
    set([subset 'A_inter_B' of the 2-dimensional manifold 'M'])
    sage: d.subsets()  # random (set output)
    set([subset 'B' of the 2-dimensional manifold 'M',
         subset 'A_union_B' of the 2-dimensional manifold 'M',
         subset 'A_inter_B' of the 2-dimensional manifold 'M',
         subset 'A' of the 2-dimensional manifold 'M'])
    sage: d._supersets  # random (set output)
    set([subset 'A_union_B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M'])

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
from sage.categories.homset import Hom
from sage.rings.infinity import Infinity
from sage.manifolds.point import TopManifoldPoint

class TopManifoldSubset(UniqueRepresentation, Parent):
    r"""
    Subset of a topological manifold.

    The class :class:`TopManifoldSubset` inherits from the generic Sage class
    :class:`~sage.structure.parent.Parent` and is declared to belong to
    the category of facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`).
    The corresponding element class is
    :class:`~sage.manifolds.point.TopManifoldPoint`. A subset acts
    as a facade for the true parent of its points, which is the whole manifold
    (see example below).

    For an open subset, use the class :class:`TopManifoldOpenSubset` instead.

    INPUT:

    - ``manifold`` -- topological manifold on which the subset is defined
    - ``name`` -- string; name (symbol) given to the subset
    - ``latex_name`` --  (default: ``None``) string: LaTeX symbol to denote the
      subset; if none is provided, it is set to ``name``

    EXAMPLES:

    A subset of a manifold::

        sage: Manifold._clear_cache_() # for doctests only
        sage: M = Manifold(2, 'M')
        sage: from sage.manifolds.subset import TopManifoldSubset
        sage: A = TopManifoldSubset(M, 'A', latex_name=r'\mathcal{A}') ; A
        subset 'A' of the 2-dimensional manifold 'M'
        sage: latex(A)
        \mathcal{A}
        sage: A.is_subset(M)
        True

    Instead of importing :class:`TopManifoldSubset` in the global namespace,
    it is recommended to use the method :meth:`subset` to create a new subset::

        sage: B = M.subset('B', latex_name=r'\mathcal{B}') ; B
        subset 'B' of the 2-dimensional manifold 'M'
        sage: M.subsets()  # random (set output)
        {subset 'A' of the 2-dimensional manifold 'M',
         subset 'B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M'}

    The manifold is itself a subset::

        sage: isinstance(M, TopManifoldSubset)
        True

    Actually, it is an instance of the subclass
    :class:`~sage.manifolds.subset.TopManifoldOpenSubset`, for it is
    (by definition) open::

        sage: isinstance(M, sage.manifolds.subset.TopManifoldOpenSubset)
        True

    Instances of :class:`TopManifoldSubset` are Sage's facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`):
    their elements are manifold points
    (class :class:`~sage.manifolds.point.TopManifoldPoint`),
    which have the manifold (and not the subset) as parent::

        sage: isinstance(A, Parent)
        True
        sage: A.category()
        Category of facade sets
        sage: A.facade_for()
        (2-dimensional manifold 'M',)
        sage: p = A.an_element() ; p
        point on 2-dimensional manifold 'M'
        sage: p.parent()
        2-dimensional manifold 'M'
        sage: p in A
        True
        sage: p in M
        True

    """

    Element = TopManifoldPoint

    def __init__(self, manifold, name, latex_name=None):
        r"""
        Construct a manifold subset.

        TESTS::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: A = M.subset('A') ; A
            subset 'A' of the 2-dimensional manifold 'M'

        """
        if not isinstance(name, str):
            raise TypeError("{} is not a string".format(name))
        # Except for the manifold itself, the subsets are facade sets:
        if self is manifold:
            Parent.__init__(self, category=Sets())
        else:
            Parent.__init__(self, category=Sets(), facade=manifold)
            for dom in manifold._subsets:
                if name == dom._name:
                    raise ValueError("The name '" + name +
                                     "' is already used for another " +
                                     "subset of the {}".format(manifold))
            manifold._subsets.add(self)
        self._manifold = manifold
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            if not isinstance(latex_name, str):
                raise TypeError("{} is not a string".format(latex_name))
            self._latex_name = latex_name
        self._supersets = set([manifold, self]) # subsets containing self
        self._subsets = set([self]) # subsets of self
        self._top_subsets = set([self]) # subsets contained in self but not
                                        # in another strict subset of self
        self._intersections = {} # dict. of intersections with other subsets
                                 # (key: subset name)
        self._unions = {} # dict. of unions with other subsets (key: subset
                          # name)

    #### Methods required for any Parent in the category of sets:
    def _element_constructor_(self, coords=None, chart=None, name=None,
                              latex_name=None, check_coords=True):
        r"""
        Construct a point on ``self`` from its coordinates in some chart.
        """
        if isinstance(coords, TopManifoldPoint):
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
        Construct some point in ``self``.

        EXAMPLES::

        """
        #!# should be improved...
        return self.element_class(self)

    #### End of methods required for any Parent in the category of sets

    def _repr_(self):
        r"""
        String representation of ``self``.
        """
        return "Subset {} of the {}".format(self._name, self._manifold)

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.
        """
        return self._latex_name

    def manifold(self):
        r"""
        Return the manifold of which ``self`` is a subset.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: A = M.subset('A')
            sage: A.manifold()
            2-dimensional manifold 'M'
            sage: B = A.subset('B')
            sage: B.manifold()
            2-dimensional manifold 'M'
            sage: M.manifold() is M
            True

        """
        return self._manifold

    def subsets(self):
        r"""
        Return the set of subsets that have been defined on ``self``.

        OUTPUT:

        - A Python set containing all the subsets that have been defined on
          ``self``.

        .. NOTE::

            To get the subsets as a list, used the method
            :meth:`list_of_subsets` instead.

        EXAMPLE:

        Subsets of a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: M.subsets()  # random (set output)
             {open subset 'U' of the 2-dimensional manifold 'M',
              subset 'V' of the 2-dimensional manifold 'M',
              2-dimensional manifold 'M'}
            sage: U in M.subsets()
            True

        The method :meth:`list_of_subsets` returns a list instead of a set::

            sage: M.list_of_subsets()
            [2-dimensional manifold 'M',
             open subset 'U' of the 2-dimensional manifold 'M',
             subset 'V' of the 2-dimensional manifold 'M']

        """
        return self._subsets

    def list_of_subsets(self):
        r"""
        Return the list of subsets that have been defined on ``self``.

        The list is sorted by the alphabetical names of the subsets.

        OUTPUT:

        - A list containing all the subsets that have been defined on
          ``self``.

        .. NOTE::

            To get the subsets as a Python set, used the method
            :meth:`subsets` instead.

        EXAMPLE:

        Subsets of a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: M.list_of_subsets()
            [2-dimensional manifold 'M',
             open subset 'U' of the 2-dimensional manifold 'M',
             subset 'V' of the 2-dimensional manifold 'M']

        The method :meth:`subsets` returns a set instead of a list::

            sage: M.subsets()  # random (set output)
            {open subset 'U' of the 2-dimensional manifold 'M',
             subset 'V' of the 2-dimensional manifold 'M',
             2-dimensional manifold 'M'}

        """
        return sorted(self._subsets, key = lambda x: x._name)

    def subset(self, name, latex_name=None, is_open=False):
        r"""
        Create a subset of ``self``.

        INPUT:

        - ``name`` -- name given to the subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``is_open`` -- (default: False) if True, the created subset is
          assumed to be open with respect to the manifold's topology

        OUTPUT:

        - the subset, as an instance of :class:`TopManifoldSubset`, or of
          :class:`TopManifoldOpenSubset` if ``is_open`` is True.

        EXAMPLES:

        Creating a subset of a manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: a = M.subset('A') ; a
            subset 'A' of the 2-dimensional manifold 'M'

        Creating a subset of A::

            sage: b = a.subset('B', latex_name=r'\mathcal{B}') ; b
            subset 'B' of the 2-dimensional manifold 'M'
            sage: latex(b)
            \mathcal{B}

        B is then a subset of A and A is a superset of B::

            sage: a._subsets  # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])
            sage: b._supersets  # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])

        Creating an open subset of A::

            sage: c = a.subset('C', is_open=True) ; c
            open subset 'C' of the 2-dimensional manifold 'M'

        """
        if is_open:
            res = TopManifoldOpenSubset(self._manifold, name, latex_name)
        else:
            res = TopManifoldSubset(self._manifold, name, latex_name)
        res._supersets.update(self._supersets)
        for sd in self._supersets:
            sd._subsets.add(res)
        self._top_subsets.add(res)
        return res

    def superset(self, name, latex_name=None, is_open=False):
        r"""
        Create a superset of ``self``.

        A *superset* is a manifold subset in which ``self`` is included.

        INPUT:

        - ``name`` -- name given to the superset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          superset; if none is provided, it is set to ``name``
        - ``is_open`` -- (default: False) if True, the created subset is
          assumed to be open with respect to the manifold's topology

        OUTPUT:

        - the superset, as an instance of :class:`TopManifoldSubset` or of
          :class:`TopManifoldOpenSubset` if ``is_open==True``.

        EXAMPLES:

        Creating some superset of a given subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: a = M.subset('A')
            sage: b = a.superset('B') ; b
            subset 'B' of the 2-dimensional manifold 'M'
            sage: b._subsets # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])
            sage: a._supersets # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])

        The superset of the manifold is itself::

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
            res = TopManifoldOpenSubset(self._manifold, name, latex_name)
        else:
            res = TopManifoldSubset(self._manifold, name, latex_name)
        res._subsets.update(self._subsets)
        for sd in self._subsets:
            sd._supersets.add(res)
        res._atlas = list(self._atlas)
        res._top_charts = list(self._top_charts)
        res._coord_changes = dict(self._coord_changes)
        res._frames = list(self._frames)
        res._top_frames = list(self._top_frames)
        res._frame_changes = dict(self._frame_changes)
        res._coframes = list(self._coframes)
        res._def_chart = self._def_chart
        res._def_frame = self._def_frame
        return res

    def intersection(self, other, name=None, latex_name=None):
        r"""
        Return the intersection of ``self`` with another subset.

        INPUT:

        - ``other`` -- another subset of the same manifold
        - ``name`` -- (default: ``None``) name given to the intersection in the
          case the latter has to be created; the default is
          ``self._name`` inter ``other._name``
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          intersection in the case the latter has to be created; the default
          is built upon the symbol `\cap`

        OUTPUT:

        - instance of :class:`TopManifoldSubset` (or :class:`TopManifoldOpenSubset`
          if both subsets are open) representing the subset that is the
          intersection of ``self`` with ``other``

        EXAMPLES:

        Intersection of two subsets::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: a = M.subset('A')
            sage: b = M.subset('B')
            sage: c = a.intersection(b) ; c
            subset 'A_inter_B' of the 2-dimensional manifold 'M'
            sage: a._subsets  # random (set output)
            set([subset 'A_inter_B' of the 2-dimensional manifold 'M',
                 subset 'A' of the 2-dimensional manifold 'M'])
            sage: b._subsets  # random (set output)
            set([subset 'B' of the 2-dimensional manifold 'M',
                 subset 'A_inter_B' of the 2-dimensional manifold 'M'])
            sage: c._supersets  # random (set output)
            set([subset 'A' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 subset 'A_inter_B' of the 2-dimensional manifold 'M',
                 subset 'B' of the 2-dimensional manifold 'M'])

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
            raise TypeError(
                "The two subsets do not belong to the same manifold.")
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
            if isinstance(self, TopManifoldOpenSubset) and isinstance(other, TopManifoldOpenSubset):
                res = self.open_subset(name, latex_name)
            else:
                res = self.subset(name, latex_name)
            res._supersets.update(other._supersets)
            for sd in other._supersets:
                sd._subsets.add(res)
            other._top_subsets.add(res)
            self._intersections[other._name] = res
            other._intersections[self._name] = res
            return res

    def union(self, other, name=None, latex_name=None):
        r"""
        Return the union of ``self`` with another subset.

        INPUT:

        - ``other`` -- another subset of the same manifold
        - ``name`` -- (default: ``None``) name given to the union in the
          case the latter has to be created; the default is
          ``self._name`` union ``other._name``
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          union in the case the latter has to be created; the default
          is built upon the symbol `\cup`

        OUTPUT:

        - instance of :class:`TopManifoldSubset` (or :class:`TopManifoldOpenSubset`
          if both subsets are open) representing the subset that is the
          union of ``self`` with ``other``

        EXAMPLES:

        Union of two subsets::

            sage: M = Manifold(2, 'M')
            sage: a = M.subset('A')
            sage: b = M.subset('B')
            sage: c = a.union(b) ; c
            subset 'A_union_B' of the 2-dimensional manifold 'M'
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
            raise TypeError(
                "The two subsets do not belong to the same manifold.")
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
            res_open = isinstance(self, TopManifoldOpenSubset) and \
                       isinstance(other, TopManifoldOpenSubset)
            res = self.superset(name, latex_name, is_open=res_open)
            res._subsets.update(other._subsets)
            res._top_subsets.add(self)
            res._top_subsets.add(other)
            for sd in other._subsets:
                sd._supersets.add(res)
            for chart in other._atlas:
                if chart not in res._atlas:
                    res._atlas.append(chart)
            for chart in other._top_charts:
                if chart not in res._top_charts:
                    res._top_charts.append(chart)
            res._coord_changes.update(other._coord_changes)
            for frame in other._frames:
                if frame not in res._frames:
                    res._frames.append(frame)
            for frame in other._top_frames:
                if frame not in res._top_frames:
                    res._top_frames.append(frame)
            res._frame_changes.update(other._frame_changes)
            for coframe in other._coframes:
                if coframe not in res._coframes:
                    res._coframes.append(coframe)
            self._unions[other._name] = res
            other._unions[self._name] = res
            return res

    def declare_union(self, dom1, dom2):
        r"""
        Declare that ``self`` is the union of two subsets, i.e.
        that

        .. MATH::

            U = U_1 \cup U_2

        where `U` is ``self``,  `U_1\subset U` and `U_2\subset U`.

        INPUT:

        - ``dom1`` -- subset `U_1`
        - ``dom2`` -- subset `U_2`

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: M.declare_union(A, B)
            sage: A.union(B)
            2-dimensional manifold 'M'

        """
        if not dom1.is_subset(self):
            raise TypeError("The " + str(dom1) + " is not a subset of " +
                            "the " + str(self) + ".")
        if not dom2.is_subset(self):
            raise TypeError("The " + str(dom2) + " is not a subset of " +
                            "the " + str(self) + ".")
        dom1._unions[dom2._name] = self
        dom2._unions[dom1._name] = self

    def is_subset(self, other):
        r"""
        Return ``True`` iff ``self`` is included in ``other``.

        EXAMPLES:

        Subsubsets on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
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
        Check whether a point is contained in ``self``.
        """
        # for efficiency, a quick test first:
        if point._subset is self:
            return True
        if point._subset.is_subset(self):
            return True
        #!# should be improved once coordinate definition have been introduced
        # in TopManifoldSubset
        return False

    def point(self, coords=None, chart=None, name=None, latex_name=None):
        r"""
        Define a point in ``self``.

        See :class:`~sage.manifolds.point.TopManifoldPoint` for a
        complete documentation.

        INPUT:

        - ``coords`` -- the point coordinates (as a tuple or a list) in the
          chart specified by ``chart``
        - ``chart`` -- (default: ``None``) chart in which the point coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the default chart of ``self``
        - ``name`` -- (default: ``None``) name given to the point
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the point;
          if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - the declared point, as an instance of
          :class:`~sage.manifolds.point.TopManifoldPoint`.

        EXAMPLES:

        Points on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: p = M.point((1,2), name='p') ; p
            point 'p' on 2-dimensional manifold 'M'
            sage: p in M
            True
            sage: a = M.open_subset('A')
            sage: c_uv.<u,v> = a.chart()
            sage: q = a.point((-1,0), name='q') ; q
            point 'q' on 2-dimensional manifold 'M'
            sage: q in a
            True
            sage: p._coordinates
            {chart (M, (x, y)): (1, 2)}
            sage: q._coordinates
            {chart (A, (u, v)): (-1, 0)}

        """
        return self.element_class(self, coords=coords, chart=chart, name=name,
                                  latex_name=latex_name)


#******************************************************************************

class TopManifoldOpenSubset(TopManifoldSubset):
    r"""
    Open subset of a differentiable manifold over `\RR`.

    The class :class:`TopManifoldOpenSubset` inherits naturally from the class
    :class:`TopManifoldSubset`. Via the latter, it inherits from the generic
    Sage class :class:`~sage.structure.parent.Parent` and is declared to belong
    to the category of facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`).
    The corresponding element class is
    :class:`~sage.manifolds.point.TopManifoldPoint`. An open subset acts
    as a facade for the true parent of its points, which is the whole manifold
    (see example below).

    INPUT:

    - ``manifold`` -- manifold on which the open subset is defined
    - ``name`` -- string; name (symbol) given to the open subset
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the open
      subset; if none is provided, it is set to ``name``

    EXAMPLES:

    An open subset of a 2-dimensional manifold::

        sage: Manifold._clear_cache_() # for doctests only
        sage: M = Manifold(2, 'M')
        sage: from sage.manifolds.subset import TopManifoldOpenSubset
        sage: A = TopManifoldOpenSubset(M, 'A', latex_name=r'\mathcal{A}') ; A
        open subset 'A' of the 2-dimensional manifold 'M'
        sage: latex(A)
        \mathcal{A}

    Instead of importing :class:`TopManifoldOpenSubset` in the global namespace,
    it is recommended to use the method :meth:`open_subset` to create a new
    open subset::

        sage: B = M.open_subset('B', latex_name=r'\mathcal{B}') ; B
        open subset 'B' of the 2-dimensional manifold 'M'
        sage: M.subsets()  # random (set output)
        {open subset 'A' of the 2-dimensional manifold 'M',
         open subset 'B' of the 2-dimensional manifold 'M',
         2-dimensional manifold 'M'}

    The manifold is itself an open subset (by definition!)::

        sage: isinstance(M, TopManifoldOpenSubset)
        True

    Instances of :class:`TopManifoldOpenSubset` are Sage's facade sets
    (see :meth:`~sage.categories.sets_cat.Sets.SubcategoryMethods.Facade`):
    their elements are manifold points
    (class :class:`~sage.manifolds.point.TopManifoldPoint`),
    which have the manifold (and not the open subset) as parent::

        sage: A.category()
        Category of facade sets
        sage: A.facade_for()
        (2-dimensional manifold 'M',)
        sage: p = A.an_element() ; p
        point on 2-dimensional manifold 'M'
        sage: p.parent()
        2-dimensional manifold 'M'

    Points can be created by providing their coordinates in some
    chart via the operator () applied to the subset::

        sage: X.<x,y> = A.chart()
        sage: p = A((-2,3)) ; p
        point on 2-dimensional manifold 'M'
        sage: p.coord()
        (-2, 3)

    Other arguments can be specified::

        sage: p = A((-2,3), chart=X, name='p') ; p
        point 'p' on 2-dimensional manifold 'M'

    It is equivalent to use the method
    :meth:`~sage.manifolds.subset.TopManifoldSubset.point`::

        sage: A((-2,3)) == A.point((-2,3))
        True

    An open subset can be defined by some coordinate conditions::

        sage: U = A.open_subset('U', coord_def={X: x<0}) ; U
        open subset 'U' of the 2-dimensional manifold 'M'
        sage: U.is_subset(A)
        True
        sage: p in U  # since p's coordinates in chart X are (x,y)=(-2,3)
        True

    The coordinate conditions are taken into account when asking for a
    generic element::

        sage: q = U.an_element() ; q
        point on 2-dimensional manifold 'M'
        sage: q in U
        True
        sage: q.coord() # coordinates in U's default chart (x,y)
        (-1, 0)

    An open subset can be called to construct a point from another one,
    provided that the latter belongs to the open subset::

        sage: p1 = U(p) ; p1
        point 'p' on 2-dimensional manifold 'M'

    The two points simply differ by the returned value of
    :meth:`~sage.manifolds.point.TopManifoldPoint.containing_set`::

        sage: p1 == p
        True
        sage: p.containing_set()
        open subset 'A' of the 2-dimensional manifold 'M'
        sage: p1.containing_set()
        open subset 'U' of the 2-dimensional manifold 'M'

    """
    def __init__(self, manifold, name, latex_name=None):
        r"""
        Construct an open subset of a manifold.

        TESTS::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x^2+y^2<1}) ; U
            open subset 'U' of the 2-dimensional manifold 'M'
            sage: U.category()
            Category of facade sets
            sage: TestSuite(U).run()

        """
        #*# from sage.manifolds.scalarfield_algebra import ScalarFieldAlgebra
        TopManifoldSubset.__init__(self, manifold, name, latex_name)
        self._atlas = []  # list of charts defined on subsets of self
        self._top_charts = []  # list of charts defined on subsets of self
                        # that are not subcharts of charts on larger subsets
        self._def_chart = None  # default chart
        self._coord_changes = {} # dictionary of transition maps
        # list of charts that individually cover self, i.e. whose
        # domains are self (if non-empty, self is a coordinate domain):
        self._covering_charts = []
        # algebra of scalar fields defined on self:
        #*# self._scalar_field_algebra = ScalarFieldAlgebra(self)
        # the zero scalar field:
        #*# self._zero_scalar_field = self._scalar_field_algebra.zero()
        # the identity map on self:
        #*# self._identity_map = Hom(self, self).one()

    def _an_element_(self):
        r"""
        Construct some (unamed) point on ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M._an_element_() ; p
            point on 2-dimensional manifold 'M'
            sage: p.coord()
            (0, 0)
            sage: U = M.open_subset('U', coord_def={X: y>1}) ; U
            open subset 'U' of the 2-dimensional manifold 'M'
            sage: p = U._an_element_() ; p
            point on 2-dimensional manifold 'M'
            sage: p in U
            True
            sage: p.coord()
            (0, 2)
            sage: V = U.open_subset('V', coord_def={X.restrict(U): x<-pi})
            sage: p = V._an_element_() ; p
            point on 2-dimensional manifold 'M'
            sage: p in V
            True
            sage: p.coord()
            (-pi - 1, 2)

        """
        if self._def_chart is None:
            return self.element_class(self)
        # Attempt to construct a point in the domain of the default chart
        chart = self._def_chart
        if self._manifold.base_field() == 'real':
            coords = []
            for coord_range in chart._bounds:
                xmin = coord_range[0][0]
                xmax = coord_range[1][0]
                if xmin == -Infinity:
                    if xmax == Infinity:
                        x = 0
                    else:
                        x = xmax - 1
                else:
                    if xmax == Infinity:
                        x = xmin + 1
                    else:
                        x = (xmin + xmax)/2
                coords.append(x)
        else:
            coords = self._manifold.dim()*[0]
        if not chart.valid_coordinates(*coords):
            # Attempt to construct a point in the domain of other charts
            if self._manifold.base_field() == 'real':
                for ch in self._atlas:
                    if ch is self._def_chart:
                        continue # since this case has already been attempted
                    coords = []
                    for coord_range in ch._bounds:
                        xmin = coord_range[0][0]
                        xmax = coord_range[1][0]
                        if xmin == -Infinity:
                            if xmax == Infinity:
                                x = 0
                            else:
                                x = xmax - 1
                        else:
                            if xmax == Infinity:
                                x = xmin + 1
                            else:
                                x = (xmin + xmax)/2
                        coords.append(x)
                    if ch.valid_coordinates(*coords):
                        chart = ch
                        break
                else:
                    # A generic element with specific coordinates could not be
                    # automatically generated, due to too complex cooordinate
                    # conditions. An element without any coordinate set is
                    # returned instead:
                    return self.element_class(self)
            else:
                # Case of manifolds over a field different from R
                for ch in self._atlas:
                    if ch is self._def_chart:
                        continue # since this case has already been attempted
                    if ch.valid_coordinates(*coords):
                        chart = ch
                        break
                else:
                    return self.element_class(self)
        # The point is constructed with check_coords=False since the check
        # has just been performed above:
        return self.element_class(self, coords=coords, chart=chart,
                                  check_coords=False)


    def _repr_(self):
        r"""
        String representation of ``self``.
        """
        return "Open subset {} of the {}".format(self._name, self._manifold)

    def __contains__(self, point):
        r"""
        Check whether a point is contained in ``self``.
        """
        # for efficiency, a quick test first:
        if point._subset is self:
            return True
        if point._subset.is_subset(self):
            return True
        for chart in self._atlas:
            if chart in point._coordinates:
                if chart.valid_coordinates( *(point._coordinates[chart]) ):
                    return True
        for chart in point._coordinates:
            for schart in chart._subcharts:
                if schart in self._atlas and schart.valid_coordinates(
                                          *(point._coordinates[chart]) ):
                    return True
        return False

    def atlas(self):
        r"""
        Return the atlas of ``self``.

        OUTPUT:

        - list of charts defined on open subsets of ``self``.

        EXAMPLES:

        Charts on subsets of `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: M.atlas()
            [chart (R^2, (x, y))]
            sage: U = M.open_subset('U', coord_def={c_cart: (y!=0,x<0)}) # U = R^2 \ half line {y=0,x>=0}
            sage: U.atlas()
            [chart (U, (x, y))]
            sage: M.atlas()
            [chart (R^2, (x, y)), chart (U, (x, y))]
            sage: c_spher.<r, ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical (polar) coordinates on U
            sage: U.atlas()
            [chart (U, (x, y)), chart (U, (r, ph))]
            sage: M.atlas()
            [chart (R^2, (x, y)), chart (U, (x, y)), chart (U, (r, ph))]

        """
        return self._atlas

    def top_charts(self):
        r"""
        Return the list of charts defined on subsets of the current set
        that are not subcharts of charts on larger subsets.

        OUTPUT:

        - list of charts defined on open subsets of ``self`` but not on
          larger subsets

        EXAMPLES:

        Charts on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: Y.<u,v> = U.chart()
            sage: M.top_charts()
            [chart (M, (x, y)), chart (U, (u, v))]

        Note that the (user) atlas contains one more chart: (U, (x,y)), which
        is not a "top" chart::

            sage: M.atlas()
            [chart (M, (x, y)), chart (U, (x, y)), chart (U, (u, v))]

        """
        return self._top_charts


    def default_chart(self):
        r"""
        Return the default chart defined on ``self``.

        Unless changed via :meth:`set_default_chart`, the *default chart*
        is the first one defined on a subset of ``self`` (possibly itself).

        OUTPUT:

        - instance of :class:`~sage.manifolds.chart.Chart`
          representing the default chart.

        EXAMPLES:

        Default chart on a 2-dimensional manifold and on some subsets::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: M.chart('x y')
            chart (M, (x, y))
            sage: M.chart('u v')
            chart (M, (u, v))
            sage: M.default_chart()
            chart (M, (x, y))
            sage: a = M.subset('A')
            sage: b = a.subset('B', is_open=True)
            sage: b.chart('t z')
            chart (B, (t, z))
            sage: a.default_chart()
            chart (B, (t, z))
            sage: b.default_chart()
            chart (B, (t, z))

        """
        return self._def_chart

    def set_default_chart(self, chart):
        r"""
        Changing the default chart on ``self``.

        INPUT:

        - ``chart`` -- a chart (must be defined on some subset ``self``)

        EXAMPLES:

        Charts on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: M.default_chart()
            chart (M, (x, y))
            sage: M.set_default_chart(c_uv)
            sage: M.default_chart()
            chart (M, (u, v))

        """
        from chart import Chart
        if not isinstance(chart, Chart):
            raise TypeError(str(chart) + " is not a chart.")
        if chart._domain is not self:
            if self.is_manifestly_coordinate_domain():
                raise TypeError("The chart domain must coincide with the " +
                                str(self) + ".")
            if chart not in self._atlas:
                raise ValueError("The chart must be defined on the " +
                                 str(self))
        self._def_chart = chart

    def coord_change(self, chart1, chart2):
        r"""
        Return a change of coordinates (transition map) defined on some
        subset of ``self``.

        INPUT:

        - ``chart1`` -- chart 1
        - ``chart2`` -- chart 2

        OUTPUT:

        - instance of :class:`~sage.manifolds.chart.CoordChange`
          representing the transition map from chart 1 to chart 2

        EXAMPLES:

        Change of coordinates on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: c_xy.transition_map(c_uv, (x+y, x-y)) # defines the coordinate change
            coordinate change from chart (M, (x, y)) to chart (M, (u, v))
            sage: M.coord_change(c_xy, c_uv) # returns the coordinate change defined above
            coordinate change from chart (M, (x, y)) to chart (M, (u, v))

        """
        if (chart1, chart2) not in self._coord_changes:
            raise TypeError("The change of coordinates from " + str(chart1) +
                            " to " + str(chart2) + " has not been " +
                            "defined on the " + str(self))
        return self._coord_changes[(chart1, chart2)]


    def coord_changes(self):
        r"""
        Return the changes of coordinates defined on ``self``.
        """
        return self._coord_changes


    def is_manifestly_coordinate_domain(self):
        r"""
        Returns True if ``self`` is known to be the domain of some coordinate
        chart and False otherwise.

        If False is returned, either ``self`` cannot be the domain of some
        coordinate chart or no such chart has been declared yet.
        """
        if not isinstance(self, TopManifoldOpenSubset):
            return False
        return not self._covering_charts == []

    def open_subset(self, name, latex_name=None, coord_def={}):
        r"""
        Create an open subset of ``self``.

        An open subset is a set that is (i) included in ``self`` and (ii)
        open with respect to the manifold's topology.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts on ``self`` and values the symbolic expressions formed by the
          coordinates to define the subset.

        OUTPUT:

        - the open subset, as an instance of :class:`TopManifoldOpenSubset`.

        EXAMPLES:

        Creating an open subset of a manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: a = M.open_subset('A') ; a
            open subset 'A' of the 2-dimensional manifold 'M'

        Creating an open subset of A::

            sage: b = a.open_subset('B') ; b
            open subset 'B' of the 2-dimensional manifold 'M'

        B is then a subset of A and A is a superset of B::

            sage: a._subsets # random (set output)
            set([open subset 'A' of the 2-dimensional manifold 'M',
                 open subset 'B' of the 2-dimensional manifold 'M'])
            sage: b._supersets # random (set output)
            set([open subset 'A' of the 2-dimensional manifold 'M',
                 2-dimensional manifold 'M',
                 open subset 'B' of the 2-dimensional manifold 'M'])

        Defining an open subset by some coordinate restrictions: the open
        unit disk in `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}) ; U
            open subset 'U' of the 2-dimensional manifold 'R^2'

        Since the argument ``coord_def`` has been set, U is automatically
        provided with a chart, which is the restriction of the Cartesian one
        to U::

            sage: U.atlas()
            [chart (U, (x, y))]

        Therefore, one can immediately check whether a point belongs to U::

            sage: M.point((0,0)) in U
            True
            sage: M.point((1/2,1/3)) in U
            True
            sage: M.point((1,2)) in U
            False

        """
        resu = self.subset(name, latex_name=latex_name, is_open=True)
        for chart, restrictions in coord_def.iteritems():
            if chart not in self._atlas:
                raise ValueError("The " + str(chart) + "does not belong to " +
                    "the atlas of " + str(self))
            chart.restrict(resu, restrictions)
        return resu

    def chart(self, coordinates='', names=None):
        r"""
        Define a chart the domain of which is the open set.

        A *chart* is a pair `(U,\varphi)`, where `U` is the current open set
        represented by ``self`` and `\varphi: U \rightarrow V \subset K^n`
        is a homeomorphism from `U` to an open subset `V` of `K^n`, `K` being
        the field on which the manifold containing the open set is defined.

        The components `(x^1,\ldots,x^n)` of `\varphi`, defined by
        `\varphi(p) = (x^1(p),\ldots,x^n(p))`, are called the *coordinates*
        of the chart `(U,\varphi)`.

        See :class:`~sage.manifolds.chart.Chart` for a complete
        documentation.

        INPUT:

        - ``coordinates`` -- single string defining the coordinate symbols and
          ranges: the coordinates are separated by ' ' (space) and each
          coordinate has at most three fields, separated by ':':

          1. The coordinate symbol (a letter or a few letters)
          2. (optional, only for manifolds over `\RR`) The interval `I`
             defining the coordinate range: if not
             provided, the coordinate is assumed to span all `\RR`; otherwise
             `I` must be provided in the form (a,b) (or equivalently ]a,b[)
             The bounds a and b can be +/-Infinity, Inf, infinity, inf or oo.
             For *singular* coordinates, non-open intervals such as [a,b] and
             (a,b] (or equivalently ]a,b]) are allowed.
             Note that the interval declaration must not contain any space
             character.
          3. (optional) The LaTeX spelling of the coordinate; if not provided the
             coordinate symbol given in the first field will be used.

          The order of the fields 2 and 3 does not matter and each of them can
          be omitted.
          If it contains any LaTeX expression, the string ``coordinates`` must
          be declared with the prefix 'r' (for "raw") to allow for a proper
          treatment of the backslash character (see examples below).
          If no interval range and no LaTeX spelling is to be provided for any
          coordinate, the argument ``coordinates`` can be omitted when the
          shortcut operator <,> is used via Sage preparser (see examples below)
        - ``names`` -- (default: ``None``) unused argument, except if
          ``coordinates`` is not provided; it must then be a tuple containing
          the coordinate symbols (this is guaranted if the shortcut operator <,>
          is used).

        OUTPUT:

        - the created chart, as an instance of
          :class:`~sage.manifolds.chart.Chart` or of the subclass
          :class:`~sage.manifolds.chart.RealChart` for manifolds over `\RR`.

        EXAMPLES:

        Chart on a 2-dimensional manifold::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: X = U.chart('x y') ; X
            chart (U, (x, y))
            sage: X[0]
            x
            sage: X[1]
            y
            sage: X[:]
            (x, y)

        The declared coordinates are not known at the global level::

            sage: y
            Traceback (most recent call last):
            ...
            NameError: name 'y' is not defined

        They can be recovered by the Chart method [:]::

            sage: (x, y) = X[:]
            sage: y
            y
            sage: type(y)
            <type 'sage.symbolic.expression.Expression'>

        But a shorter way to proceed is to use the operator <,> in the
        left-hand side of the chart declaration (there is then no need to
        pass the string 'x y' to chart())::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: X.<x,y> = U.chart() ; X
            chart (U, (x, y))

        Indeed, the declared coordinates are then known at the global level::

            sage: y
            y
            sage: (x,y) == X[:]
            True

        Actually the instruction ``X.<x,y> = U.chart()`` is
        equivalent to the two instructions ``X = U.chart('x y')``
        and ``(x,y) = X[:]``.

        See the documentation of class
        :class:`~sage.manifolds.chart.Chart` for more examples,
        especially regarding the coordinates ranges and restrictions.

        """
        from sage.manifolds.chart import Chart, RealChart
        if self._manifold.base_field() == 'real':
            return RealChart(self, coordinates=coordinates, names=names)
        return Chart(self, coordinates=coordinates, names=names)

    def scalar_field_algebra(self):
        r"""
        Returns the algebra of scalar fields defined on ``self``.

        See :class:`~sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra`
        for a complete documentation.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra`
          representing the algebra `C^\infty(U)` of all scalar fields defined
          on `U` = ``self``.

        EXAMPLE:

        Scalar algebra of a 3-dimensional open subset::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: CU = U.scalar_field_algebra() ; CU
            algebra of scalar fields on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: CU.category()
            Category of commutative algebras over Symbolic Ring
            sage: CU.zero()
            scalar field 'zero' on the open subset 'U' of the 3-dimensional manifold 'M'

        """
        return self._scalar_field_algebra

    def _Hom_(self, other, category=None):
        r"""
        Construct the set of morphisms (i.e. continuous maps)
        ``self`` --> ``other``.

        INPUT:

        - ``other`` -- an open subset of some manifold
        - ``category`` -- (default: ``None``) not used here (to ensure
          compatibility with generic hook ``_Hom_``)

        OUTPUT:

        - the homset Hom(U,V), where U is ``self`` and V is ``other``

        See class
        :class:`~sage.manifolds.manifold_homset.ManifoldHomset`
        for more documentation.

        """
        from sage.manifolds.manifold_homset import ManifoldHomset
        return ManifoldHomset(self, other)

    def scalar_field(self, coord_expression=None, chart=None, name=None,
                     latex_name=None):
        r"""
        Define a scalar field on the open set.

        See :class:`~sage.manifolds.scalarfield.ScalarField` for a
        complete documentation.

        INPUT:

        - ``coord_expression`` -- (default: ``None``) coordinate expression(s)
          of the scalar field; this can be either

          - a single coordinate expression; if the argument ``chart`` is
            ``'all'``, this expression is set to all the charts defined
            on the open set; otherwise, the expression is set in the
            specific chart provided by the argument ``chart``
          - a dictionary of coordinate expressions, with the charts as keys.

          If ``coord_expression`` is ``None`` or does not fully specified the
          scalar field, other coordinate expressions can be added subsequently
          by means of the methods
          :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr`,
          :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr_by_continuation`,
          or :meth:`~sage.manifolds.scalarfield.ScalarField.set_expr`
        - ``chart`` -- (default: ``None``) chart defining the coordinates used
          in ``coord_expression`` when the latter is a single coordinate
          expression; if none is provided (default), the default chart of the
          open set is assumed. If ``chart=='all'``, ``coord_expression`` is
          assumed to be independent of the chart (constant scalar field).
        - ``name`` -- (default: ``None``) name given to the scalar field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the scalar
          field; if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of :class:`~sage.manifolds.scalarfield.ScalarField`
          representing the defined scalar field.

        EXAMPLES:

        A scalar field defined by its coordinate expression in the open
        set's default chart::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: c_xyz.<x,y,z> = U.chart()
            sage: f = U.scalar_field(sin(x)*cos(y) + z, name='F'); f
            scalar field 'F' on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: f.display()
            F: U --> R
               (x, y, z) |--> cos(y)*sin(x) + z
            sage: f.parent()
            algebra of scalar fields on the open subset 'U' of the 3-dimensional manifold 'M'
            sage: f in U.scalar_field_algebra()
            True

        Equivalent definition with the chart specified::

            sage: f = U.scalar_field(sin(x)*cos(y) + z, chart=c_xyz, name='F')
            sage: f.display()
            F: U --> R
               (x, y, z) |--> cos(y)*sin(x) + z

        Equivalent definition with a dictionary of coordinate expression(s)::

            sage: f = U.scalar_field({c_xyz: sin(x)*cos(y) + z}, name='F')
            sage: f.display()
            F: U --> R
               (x, y, z) |--> cos(y)*sin(x) + z

        See the documentation of class
        :class:`~sage.manifolds.scalarfield.ScalarField` for more
        examples.

        .. SEEALSO::

            :meth:`constant_scalar_field`, :meth:`zero_scalar_field`

        """
        if isinstance(coord_expression, dict):
            # check validity of entry
            for chart in coord_expression:
                if not chart._domain.is_subset(self):
                    raise ValueError("The " + str(chart) + " is not defined " +
                                     "on some subset of the " + str(self))
        elif coord_expression is not None and chart != 'all':
            # coord_expression is valid only in a specific chart
            if chart is None:
                chart = self._def_chart
            coord_expression = {chart: coord_expression}
        return self.scalar_field_algebra()._element_constructor_(
                                            coord_expression=coord_expression,
                                            name=name, latex_name=latex_name)

    def constant_scalar_field(self, value, name=None, latex_name=None):
        r"""
        Define a constant scalar field on the open set.

        INPUT:

        - ``value`` -- constant value of the scalar field, either a numerical
          value or a symbolic expression not involving any chart coordinates
        - ``name`` -- (default: ``None``) name given to the scalar field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the scalar
          field; if none is provided, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of :class:`~sage.manifolds.scalarfield.ScalarField`
          representing the scalar field whose constant value is ``value``

        EXAMPLES:

        A constant scalar field on the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W',
            ....:                                restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: f = M.constant_scalar_field(1) ; f
            scalar field on the 2-dimensional manifold 'M'
            sage: f.display()
            M --> R
            on U: (x, y) |--> 1
            on V: (u, v) |--> 1

        We have::

            sage: f.restrict(U) == U.constant_scalar_field(1)
            True
            sage: M.constant_scalar_field(0) is M.zero_scalar_field()
            True

        .. SEEALSO::

            :meth:`zero_scalar_field`
        """
        return self.scalar_field_algebra()._element_constructor_(
                                              coord_expression=value,
                                              name=name, latex_name=latex_name)

    def zero_scalar_field(self):
        r"""
        Return the zero scalar field defined on the open set.

        EXAMPLE::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.zero_scalar_field() ; f
            scalar field 'zero' on the 2-dimensional manifold 'M'
            sage: f.display()
            zero: M --> R
               (x, y) |--> 0
            sage: f.parent()
            algebra of scalar fields on the 2-dimensional manifold 'M'
            sage: f is M.scalar_field_algebra().zero()
            True

        """
        return self._zero_scalar_field


    def curve(self, coord_expression, param, chart=None, name=None,
              latex_name=None):
        r"""
        Define a curve in ``self``.

        See :class:`~sage.manifolds.curve.ManifoldCurve` for details.

        INPUT:

        - ``coord_expression`` -- either

          - (i) a dictionary whose keys are charts on ``self`` and values
            the coordinate expressions (as lists or tuples) of the curve in
            the given chart
          - (ii) a single coordinate expression in a given chart on ``self``,
            the latter being provided by the argument ``chart``

          In both cases, if the dimension of the arrival manifold is 1,
          a single coordinate expression can be passed instead of a tuple with
          a single element
        - ``param`` -- a tuple of the type ``(t, t_min, t_max)``, where ``t``
          is the curve parameter used in ``coord_expression``, ``t_min`` is its
          minimal value and ``t_max`` its maximal value; if ``t_min=-Infinity``
          and ``t_max=+Infinity``, they can be omitted and ``t`` can be passed
          for ``param``, instead of the tuple ``(t, t_min, t_max)``
        - ``chart`` -- (default: ``None``) chart on ``self`` used for case (ii)
          above; if ``None`` the default chart of ``self`` is assumed.
        - ``name`` -- (default: ``None``) string; symbol given to the curve
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
          the curve; if none is provided, ``name`` will be used

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.curve.ManifoldCurve`

        EXAMPLES:

        The lemniscate of Gerono in the 2-dimensional Euclidean plane::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
            sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c') ; c
            Curve 'c' in the 2-dimensional manifold 'M'

        The same definition with the coordinate expression passed as a
        dictionary::

            sage: c = M.curve({X: [sin(t), sin(2*t)/2]}, (t, 0, 2*pi), name='c') ; c
            Curve 'c' in the 2-dimensional manifold 'M'

        An example of definition with ``t_min`` and ``t_max`` omitted: a helix
        in `\RR^3`::

            sage: R3 = Manifold(3, 'R^3')
            sage: X.<x,y,z> = R3.chart()
            sage: c = R3.curve([cos(t), sin(t), t], t, name='c') ; c
            Curve 'c' in the 3-dimensional manifold 'R^3'
            sage: c.domain() # check that t is unbounded
            field R of real numbers

        See the documentation of
        :class:`~sage.manifolds.curve.ManifoldCurve` for more
        examples.

        """
        from sage.rings.infinity import Infinity
        from sage.manifolds.smooth.manifold import RealLine
        if not isinstance(param, (tuple, list)):
            param = (param, -Infinity, Infinity)
        elif len(param) != 3:
            raise TypeError("the argument 'param' must be of the type " +
                            "(t, t_min, t_max)")
        t = param[0]
        t_min = param[1]
        t_max = param[2]
        real_field = RealLine(names=(repr(t),))
        interval = real_field.open_interval(t_min, t_max)
        curve_set = Hom(interval, self)
        if not isinstance(coord_expression, dict):
            # Turn coord_expression into a dictionary:
            if chart is None:
                chart = self._def_chart
            elif chart not in self._atlas:
                raise ValueError("the {} has not been".format(chart) +
                                     " defined on the {}".format(self))
            if isinstance(coord_expression, (tuple, list)):
                coord_expression = {chart: coord_expression}
            else:
                # case self.dim()=1
                coord_expression = {chart: (coord_expression,)}
        return curve_set(coord_expression, name=name, latex_name=latex_name)

    def continuous_mapping(self, codomain, coord_functions=None, chart1=None,
                     chart2=None, name=None, latex_name=None):
        r"""
        Define a continuous mapping between ``self`` and another
        subset (possibly on another manifold).

        See :class:`~sage.manifolds.diffmapping.DiffMapping` for a
        complete documentation.

        INPUT:

        - ``codomain`` -- mapping's codomain (the arrival manifold or some
          subset of it)
        - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
          either

          - (i) a dictionary of
            the coordinate expressions (as lists (or tuples) of the
            coordinates of the image expressed in terms of the coordinates of
            the considered point) with the pairs of charts (chart1, chart2)
            as keys (chart1 being a chart on ``self`` and chart2 a chart on
            ``codomain``)
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``

          In both cases, if the dimension of the arrival manifold is 1,
          a single coordinate expression can be passed instead of a tuple with
          a single element
        - ``chart1`` -- (default: ``None``; used only in case (ii) above) chart
          on ``self`` defining the start coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``self``
        - ``chart2`` -- (default: ``None``; used only in case (ii) above) chart
          on ``codomain`` defining the arrival coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``codomain``
        - ``name`` -- (default: ``None``) name given to the differentiable
          mapping
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          differentiable mapping; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - the differentiable mapping, as an instance of
          :class:`~sage.manifolds.diffmapping.DiffMapping`

        EXAMPLES:

        A mapping between an open subset of `S^2` covered by regular spherical
        coordinates and `\RR^3`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'S^2')
            sage: U = M.open_subset('U')
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: N = Manifold(3, 'R^3', r'\RR^3')
            sage: c_cart.<x,y,z> = N.chart()  # Cartesian coord. on R^3
            sage: Phi = U.diff_mapping(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)),
            ....:                      name='Phi', latex_name=r'\Phi') ; Phi
            differentiable mapping 'Phi' from the open subset 'U' of the
             2-dimensional manifold 'S^2' to the 3-dimensional manifold 'R^3'

        The same definition, but with a dictionary with pairs of charts as
        keys (case (i) above)::

            sage: Phi1 = U.diff_mapping(N,
            ....:        {(c_spher, c_cart): (sin(th)*cos(ph), sin(th)*sin(ph), cos(th))},
            ....:        name='Phi', latex_name=r'\Phi')
            sage: Phi1 == Phi
            True

        The differentiable mapping acting on a point::

            sage: p = U.point((pi/2, pi)) ; p
            point on 2-dimensional manifold 'S^2'
            sage: Phi(p)
            point on 3-dimensional manifold 'R^3'
            sage: Phi(p).coord(c_cart)
            (-1, 0, 0)
            sage: Phi1(p) == Phi(p)
            True

        See the documentation of class
        :class:`~sage.manifolds.diffmapping.DiffMapping` for more
        examples.

        """
        homset = Hom(self, codomain)
        if coord_functions is None:
            coord_functions = {}
        if not isinstance(coord_functions, dict):
            # Turn coord_functions into a dictionary:
            if chart1 is None:
                chart1 = self._def_chart
            elif chart1 not in self._atlas:
                raise ValueError("{} is not a chart ".format(chart1) +
                                 "defined on the {}".format(self))
            if chart2 is None:
                chart2 = codomain._def_chart
            elif chart2 not in codomain._atlas:
                raise ValueError("{} is not a chart ".format(chart2) +
                                 " defined on the {}".format(codomain))
            coord_functions = {(chart1, chart2): coord_functions}
        return homset(coord_functions, name=name, latex_name=latex_name)

    def homeomorphism(self, codomain, coord_functions=None, chart1=None,
                       chart2=None, name=None, latex_name=None):
        r"""
        Define a homeomorphism between ``self`` and another open subset
        (possibly on another manifold).

        See :class:`~sage.manifolds.diffmapping.DiffMapping` for a
        complete documentation.

        INPUT:

        - ``codomain`` -- mapping's codomain (the arrival manifold or some
          subset of it)
        - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
          either

          - (i) a dictionary of
            the coordinate expressions (as lists (or tuples) of the
            coordinates of the image expressed in terms of the coordinates of
            the considered point) with the pairs of charts (chart1, chart2)
            as keys (chart1 being a chart on ``self`` and chart2 a chart on
            ``codomain``)
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``

          In both cases, if the dimension of the arrival manifold is 1,
          a single coordinate expression can be passed instead of a tuple with
          a single element
        - ``chart1`` -- (default: ``None``; used only in case (ii) above) chart
          on ``self`` defining the start coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``self``
        - ``chart2`` -- (default: ``None``; used only in case (ii) above) chart
          on ``codomain`` defining the arrival coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``codomain``
        - ``name`` -- (default: ``None``) name given to the diffeomorphism
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          diffeomorphism; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - the diffeomorphism, as an instance of
          :class:`~sage.manifolds.diffmapping.DiffMapping`

        EXAMPLE:

        A diffeomorphism between two 2-dimensional subsets::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M', r'{\cal M}')
            sage: U = M.open_subset('U')
            sage: c_xv.<x,y> = U.chart(r'x:(-pi/2,+oo) y:(-pi/2,+oo)')
            sage: N = Manifold(2, 'N', r'{\cal N}')
            sage: V = N.open_subset('V')
            sage: c_zt.<z,t> = V.chart(r'z t')
            sage: Phi = U.diffeomorphism(V, (arctan(x), arctan(y)), name='Phi',
            ....:                        latex_name=r'\Phi')

        See the documentation of class
        :class:`~sage.manifolds.diffmapping.DiffMapping` for more
        examples.

        """
        homset = Hom(self, codomain)
        if coord_functions is None:
            coord_functions = {}
        if not isinstance(coord_functions, dict):
            # Turn coord_functions into a dictionary:
            if chart1 is None:
                chart1 = self._def_chart
            elif chart1 not in self._atlas:
                raise ValueError("{} is not a chart ".format(chart1) +
                                 "defined on the {}".format(self))
            if chart2 is None:
                chart2 = codomain._def_chart
            elif chart2 not in codomain._atlas:
                raise ValueError("{} is not a chart ".format(chart2) +
                                 " defined on the {}".format(codomain))
            coord_functions = {(chart1, chart2): coord_functions}
        return homset(coord_functions, name=name, latex_name=latex_name,
                      is_diffeomorphism=True)

    def identity_map(self):
        r"""
        Identity map on ``self``.

        See :class:`~sage.manifolds.diffmapping.DiffMapping` for a
        complete documentation.

        OUTPUT:

        - the identity map, as an instance of
          :class:`~sage.manifolds.diffmapping.DiffMapping`

        """
        return self._identity_map

