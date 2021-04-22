r"""
Immutable Sets of Subsets of Topological Manifolds

The class :class:`ManifoldSubsetSet` is a subclass of the built-in
``frozenset`` class that provides specialized ``__repr__`` and
``_latex_`` methods.

:class:`ManifoldSubsetSet` instances are totally ordered according
to their lexicographically ordered element (subset) names.

"""
#*****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from functools import total_ordering

@total_ordering
class ManifoldSubsetSet(frozenset):

    def __new__(cls, *args):
        self = super().__new__(cls, *args)
        names_and_latex_names = sorted((subset._name, subset._latex_name)
                                       for subset in self)
        self._names = tuple(name for name, latex_name in names_and_latex_names)
        self._name = '{' + ', '.join(self._names) + '}'
        latex_names = (latex_name for name, latex_name in names_and_latex_names)
        self._latex_name = r'\{' + ', '.join(latex_names) + r'\}'
        try:
            subset_iter = iter(self)
            self._manifold = next(subset_iter)._manifold
        except StopIteration:
            pass
        else:
            if not all(subset._manifold == self._manifold for subset in subset_iter):
                raise TypeError('all elements must be subsets of the same manifold')
        return self

    def __repr__(self):
        r"""
        String representation of the object.

        TESTS::

            sage: from sage.manifolds.subset import ManifoldSubsetSet
            sage: ManifoldSubsetSet().__repr__()
            '{}'
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: ManifoldSubsetSet([A, B]).__repr__()
            'Set {A, B} of subsets of the 2-dimensional topological manifold M'

        """
        if self:
            return "Set {} of subsets of the {}".format(self._name, self._manifold)
        else:
            return "{}"

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        TESTS::

            sage: from sage.manifolds.subset import ManifoldSubsetSet
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: ManifoldSubsetSet([B, A])._latex_()
            '\\{A, B\\}'
        """
        return self._latex_name

    def __lt__(self, other):
        r"""
        Implement the total order on instances of :class:`ManifoldSubsetSet`.

        TESTS::

            sage: from sage.manifolds.subset import ManifoldSubsetSet
            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: sorted([ManifoldSubsetSet([A, B]), ManifoldSubsetSet([]),
            ....:         ManifoldSubsetSet([B]), ManifoldSubsetSet([A])])
            [{},
             Set {A} of subsets of the 2-dimensional topological manifold M,
             Set {A, B} of subsets of the 2-dimensional topological manifold M,
             Set {B} of subsets of the 2-dimensional topological manifold M]
        """
        return self._names < other._names
