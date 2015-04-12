r"""
Topological Spaces
"""
#*****************************************************************************
#  Copyright (C) 2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.categories.category import Category
from sage.categories.covariant_functorial_construction import RegressiveCovariantConstructionCategory

class MetricSpacesCategory(RegressiveCovariantConstructionCategory):

    _functor_category = "Metric"

    @classmethod
    def default_super_categories(cls, category):
        """
        Return the default super categories of ``category.Metric()``.

        Mathematical meaning: if `A` is a metric space in the
        category `C`, then `A` is also a topological space.

        INPUT:

         - ``cls`` -- the class ``MetricSpaces``
         - ``category`` -- a category `Cat`

        OUTPUT:

        A (join) category

        In practice, this returns ``category.Subquotients()``, joined
        together with the result of the method
        :meth:`RegressiveCovariantConstructionCategory.default_super_categories()
        <sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory.default_super_categories>`
        (that is the join of ``category`` and ``cat.Metric()`` for
        each ``cat`` in the super categories of ``category``).

        EXAMPLES:

        Consider ``category=Groups()``, which has ``cat=Monoids()`` as
        super category. Then, a subgroup of a group `G` is
        simultaneously a subquotient of `G`, a group by itself, and a
        submonoid of `G`::

            sage: Groups().Subobjects().super_categories()
            [Category of groups, Category of subquotients of monoids, Category of subobjects of sets]

        Mind the last item above: there is indeed currently nothing
        implemented about submonoids.

        This resulted from the following call::

            sage: sage.categories.subobjects.SubobjectsCategory.default_super_categories(Groups())
            Join of Category of groups and Category of subquotients of monoids and Category of subobjects of sets
        """
        return Category.join([category.Topological(),
                              super(MetricSpaces, cls).default_super_categories(category)])

    # We currently don't have a use for this, but we probably will
    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Groups().Metric()  # indirect doctest
            Category of metric groups
        """
        return "metric {}".format(self.base_category()._repr_object_names())

class MetricSpaces(MetricSpacesCategory):
    r"""
    The category of metric spaces.

    A *metric* on a set `S` is a function `d : S \times S \to \RR`
    such that:

    - `d(a, b) \geq 0`,
    - `d(a, b) = 0` if and only if `a = b`.

    A metric space is a set `S` with a distinguished metric.
    """
    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Sets().Metric()  # indirect doctest
            Category of metric spaces
        """
        return "metric spaces"

    class ParentMethods:
        def _test_metric(self, **options):
            r"""
            Test that this metric space has a properly implemented metric.

            INPUT::

            - ``options`` -- any keyword arguments accepted
              by :meth:`_tester`
            """
            tester = self._tester(**options)
            S = tester.some_elements()
            dist = self.metric()
            for a in S:
                for b in S:
                    d = dist(a, b)
                    if a != b:
                        tester.assertGreater(d, 0)
                    else:
                        tester.assertEqual(d, 0)

        def metric(self):
            """
            Return the metric of ``self``.
            """
            return lambda a,b: self.dist(a, b)

        @abstract_method
        def dist(self, a, b):
            """
            Return the distance between ``a`` and ``b`` in ``self``.
            """

    class ElementMethods:
        def dist(self, b):
            """
            Return the distance between ``self`` and ``b``.
            """
            return self.parent().dist(self, b)

