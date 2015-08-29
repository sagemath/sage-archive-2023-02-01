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
from sage.categories.with_realizations import WithRealizationsCategory

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

        In practice, this returns ``category.Metric()``, joined
        together with the result of the method
        :meth:`RegressiveCovariantConstructionCategory.default_super_categories()
        <sage.categories.covariant_functorial_construction.RegressiveCovariantConstructionCategory.default_super_categories>`
        (that is the join of ``category`` and ``cat.Metric()`` for
        each ``cat`` in the super categories of ``category``).

        EXAMPLES:

        Consider ``category=Groups()``. Then, a group `G` with a metric
        is simultaneously a topological group by itself, and a
        topogological space::

            sage: Groups().Metric().super_categories()
            [Category of topological groups, Category of metric spaces]

        This resulted from the following call::

            sage: sage.categories.metric_spaces.MetricSpacesCategory.default_super_categories(Groups())
            Join of Category of topological groups and Category of metric spaces
        """
        return Category.join([category.Topological(),
                              super(MetricSpacesCategory, cls).default_super_categories(category)])

    # We currently don't have a use for this, but we probably will
    def _repr_object_names(self):
        """
        EXAMPLES::

            sage: Groups().Metric()  # indirect doctest
            Join of Category of topological groups and Category of metric spaces
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

            INPUT:

            - ``options`` -- any keyword arguments accepted
              by :meth:`_tester`

            EXAMPLES::

                sage: UHP = HyperbolicPlane().UHP()
                sage: UHP._test_metric()
                sage: elts = [UHP.random_element() for i in range(5)]
                sage: UHP._test_metric(some_elements=elts)
            """
            tester = self._tester(**options)
            S = tester.some_elements()
            dist = self.metric()
            for a in S:
                for b in S:
                    d = dist(a, b)
                    if a is not b:
                        tester.assertGreater(d, 0)
                    else:
                        tester.assertEqual(d, 0)

        def metric(self):
            """
            Return the metric of ``self``.

            EXAMPLES::

                sage: UHP = HyperbolicPlane().UHP()
                sage: m = UHP.metric()
                sage: p1 = UHP.get_point(5 + 7*I)
                sage: p2 = UHP.get_point(1.0 + I)
                sage: m(p1, p2)
                2.23230104635820
            """
            return lambda a,b: self.dist(a, b)

        @abstract_method
        def dist(self, a, b):
            """
            Return the distance between ``a`` and ``b`` in ``self``.

            EXAMPLES::

                sage: UHP = HyperbolicPlane().UHP()
                sage: p1 = UHP.get_point(5 + 7*I)
                sage: p2 = UHP.get_point(1.0 + I)
                sage: UHP.dist(p1, p2)
                2.23230104635820

                sage: PD = HyperbolicPlane().PD()
                sage: PD.dist(PD.get_point(0), PD.get_point(I/2))
                arccosh(5/3)
            """

    class ElementMethods:
        def dist(self, b):
            """
            Return the distance between ``self`` and ``other``.

            EXAMPLES::

                sage: UHP = HyperbolicPlane().UHP()
                sage: p1 = UHP.get_point(5 + 7*I)
                sage: p2 = UHP.get_point(1 + I)
                sage: p1.dist(p2)
                arccosh(33/7)
            """
            return self.parent().dist(self, b)

    class WithRealizations(WithRealizationsCategory):
        class ParentMethods:
            def dist(self, a, b):
                """
                Return the distance between ``a`` and ``b`` by converting them
                to a realization of ``self`` and doing the computation.

                EXAMPLES::

                    sage: H = HyperbolicPlane()
                    sage: PD = H.PD()
                    sage: p1 = PD.get_point(0)
                    sage: p2 = PD.get_point(I/2)
                    sage: H.dist(p1, p2)
                    arccosh(5/3)
                """
                R = self.a_realization()
                return R.dist(R(a), R(b))

