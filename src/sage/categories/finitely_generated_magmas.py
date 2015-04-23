r"""
Finitely generated magmas
"""
#*****************************************************************************
#  Copyright (C) 2014 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.magmas import Magmas

class FinitelyGeneratedMagmas(CategoryWithAxiom):
    r"""
    The category of finitely generated (multiplicative) magmas.

    See :meth:`Magmas.SubcategoryMethods.FinitelyGeneratedAsMagma` for
    details.

    EXAMPLES::

        sage: C = Magmas().FinitelyGeneratedAsMagma(); C
        Category of finitely generated magmas
        sage: C.super_categories()
        [Category of magmas]
        sage: sorted(C.axioms())
        ['FinitelyGeneratedAsMagma']

    TESTS::

        sage: TestSuite(C).run()
    """

    _base_category_class_and_axiom = (Magmas, "FinitelyGeneratedAsMagma")

    class ParentMethods:

        @abstract_method
        def magma_generators(self):
            """
            Return distinguished magma generators for ``self``.

            OUTPUT: a finite family

            This method should be implemented by all
            :class:`finitely generated magmas <FinitelyGeneratedMagmas>`_.

            EXAMPLES::

                sage: S = FiniteSemigroups().example()
                sage: S.magma_generators()
                Family ('a', 'b', 'c', 'd')
            """

    class Unital(CategoryWithAxiom):
        class ParentMethods:
            def is_empty(self):
                r"""
                Return whether ``self`` is empty.

                Since this set is a unital magma it is not empty and this method
                always return ``False``.

                EXAMPLES::

                    sage: SymmetricGroup(2).is_empty()
                    False

                TESTS::

                    sage: SymmetricGroup(2).is_empty.__module__
                    'sage.categories.finitely_generated_magmas'

                .. TODO::

                    This method belong to ``Magmas.Unital.ParentMethods``. But due to
                    another implementation of ``is_empty`` in ``EnumeratedSets``
                    and the fact that ``Enumerated`` is not an axiom we get that
                    for an enumerated unital magmas the ``is_empty`` which is
                    used is the one from enumerated sets. But it should be this
                    one!
                """
                return False
