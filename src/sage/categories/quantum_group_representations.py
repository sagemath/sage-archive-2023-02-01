r"""
Quantum Group Representations

AUTHORS:

- Travis Scrimshaw (2018): initial version
"""

#*****************************************************************************
#       Copyright (C) 2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.categories.modules import Modules
from sage.categories.category_types import Category_module
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring

class QuantumGroupRepresentations(Category_module):
    """
    The category of quantum group representations.
    """
    @cached_method
    def super_categories(self):
        """
        Return the super categories of ``self``.

        EXAMPLES::

            sage: from sage.categories.quantum_group_representations import QuantumGroupRepresentations
            sage: QuantumGroupRepresentations(ZZ['q'].fraction_field()).super_categories()
            [Category of vector spaces over
             Fraction Field of Univariate Polynomial Ring in q over Integer Ring]
        """
        return [Modules(self.base_ring())]

    def example(self):
        """
        Return an example of a quantum group representation as per
        :meth:`Category.example <sage.categories.category.Category.example>`.

        EXAMPLES::

            sage: from sage.categories.quantum_group_representations import QuantumGroupRepresentations
            sage: Cat = QuantumGroupRepresentations(ZZ['q'].fraction_field())
            sage: Cat.example()
            V((2, 1, 0))
        """
        from sage.algebras.quantum_groups.representations import AdjointRepresentation
        from sage.combinat.crystals.tensor_product import CrystalOfTableaux
        T = CrystalOfTableaux(['A',2], shape=[2,1])
        return AdjointRepresentation(self.base_ring(), T)

    class WithBasis(CategoryWithAxiom_over_base_ring):
        """
        The category of quantum group representations with a
        distinguished basis.
        """
        class ElementMethods:
            def e(self, i):
                r"""
                Return the action of `e_i` on ``self``.

                INPUT:

                - ``i`` -- an element of the index set

                EXAMPLES::

                    sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
                    sage: C = crystals.Tableaux(['G',2], shape=[1,1])
                    sage: R = ZZ['q'].fraction_field()
                    sage: V = AdjointRepresentation(R, C)
                    sage: v = V.an_element(); v
                    2*B[[[1], [2]]] + 2*B[[[1], [3]]] + 3*B[[[2], [3]]]
                    sage: v.e(1)
                    ((3*q^4+3*q^2+3)/q^2)*B[[[1], [3]]]
                    sage: v.e(2)
                    2*B[[[1], [2]]]
                """
                F = self.parent()
                mc = self.monomial_coefficients(copy=False)
                return F.linear_combination( (F.e_on_basis(i, m), c)
                                             for m,c in mc.iteritems() )

            def f(self, i):
                r"""
                Return the action of `f_i` on ``self``.

                INPUT:

                - ``i`` -- an element of the index set

                EXAMPLES::

                    sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
                    sage: K = crystals.KirillovReshetikhin(['D',4,1], 2,1)
                    sage: R = ZZ['q'].fraction_field()
                    sage: V = AdjointRepresentation(R, K)
                    sage: v = V.an_element(); v
                    2*B[[]] + 2*B[[[1], [2]]] + 3*B[[[1], [3]]]
                    sage: v.f(0)
                    ((2*q^2+2)/q)*B[[[1], [2]]]
                    sage: v.f(1)
                    3*B[[[2], [3]]]
                    sage: v.f(2)
                    2*B[[[1], [3]]]
                    sage: v.f(3)
                    3*B[[[1], [4]]]
                    sage: v.f(4)
                    3*B[[[1], [-4]]]
                """
                F = self.parent()
                mc = self.monomial_coefficients(copy=False)
                return F.linear_combination( (F.f_on_basis(i, m), c)
                                             for m,c in mc.iteritems() )

            def K(self, i, power=1):
                r"""
                Return the action of `K_i` on ``self`` to the power ``power``.

                INPUT:

                - ``i`` -- an element of the index set
                - ``power`` -- (default: 1) the power of `K_i`

                EXAMPLES::

                    sage: from sage.algebras.quantum_groups.representations import AdjointRepresentation
                    sage: K = crystals.KirillovReshetikhin(['D',4,2], 1,1)
                    sage: R = ZZ['q'].fraction_field()
                    sage: V = AdjointRepresentation(R, K)
                    sage: v = V.an_element(); v
                    2*B[[]] + 2*B[[[1]]] + 3*B[[[2]]]
                    sage: v.K(0)
                    2*B[[]] + 2/q^2*B[[[1]]] + 3*B[[[2]]]
                    sage: v.K(1)
                    2*B[[]] + 2*q^2*B[[[1]]] + 3/q^2*B[[[2]]]
                    sage: v.K(1, 2)
                    2*B[[]] + 2*q^4*B[[[1]]] + 3/q^4*B[[[2]]]
                    sage: v.K(1, -1)
                    2*B[[]] + 2/q^2*B[[[1]]] + 3*q^2*B[[[2]]]
                """
                F = self.parent()
                mc = self.monomial_coefficients(copy=False)
                return F.linear_combination( (F.K_on_basis(i, m, power), c)
                                             for m,c in mc.iteritems() )

    class ParentMethods:
        def _test_representation(self, tester=None, **options):
            """
            Test the quantum group relations on ``self``.

            .. SEEALSO:: :class:`TestSuite`

            EXAMPLES::

                sage: from sage.algebras.quantum_groups.representations import \
                ....:  MinusculeRepresentation, AdjointRepresentation
                sage: C = crystals.Tableaux(['G',2], shape=[1,1])
                sage: R = ZZ['q'].fraction_field()
                sage: V = AdjointRepresentation(R, C)
                sage: V._test_representation()

            We verify that ``C`` does not define a minuscule
            representation::

                sage: M = MinusculeRepresentation(R, C)
                sage: M._test_representation()
                Traceback (most recent call last):
                ...
                AssertionError: [e,f] = (K-K^-1)/(q_i-q_i^-1) -- i: 1 j: 1
            """
            tester = self._tester(**options)
            ct = self.cartan_type()
            d = ct.symmetrizer()
            I = ct.index_set()
            A = ct.cartan_matrix()
            al = ct.root_system().weight_lattice().simple_roots()
            ac = ct.root_system().weight_lattice().simple_coroots()
            q = self.q()
            from sage.algebras.quantum_groups.q_numbers import q_factorial
            def apply_e(d, elt):
                for i in d:
                    elt = elt.e(i)
                return elt

            def apply_f(d, elt):
                for i in d:
                    elt = elt.f(i)
                return elt

            count = 0
            for x in self.basis():
                for i in I:
                    for j in I:
                        tester.assertEqual(x.K(j,-1).f(i).K(j,1),
                                           q**-(al[i].scalar(ac[j]) * d[j]) * x.f(i),
                                           "KfK^-1 -- i: {}, j: {}".format(i,j))
                        tester.assertEqual(x.K(j,-1).e(i).K(j,1),
                                           q**(al[i].scalar(ac[j]) * d[j]) * x.e(i),
                                           "KeK^-1 -- i: {}, j: {}".format(i,j))
                        if i == j:
                            tester.assertEqual(x.f(i).e(i) - x.e(i).f(i),
                                               (x.K(i,1) - x.K(i,-1)) / (q**d[i] - q**(-d[i])),
                                               "[e,f] = (K-K^-1)/(q_i-q_i^-1) -- i: {} j: {}".format(i, j))
                            continue
                        tester.assertEqual(x.f(j).e(i) - x.e(i).f(j), 0,
                                           "[e,f] = 0 -- i: {} j: {}".format(i, j))
                        # Check quantum Serre
                        aij = A[I.index(i),I.index(j)]
                        tester.assertEqual(0,
                                           sum((-1)**n
                                               * q_factorial(1-aij, q**d[i])
                                               / (q_factorial(n, q**d[i])
                                                  * q_factorial(1-aij-n, q**d[i]))
                                               * apply_e([i]*(1-aij-n) + [j] + [i]*n, x)
                                               for n in range(1-aij+1)),
                                           "quantum Serre e -- i: {}, j: {}".format(i,j))
                        tester.assertEqual(0,
                                           sum((-1)**n
                                               * q_factorial(1-aij, q**d[i])
                                               / (q_factorial(n, q**d[i])
                                                  * q_factorial(1-aij-n, q**d[i]))
                                               * apply_f([i]*(1-aij-n) + [j] + [i]*n, x)
                                               for n in range(1-aij+1)),
                                           "quantum Serre f -- i: {}, j: {}".format(i,j))
                count += 1
                if count > tester._max_runs:
                    return

        @abstract_method
        def cartan_type(self):
            """
            Return the Cartan type of ``self``.

            EXAMPLES::

                sage: from sage.algebras.quantum_groups.representations import MinusculeRepresentation
                sage: C = crystals.Tableaux(['C',4], shape=[1])
                sage: R = ZZ['q'].fraction_field()
                sage: V = MinusculeRepresentation(R, C)
                sage: V.cartan_type()
                ['C', 4]
            """

        @cached_method
        def index_set(self):
            """
            Return the index set of ``self``.

            EXAMPLES::

                sage: from sage.algebras.quantum_groups.representations import MinusculeRepresentation
                sage: C = crystals.Tableaux(['C',4], shape=[1])
                sage: R = ZZ['q'].fraction_field()
                sage: V = MinusculeRepresentation(R, C)
                sage: V.index_set()
                (1, 2, 3, 4)
            """
            return self.cartan_type().index_set()

        def q(self):
            r"""
            Return the quantum parameter `q` of ``self``.

            EXAMPLES::

                sage: from sage.algebras.quantum_groups.representations import MinusculeRepresentation
                sage: C = crystals.Tableaux(['C',4], shape=[1])
                sage: R = ZZ['q'].fraction_field()
                sage: V = MinusculeRepresentation(R, C)
                sage: V.q()
                q
            """
            return self._q

