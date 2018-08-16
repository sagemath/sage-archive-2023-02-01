r"""
Nilpotent Lie algebras

AUTHORS:

- Eero Hakavuori (2018-08-16): initial version
"""

# ****************************************************************************
#       Copyright (C) 2018 Eero Hakavuori <eero.hakavuori@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.lie_algebras import LieAlgebras
from sage.misc.abstract_method import abstract_method


class NilpotentLieAlgebras(CategoryWithAxiom_over_base_ring):
    r"""
    Category of nilpotent Lie algebras.

    TESTS::

        sage: C = LieAlgebras(QQ).Nilpotent()
        sage: TestSuite(C).run()

    """

    class ParentMethods:
        @abstract_method
        def step(self):
            r"""
            Returns the nilpotency step of the Lie algebra.
            """


class FiniteDimensionalNilpotentLieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    """
    Category of finite dimensional nilpotent Lie algebras with basis.

    TESTS::

        sage: C1 = LieAlgebras(QQ).FiniteDimensional().WithBasis().Nilpotent()
        sage: C2 = LieAlgebras(QQ).Nilpotent().FiniteDimensional().WithBasis()
        sage: C1 is C2
        True
        sage: TestSuite(C1).run()

    """
    _base_category_class_and_axiom = [LieAlgebras.FiniteDimensional.WithBasis,
                                      "Nilpotent"]

    class ParentMethods:

        def step(self):
            r"""
            Returns the nilpotency step of the Lie algebra.
            
            EXAMPLES::
            
                sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                sage: NilpotentLieAlgebra(QQ, {('X','Y'): {'Z': 1}}).step()
                2
                sage: sc = {('X','Y'): {'Z': 1}, ('X','Z'): {'W': 1}}
                sage: NilpotentLieAlgebra(QQ, sc).step()
                3
            """
            if not hasattr(self, '_step'):
                self._step = len(self.lower_central_series(submodule=True)) - 1
            return self._step

        def _test_nilpotency(self, **options):
            r"""
            Tests that the Lie algebra is nilpotent and has the correct step.

            INPUT:

            - ``options`` -- any keyword arguments accepted by :meth:`_tester`.

            EXAMPLES::
            
                sage: from sage.algebras.lie_algebras.nilpotent_lie_algebra import NilpotentLieAlgebra
                sage: L = NilpotentLieAlgebra(QQ, {('X','Y'): {'Z': 1}})
                sage: L._test_nilpotency()
                sage: L = NilpotentLieAlgebra(QQ, {('X','Y'): {'Z': 1}}, step = 3)
                sage: L._test_nilpotency()
                Traceback (most recent call last):
                ...
                AssertionError: claimed nilpotency step 3 does not match the actual nilpotency step 2
                sage: L = NilpotentLieAlgebra(QQ, {('X','Y'): {'X': 1}})
                sage: L._test_nilpotency()
                Traceback (most recent call last):
                ...
                AssertionError: final term of lower central series is non-zero

            See the documentation for :class:`TestSuite` for more information.
            """
            tester = self._tester(**options)

            lcs = self.lower_central_series(submodule=True)
            tester.assertEqual(lcs[-1].dimension(), 0,
                msg="final term of lower central series is non-zero")

            if hasattr(self, '_step'):
                tester.assertEqual(len(lcs) - 1, self._step,
                    msg="claimed nilpotency step %d does not match the "
                    "actual nilpotency step %d" % (self._step, len(lcs) - 1))
