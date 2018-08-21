r"""
Finite Dimensional Nilpotent Lie Algebras With Basis

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
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.categories.lie_algebras import LieAlgebras

class FiniteDimensionalNilpotentLieAlgebrasWithBasis(CategoryWithAxiom_over_base_ring):
    r"""
    Category of finite dimensional nilpotent Lie algebras with basis.

    TESTS::

        sage: C1 = LieAlgebras(QQ).FiniteDimensional().WithBasis().Nilpotent()
        sage: C2 = LieAlgebras(QQ).FiniteDimensional().Nilpotent().WithBasis()
        sage: C3 = LieAlgebras(QQ).Nilpotent().FiniteDimensional().WithBasis()
        sage: C4 = LieAlgebras(QQ).Nilpotent().WithBasis().FiniteDimensional()
        sage: C5 = LieAlgebras(QQ).WithBasis().Nilpotent().FiniteDimensional()
        sage: C6 = LieAlgebras(QQ).WithBasis().FiniteDimensional().Nilpotent()
        sage: C1 is C2
        True
        sage: C2 is C3
        True
        sage: C3 is C4
        True
        sage: C4 is C5
        True
        sage: C5 is C6
        True
        sage: TestSuite(C1).run()
    """
    _base_category_class_and_axiom = [LieAlgebras.FiniteDimensional.WithBasis, "Nilpotent"]

    class ParentMethods:
        def _test_nilpotency(self, **options):
            r"""
            Tests that the Lie algebra is nilpotent and has
            the correct step.

            INPUT:

            - ``options`` -- any keyword arguments accepted by
              :meth:`_tester`.

            EXAMPLES::
            
                sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
                sage: L._test_nilpotency()
                sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1}},
                ....:                nilpotent=True, step = 3)
                sage: L._test_nilpotency()
                Traceback (most recent call last):
                ...
                AssertionError: claimed nilpotency step 3 does not match the actual nilpotency step 2
                sage: L = LieAlgebra(QQ, {('X','Y'): {'X': 1}}, nilpotent=True)
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

            step = self.step()
            tester.assertEqual(len(lcs) - 1, step,
                msg="claimed nilpotency step %d does not match the "
                "actual nilpotency step %d" % (step, len(lcs) - 1))

        def step(self):
            r"""
            Return the nilpotency step of the Lie algebra.
            
            EXAMPLES::
            
                sage: L = LieAlgebra(QQ, {('X','Y'): {'Z': 1}}, nilpotent=True)
                sage: L.step()
                2
                sage: sc = {('X','Y'): {'Z': 1}, ('X','Z'): {'W': 1}}
                sage: LieAlgebra(QQ, sc, nilpotent=True).step()
                3
            """
            if not hasattr(self, '_step'):
                self._step = len(self.lower_central_series(submodule=True)) - 1
            return self._step

        def is_nilpotent(self):
            r"""
            Return ``True`` since ``self`` is nilpotent.

            EXAMPLES::

                sage: L = LieAlgebra(QQ, {('x','y'): {'z': 1}}, nilpotent=True)
                sage: L.is_nilpotent()
                True
            """
            return True

