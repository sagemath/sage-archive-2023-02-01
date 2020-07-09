r"""
Finitely Generated Lie Conformal Algebra

AUTHORS:

- Reimundo Heluani (2020-07-09): Initial implementation.
"""


#******************************************************************************
#       Copyright (C) 2020 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .lie_conformal_algebra import LieConformalAlgebra

class FinitelyGeneratedLieConformalAlgebra(LieConformalAlgebra):
    """
    Base class for finitely generated (super) Lie conformal algebras.
    """
    def gen(self,i):
        r"""
        The ``i``-th generator of this Lie conformal algebra.

        EXAMPLES::

            sage: V = lie_conformal_algebras.Affine(QQ, 'A1')
            sage: V.gens()
            (B[alpha[1]], B[alphacheck[1]], B[-alpha[1]], B['K'])
            sage: V.gen(0)
            B[alpha[1]]
            sage: V.1
            B[alphacheck[1]]
        """
        return self.gens()[i]

    def some_elements(self):
        """
        Some elements of this Lie conformal algebra.

        This method returns a list with elements containing at
        least the generators.

        EXAMPLES::

            sage: V = lie_conformal_algebras.Affine(QQ, 'A1', names=('e', 'h', 'f'))
            sage: V.some_elements()
            [e, h, f, K, Th + 4*T^(2)e, 4*T^(2)h, Te + 4*T^(2)e, Te + 4*T^(2)h]
        """
        S = list(self.gens())
        from sage.misc.misc import some_tuples
        for x,y in some_tuples(S, 2, 0, max_samples=self.ngens()):
            S.append(x.T() + 2*y.T(2))
        return S

    def ngens(self):
        r"""
        The number of generators of this Lie conformal algebra.

        EXAMPLES::

            sage: Vir = lie_conformal_algebras.Virasoro(QQ)
            sage: Vir.ngens()
            2

            sage: V = lie_conformal_algebras.Affine(QQ, 'A2')
            sage: V.ngens()
            9
        """
        return len(self.gens())
