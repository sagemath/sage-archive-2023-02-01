"""
Morphisms between finite algebras
"""

#*****************************************************************************
#  Copyright (C) 2011 Johan Bosman <johan.g.bosman@gmail.com>
#  Copyright (C) 2011, 2013 Peter Bruin <peter.bruin@math.uzh.ch>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.homset import Hom
from sage.rings.morphism import RingHomomorphism_im_gens


class FiniteAlgebraMorphism(RingHomomorphism_im_gens):
    """
    Create a morphism between FiniteAlgebras.

    INPUT:

    - ``A`` -- FiniteAlgebra

    - ``B`` -- FiniteAlgebra

    - ``f`` -- matrix of the underlying k-linear map

    - ``unitary`` -- boolean (default: True) - if True and ``check``
      is also True, raise a ValueError unless `A` and `B` are unitary
      and `f` respects unit elements

    - ``check`` -- boolean (default: True) - check whether the given
      k-linear map really defines a (not necessarily unitary)
      k-algebra homomorphism

    The algebras `A` and `B` must be defined over the same base field.

    EXAMPLES:

        sage: from sage.algebras.finite_algebras.finite_algebra_morphism import FiniteAlgebraMorphism
        sage: A = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
        sage: B = FiniteAlgebra(QQ, [Matrix([1])])
        sage: f = FiniteAlgebraMorphism(A, B, Matrix([[1], [0]]))
        sage: f.domain() is A
        True
        sage: f.codomain() is B
        True
        sage: f(A.basis()[0])
        e
        sage: f(A.basis()[1])
        0

    TODO: example illustrating unitary flag
    """

    def __init__(self, A, B, f, unitary=True, check=True):
        """
        TEST::

            sage: from sage.algebras.finite_algebras.finite_algebra_morphism import FiniteAlgebraMorphism
            sage: A = FiniteAlgebra(QQ, [Matrix([1])])
            sage: B = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: FiniteAlgebraMorphism(A, B, Matrix([[1, 0]]))
            Morphism from Finite algebra of degree 1 over Rational Field to Finite algebra of degree 2 over Rational Field given by matrix
            [1 0]
        """
        RingHomomorphism_im_gens.__init__(self, parent=Hom(A, B), im_gens=f.rows(), check=check)
        self._matrix = f
        if not (unitary and check):
            return
        if (not A.is_unitary()
            or not B.is_unitary()
            or self(A.one()) != B.one()):
            raise ValueError("homomorphism does not respect unit elements")

    def __call__(self, x):
        """
        TEST::

            sage: A = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.maximal_ideal()
            sage: q = A.quotient_map(I)
            sage: q(0) == 0 and q(1) == 1
            True
        """
        x = self.domain()(x)
        B = self.codomain()
        return B.element_class(B, x.vector() * self.matrix())

    def _repr_(self):
        """
        TEST::

            sage: A = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.maximal_ideal()
            sage: q = A.quotient_map(I)
            sage: q._repr_()
            'Morphism from Finite algebra of degree 2 over Rational Field to Finite algebra of degree 1 over Rational Field given by matrix\n[1]\n[0]'
        """
        return "Morphism from %s to %s given by matrix\n%s" \
            % (self.domain(), self.codomain(), self.matrix())

    def matrix(self):
        """
        Return the matrix of ``self``.

        EXAMPLE::

            sage: from sage.algebras.finite_algebras.finite_algebra_morphism import FiniteAlgebraMorphism
            sage: A = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: B = FiniteAlgebra(QQ, [Matrix([1])])
            sage: M = Matrix([[1], [0]])
            sage: f = FiniteAlgebraMorphism(A, B, M)
            sage: f.matrix() == M
            True
        """
        return self._matrix

    def inverse_image(self, I):
        """
        Return the inverse image of `I` under ``self``.

        INPUT:

        - ``I`` -- ``FiniteAlgebraIdeal``, an ideal of ``self.codomain()``

        OUTPUT:

        -- ``FiniteAlgebraIdeal``, the inverse image of `I` under ``self``.

        EXAMPLE::

            sage: A = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.maximal_ideal()
            sage: q = A.quotient_map(I)
            sage: B = q.codomain()
            sage: q.inverse_image(B.zero_ideal()) == I
            True
        """
        coker_I = I.basis_matrix().transpose().kernel().basis_matrix().transpose()
        return self.domain().ideal((self.matrix() * coker_I).kernel().basis_matrix(), given_by_matrix=True)
