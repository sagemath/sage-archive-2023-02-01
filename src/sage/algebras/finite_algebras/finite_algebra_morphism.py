"""
Morphisms Between Finite Algebras
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
from sage.rings.homset import RingHomset_generic
from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix

class FiniteAlgebraMorphism(RingHomomorphism_im_gens):
    """
    Create a morphism between two :class:`finite algebras <FiniteAlgebra>`.

    INPUT:

    - ``A``, ``B`` -- finite algebras

    - ``f`` -- matrix of the underlying `k`-linear map

    - ``unitary`` -- boolean (default: ``True``); if ``True`` and ``check``
      is also ``True``, raise a ``ValueError`` unless ``A`` and ``B`` are
      unitary and ``f`` respects unit elements

    - ``check`` -- boolean (default: ``True``); check whether the given
      k-linear map really defines a (not necessarily unitary)
      k-algebra homomorphism

    The algebras ``A`` and ``B`` must be defined over the same base field.

    EXAMPLES::

        sage: from sage.algebras.finite_algebras.finite_algebra_morphism import FiniteAlgebraMorphism
        sage: A = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
        sage: B = FiniteAlgebra(QQ, [Matrix([1])])
        sage: H = Hom(A, B)
        sage: f = FiniteAlgebraMorphism(H, Matrix([[1], [0]]))
        sage: f.domain() is A
        True
        sage: f.codomain() is B
        True
        sage: f(A.basis()[0])
        e
        sage: f(A.basis()[1])
        0

    .. TODO:: An example illustrating unitary flag.
    """
    def __init__(self, parent, im_gens, check=True, unitary=True, f=None):
        """
        TEST::

            sage: from sage.algebras.finite_algebras.finite_algebra_morphism import FiniteAlgebraMorphism
            sage: A = FiniteAlgebra(QQ, [Matrix([1])])
            sage: B = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: H = Hom(A, B)
            sage: phi = FiniteAlgebraMorphism(H, Matrix([[1, 0]]))
            sage: TestSuite(phi).run(skip="_test_category") # Currently ring morphisms are not using the category framework
        """
        A = parent.domain()
        B = parent.codomain()
        if not (unitary and check) and (not A.is_unitary()
                                        or not B.is_unitary()
                                        or self(A.one()) != B.one()):
            raise ValueError("homomorphism does not respect unit elements")

        if f is None:
            if is_Matrix(im_gens):
                f = im_gens
                im_gens = f.rows()
            else:
                f = matrix(im_gens)

        RingHomomorphism_im_gens.__init__(self, parent=parent, im_gens=im_gens, check=check)
        self._matrix = f

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
        return B.element_class(B, x.vector() * self._matrix)

    def _repr_(self):
        """
        TEST::

            sage: A = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: I = A.maximal_ideal()
            sage: q = A.quotient_map(I)
            sage: q._repr_()
            'Morphism from Finite algebra of degree 2 over Rational Field to Finite algebra of degree 1 over Rational Field given by matrix\n[1]\n[0]'
        """
        return "Morphism from {} to {} given by matrix\n{}".format(
                self.domain(), self.codomain(), self._matrix)

    def matrix(self):
        """
        Return the matrix of ``self``.

        EXAMPLE::

            sage: from sage.algebras.finite_algebras.finite_algebra_morphism import FiniteAlgebraMorphism
            sage: A = FiniteAlgebra(QQ, [Matrix([[1, 0], [0, 1]]), Matrix([[0, 1], [0, 0]])])
            sage: B = FiniteAlgebra(QQ, [Matrix([1])])
            sage: M = Matrix([[1], [0]])
            sage: H = Hom(A, B)
            sage: f = FiniteAlgebraMorphism(H, M)
            sage: f.matrix() == M
            True
        """
        return self._matrix

    def inverse_image(self, I):
        """
        Return the inverse image of ``I`` under ``self``.

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
        return self.domain().ideal((self._matrix * coker_I).kernel().basis_matrix(), given_by_matrix=True)

class FiniteAlgebraHomset(RingHomset_generic):
    """
    Set of morphisms between two finite algebras.
    """
    def __call__(self, im_gens, check=True, unitary=True):
        """
        Construct a homomorphism.
        """
        if is_Matrix(im_gens):
            return FiniteAlgebraMorphism(self, im_gens.rows(), check=check, unitary=unitary, f=im_gens)
        if isinstance(im_gens, FiniteAlgebraMorphism):
            if im_gens.parent() == self:
                return FiniteAlgebraMorphism(self, im_gens.im_gens(), check=check, unitary=unitary, f=im_gens._matrix)
        try:
            return FiniteAlgebraMorphism(self, im_gens.rows(), check=check, unitary=unitary, f=im_gens)
        except (NotImplementedError, ValueError):
            return RingHomset_generic.__call__(self, im_gens, check)

