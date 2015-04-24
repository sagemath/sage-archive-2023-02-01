r"""
Semisimple Algebras
"""
#*****************************************************************************
#  Copyright (C) 2011-2015 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.bindable_class import BoundClass
from category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.misc.cachefunc import cached_method
from algebras import Algebras


class SemisimpleAlgebras(Category_over_base_ring):
    """
    The category of semisimple algebras over a given base ring.

    EXAMPLES::

        sage: from sage.categories.semisimple_algebras import SemisimpleAlgebras
        sage: C = SemisimpleAlgebras(QQ); C
        Category of semisimple algebras over Rational Field

    This category is best constructed as::

        sage: D = Algebras(QQ).Semisimple(); D
        Category of semisimple algebras over Rational Field
        sage: D is C
        True

        sage: C.super_categories()
        [Category of algebras over Rational Field]

    Typically, finite group algebras are semisimple::

        sage: DihedralGroup(5).algebra(QQ) in SemisimpleAlgebras
        True

    Unless the characteristic of the field divides the order of the group::

        sage: DihedralGroup(5).algebra(IntegerModRing(5)) in SemisimpleAlgebras
        False

        sage: DihedralGroup(5).algebra(IntegerModRing(7)) in SemisimpleAlgebras
        True

    .. SEEALSO:: `<http://en.wikipedia.org/wiki/Semisimple_algebra>`_

    TESTS::

        sage: TestSuite(C).run()
    """
    @staticmethod
    def __classget__(cls, base_category, base_category_class):
        """
        Implement the shorthand ``Algebras(K).Semisimple()`` for ``SemisimpleAlgebras(K)``.

        This magic mimics the syntax of axioms for a smooth transition
        if ``Semisimple`` becomes one.

        EXAMPLES::

            sage: Algebras(QQ).Semisimple()
            Category of semisimple algebras over Rational Field
            sage: Algebras.Semisimple
            <class 'sage.categories.semisimple_algebras.SemisimpleAlgebras'>
        """
        if base_category is None:
            return cls
        return BoundClass(cls, base_category.base_ring())

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: Algebras(QQ).Semisimple().super_categories()
            [Category of algebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R)]

    class ParentMethods:

        def radical_basis(self, **keywords):
            r"""
            Return a basis of the Jacobson radical of this algebra.

            - ``keywords`` -- for compatibility; ignored.

            OUTPUT: the empty list since this algebra is semisimple.

            EXAMPLES::

                sage: A = SymmetricGroup(4).algebra(QQ)
                sage: A.radical_basis()
                []

            TESTS::

                sage: A.radical_basis.__module__
                'sage.categories.finite_dimensional_algebras_with_basis'
            """
            return []

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):

        class WithBasis(CategoryWithAxiom_over_base_ring):

            class ParentMethods:

                @cached_method
                def central_orthogonal_idempotents(self):
                    r"""
                    Return a maximal list of central orthogonal
                    idempotents of ``self``.

                    *Central orthogonal idempotents* of an algebra `A`
                    are idempotents `(e_1, \dots, e_n)` in the center
                    of `A` such that `e_i e_j = 0` whenever `i \neq j`.

                    With the maximality condition, they sum up to `1`
                    and are uniquely determined (up to order).

                    INPUT:

                    - ``self`` -- a semisimple algebra.

                    EXAMPLES:

                    For the algebra of the symmetric group `S_3`, we
                    recover the sum and alternating sum of all
                    permutations, together with a third idempotent::

                        sage: A3 = SymmetricGroup(3).algebra(QQ)
                        sage: idempotents = A3.central_orthogonal_idempotents()
                        sage: idempotents
                        [2/3*() - 1/3*(1,2,3) - 1/3*(1,3,2),
                         1/6*() + 1/6*(2,3) + 1/6*(1,2) + 1/6*(1,2,3) + 1/6*(1,3,2) + 1/6*(1,3),
                         1/6*() - 1/6*(2,3) - 1/6*(1,2) + 1/6*(1,2,3) + 1/6*(1,3,2) - 1/6*(1,3)]
                        sage: A3.is_identity_decomposition_into_orthogonal_idempotents(idempotents)
                        True

                    For the semisimple quotient of a quiver algebra,
                    we recover the vertices of the quiver::

                        sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                        An example of a finite dimensional algebra with basis:
                        the path algebra of the Kronecker quiver (containing
                        the arrows a:x->y and b:x->y) over Rational Field
                        sage: Aquo = A.semisimple_quotient()
                        sage: Aquo.central_orthogonal_idempotents()
                        [B['y'], B['x']]
                    """
                    return [x.lift()
                            for x in self.center().central_orthogonal_idempotents()]


            class Commutative(CategoryWithAxiom_over_base_ring):

                class ParentMethods:

                    @cached_method
                    def _orthogonal_decomposition(self, generators=None):
                        r"""
                        Return a list of orthogonal quasi-idempotents of this
                        semisimple commutative finite dimensional algebra.

                        INPUT:

                        - ``generators`` -- a list of generators of
                          ``self`` (default: the basis of ``self``)

                        OUTPUT:

                        - list of elements of ``self`` each generating a one
                          dimensional simple submodule summand of ``self``. The
                          list is maximal in the sense that no quasi-idempotent
                          `e` can be decomposed as a sum `e = e_1 + e_2` of
                          quasi-idempotents elements.

                        ALGORITHM:

                        Thanks to Schur's Lemma, a commutative
                        semisimple algebra `A` is a direct sum of
                        dimension 1 subalgebras. The algorithm is
                        recursive and proceeds as follows:

                        0. If `A` is of dimension 1, return a non zero
                           element.

                        1. Otherwise: find one of the generators such
                           that the morphism `x \mapsto ax` has at
                           least two (right) eigenspaces.

                        2. Decompose both eigenspaces recursively.

                        EXAMPLES::

                            sage: G5 = SymmetricGroup(5)
                            sage: A5 = G5.algebra(QQ)
                            sage: Z5 = A5.center()
                            sage: Z5._orthogonal_decomposition()
                            [B[0] - 1/3*B[2] + 1/6*B[6], B[0] + B[1] + B[2] +
                            B[3] + B[4] + B[5] + B[6], B[0] - B[1] + B[2] + B[3]
                            - B[4] - B[5] + B[6], B[0] + 1/5*B[1] + 1/5*B[2] -
                            1/5*B[3] + 1/5*B[4] - 1/5*B[5], B[0] - 1/5*B[1] +
                            1/5*B[2] - 1/5*B[3] - 1/5*B[4] + 1/5*B[5], B[0] +
                            1/2*B[1] + 1/4*B[3] - 1/4*B[4] - 1/4*B[6], B[0] -
                            1/2*B[1] + 1/4*B[3] + 1/4*B[4] - 1/4*B[6]]

                        .. TODO::

                            Improve speed by using matrix operations
                            only, if possible.
                        """
                        category = Algebras(self.base_ring()).Semisimple().WithBasis().FiniteDimensional().Commutative().Subobjects()
                        if generators is None:
                            generators = self.basis().list()
                        if self.dimension() == 1:
                            return self.basis().list()

                        def rec(space, generators=None):
                            if generators is None:
                                generators = space.basis().list()
                            generators = list(generators)

                            if space.dimension() == 1:
                                return space.basis().list()
                            eigenspaces = []

                            # Searching for a good generator ...
                            while len(eigenspaces) < 2:
                                if generators == []:
                                    raise Exception("Unable to fully decompose...")
                                gen = generators.pop()
                                phi = space.module_morphism(on_basis=lambda i:
                                        gen*space.basis()[i],
                                        codomain=space,
                                        triangular='lower')
                                eigenspaces = phi.matrix().eigenspaces_right()

                            # Gotcha! Let's split the algebra ...
                            subalgebras = \
                            [space.submodule(map(space.from_vector, eigenspace.basis()),
                                             category=category)
                             for eigenvalue, eigenspace in eigenspaces]

                            # Decompose recursively each eigenspace
                            return [idempotent.lift()
                                    for subalgebra in subalgebras
                                    for idempotent in rec(subalgebra)]

                        return rec(self, generators)

                    @cached_method
                    def central_orthogonal_idempotents(self):
                        r"""
                        Return the central orthogonal idempotents of
                        this semisimple commutative algebra.

                        Those idempotents form a maximal decomposition
                        of the identity into primitive orthogonal
                        idempotents.

                        OUTPUT:

                        - a list of orthogonal idempotents of ``self``.

                        EXAMPLES::

                            sage: A4 = SymmetricGroup(4).algebra(QQ)
                            sage: Z4 = A4.center()
                            sage: idempotents = Z4.central_orthogonal_idempotents()
                            sage: idempotents
                            [1/24*B[0] + 1/24*B[1] + 1/24*B[2] + 1/24*B[3] + 1/24*B[4],
                             3/8*B[0] - 1/8*B[1] - 1/8*B[2] + 1/8*B[4],
                             1/6*B[0] + 1/6*B[2] - 1/12*B[3],
                             3/8*B[0] + 1/8*B[1] - 1/8*B[2] - 1/8*B[4],
                             1/24*B[0] - 1/24*B[1] + 1/24*B[2] + 1/24*B[3] - 1/24*B[4]]

                        Lifting those idempotents from the center, we
                        recognize among them the sum and alternating
                        sum of all permutations::

                            sage: [e.lift() for e in idempotents]
                            [1/24*() + 1/24*(3,4) + 1/24*(2,3) + 1/24*(2,3,4) + 1/24*(2,4,3)
                             + 1/24*(2,4) + 1/24*(1,2) + 1/24*(1,2)(3,4) + 1/24*(1,2,3)
                             + 1/24*(1,2,3,4) + 1/24*(1,2,4,3) + 1/24*(1,2,4) + 1/24*(1,3,2)
                             + 1/24*(1,3,4,2) + 1/24*(1,3) + 1/24*(1,3,4) + 1/24*(1,3)(2,4)
                             + 1/24*(1,3,2,4) + 1/24*(1,4,3,2) + 1/24*(1,4,2) + 1/24*(1,4,3)
                             + 1/24*(1,4) + 1/24*(1,4,2,3) + 1/24*(1,4)(2,3),
                             ...,
                             1/24*() - 1/24*(3,4) - 1/24*(2,3) + 1/24*(2,3,4) + 1/24*(2,4,3)
                             - 1/24*(2,4) - 1/24*(1,2) + 1/24*(1,2)(3,4) + 1/24*(1,2,3)
                             - 1/24*(1,2,3,4) - 1/24*(1,2,4,3) + 1/24*(1,2,4) + 1/24*(1,3,2)
                             - 1/24*(1,3,4,2) - 1/24*(1,3) + 1/24*(1,3,4) + 1/24*(1,3)(2,4)
                             - 1/24*(1,3,2,4) - 1/24*(1,4,3,2) + 1/24*(1,4,2) + 1/24*(1,4,3)
                             - 1/24*(1,4) - 1/24*(1,4,2,3) + 1/24*(1,4)(2,3)]

                        We check that they indeed form a decomposition
                        of the identity of `Z_4` into orthogonal idempotents::

                            sage: Z4.is_identity_decomposition_into_orthogonal_idempotents(idempotents)
                            True
                        """
                        return [(e.leading_coefficient()/(e*e).leading_coefficient())*e
                                for e in self._orthogonal_decomposition()]
