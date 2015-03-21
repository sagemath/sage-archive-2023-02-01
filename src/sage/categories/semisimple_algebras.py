r"""
Semisimple Algebras
"""
#*****************************************************************************
#  Copyright (C) 2011 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.category_with_axiom import CategoryWithAxiom_over_base_ring
from sage.misc.cachefunc import cached_method
from algebras import Algebras


class SemisimpleAlgebras(Category_over_base_ring):
    """
    The category of semisimple algebras over a given base ring.

    EXAMPLES::

        from sage.categories.semisimple_algebras import SemisimpleAlgebras
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
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: SemisimpleAlgebras(QQ).super_categories()
            [Category of algebras over Rational Field]
        """
        R = self.base_ring()
        return [Algebras(R)]

    class FiniteDimensional(CategoryWithAxiom_over_base_ring):

        class WithBasis(CategoryWithAxiom_over_base_ring):

            class ParentMethods:

                @cached_method
                def orthogonal_idempotents(self):
                    r"""
                    Return a maximal list of orthogonal idempotents of
                    ``self``.

                    INPUT:

                    - ``self`` -- semisimple algebra

                    EXAMPLES::

                        sage: A3 = SymmetricGroup(3).algebra(QQ)
                        sage: A3.orthogonal_idempotents()
                        [2/3*() - 1/3*(1,2,3) - 1/3*(1,3,2), 1/6*() + 1/6*(2,3)
                        + 1/6*(1,2) + 1/6*(1,2,3) + 1/6*(1,3,2) + 1/6*(1,3),
                        1/6*() - 1/6*(2,3) - 1/6*(1,2) + 1/6*(1,2,3) +
                        1/6*(1,3,2) - 1/6*(1,3)]


                    ::

                        sage: A = FiniteDimensionalAlgebrasWithBasis(QQ).example(); A
                        An example of a finite dimensional algebra with basis:
                        the path algebra of the Kronecker quiver (containing
                        the arrows a:x->y and b:x->y) over Rational Field 
                        sage: Aquo = A.semisimple_quotient()
                        sage: Aquo.orthogonal_idempotents()
                        [B['y'], B['x']]
                    """
                    return [x.lift()
                            for x in self.center().orthogonal_idempotents()]

            class Commutative(CategoryWithAxiom_over_base_ring):

                class ParentMethods:

                    @cached_method
                    def _orthogonal_decomposition(self, generators=None):
                        r"""
                        Return a list of orthogonal idempotents of a semisimple
                        commutative finite dimensional algebra ``self``.

                        INPUT:

                        - ``self`` a finite dimensional semisimple commutative
                          algebra.
                        - ``generators`` a list of generators of ``self``. By
                          default it will be the basis of ``self``.

                        OUTPUT:

                        - list of elements of ``self`` each generating a one
                          dimensional simple submodule of ``self`` in direct
                          sum with the others. The list is maximal.

                        ALGORITHM:

                            A commutative semisimple algebra `A` is a direct
                            sum of dimension 1 sub-algebras thanks to Schur
                            Lemma. The algorithm is recursive a proceed in two
                            steps:

                            0. If `A` is of dimension 1, return a non zero
                            element.
                            1. Find a generator `a \in A` such that the
                            morphism `x \mapsto ax` has at least two (right)
                            eigenspaces.
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

                            Improve speed by only using matrix operations, if
                            possible.
                        """
                        # Main program
                        if generators is None:
                            generators = self.basis().list()
                        if self.dimension() == 1:
                            return self.basis().list()

                        def rec(space, generators):
                            if space.dimension() == 1:
                                return space.basis().list()
                            eigenspaces = []

                            # Searching a good generator...
                            while len(eigenspaces) < 2:
                                if generators == []:
                                    raise Exception("Unable to fully decompose...")
                                gen = generators.pop()
                                phi = space.module_morphism(on_basis=lambda i:
                                        gen*space.basis()[i],
                                        codomain=space,
                                        triangular=True)
                                eigenspaces = phi.matrix(space.base_ring()).eigenspaces_right()

                            # Gotcha! Let's settle the algebra...
                            eigenspaces = [vector_space for _, vector_space in eigenspaces]
                            eigenspaces = [[space.from_vector(vector)
                                for vector in eigenspace.basis()]
                                for eigenspace in eigenspaces]
                            decomp = [space.submodule(v,
                                category=SemisimpleAlgebras(space.base_ring()).WithBasis().FiniteDimensional().Commutative().Subobjects())
                                for v in eigenspaces]

                            # Decompose recursively each eigenspace
                            return map(lambda x: x.lift(), [x for eigenspace in
                                decomp for x in rec(eigenspace, eigenspace.basis().list())])

                        return rec(self, generators)

                    @cached_method
                    def orthogonal_idempotents(self):
                        r"""
                        Return the minimum set of central orthogonal idempotents of ``self``.

                        INPUT:

                        - ``self`` a commutative semisimple algebra

                        OUTPUT:

                        - list of idempotents of ``self``

                        EXAMPLES::

                            sage: A5 = SymmetricGroup(5).algebra(QQ)
                            sage: Z5 = A5.center()
                            sage: orth = Z5.orthogonal_idempotents()
                            sage: orth
                            [3/10*B[0] - 1/10*B[2] + 1/20*B[6], 1/120*B[0] +
                            1/120*B[1] + 1/120*B[2] + 1/120*B[3] + 1/120*B[4] +
                            1/120*B[5] + 1/120*B[6], 1/120*B[0] - 1/120*B[1] +
                            1/120*B[2] + 1/120*B[3] - 1/120*B[4] - 1/120*B[5] +
                            1/120*B[6], 5/24*B[0] + 1/24*B[1] + 1/24*B[2] -
                            1/24*B[3] + 1/24*B[4] - 1/24*B[5], 5/24*B[0] -
                            1/24*B[1] + 1/24*B[2] - 1/24*B[3] - 1/24*B[4] +
                            1/24*B[5], 2/15*B[0] + 1/15*B[1] + 1/30*B[3] -
                            1/30*B[4] - 1/30*B[6], 2/15*B[0] - 1/15*B[1] +
                            1/30*B[3] + 1/30*B[4] - 1/30*B[6]]
                            sage: orth[2] * orth[4]
                            0
                            sage: orth[1] ** 2 == orth[1]
                            True
                        """
                        return [(e.leading_coefficient()/(e*e).leading_coefficient())*e for
                            e in self._orthogonal_decomposition()]
