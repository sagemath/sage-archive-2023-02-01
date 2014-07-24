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
from algebras import Algebras
from sage.misc.cachefunc import cached_method
import operator

from sage.rings.integer_ring import ZZ

class SemisimpleAlgebras(Category_over_base_ring):
    """
    The category of semisimple algebras over a given base ring.

    EXAMPLES::

        sage: SemisimpleAlgebras(QQ)
        Category of semisimple algebras over Rational Field
        sage: SemisimpleAlgebras(QQ).super_categories()
        [Category of algebras over Rational Field]

    Typically, finite group algebras are semisimple::

        sage: DihedralGroup(5).algebra(QQ) in SemisimpleAlgebras
        True

    Unless the characteristic of the field divides the order of the group::

        sage: DihedralGroup(5).algebra(IntegerModRing(5)) in SemisimpleAlgebras
        False

        sage: DihedralGroup(5).algebra(IntegerModRing(7)) in SemisimpleAlgebras # todo: not implemented
        True

    .. seealso:: `<http://en.wikipedia.org/wiki/Semisimple_algebra>`_

    TESTS::

        sage: TestSuite(SemisimpleAlgebras(QQ)).run()
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

    class Commutative(CategoryWithAxiom_over_base_ring):

        from sage.rings.integer_ring import ZZ
        from sage.categories.associative_algebras import AssociativeAlgebras

        class ParentMethods:

            def semi_simple_commutative_decomposition(self, listGen=None, topLevel=True):
                r"""
                Decompose a commutative semi-simple algebra ``A`` into a direct sum of simple A-modules.

                Return a list of generators of simple A-modules.

                EXAMPLES:

                    sage: A5 = SymmetricGroupAlgebra(QQ,5)
                    sage: Z5 = A5.center()                          
                    sage: semi_simple_commutative_decomposition(Z5)
                    [B[0] - 1/3*B[3] + 1/6*B[6],
                    B[0] + B[1] + B[2] + B[3] + B[4] + B[5] + B[6],
                    B[0] - B[1] + B[2] + B[3] - B[4] - B[5] + B[6],
                    B[0] + 1/2*B[1] + 1/4*B[2] - 1/4*B[5] - 1/4*B[6],
                    B[0] - 1/2*B[1] + 1/4*B[2] + 1/4*B[5] - 1/4*B[6],
                    B[0] + 1/5*B[1] - 1/5*B[2] + 1/5*B[3] - 1/5*B[4] + 1/5*B[5],
                    B[0] - 1/5*B[1] - 1/5*B[2] + 1/5*B[3] + 1/5*B[4] - 1/5*B[5]]
                """
                #Terminal case and stuffs
                if listGen==None:
                    listGen = self.basis().list()
                if self.dimension() == 1:
                    if topLevel:
                        return self.basis().list()
                    else:
                        return [x.lift() for x in self.basis()]

                #Searching for a good generator...
                res = []
                B = self.basis()
                while len(res)<2:
                    if listGen==[]:
                        raise Exception("Unable to fully decompose...")
                    curGen = listGen[ZZ.random_element(len(listGen))]
                    listGen.remove(curGen)
                    phi = self.module_morphism(on_basis=lambda i:
                            curGen*B[i],
                            codomain=self,
                            triangular=True)
                    # phi = self.module_morphism(on_basis=lambda i:
                    #         curGen*self.monomial(i),
                    #         codomain=self,
                    #         triangular=True)
                    return phi
                    aMat = phi.matrix(self.base_ring())
                    res = aMat.eigenspaces_right()

                #Gotcha! Let's settle the algebra...
                
                res = [vector_space for _,vector_space in res]
                res = [[self.from_vector(vector) for vector in eigenspace.basis()]
                        for eigenspace in res]
                
                decomp = [self.submodule(v, category=AssociativeAlgebras(self.base_ring()).WithBasis().FiniteDimensional().Subobjects()) for v in res]

                #Recursive on decomp
                res = [x for space in decomp for x in space.semi_simple_commutative_decomposition(topLevel=False)]
                if topLevel:
                    return res
                else:
                    return map( lambda x: x.lift(), res)
