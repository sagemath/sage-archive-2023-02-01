r"""
FiniteDimensionalAlgebrasWithBasis
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from category_types import Category_over_base_ring
from sage.categories.all import FiniteDimensionalModulesWithBasis, AlgebrasWithBasis
from sage.misc.cachefunc import cached_method

class FiniteDimensionalAlgebrasWithBasis(Category_over_base_ring):
    """
    The category of finite dimensional algebras with a distinguished basis

    EXAMPLES::

      sage: FiniteDimensionalAlgebrasWithBasis(QQ)
      Category of finite dimensional algebras with basis over Rational Field
      sage: FiniteDimensionalAlgebrasWithBasis(QQ).super_categories()
      [Category of finite dimensional modules with basis over Rational Field, Category of algebras with basis over Rational Field]

    TESTS::

        sage: TestSuite(FiniteDimensionalAlgebrasWithBasis(ZZ)).run()
    """

    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: FiniteDimensionalAlgebrasWithBasis(QQ).super_categories()
            [Category of finite dimensional modules with basis over Rational Field, Category of algebras with basis over Rational Field]
        """
        R = self.base_ring()
        return [FiniteDimensionalModulesWithBasis(R), AlgebrasWithBasis(R)]

    class ParentMethods:
        pass

    class ElementMethods:
        def on_left_matrix(self, new_BR = None):
            """
            Returns the matrix of the action of self on the algebra my
            multiplication on the left

            If new_BR is specified, then the matrix will be over new_BR.

            TODO: split into to parts
             - build the endomorphism of multiplication on the left
             - build the matrix of an endomorphism

            EXAMPLES::

                sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
                sage: a = QS3([2,1,3])
                sage: a.on_left_matrix()
                [0 0 1 0 0 0]
                [0 0 0 0 1 0]
                [1 0 0 0 0 0]
                [0 0 0 0 0 1]
                [0 1 0 0 0 0]
                [0 0 0 1 0 0]
                sage: a.on_left_matrix(RDF)
                [0.0 0.0 1.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 0.0 1.0 0.0]
                [1.0 0.0 0.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 0.0 0.0 1.0]
                [0.0 1.0 0.0 0.0 0.0 0.0]
                [0.0 0.0 0.0 1.0 0.0 0.0]

            AUTHOR: Mike Hansen
            """
            parent = self.parent()

            if parent.get_order() is None:
                cc = parent._combinatorial_class
            else:
                cc = parent.get_order()

            BR = parent.base_ring()
            if new_BR is None:
                new_BR = BR

            from sage.matrix.all import MatrixSpace
            MS = MatrixSpace(new_BR, parent.dimension(), parent.dimension())
            l = [ (self*parent(m)).to_vector() for m in cc ]
            return MS(l).transpose()

        _matrix_ = on_left_matrix  # For temporary backward compatibility
        to_matrix = on_left_matrix
