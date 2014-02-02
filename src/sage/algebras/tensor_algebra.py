r"""
Tensor Algebras

AUTHORS:

- Travis Scrimshaw (2014-01-24): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2014 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.modules.tensor_module import TensorModule
from sage.categories.graded_hopf_algebras_with_basis import GradedHopfAlgebrasWithBasis
from sage.misc.cachefunc import cached_method

class TensorAlgebra(TensorModule):
    r"""
    The tensor algebra `T(M)` of a module `M`.

    Let `\{ b_i \}` be a basis of `M`, then the tensor algebra has a basis of
    `\{ b_{i_1} \otimes b_{i_2} \otimes \cdots \otimes b_{i_n} \}` and a
    product given by:

    .. MATH::

        (b_{i_1} \otimes \cdots \otimes b_{i_m}) \cdot (b_{j_1} \otimes
        \cdots \otimes b_{j_n}) = b_{i_1} \otimes \cdots \otimes b_{i_m}
        \otimesb_{j_1} \otimes \cdots \otimes b_{j_n}.
    """
    def __init__(self, M, prefix='T', category=None, **options):
        r"""
        Initialize ``self``.

        EXAMPLES::
        """
        R = M.base_ring()
        category = GradedHopfAlgebrasWithBasis(R).or_subcategory(category)
        TensorModule.__init__(self, M, prefix, category, **options)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

        """
        return "Tensor Algebra of {}".format(self._base_module)

    @cached_method
    def one_basis(self):
        r"""
        Return the empty word, which index of `1` of this module.

        EXAMPLES::

        """
        return self._indices.one()

    def product_on_basis(self, a, b):
        r"""
        Return the product of simple tensor elements ``a`` and ``b``, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis()`.

        INPUT:

        - ``a``, ``b`` -- simple tensors (i.e. basis elements)

        EXAMPLES::
        """
        return self.monomial(a * b)

    def counit(self, x):
        """
        Return the counit of ``x``.

        INPUT:

        - ``x`` -- an element of ``self``

        EXAMPLES::
        """
        #dic = x.monomial_coefficients()
        #if Word() not in dic:
        #    return 0
        #return dic[Word()]

    def coproduct_on_basis(self, m):
        """
        Return the coproduct of the element of the basis indexed
        by the simple tensor ``m``.

        INPUT:

        - ``m`` -- a simple tensor

        EXAMPLES::
        """
        if len(w) == 0:
            return self.tensor_square().one()

        if len(w) == 1:
            S = self.tensor_square()
            #A = self._alphabet
            #W = self.basis().keys()._construct_word
            #i = self._alphabet.index(w[0]) + 1
            #l = [(W([A[j-1]]), W([A[i-j-1]])) for j in range(1, i)] + [(w, W([])), (W([]), w)]
            #return S.sum_of_monomials(l)

        #B = self.basis()
        #result = self.coproduct_on_basis(Word([w[0]]))
        #for i in w[1:]:
        #    temp1 = self.coproduct_on_basis(Word([i]))
        #    temp2 = 0
        #    for ((u1, u2), coeff1) in list(temp1):
        #        for ((v1, v2), coeff2) in list(result):
        #            temp2 += coeff1 * coeff2 * tensor((St[Word(v1)*Word(u1)], St[Word(v2)*Word(u2)]))
        #    result = temp2
        #return result

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        EXAMPLES::
        """
        return Family(self._indices.indices(),
                      lambda i: self.monomial(self._indices(i)),
                      name='generator')

    gens = algebra_generators

