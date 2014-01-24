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
#from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.categories.tensor import tensor
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.rings.infinity import infinity
from sage.rings.all import NN

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
    def __init__(self, M, prefix="St", category=None):
        r"""
        Initialize ``self``.

        EXAMPLES::
        """
        if category is None:
            category = [HopfAlgebrasWithBasis(R)]
            if alg in CommutativeAlgebras(R):
                category.append(CommutativeAlgebras(R))
            #if alg in GradedAlgebras(R):
            #    cat.append(GradedAlgebras(R))
            category = tuple(category)

        TensorModule.__init__(self, M, prefix, category)

    def _repr_(self):
        r"""
        Text representation of this stuffle algebra.

        EXAMPLES::

            sage: StuffleAlgebra(QQ)
            Stuffle Algebra with variables y over Rational Field

            sage: StuffleAlgebra(ZZ)
            Stuffle Algebra with variables y over Integer Ring
        """
        return "Stuffle Algebra of {}".format(self._base_algebra)

    @cached_method
    def one_basis(self):
        r"""
        Return the empty word, which index of `1` of this module.

        EXAMPLES::

            sage: A = StuffleAlgebra(QQ)
            sage: A.one_basis()
            word:
            sage: A.one()
            St[word: ]
        """
        return self._indices.one()

    def product_on_basis(self, w1, w2):
        r"""
        Return the product of basis elements ``w1`` and ``w2``, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis()`.

        INPUT:

        - ``w1``, ``w2`` -- basis elements

        EXAMPLES::

            sage: A = StuffleAlgebra(QQ)
            sage: W = A.basis().keys()
            sage: A.product_on_basis(W([W.alphabet()[0]]),W([W.alphabet()[1]]))
            St[word: y1,y2] + St[word: y2,y1] + St[word: y3]
        """
        return self._stuffle(w1, w2)

    def _stuffle(self, w1, w2):
        r"""
        Return a element of ``self`` obtained by the stuffle
        (quasi-shuffle) product of words ``w1`` and ``w2``.

        EXAMPLES::

            sage: from sage.algebras.stuffle_algebra import stuffle
            sage: A = build_alphabet(10, 'y')
            sage: W = Words(A)
            sage: w1 = W([W.alphabet()[0], W.alphabet()[0]])
            sage: w2 = W([W.alphabet()[1]])
            sage: stuffle(w1, w2)
            [word: y1,y1,y2, word: y1,y2,y1, word: y1,y3, word: y2,y1,y1, word: y3,y1]

            sage: A = build_alphabet(10,'g')
            sage: W = Words(A)
            sage: w1 = W([W.alphabet()[0],W.alphabet()[0]])
            sage: w2 = W([W.alphabet()[1]])
            sage: stuffle(w1, w2, A)
            [word: g1,g1,g2, word: g1,g2,g1, word: g1,g3, word: g2,g1,g1, word: g3,g1]
        """
        if len(w1) == 0:
            return self.monomial(w2)
        if len(w2) == 0:
            return self.monomial(w1)

        # Split off the first letter
        u = list(w1)
        v = list(w2)
        lu = u[0]
        lv = v[0]
        yi = lu[0]
        yj = lv[0]
        if u[0][1] == 1:
            u.pop(0)
        else:
            u[0] = (lu[0], lu[1] - 1)
        if v[0][1] == 1:
            v.pop(0)
        else:
            v[0] = (lv[0], lv[1] - 1)
        u = self._indices(u)
        v = self._indices(v)

        ret = self.sum_of_terms((yi*w, c) for w,c in self._stuffle(u, w2))
        ret += self.sum_of_terms((yj*w, c) for w,c in self._stuffle(w1, v))
        s = self._stuffle(u, v)

        yi = yi.leading_support()
        yj = yj.leading_support()
        p = self._base_algebra.monomial(yi) * self._base_algebra.monomial(yj)
        return ret + self.sum_of_terms((b*w, cs*cp) for w,cs in s for b,cp in p)

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

    def coproduct_on_basis(self, w):
        """
        Return the coproduct of the element of the basis indexed
        by the word ``w``.

        INPUT:

        - ``w`` -- a word

        EXAMPLES::

            sage: F = StuffleAlgebra(QQ)
            sage: F.coproduct_on_basis(Word(['y3']))
            St[word: ] # St[word: y3] + St[word: y1] # St[word: y2] + St[word: y2] # St[word: y1] + St[word: y3] # St[word: ]
            sage: F.coproduct_on_basis(Word(['y1','y2']))
            St[word: ] # St[word: y1,y2] + St[word: y1] # St[word: y1,y1] + St[word: y1] # St[word: y2] + St[word: y1,y1] # St[word: y1] + St[word: y1,y2] # St[word: ] + St[word: y2] # St[word: y1]
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

            sage: A = StuffleAlgebra(QQ); A
            Stuffle Algebra over Rational Field with variables y
            sage: A.algebra_generators()
            Lazy family (generator(i))_{i in Free monoid on 1 generators (y,)}
        """
        return Family(self._indices.indices(),
                      lambda i: self.monomial(self._indices(i)),
                      name='generator')

    gens = algebra_generators

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = StuffleAlgebra(QQ)
            sage: x, y = R.gens()[0:2]
            sage: x,y
            (St[word: y1], St[word: y2])
            sage: R(3)
            3*St[word: ]
            sage: R(x)
            St[word: y1]
        """
        P = x.parent()
        if isinstance(P, StuffleAlgebra):
            if P is self:
                return x
            if P is not self.base_ring():
                return self.element_class(self, x.monomial_coefficients())

        # ok, not a stuffle algebra element (or should not be viewed as one).
        if isinstance(x, basestring):
            from sage.misc.sage_eval import sage_eval
            return sage_eval(x, locals=self.gens_dict())

        R = self.base_ring()
        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.zero()
        return self.from_base_ring_from_one_basis(x)

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of ``x`` into ``self``.

        Here is what canonically coerces to ``self``:

        - this stuffle algebra,

        - anything that coerces to the base ring of this stuffle algebra,

        - any stuffle algebra on the same variables, whose base ring
          coerces to the base ring of this stuffle algebra.

        EXAMPLES::

            sage: F = StuffleAlgebra(GF(7)); F
            Stuffle Algebra over Finite Field of size 7 with variables y

        Elements of the stuffle algebra canonically coerce in::

            sage: x, y, z = F.gens()[0:3]
            sage: F.coerce(x*y) # indirect doctest
            St[word: y1,y2] + St[word: y2,y1] + St[word: y3]

        Elements of the integers coerce in, since there is a coerce map
        from `\ZZ` to GF(7)::

            sage: F.coerce(1)       # indirect doctest
            St[word: ]

        There is no coerce map from `\QQ` to `\GF{7}`::

            sage: F.coerce(2/3)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Stuffle Algebra over Finite Field of size 7 with variables y

        Elements of the base ring coerce in::

            sage: F.coerce(GF(7)(5))
            5*St[word: ]

        The stuffle algebra over `\ZZ` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = StuffleAlgebra(ZZ)
            sage: Gx,Gy,Gz = G.gens()[0:3]
            sage: z = F.coerce(Gx**2 * Gy);z
            2*St[word: y1,y1,y2] + 2*St[word: y1,y2,y1] + 2*St[word: y1,y3] + 2*St[word: y2,y1,y1] + 2*St[word: y2,y2] + 2*St[word: y3,y1] + St[word: y4]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the stuffle
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(x^3*y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Stuffle Algebra over Finite Field of size 7 with variables y to Stuffle Algebra over Integer Ring with variables y
        """
        try:
            R = x.parent()

            # stuffle algebras in the same variables over any base
            # that coerces in:
            if isinstance(R, StuffleAlgebra) and R.indices() == self.indices():
                if self.has_coerce_map_from(R.base_ring()):
                    return self(x)

                raise TypeError("no natural map between bases of stuffle algebras")

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this stuffle algebra.
        return self._coerce_try(x, [self.base_ring()])

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are:

        - Stuffle Algebras in the same variables over a base with a coercion
          map into ``self.base_ring()``.

        - Anything with a coercion into ``self.base_ring()``.

        TESTS::

            sage: F = StuffleAlgebra(ZZ)
            sage: G = StuffleAlgebra(QQ)
            sage: H = StuffleAlgebra(ZZ, 'g')
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            False
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            True
        """
        # stuffle algebras in the same variable over any base that coerces in:
        if isinstance(R, StuffleAlgebra) and R.indices() == self.indices():
            return self.base_ring().has_coerce_map_from(R.base_ring())

        return self.base_ring().has_coerce_map_from(R)

