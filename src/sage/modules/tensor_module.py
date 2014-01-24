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

#from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.hopf_algebras_with_basis import HopfAlgebrasWithBasis
from sage.categories.tensor import tensor
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.algebras.free_algebra import FreeAlgebra
from sage.combinat.free_module import CombinatorialFreeModule
from sage.monoids.indexed_monoid import IndexedFreeMonoid
from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.rings.infinity import infinity
from sage.rings.all import NN

class TensorAlgebra(CombinatorialFreeModule):
    r"""
    The tensor module `T(M)` of a module `M`.

    Let `\{ b_i \}` be a basis of `M`, then the tensor module is the span of
    `\{ b_{i_1} \otimes b_{i_2} \otimes \cdots \otimes b_{i_k} \}`.
    """
    def __init__(self, M, prefix="T", category=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: F = StuffleAlgebra(QQ)
            sage: TestSuite(F).run()
        """
        self._base_module = M
        R = M.base_ring()

        if category is None:
            category = [HopfAlgebrasWithBasis(R)]
            if alg in CommutativeAlgebras(R):
                category.append(CommutativeAlgebras(R))
            #if alg in GradedAlgebras(R):
            #    cat.append(GradedAlgebras(R))
            category = tuple(category)

        CombinatorialFreeModule.__init__(self, R, IndexedFreeMonoid(M.indices()),
            prefix=prefix, latex_prefix="", category=category)

    def _repr_(self):
        r"""
        Text representation of ``self``.

        EXAMPLES::
        """
        return "Tensor Module of {}".format(self._base_module)

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

