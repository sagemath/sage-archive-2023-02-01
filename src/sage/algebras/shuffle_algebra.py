# -*- coding: utf-8 -*-
r"""
Shuffle algebras

AUTHORS:

- Frédéric Chapoton (2013-03)
"""

#*****************************************************************************
#  Copyright (C) 2013 Frédéric Chapoton <chapoton-math-univ-lyon1-fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.rings import Rings
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.words.words import Words
from sage.misc.cachefunc import cached_method
from sage.sets.family import Family

class ShuffleAlgebra(CombinatorialFreeModule):
    r"""
    The shuffle algebra on some generators over a base ring.

    Shuffle algebras are commutative and associative algebras, with a
    basis indexed by words. The product of two words `w_1 \cdot w_2` is given
    by the sum over the shuffle product of `w_1` and `w_2`.

    .. SEEALSO::

        For more on shuffle products, see
        :mod:`~sage.combinat.words.shuffle_product` and
        :meth:`~sage.combinat.words.finite_word.FiniteWord_class.shuffle()`.

    REFERENCES:

    - :wikipedia:`Shuffle algebra`

    INPUT:

    -  ``R`` -- ring

    -  ``names`` -- generator names (string)

    EXAMPLES::

        sage: F = ShuffleAlgebra(QQ, 'xyz'); F
        Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Rational Field

        sage: mul(F.gens())
        B[word: xyz] + B[word: xzy] + B[word: yxz] + B[word: yzx] + B[word: zxy] + B[word: zyx]

        sage: mul([ F.gen(i) for i in range(2) ]) + mul([ F.gen(i+1) for i in range(2) ])
        B[word: xy] + B[word: yx] + B[word: yz] + B[word: zy]

        sage: S = ShuffleAlgebra(ZZ, 'abcabc'); S
        Shuffle Algebra on 3 generators ['a', 'b', 'c'] over Integer Ring
        sage: S.base_ring()
        Integer Ring

        sage: G = ShuffleAlgebra(S, 'mn'); G
        Shuffle Algebra on 2 generators ['m', 'n'] over Shuffle Algebra on 3 generators ['a', 'b', 'c'] over Integer Ring
        sage: G.base_ring()
        Shuffle Algebra on 3 generators ['a', 'b', 'c'] over Integer Ring

    Shuffle algebras commute with their base ring::

        sage: K = ShuffleAlgebra(QQ,'ab')
        sage: a,b = K.gens()
        sage: K.is_commutative()
        True
        sage: L = ShuffleAlgebra(K,'cd')
        sage: c,d = L.gens()
        sage: L.is_commutative()
        True
        sage: s = a*b^2 * c^3; s
        (12*B[word:abb]+12*B[word:bab]+12*B[word:bba])*B[word: ccc]
        sage: parent(s)
        Shuffle Algebra on 2 generators ['c', 'd'] over Shuffle Algebra on 2 generators ['a', 'b'] over Rational Field
        sage: c^3 * a * b^2
        (12*B[word:abb]+12*B[word:bab]+12*B[word:bba])*B[word: ccc]

    Shuffle algebras are commutative::

        sage: c^3 * b * a * b == c * a * c * b^2 * c
        True

    We can also manipulate elements in the basis and coerce elements from our
    base field::

        sage: F = ShuffleAlgebra(QQ, 'abc')
        sage: B = F.basis()
        sage: B[Word('bb')] * B[Word('ca')]
        B[word: bbca] + B[word: bcab] + B[word: bcba] + B[word: cabb] + B[word: cbab] + B[word: cbba]
        sage: 1 - B[Word('bb')] * B[Word('ca')] / 2
        B[word: ] - 1/2*B[word: bbca] - 1/2*B[word: bcab] - 1/2*B[word: bcba] - 1/2*B[word: cabb] - 1/2*B[word: cbab] - 1/2*B[word: cbba]
    """
    def __init__(self, R, names):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: F = ShuffleAlgebra(QQ, 'xyz'); F
            Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Rational Field
            sage: TestSuite(F).run()

        TESTS::

            sage: ShuffleAlgebra(24, 'toto')
            Traceback (most recent call last):
            ...
            TypeError: argument R must be a ring
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        self._alphabet = Alphabet(names)
        self.__ngens = self._alphabet.cardinality()
        CombinatorialFreeModule.__init__(self, R, Words(names),
            latex_prefix = "",
            category = (AlgebrasWithBasis(R), CommutativeAlgebras(R)))

    def variable_names(self):
        r"""
        Return the names of the variables.

        EXAMPLES::

            sage: R = ShuffleAlgebra(QQ,'xy')
            sage: R.variable_names()
            {'x', 'y'}
        """
        return self._alphabet

    def is_commutative(self):
        r"""
        Return ``True`` as the shuffle algebra is commutative.

        EXAMPLES::

            sage: R = ShuffleAlgebra(QQ,'x')
            sage: R.is_commutative()
            True
            sage: R = ShuffleAlgebra(QQ,'xy')
            sage: R.is_commutative()
            True
        """
        return True

    def _repr_(self):
        r"""
        Text representation of this shuffle algebra.

        EXAMPLES::

            sage: F = ShuffleAlgebra(QQ,'xyz')
            sage: F  # indirect doctest
            Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Rational Field

            sage: ShuffleAlgebra(ZZ,'a')
            Shuffle Algebra on one generator ['a'] over Integer Ring
        """
        if self.__ngens == 1:
            gen = "one generator"
        else:
            gen = "%s generators" %self.__ngens
        return "Shuffle Algebra on "+ gen +" %s over %s"%(
            self._alphabet.list(), self.base_ring())

    @cached_method
    def one_basis(self):
        r"""
        Return the empty word, which index of `1` of this algebra,
        as per :meth:`AlgebrasWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = ShuffleAlgebra(QQ,'a')
            sage: A.one_basis()
            word:
            sage: A.one()
            B[word: ]
        """
        return self.basis().keys()([])

    def product_on_basis(self, w1, w2):
        r"""
        Return the product of basis elements ``w1`` and ``w2``, as per
        :meth:`AlgebrasWithBasis.ParentMethods.product_on_basis()`.

        INPUT:

        - ``w1``, ``w2`` -- Basis elements

        EXAMPLES::

            sage: A = ShuffleAlgebra(QQ,'abc')
            sage: W = A.basis().keys()
            sage: A.product_on_basis(W("acb"), W("cba"))
            B[word: acbacb] + B[word: acbcab] + 2*B[word: acbcba] + 2*B[word: accbab] + 4*B[word: accbba] + B[word: cabacb] + B[word: cabcab] + B[word: cabcba] + B[word: cacbab] + 2*B[word: cacbba] + 2*B[word: cbaacb] + B[word: cbacab] + B[word: cbacba]

            sage: (a,b,c) = A.algebra_generators()
            sage: a * (1-b)^2 * c
            2*B[word: abbc] - 2*B[word: abc] + 2*B[word: abcb] + B[word: ac] - 2*B[word: acb] + 2*B[word: acbb] + 2*B[word: babc] - 2*B[word: bac] + 2*B[word: bacb] + 2*B[word: bbac] + 2*B[word: bbca] - 2*B[word: bca] + 2*B[word: bcab] + 2*B[word: bcba] + B[word: ca] - 2*B[word: cab] + 2*B[word: cabb] - 2*B[word: cba] + 2*B[word: cbab] + 2*B[word: cbba]
        """
        return sum(self.basis()[u] for u in w1.shuffle(w2))

    def gen(self,i):
        r"""
        The ``i``-th generator of the algebra.

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: F = ShuffleAlgebra(ZZ,'xyz')
            sage: F.gen(0)
            B[word: x]

            sage: F.gen(4)
            Traceback (most recent call last):
            ...
            IndexError: argument i (= 4) must be between 0 and 2
        """
        n = self.__ngens
        if i < 0 or not i < n:
            raise IndexError("argument i (= %s) must be between 0 and %s"%(i, n-1))
        return self.algebra_generators()[i]

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        EXAMPLES::

            sage: A = ShuffleAlgebra(ZZ,'fgh'); A
            Shuffle Algebra on 3 generators ['f', 'g', 'h'] over Integer Ring
            sage: A.algebra_generators()
            Family (B[word: f], B[word: g], B[word: h])
        """
        Words = self.basis().keys()
        return Family( [self.monomial(Words(a)) for a in self._alphabet] )
        # FIXME: use this once the keys argument of FiniteFamily will be honoured
        # for the specifying the order of the elements in the family
        #return Family(self._alphabet, lambda a: self.term(self.basis().keys()(a)))

    gens = algebra_generators

    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = ShuffleAlgebra(QQ,'xy')
            sage: x, y = R.gens()
            sage: R(3) # indirect doctest
            3*B[word: ]
            sage: R(x)
            B[word: x]
        """
        P = x.parent()
        if isinstance(P, ShuffleAlgebra):
            if P is self:
                return x
            if not (P is self.base_ring()):
                return self.element_class(self, x.monomial_coefficients())
        # ok, not a shuffle algebra element (or should not be viewed as one).
        if isinstance(x, basestring):
            from sage.misc.sage_eval import sage_eval
            return sage_eval(x,locals=self.gens_dict())
        R = self.base_ring()
        # coercion via base ring
        x = R(x)
        if x == 0:
            return self.element_class(self,{})
        else:
            return self.from_base_ring_from_one_basis(x)

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of ``x`` into ``self``.

        Here is what canonically coerces to ``self``:

        - this shuffle algebra,

        - anything that coerces to the base ring of this shuffle algebra,

        - any shuffle algebra on the same variables, whose base ring
          coerces to the base ring of this shuffle algebra.

        EXAMPLES::

            sage: F = ShuffleAlgebra(GF(7), 'xyz'); F
            Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Finite Field of size 7

        Elements of the shuffle algebra canonically coerce in::

            sage: x, y, z = F.gens()
            sage: F.coerce(x*y) # indirect doctest
            B[word: xy] + B[word: yx]

        Elements of the integers coerce in, since there is a coerce map
        from `\ZZ` to GF(7)::

            sage: F.coerce(1)       # indirect doctest
            B[word: ]

        There is no coerce map from `\QQ` to `\GF{7}`::

            sage: F.coerce(2/3)  # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Rational Field to Shuffle Algebra on 3 generators ['x', 'y', 'z'] over Finite Field of size 7

        Elements of the base ring coerce in::

            sage: F.coerce(GF(7)(5))
            5*B[word: ]

        The shuffle algebra over `\ZZ` on `x, y, z` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = ShuffleAlgebra(ZZ,'xyz')
            sage: Gx,Gy,Gz = G.gens()
            sage: z = F.coerce(Gx**2 * Gy);z
            2*B[word: xxy] + 2*B[word: xyx] + 2*B[word: yxx]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the shuffle
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(x^3*y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Shuffle Algebra on 3 generators
            ['x', 'y', 'z'] over Finite Field of size 7 to Shuffle Algebra on 3
            generators ['x', 'y', 'z'] over Integer Ring
        """
        try:
            R = x.parent()

            # shuffle algebras in the same variables over any base
            # that coerces in:
            if isinstance(R,ShuffleAlgebra):
                if R.variable_names() == self.variable_names():
                    if self.has_coerce_map_from(R.base_ring()):
                        return self(x)
                    else:
                        raise TypeError("no natural map between bases of shuffle algebras")

        except AttributeError:
            pass

        # any ring that coerces to the base ring of this shuffle algebra.
        return self._coerce_try(x, [self.base_ring()])

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - Shuffle Algebras in the same variables over a base with a coercion
          map into ``self.base_ring()``.

        - Anything with a coercion into ``self.base_ring()``.

        TESTS::

            sage: F = ShuffleAlgebra(ZZ, 'xyz')
            sage: G = ShuffleAlgebra(QQ, 'xyz')
            sage: H = ShuffleAlgebra(ZZ, 'y')
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
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        # shuffle algebras in the same variable over any base that coerces in:
        if isinstance(R, ShuffleAlgebra):
            if R.variable_names() == self.variable_names():
                if self.base_ring().has_coerce_map_from(R.base_ring()):
                    return True
                else:
                    return False

        return self.base_ring().has_coerce_map_from(R)

