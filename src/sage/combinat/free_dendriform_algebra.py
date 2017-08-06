# -*- coding: utf-8 -*-
r"""
Free Dendriform Algebras

AUTHORS:

Frédéric Chapoton (2017)
"""
# ****************************************************************************
#       Copyright (C) 2010-2015 Frédéric Chapoton <chapoton@unistra.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.algebras import Algebras
from sage.categories.magmatic_algebras import MagmaticAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.binary_tree import (BinaryTrees, BinaryTree,
                                       LabelledBinaryTrees,
                                       LabelledBinaryTree)
from sage.combinat.ordered_tree import OrderedTree, LabelledOrderedTree
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.categories.rings import Rings
from sage.sets.family import Family


class FreeDendriformAlgebra(CombinatorialFreeModule):
    r"""
    The free dendriform algebra.

    Dendriform algebras are associative algebras, where the associative
    product `*` is decomposed as a sum of two binary operations

    .. MATH::

        x * y = x \succ y + x \prec y

    that satisfy the axioms:

    .. MATH::

        (x \succ y) \prec z = x \succ (y \prec z),

    .. MATH::

        (x \prec y) \prec z = x \prec (y * z).

    .. MATH::

        (x * y) \succ z = x \succ (y \succ z).

    The free Dendriform algebra on a given set `E` has an explicit
    description using (planar) binary trees, just as the free
    associative algebra can be described using words. The underlying
    vector space has a basis indexed by finite binary trees endowed
    with a map from their vertices to `E`. In this basis, the
    associative product of two (decorated) binary trees `S * T` is the
    sum over all possible ways of identifying (glueing) the rightmost path in
    `S` and the leftmost path in `T`.

    The decomposition of the associative product as the sum of two
    binary operations `\succ` and
    `\prec` is made by separating the terms according to the origin of
    the root vertex. For `x \succ y`, one keeps the terms where the root
    vertex comes from `y`, whereas for `x \prec y` one keeps the terms
    where the root vertex comes from `x`.

    The free dendriform algebra can also be considered as the free
    algebra over the Dendriform operad.

    .. NOTE::

        The usual binary operator `*` is used for the
        associative product.

    EXAMPLES::

        sage: F = algebras.FreeDendriform(ZZ, 'xyz')
        sage: x,y,z = F.gens()
        sage: (x * y) * z
        B[x[., y[., z[., .]]]] + B[x[., z[y[., .], .]]] + B[y[x[., .], z[., .]]] + B[z[x[., y[., .]], .]] + B[z[y[x[., .], .], .]]

    The free dendriform algebra is associative::

        sage: x * (y * z) == (x * y) * z
        True

    The associative product decomposes in two parts::

        sage: x * y == F.prec(x, y) + F.succ(x, y)
        True

    The axioms hold::

        sage: F.prec(F.succ(x, y), z) == F.succ(x, F.prec(y, z))
        True
        sage: F.prec(F.prec(x, y), z) == F.prec(x, y * z)
        True
        sage: F.succ(x * y, z) == F.succ(x, F.succ(y, z))
        True

    When there is only one generator, unlabelled trees are used instead::

        sage: F1 = algebras.FreeDendriform(QQ, 'w')
        sage: w = F1.gen(0); w
        B[[., .]]
        sage: w * w * w
        B[[., [., [., .]]]] + B[[., [[., .], .]]] + B[[[., .], [., .]]] + B[[[., [., .]], .]] + B[[[[., .], .], .]]

    REFERENCES:

    - [LodayRonco]_
    """
    @staticmethod
    def __classcall_private__(cls, R, names):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: F1 = algebras.FreeDendriform(QQ, 'xyz')
            sage: F2 = algebras.FreeDendriform(QQ, ['x','y','z'])
            sage: F3 = algebras.FreeDendriform(QQ, Alphabet('xyz'))
            sage: F1 is F2 and F1 is F3
            True
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        return super(FreeDendriformAlgebra, cls).__classcall__(cls, R,
                                                               Alphabet(names))

    def __init__(self, R, names=None):
        """
        Initialize ``self``.

        TESTS::

            sage: A = algebras.FreeDendriform(QQ, '@'); A
            Free Dendriform algebra on one generator ['@'] over Rational Field
            sage: TestSuite(A).run()

            sage: F = algebras.FreeDendriform(QQ, 'xy')
            sage: TestSuite(F).run() # long time
        """
        if names.cardinality() == 1:
            Trees = BinaryTrees()
            key = BinaryTree._sort_key
        else:
            Trees = LabelledBinaryTrees()
            key = LabelledBinaryTree._sort_key
        # Here one would need LabelledBinaryTrees(names)
        # so that one can restrict the labels to some fixed set
        self._alphabet = names
        cat = Algebras(R).WithBasis().Graded()
        CombinatorialFreeModule.__init__(self, R, Trees,
                                         latex_prefix="",
                                         sorting_key=key,
                                         category=cat)

    def variable_names(self):
        r"""
        Return the names of the variables.

        EXAMPLES::

            sage: R = algebras.FreeDendriform(QQ, 'xy')
            sage: R.variable_names()
            {'x', 'y'}
        """
        return self._alphabet

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: algebras.FreeDendriform(QQ, '@')  # indirect doctest
            Free Dendriform algebra on one generator ['@'] over Rational Field
        """
        n = self.algebra_generators().cardinality()
        if n == 1:
            gen = "one generator"
        else:
            gen = "{} generators".format(n)
        s = "Free Dendriform algebra on {} {} over {}"
        try:
            return s.format(gen, self._alphabet.list(), self.base_ring())
        except NotImplementedError:
            return s.format(gen, self._alphabet, self.base_ring())

    def gen(self, i):
        r"""
        Return the ``i``-th generator of the algebra.

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: F = algebras.FreeDendriform(ZZ, 'xyz')
            sage: F.gen(0)
            B[x[., .]]

            sage: F.gen(4)
            Traceback (most recent call last):
            ...
            IndexError: argument i (= 4) must be between 0 and 2
        """
        G = self.algebra_generators()
        n = G.cardinality()
        if i < 0 or not i < n:
            m = "argument i (= {}) must be between 0 and {}".format(i, n - 1)
            raise IndexError(m)
        return G[G.keys().unrank(i)]

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        These are the binary trees with just one vertex.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(ZZ, 'fgh'); A
            Free Dendriform algebra on 3 generators ['f', 'g', 'h']
             over Integer Ring
            sage: list(A.algebra_generators())
            [B[f[., .]], B[g[., .]], B[h[., .]]]

            sage: A = algebras.FreeDendriform(QQ, ['x1','x2'])
            sage: list(A.algebra_generators())
            [B[x1[., .]], B[x2[., .]]]
        """
        Trees = self.basis().keys()
        return Family(self._alphabet, lambda a: self.monomial(Trees([], a)))

    def gens(self):
        """
        Return the generators of ``self`` (as an algebra).

        EXAMPLES::

            sage: A = algebras.FreeDendriform(ZZ, 'fgh')
            sage: A.gens()
            (B[f[., .]], B[g[., .]], B[h[., .]])
        """
        return tuple(self.algebra_generators())

    def degree_on_basis(self, t):
        """
        Return the degree of a binary tree in the free Dendriform algebra.

        This is the number of vertices.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ,'@')
            sage: RT = A.basis().keys()
            sage: u = RT([], '@')
            sage: A.degree_on_basis(u.over(u))
            2
        """
        return t.node_number()

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, 'xy')
            sage: A.an_element()
            B[x[., .]] + 2*B[x[., x[., .]]] + 2*B[x[x[., .], .]]
        """
        o = self.gen(0)
        return o + 2 * o * o

    def some_elements(self):
        """
        Return some elements of the free dendriform algebra.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ,'@')
            sage: A.some_elements()
            [B[.],
             B[[., .]],
             B[[., [., .]]] + B[[[., .], .]],
             B[.] + B[[., [., .]]] + B[[[., .], .]]]

        With several generators::

            sage: A = algebras.FreeDendriform(QQ, 'xy')
            sage: A.some_elements()
            [B[.],
             B[x[., .]],
             B[x[., x[., .]]] + B[x[x[., .], .]],
             B[.] + B[x[., x[., .]]] + B[x[x[., .], .]]]
        """
        u = self.one()
        o = self.gen(0)
        x = o * o
        y = u + x
        return [u, o, x, y]

    @cached_method
    def one(self):
        r"""
        Return the element `1` of ``self``.

        This is the unit for the associative dendriform product `*`.

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, '@')
            sage: A.one()
            B[.]
            sage: A = algebras.FreeDendriform(QQ, 'xy')
            sage: A.one()
            B[.]
        """
        Trees = self.basis().keys()
        return self._monomial(Trees(None))

    def product_on_basis(self, x, y):
        r"""
        Return the `*` associative dendriform product of two trees.

        This is the sum over all possible ways of identifying the
        rightmost path in `x` and the leftmost path in `y`. Every term
        corresponds to a shuffle of the vertices on the rightmost path
        in `x` and the vertices on the leftmost path in `y`.

        .. SEEALSO::

            - :meth:`succ_product_on_basis`, :meth:`prec_product_on_basis`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, '@')
            sage: RT = A.basis().keys()
            sage: x = RT([], '@')
            sage: A.product_on_basis(x, x)
            B[[., [., .]]] + B[[[., .], .]]
        """
        return self.sum(self.basis()[u] for u in x.dendriform_shuffle(y))

    def succ_product_on_basis(self, x, y):
        r"""
        Return the `\succ` dendriform product of two trees.

        This is the sum over all possible ways to identify the rightmost path
        in `x` and the leftmost path in `y`, with the additional condition
        that the root vertex of the result comes from `y`.

        The usual symbol for this operation is `\succ`.

        .. SEEALSO::

            - :meth:`product_on_basis`, :meth:`prec_product_on_basis`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, '@')
            sage: RT = A.basis().keys()
            sage: x = RT([], '@')
            sage: A.succ_product_on_basis(x, x)
            B[[[., .], .]]

        TESTS::

            sage: u = A.one().support()[0]
            sage: A.succ_product_on_basis(u, u)
            Traceback (most recent call last):
            ...
            ValueError: dendriform products | < | and | > | are not defined
        """
        if y.is_empty():
            if x.is_empty():
                raise ValueError("dendriform products | < | and | > | are "
                                 "not defined")
            else:
                return []
        if x.is_empty():
            return [y]
        K = self.basis().keys()
        if hasattr(y, 'label'):
            return self.sum(self.basis()[K([u, y[1]], y.label())]
                            for u in x.dendriform_shuffle(y[0]))

        return self.sum(self.basis()[K([u, y[1]])]
                        for u in x.dendriform_shuffle(y[0]))

    @lazy_attribute
    def succ(self):
        r"""
        Return the `\succ` dendriform product.

        This is the sum over all possible ways of identifying the
        rightmost path in `x` and the leftmost path in `y`, with the
        additional condition that the root vertex of the result comes
        from `y`.

        The usual symbol for this operation is `\succ`.

        .. SEEALSO::

            :meth:`product`, :meth:`prec`, :meth:`over`, :meth:`under`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, '@')
            sage: RT = A.basis().keys()
            sage: x = A.gen(0)
            sage: A.succ(x, x)
            B[[[., .], .]]
        """
        suc = self.succ_product_on_basis
        return self._module_morphism(self._module_morphism(suc, position=0,
                                                           codomain=self),
                                     position=1)

    def prec_product_on_basis(self, x, y):
        r"""
        Return the `\prec` dendriform product of two trees.

        This is the sum over all possible ways of identifying the
        rightmost path in `x` and the leftmost path in `y`, with the
        additional condition that the root vertex of the result comes
        from `x`.

        The usual symbol for this operation is `\prec`.

        .. SEEALSO::

            - :meth:`product_on_basis`, :meth:`succ_product_on_basis`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, '@')
            sage: RT = A.basis().keys()
            sage: x = RT([], '@')
            sage: A.prec_product_on_basis(x, x)
            B[[., [., .]]]

        TESTS::

            sage: u = A.one().support()[0]
            sage: A.prec_product_on_basis(u, u)
            Traceback (most recent call last):
            ...
            ValueError: dendriform products | < | and | > | are not defined
        """
        if x.is_empty() and y.is_empty():
            raise ValueError("dendriform products | < | and | > | are "
                             "not defined")
        if x.is_empty():
            return []
        if y.is_empty():
            return [x]
        K = self.basis().keys()
        if hasattr(y, 'label'):
            return self.sum(self.basis()[K([x[0], u], x.label())]
                            for u in x[1].dendriform_shuffle(y))

        return self.sum(self.basis()[K([x[0], u])]
                        for u in x[1].dendriform_shuffle(y))

    @lazy_attribute
    def prec(self):
        r"""
        Return the `\prec` dendriform product.

        This is the sum over all possible ways to identify the rightmost path
        in `x` and the leftmost path in `y`, with the additional condition
        that the root vertex of the result comes from `x`.

        The usual symbol for this operation is `\prec`.

        .. SEEALSO::

            :meth:`product`, :meth:`succ`, :meth:`over`, :meth:`under`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, '@')
            sage: RT = A.basis().keys()
            sage: x = A.gen(0)
            sage: A.prec(x, x)
            B[[., [., .]]]
        """
        pre = self.prec_product_on_basis
        return self._module_morphism(self._module_morphism(pre, position=0,
                                                           codomain=self),
                                     position=1)

    @lazy_attribute
    def over(self):
        r"""
        Return the over product.

        The over product `x/y` is the binary tree obtained by
        grafting the root of `y` at the rightmost leaf of `x`.

        The usual symbol for this operation is `/`.

        .. SEEALSO::

            :meth:`product`, :meth:`succ`, :meth:`prec`, :meth:`under`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, '@')
            sage: RT = A.basis().keys()
            sage: x = A.gen(0)
            sage: A.over(x, x)
            B[[., [., .]]]
        """
        ov = lambda x, y: self._monomial(x.over(y))
        return self._module_morphism(self._module_morphism(ov, position=0,
                                                           codomain=self),
                                     position=1)

    @lazy_attribute
    def under(self):
        r"""
        Return the under product.

        The over product `x \backslash y` is the binary tree obtained by
        grafting the root of `x` at the leftmost leaf of `y`.

        The usual symbol for this operation is `\backslash`.

        .. SEEALSO::

            :meth:`product`, :meth:`succ`, :meth:`prec`, :meth:`over`

        EXAMPLES::

            sage: A = algebras.FreeDendriform(QQ, '@')
            sage: RT = A.basis().keys()
            sage: x = A.gen(0)
            sage: A.under(x, x)
            B[[[., .], .]]
        """
        und = lambda x, y: self._monomial(x.under(y))
        return self._module_morphism(self._module_morphism(und, position=0,
                                                           codomain=self),
                                     position=1)

    # after this line : coercion
    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = algebras.FreeDendriform(QQ, 'xy')
            sage: x, y = R.gens()
            sage: R(x)
            B[x[., .]]
            sage: R(x+4*y)
            B[x[., .]] + 4*B[y[., .]]

            sage: Trees = R.basis().keys()
            sage: R(Trees([],'x'))
            B[x[., .]]
            sage: D = algebras.FreeDendriform(ZZ, 'xy')
            sage: X, Y = D.gens()
            sage: R(X-Y).parent()
            Free Dendriform algebra on 2 generators ['x', 'y'] over Rational Field
        """
        if x in self.basis().keys():
            return self.monomial(x)
        try:
            P = x.parent()
            if isinstance(P, FreeDendriformAlgebra):
                if P is self:
                    return x
                return self.element_class(self, x.monomial_coefficients())
        except AttributeError:
            raise TypeError('not able to coerce this in this algebra')
        # Ok, not a dendriform algebra element (or should not be viewed as one).

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - free dendriform algebras in the same variables over a base with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: F = algebras.FreeDendriform(GF(7), 'xyz'); F
            Free Dendriform algebra on 3 generators ['x', 'y', 'z']
             over Finite Field of size 7

        Elements of the free dendriform algebra canonically coerce in::

            sage: x, y, z = F.gens()
            sage: F.coerce(x+y) == x+y
            True

        The free dendriform algebra over `\ZZ` on `x, y, z` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = algebras.FreeDendriform(ZZ, 'xyz')
            sage: Gx,Gy,Gz = G.gens()
            sage: z = F.coerce(Gx+Gy); z
            B[x[., .]] + B[y[., .]]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the free dendriform
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Free Dendriform algebra
             on 3 generators ['x', 'y', 'z'] over Finite Field of size
             7 to Free Dendriform algebra on 3 generators ['x', 'y', 'z']
             over Integer Ring

        TESTS::

            sage: F = algebras.FreeDendriform(ZZ, 'xyz')
            sage: G = algebras.FreeDendriform(QQ, 'xyz')
            sage: H = algebras.FreeDendriform(ZZ, 'y')
            sage: F._coerce_map_from_(G)
            False
            sage: G._coerce_map_from_(F)
            True
            sage: F._coerce_map_from_(H)
            False
            sage: F._coerce_map_from_(QQ)
            False
            sage: G._coerce_map_from_(QQ)
            False
            sage: F.has_coerce_map_from(PolynomialRing(ZZ, 3, 'x,y,z'))
            False
        """
        # free prelie algebras in the same variables
        # over any base that coerces in:
        if isinstance(R, FreeDendriformAlgebra):
            if R.variable_names() == self.variable_names():
                if self.base_ring().has_coerce_map_from(R.base_ring()):
                    return True
        return False
