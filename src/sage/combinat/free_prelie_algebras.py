# -*- coding: utf-8 -*-
r"""
The free Pre-Lie algebras

AUTHORS:

- Florent Hivert, Frédéric Chapoton (2011)
"""

#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@lri.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.all import GradedModulesWithBasis
from sage.combinat.free_module import (CombinatorialFreeModule,
                                       CombinatorialFreeModuleElement)
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.rooted_tree import RootedTrees, LabelledRootedTrees
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from sage.categories.rings import Rings
from sage.sets.family import Family


class FreePreLieAlgebra(CombinatorialFreeModule):
    r"""
    An example of an algebra over an operad with basis: the free
    algebras over the PreLie operad

    EXAMPLES::

        sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
        sage: P = FreePreLieAlgebra(ZZ, 'abc'); P
        Free PreLie algebra on 3 generators ['a', 'b', 'c'] over Integer Ring
    """
    @staticmethod
    def __classcall_private__(cls, R, names):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: F1 = FreePreLieAlgebra(QQ, 'xyz')
            sage: F2 = FreePreLieAlgebra(QQ, ['x','y','z'])
            sage: F3 = FreePreLieAlgebra(QQ, Alphabet('xyz'))
            sage: F1 is F2 and F1 is F3
            True
        """
        return super(FreePreLieAlgebra, cls).__classcall__(cls, R,
                                                           Alphabet(names))

    def __init__(self, R, names=None):
        """
        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: A = FreePreLieAlgebra(QQ, '@'); A
            Free PreLie algebra on one generator ['@'] over Rational Field
            sage: TestSuite(A).run()
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        if len(names) == 1:
            Trees = RootedTrees()
        else:
            Trees = LabelledRootedTrees()
        # on aurait besoin ici de LabelledRootedTrees(names)
        # pour restreindre les etiquettes aux valeurs autorisees

        self._alphabet = names
        self.__ngens = self._alphabet.cardinality()
        CombinatorialFreeModule.__init__(self, R, Trees,
                                         latex_prefix="",
                                         category=GradedModulesWithBasis(R))

    def variable_names(self):
        r"""
        Return the names of the variables.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: R = FreePreLieAlgebra(QQ, 'xy')
            sage: R.variable_names()
            {'x', 'y'}
        """
        return self._alphabet

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: FreePreLieAlgebra(QQ, '@')  # indirect doctest
            Free PreLie algebra on one generator ['@'] over Rational Field
        """
        if self.__ngens == 1:
            gen = "one generator"
        else:
            gen = "{} generators".format(self.__ngens)
        s = "Free PreLie algebra on {} {} over {}"
        return s.format(gen, self._alphabet.list(), self.base_ring())

    def gen(self, i):
        r"""
        Return the ``i``-th generator of the algebra.

        INPUT:

        - ``i`` -- an integer

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: A = FreePreLieAlgebra(QQ,'@')
            sage: F = FreePreLieAlgebra(ZZ,'xyz')
            sage: F.gen(0)
            B[x[]]

            sage: F.gen(4)
            Traceback (most recent call last):
            ...
            IndexError: argument i (= 4) must be between 0 and 2
        """
        n = self.__ngens
        if i < 0 or not i < n:
            m = "argument i (= {}) must be between 0 and {}".format(i, n - 1)
            raise IndexError(m)
        return self.algebra_generators()[i]

    @cached_method
    def algebra_generators(self):
        r"""
        Return the generators of this algebra.

        These are the rooted trees with just one vertex.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: A = FreePreLieAlgebra(ZZ,'fgh'); A
            Free PreLie algebra on 3 generators ['f', 'g', 'h'] over
            Integer Ring
            sage: A.algebra_generators()
            Family (B[f[]], B[g[]], B[h[]])

            sage: A = FreePreLieAlgebra(QQ, ['x1','x2'])
            sage: A.algebra_generators()
            Family (B[x1[]], B[x2[]])
        """
        Trees = self.basis().keys()
        return Family([self.monomial(Trees([], a)) for a in self._alphabet])
        # FIXME: use this once the keys argument of FiniteFamily will
        # be honoured for the specifying the order of the elements in
        # the family
        #return Family(self._alphabet, lambda a:
        #self.term(self.basis().keys()(a)))

    gens = algebra_generators

    def degree_on_basis(self, t):
        """
        Return the degree of a rooted tree in the free Pre-Lie algebra.

        This is the number of vertices.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: A = FreePreLieAlgebra(QQ,'@')
            sage: RT = A.basis().keys()
            sage: A.degree_on_basis(RT([RT([])]))
            2
        """
        return t.node_number()

    def some_elements(self):
        """
        Return some elements of the free pre-Lie algebra.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: A = FreePreLieAlgebra(QQ,'@')
            sage: A.some_elements()
            [B[[]], B[[[]]], B[[[], [[]]]] + B[[[[[]]]]],
            B[[[], []]] + B[[[[]]]], B[[]]]

        With several generators::

            sage: A = FreePreLieAlgebra(QQ,'xy')
            sage: A.some_elements()
            [B[x[]],
            B[x[x[]]],
            B[x[x[], x[x[]]]] + B[x[x[x[x[]]]]],
            B[x[x[], x[]]] + B[x[x[x[]]]],
            B[y[]]]
        """
        o = self.gen(0)
        x = o < o
        y = o
        for w in self.gens():
            y = (y < w)
        return [o, x, x < x, x < o, w]

    def pre_Lie_product_on_basis(self, x, y):
        """
        Return the pre-Lie product of two trees.

        This is the sum over all graftings of the root of `y` over a vertex
        of `x`. The root of the resulting trees is the root of `x`.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: A = FreePreLieAlgebra(QQ,'@')
            sage: RT = A.basis().keys()
            sage: x = RT([RT([])])
            sage: A.pre_Lie_product_on_basis(x, x)
            B[[[], [[]]]] + B[[[[[]]]]]
        """
        return sum(self.basis()[u] for u in x.graft_list(y))

    @lazy_attribute
    def pre_Lie_product(self):
        """
        Return the pre-Lie product.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: A = FreePreLieAlgebra(QQ,'@')
            sage: RT = A.basis().keys()
            sage: x = A(RT([RT([])]))
            sage: A.pre_Lie_product(x, x)
            B[[[], [[]]]] + B[[[[[]]]]]
        """
        plb = self.pre_Lie_product_on_basis
        return self._module_morphism(self._module_morphism(plb, position=0,
                                                           codomain=self),
                                     position=1)

    def nap_product_on_basis(self, x, y):
        """
        Return the nap product of two trees.

        This is the grafting of the root of `y` over the root
        of `x`. The root of the resulting tree is the root of `x`.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: A = FreePreLieAlgebra(QQ,'@')
            sage: RT = A.basis().keys()
            sage: x = RT([RT([])])
            sage: A.nap_product_on_basis(x, x)
            B[[[], [[]]]]
        """
        return self.basis()[x.graft_on_root(y)]

    @lazy_attribute
    def nap_product(self):
        """
        Return the nap product.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: A = FreePreLieAlgebra(QQ,'@')
            sage: RT = A.basis().keys()
            sage: x = A(RT([RT([])]))
            sage: A.nap_product(x, x)
            B[[[], [[]]]]
        """
        npb = self.nap_product_on_basis
        return self._module_morphism(self._module_morphism(npb,
                                                           position=0,
                                                           codomain=self),
                                     position=1)

    # after this line : coercion
    # def _element_constructor_(self, x):
    #     r"""
    #     Convert ``x`` into ``self``.

    #     EXAMPLES::

    #         sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
    #         sage: R = FreePreLieAlgebra(QQ, 'xy')
    #         sage: x, y = R.gens()
    #         sage: R(x)
    #         B[[x]]
    #         sage: R(x+4*y)
    #         4*B[y[]] + B[x[]]

    #         sage: Trees = R.basis().keys()
    #         sage: R(Trees([],'x')

    #         sage: D = FreePreLieAlgebra(ZZ, 'xy')
    #         sage: X, Y = D.gens()
    #         sage: R(X-Y).parent()
    #         Free PreLie algebra on 2 generators ['x', 'y'] over Rational Field
    #     """
    #     if x in self.basis().keys():
    #         return self.monomial(x)
    #     try:
    #         P = x.parent()
    #         if isinstance(P, FreePreLieAlgebra):
    #             if P is self:
    #                 return x
    #             return self.element_class(self, x.monomial_coefficients())
    #     except AttributeError:
    #         raise TypeError('not able to coerce this in this algebra')

    #     # ok, not a pre-Lie algebra element (or should not be viewed as one).

    def _coerce_impl(self, x):
        r"""
        Canonical coercion of ``x`` into ``self``.

        Here is what canonically coerces to ``self``:

        - this free prelie algebra,

        - any free prelie algebra on the same variables, whose base ring
          coerces to the base ring of this free prelie algebra.

        EXAMPLES::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: F = FreePreLieAlgebra(GF(7), 'xyz'); F
            Free PreLie algebra on 3 generators ['x', 'y', 'z'] over
            Finite Field of size 7

        Elements of the free PreLie algebra canonically coerce in::

            sage: x, y, z = F.gens()
            sage: F.coerce(x+y) == x+y
            True

        The free prelie algebra over `\ZZ` on `x, y, z` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = FreePreLieAlgebra(ZZ, 'xyz')
            sage: Gx,Gy,Gz = G.gens()
            sage: z = F.coerce(Gx+Gy); z
            ??<>??
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the free prelie
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Free PreLie algebra
            on 3 generators ['x', 'y', 'z'] over Finite Field of size
            7 to Free PreLie algebra on 3 generators ['x', 'y', 'z']
            over Integer Ring
        """
        try:
            R = x.parent()

            # prelie algebras in the same variables over any base
            # that coerces in:
            if isinstance(R, FreePreLieAlgebra):
                if R.variable_names() == self.variable_names():
                    if self.base_ring().has_coerce_map_from(R.base_ring()):
                        return self(x)
                    raise TypeError("no natural map between bases of "
                                    "Free PreLie algebras")

        except AttributeError:
            pass

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - PreLie Algebras in the same variables over a base with a coercion
          map into ``self.base_ring()``

        TESTS::

            sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
            sage: F = FreePreLieAlgebra(ZZ, 'xyz')
            sage: G = FreePreLieAlgebra(QQ, 'xyz')
            sage: H = FreePreLieAlgebra(ZZ, 'y')
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
        if isinstance(R, FreePreLieAlgebra):
            if R.variable_names() == self.variable_names():
                if self.base_ring().has_coerce_map_from(R.base_ring()):
                    return True
        return False

    class Element(CombinatorialFreeModuleElement):
        def __lt__(self, other):
            r"""
            Shortcut for the prelie product using the symbol ``<``.

            EXAMPLES::

                sage: from sage.combinat.free_prelie_algebras import FreePreLieAlgebra
                sage: A = FreePreLieAlgebra(QQ,'@')
                sage: a = A.gen(0)
                sage: a < a
                B[[[]]]

            .. WARNING::

                Due to priority rules for operators, term must be put
                within parentheses inside sum, product... For example you must
                write::

                    sage: a = A.gen(0)
                    sage: (a<a) + a
                    B[[]] + B[[[]]]

                Indeed ``a<a + a`` is understood as ``a< (a + a)``

                    sage: (a<a + a) - (a < (a + a))
                    0
            """
            parent = self.parent()
            assert(parent == other.parent())
            return parent.pre_Lie_product(self, other)
