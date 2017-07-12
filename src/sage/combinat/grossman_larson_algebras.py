# -*- coding: utf-8 -*-
r"""
Grossman-Larson Hopf Algebras

AUTHORS:

- Frédéric Chapoton (2017)
"""

#*****************************************************************************
#       Copyright (C) 2017 Frédéric Chapoton
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.hopf_algebras import HopfAlgebras
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.alphabet import Alphabet
from sage.combinat.rooted_tree import (RootedTrees, RootedTree,
                                       LabelledRootedTrees,
                                       LabelledRootedTree)
from sage.misc.cachefunc import cached_method
from sage.categories.rings import Rings
from sage.sets.family import Family
from sage.misc.misc import powerset
from sage.categories.cartesian_product import cartesian_product

# we use a special symbol for the fake root
ROOT = '#'


class GrossmanLarsonAlgebra(CombinatorialFreeModule):
    r"""
    The Grossman-Larson Hopf Algebra.

    The Grossman-Larson Hopf Algebras are Hopf algebras with a basis
    indexed by forests of decorated rooted trees. They are the
    universal enveloping algebras of free pre-Lie algebras, seen
    as Lie algebras.

    The Grossman-Larson Hopf algebra on a given set `E` has an
    explicit description using rooted forests. The underlying vector
    space has a basis indexed by finite rooted forests endowed with a
    map from their vertices to `E`. In this basis, the product of two
    (decorated) rooted forests `S * T` is the sum over all maps from
    the set of roots of `T` to the union of a singleton `\{\#\}` and
    the set of vertices of `S`. Given such a map, one defines a new
    forest as follows. Starting from the disjoint union of all trees
    of `S` and `T`, one adds an edge from every root of `T` to its
    image when this image is not the fake vertex labelled ``#``.

    .. NOTE::

        Instead of using forests as labels for the basis, we use rooted trees
        and ignore their root (labelled ``#`` when labels are present).

    EXAMPLES::

        sage: G = algebras.GrossmanLarson(QQ, 'xy')
        sage: x, y = G.gens()
        sage: unicode_art(x*y)
        B  + B
         #    ╭#╮
         │    │ │
         y    x y
         │
         x

        sage: unicode_art(x*x*x)
        B  + B    + 3*B    + B
         #     #       ╭#╮    ╭─#─╮
         │     │       │ │    │ │ │
         x    ╭x╮      x x    x x x
         │    │ │        │
         x    x x        x
         │
         x

    The Grossman-Larson algebra is associative::

        sage: z = x * y
        sage: x * (y * z) == (x * y) * z
        True

    When there is only one generator, unlabelled trees are used instead::

        sage: G = algebras.GrossmanLarson(QQ, 1)
        sage: x = G.gens()[0]
        sage: unicode_art(x*x)
        B  + B
         o    ╭o╮
         │    │ │
         o    o o
         │
         o

    REFERENCES:

    .. [Pana2002] \F. Panaite, *Relating the Connes-Kreimer and Grossman-Larson  Hopf algebras built on rooted trees*, Lett. Math. Phys. 51 (2000), no. 3, pages 211-219.

    .. [GroLar1] \R. Grossman and R. G. Larson, *Hopf algebraic structure of families of trees*, J. Algebra 126 (1) (1989), 184-210.
    """
    @staticmethod
    def __classcall_private__(cls, R, names):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: F1 = algebras.GrossmanLarson(QQ, 'xyz')
            sage: F2 = algebras.GrossmanLarson(QQ, ['x','y','z'])
            sage: F3 = algebras.GrossmanLarson(QQ, Alphabet('xyz'))
            sage: F1 is F2 and F1 is F3
            True
        """
        if R not in Rings():
            raise TypeError("argument R must be a ring")
        return super(GrossmanLarsonAlgebra, cls).__classcall__(cls, R,
                                                               Alphabet(names))

    def __init__(self, R, names=None):
        """
        Initialize ``self``.

        TESTS::

            sage: A = algebras.GrossmanLarson(QQ, '@'); A
            Grossman-Larson Hopf algebra on one generator ['@']
            over Rational Field
            sage: TestSuite(A).run(skip='_test_antipode')

            sage: F = algebras.GrossmanLarson(QQ, 'xy')
            sage: TestSuite(F).run(skip='_test_antipode') # long time
        """
        if names.cardinality() == 1:
            Trees = RootedTrees()
            key = RootedTree.sort_key
        else:
            Trees = LabelledRootedTrees()
            key = LabelledRootedTree.sort_key
        # Here one would need LabelledRootedTrees(names)
        # so that one can restrict the labels to some fixed set

        self._alphabet = names
        cat = HopfAlgebras(R).WithBasis().Graded()
        CombinatorialFreeModule.__init__(self, R, Trees,
                                         latex_prefix="",
                                         sorting_key=key,
                                         category=cat)

    def variable_names(self):
        r"""
        Return the names of the variables.

        EXAMPLES::

            sage: R = algebras.GrossmanLarson(QQ, 'xy')
            sage: R.variable_names()
            {'x', 'y'}
        """
        return self._alphabet

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: algebras.GrossmanLarson(QQ, '@')  # indirect doctest
            Grossman-Larson Hopf algebra on one generator ['@'] over Rational Field
        """
        n = self.algebra_generators().cardinality()
        if n == 1:
            gen = "one generator"
        else:
            gen = "{} generators".format(n)
        s = "Grossman-Larson Hopf algebra on {} {} over {}"
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

            sage: F = algebras.GrossmanLarson(ZZ, 'xyz')
            sage: F.gen(0)
            B[#[x[]]]

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

        These are the rooted forests with just one vertex.
        In our representation, this means a rooted tree with 2 vertices,
        including the fake root ``#``.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(ZZ, 'fgh'); A
            Grossman-Larson Hopf algebra on 3 generators ['f', 'g', 'h']
             over Integer Ring
            sage: list(A.algebra_generators())
            [B[#[f[]]], B[#[g[]]], B[#[h[]]]]

            sage: A = algebras.GrossmanLarson(QQ, ['x1','x2'])
            sage: list(A.algebra_generators())
            [B[#[x1[]]], B[#[x2[]]]]
        """
        Trees = self.basis().keys()
        return Family(self._alphabet,
                      lambda a: self.monomial(Trees([Trees([], a)], ROOT)))

    def gens(self):
        """
        Return the generators of ``self`` (as an algebra).

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(ZZ, 'fgh')
            sage: A.gens()
            (B[#[f[]]], B[#[g[]]], B[#[h[]]])
        """
        return tuple(self.algebra_generators())

    def degree_on_basis(self, t):
        """
        Return the degree of a rooted tree in the free GL algebra.

        This is the number of vertices minus one.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ,'@')
            sage: RT = A.basis().keys()
            sage: A.degree_on_basis(RT([RT([])]))
            1
        """
        return t.node_number() - 1

    @cached_method
    def an_element(self):
        """
        Return an element of ``self``.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, 'xy')
            sage: A.an_element()
            B[#[x[]]] + 2*B[#[x[x[]]]] + 2*B[#[x[], x[]]]
        """
        o = self.gen(0)
        return o + 2 * o * o

    def some_elements(self):
        """
        Return some elements of the Grossman-Larson Hopf algebra.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ,'@')
            sage: A.some_elements()
            [B[[[]]],
             B[[[[]]]] + B[[[], []]],
             4*B[[[[]]]] + 4*B[[[], []]],
             B[[]] + B[[[[]]]] + B[[[], []]]]

        With several generators::

            sage: A = algebras.GrossmanLarson(QQ, 'xy')
            sage: A.some_elements()
            [B[#[x[]]],
             B[#[x[x[]]]] + B[#[x[], x[]]],
             B[#[x[x[]]]] + 3*B[#[y[x[]]]] + B[#[x[], x[]]] + 3*B[#[x[], y[]]],
             B[#[]] + B[#[x[x[]]]] + B[#[x[], x[]]]]
        """
        o = self.gen(0)
        o1 = list(self.algebra_generators())[-1]
        x = o * o
        y = o * o1
        return [o, x, x + 3 * y, 1 + x]

    def product_on_basis(self, x, y):
        """
        Return the GL product of two trees.

        This is the sum over all graftings of the root of `y` over a vertex
        of `x`. The root of the resulting trees is the root of `x`.

        Sum over the list of rooted trees obtained from a rooted tree `x`
        and a rooted tree `y` by the composition `x o_i y` where `i` is the root
        of `x`.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, '@')
            sage: RT = A.basis().keys()
            sage: x = RT([RT([])])
            sage: A.product_on_basis(x, x)
            B[[[[]]]] + B[[[], []]]
        """
        return self.sum(self.basis()[y.single_graft(x, graftingFunction)]
                        for graftingFunction in
                        cartesian_product([list(y.paths())] * len(x)))

    def one_basis(self):
        """
        Return the empty rooted forest.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, 'ab')
            sage: A.one_basis()
            #[]

            sage: A = algebras.GrossmanLarson(QQ, 1)
            sage: A.one_basis()
            []
        """
        Trees = self.basis().keys()
        return Trees([], ROOT)

    def coproduct_on_basis(self, x):
        """
        Return the coproduct of a tree.

        EXAMPLES::

            sage: G = algebras.GrossmanLarson(QQ,2)
            sage: x, y = G.gens()
            sage: unicode_art(G.coproduct(x))  # indirect doctest
            1 # B  + B  # 1
                 #    #
                 │    │
                 0    0

            sage: unicode_art(G.coproduct(x*y))  # indirect doctest
            1 # B    + 1 # B  + B  # B  + B    # 1 + B  # B  + B  # 1
                 ╭#╮        #    #    #    ╭#╮        #    #    #
                 │ │        │    │    │    │ │        │    │    │
                 0 1        1    0    1    0 1        1    0    1
                            │                                   │
                            0                                   0
        """
        B = self.basis()
        Trees = B.keys()
        subtrees = list(x)
        indx = list(range(len(subtrees)))
        return sum(B[Trees([subtrees[i] for i in S], ROOT)].tensor(
                   B[Trees([subtrees[i] for i in indx if i not in S], ROOT)])
                   for S in powerset(indx))

    def counit_on_basis(self, x):
        """
        Return the counit on a basis element.

        This is zero unless the forest `x` is empty.

        EXAMPLES::

            sage: A = algebras.GrossmanLarson(QQ, 'xy')
            sage: RT = A.basis().keys()
            sage: x = RT([RT([],'x')],'#')
            sage: A.counit_on_basis(x)
            0
            sage: A.counit_on_basis(RT([],'#'))
            1
        """
        if x.node_number() == 1:
            return self.base_ring().one()
        return self.base_ring().zero()

    # after this line : coercion
    def _element_constructor_(self, x):
        r"""
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R = algebras.GrossmanLarson(QQ, 'xy')
            sage: x, y = R.gens()
            sage: R(x)
            B[#[x[]]]
            sage: R(x+4*y)
            B[#[x[]]] + 4*B[#[y[]]]

            sage: Trees = R.basis().keys()
            sage: R(Trees([],'#'))
            B[#[]]

            sage: D = algebras.GrossmanLarson(ZZ, 'xy')
            sage: X, Y = D.gens()
            sage: R(X-Y).parent()
            Grossman-Larson Hopf algebra on 2 generators ['x', 'y'] over Rational Field

        TESTS::

            sage: Trees = R.basis().keys()
            sage: R(Trees([],'x'))
            Traceback (most recent call last):
            ...
            ValueError: incorrect root label
        """
        if x in self.basis().keys():
            if hasattr(x, 'label') and x.label() != ROOT:
                raise ValueError('incorrect root label')
            return self.monomial(x)
        try:
            P = x.parent()
            if isinstance(P, GrossmanLarsonAlgebra):
                if P is self:
                    return x
                return self.element_class(self, x.monomial_coefficients())
        except AttributeError:
            raise TypeError('not able to coerce this in this algebra')
        # Ok, not an element (or should not be viewed as one).

    def _coerce_map_from_(self, R):
        r"""
        Return ``True`` if there is a coercion from ``R`` into ``self``
        and ``False`` otherwise.

        The things that coerce into ``self`` are

        - Grossman-Larson Hopf algebras in the same variables over a base with
          a coercion map into ``self.base_ring()``

        EXAMPLES::

            sage: F = algebras.GrossmanLarson(GF(7), 'xyz'); F
            Grossman-Larson Hopf algebra on 3 generators ['x', 'y', 'z']
             over Finite Field of size 7

        Elements of the Grossman-Larson Hopf algebra canonically coerce in::

            sage: x, y, z = F.gens()
            sage: F.coerce(x+y) == x+y
            True

        The Grossman-Larson Hopf algebra over `\ZZ` on `x, y, z` coerces in, since
        `\ZZ` coerces to `\GF{7}`::

            sage: G = algebras.GrossmanLarson(ZZ, 'xyz')
            sage: Gx,Gy,Gz = G.gens()
            sage: z = F.coerce(Gx+Gy); z
            B[#[x[]]] + B[#[y[]]]
            sage: z.parent() is F
            True

        However, `\GF{7}` does not coerce to `\ZZ`, so the Grossman-Larson
        algebra over `\GF{7}` does not coerce to the one over `\ZZ`::

            sage: G.coerce(y)
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Grossman-Larson Hopf algebra
             on 3 generators ['x', 'y', 'z'] over Finite Field of size
             7 to Grossman-Larson Hopf algebra on 3 generators ['x', 'y', 'z']
             over Integer Ring

        TESTS::

            sage: F = algebras.GrossmanLarson(ZZ, 'xyz')
            sage: G = algebras.GrossmanLarson(QQ, 'xyz')
            sage: H = algebras.GrossmanLarson(ZZ, 'y')
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
        # Grossman-Larson algebras in the same variables
        # over any base that coerces in:
        if isinstance(R, GrossmanLarsonAlgebra):
            if R.variable_names() == self.variable_names():
                if self.base_ring().has_coerce_map_from(R.base_ring()):
                    return True
        return False
